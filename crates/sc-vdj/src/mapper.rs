use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::reference::VdjReference;
use crate::sequence::{encode_kmer, encode_reference_kmers, reverse_complement};
use crate::types::{Chain, Orientation, SegmentHit, SegmentKind, VdjCandidate, VdjSegment};

#[derive(Debug, Clone)]
pub struct VdjMapperConfig {
    pub k: usize,
    pub min_v_hits: u32,
    pub min_j_hits: u32,
    pub min_c_hits: u32,
}

impl Default for VdjMapperConfig {
    fn default() -> Self {
        Self {
            k: 9,
            min_v_hits: 4,
            min_j_hits: 3,
            min_c_hits: 3,
        }
    }
}

#[derive(Debug, Clone)]
struct SeedOccurrence {
    segment_index: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CoverageGroup {
    chain: Chain,
    kind: SegmentKind,
    /// Relative-position bins in which this kmer occurs within this family.
    /// Only a single-bit mask is trusted as positional evidence downstream.
    bin_mask: u8,
}

#[derive(Debug, Clone)]
struct CoverageSeedEntry {
    segments: Vec<SeedOccurrence>,
    groups: Vec<CoverageGroup>,
    family_count: u8,
    unique_segment_count: u16,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct IdentityOccurrence {
    segment_index: usize,
    /// Relative-position bins for this kmer within this one germline segment.
    bin_mask: u8,
    /// Number of biological segments in the same chain/kind family carrying
    /// this identity kmer. Used to downweight within-family repetition.
    family_multiplicity: u16,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CoverageFamilyScore {
    pub chain: Chain,
    pub kind: SegmentKind,
    pub distinct_kmers: u32,
    pub weighted_score: u32,
    pub bin_mask: u8,
}

impl CoverageFamilyScore {
    pub fn covered_bins(self) -> u32 {
        self.bin_mask.count_ones()
    }

    pub fn has_five_prime(self) -> bool {
        self.bin_mask & 1 != 0
    }

    pub fn has_three_prime(self, bins: u8) -> bool {
        bins > 0 && self.bin_mask & (1u8 << (bins - 1)) != 0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IdentitySeedScore {
    pub segment_index: usize,
    pub distinct_kmers: u32,
    pub weighted_score: u32,
    pub bin_mask: u8,
}

impl IdentitySeedScore {
    pub fn covered_bins(self) -> u32 {
        self.bin_mask.count_ones()
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SeedIndexStats {
    pub k: usize,
    pub concrete_kmers: usize,
    pub unique_kmers: usize,
    pub ambiguous_kmers: usize,
    pub max_segments_per_kmer: usize,
}

#[derive(Debug, Clone, Default)]
struct HitAccumulator {
    score: u32,
    first_query_pos: usize,
    last_query_pos: usize,
}

#[derive(Debug, Clone)]
pub struct VdjMapper {
    reference: VdjReference,
    config: VdjMapperConfig,
    coverage_bins: u8,
    coverage_seeds: HashMap<u64, CoverageSeedEntry>,
    identity_k: usize,
    identity_seeds: HashMap<u64, Vec<IdentityOccurrence>>,
    terminal_k: usize,
    v_terminal_seeds: HashMap<u64, Vec<SeedOccurrence>>,
    j_terminal_seeds: HashMap<u64, Vec<SeedOccurrence>>,
}

impl VdjMapper {
    pub fn new(reference: VdjReference, config: VdjMapperConfig) -> Self {
        assert!((5..=31).contains(&config.k), "VDJ k must be 5..=31");
        let coverage_bins = 6;
        let coverage_seeds = build_coverage_seed_index(&reference, config.k, coverage_bins);
        let stats = coverage_seed_index_stats(&coverage_seeds, config.k);
        let identity_k = 13;
        let identity_seeds = build_identity_seed_index(&reference, identity_k, coverage_bins);
        let terminal_k = 13;
        let v_terminal_seeds = build_terminal_seed_index(&reference, SegmentKind::V, terminal_k, 120);
        let j_terminal_seeds = build_terminal_seed_index(&reference, SegmentKind::J, terminal_k, 120);
        eprintln!(
            "Seed indices:\n  coverage k = {}\n  coverage bins = {}\n  concrete kmers = {}\n  uniquely mapping kmers = {}\n  ambiguous kmers = {}\n  max segments/kmer = {}\n  identity k = {}\n  identity kmers = {}\n  terminal k = {}\n  V-terminal kmers = {}\n  J-terminal kmers = {}",
            stats.k,
            coverage_bins,
            stats.concrete_kmers,
            stats.unique_kmers,
            stats.ambiguous_kmers,
            stats.max_segments_per_kmer,
            identity_k,
            identity_seeds.len(),
            terminal_k,
            v_terminal_seeds.len(),
            j_terminal_seeds.len(),
        );
        Self {
            reference,
            config,
            coverage_bins,
            coverage_seeds,
            identity_k,
            identity_seeds,
            terminal_k,
            v_terminal_seeds,
            j_terminal_seeds,
        }
    }

    pub fn reference(&self) -> &VdjReference {
        &self.reference
    }

    /// Persist the biological reference together with two complementary search
    /// structures: a sensitive family-coverage index and a discriminative
    /// within-family identity index.
    pub fn save_index<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path = path.as_ref();
        let mut w = BufWriter::new(
            File::create(path).with_context(|| format!("creating {}", path.display()))?,
        );
        w.write_all(b"LVDJIDX4")?;
        write_u32(&mut w, 4)?;
        write_u32(&mut w, self.config.k as u32)?;
        write_u32(&mut w, self.config.min_v_hits)?;
        write_u32(&mut w, self.config.min_j_hits)?;
        write_u32(&mut w, self.config.min_c_hits)?;
        write_u32(&mut w, self.reference.segments.len() as u32)?;
        for segment in &self.reference.segments {
            write_segment(&mut w, segment)?;
        }
        write_u8(&mut w, self.coverage_bins)?;
        write_coverage_seed_map(&mut w, &self.coverage_seeds)?;
        write_u32(&mut w, self.identity_k as u32)?;
        write_identity_seed_map(&mut w, &self.identity_seeds)?;
        write_u32(&mut w, self.terminal_k as u32)?;
        write_seed_occurrence_map(&mut w, &self.v_terminal_seeds)?;
        write_seed_occurrence_map(&mut w, &self.j_terminal_seeds)?;
        w.flush()?;
        Ok(())
    }

    /// Load a versioned compiled VDJ index without reparsing GTF/FASTA or
    /// rebuilding either search structure.
    pub fn load_index<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let mut r = BufReader::new(
            File::open(path).with_context(|| format!("opening {}", path.display()))?,
        );
        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)?;
        if &magic == b"LVDJIDX1" || &magic == b"LVDJIDX2" || &magic == b"LVDJIDX3" {
            bail!(
                "{} is an older VDJ index; regenerate it with the current vdj-index so coverage/identity/terminal metadata is available",
                path.display()
            );
        }
        if &magic != b"LVDJIDX4" {
            bail!("{} is not a lumrik VDJ index v4 (bad magic)", path.display());
        }
        let version = read_u32(&mut r)?;
        if version != 4 {
            bail!("unsupported VDJ index version {version} in {}", path.display());
        }
        let config = VdjMapperConfig {
            k: read_u32(&mut r)? as usize,
            min_v_hits: read_u32(&mut r)?,
            min_j_hits: read_u32(&mut r)?,
            min_c_hits: read_u32(&mut r)?,
        };
        if !(5..=31).contains(&config.k) {
            bail!("invalid k={} in {}", config.k, path.display());
        }
        let n_segments = read_u32(&mut r)? as usize;
        let mut segments = Vec::with_capacity(n_segments);
        for _ in 0..n_segments {
            segments.push(read_segment(&mut r)?);
        }
        let coverage_bins = read_u8(&mut r)?;
        if coverage_bins == 0 || coverage_bins > 8 {
            bail!("invalid coverage bin count {coverage_bins} in {}", path.display());
        }
        let coverage_seeds = read_coverage_seed_map(&mut r, n_segments, path)?;
        let identity_k = read_u32(&mut r)? as usize;
        if !(5..=31).contains(&identity_k) {
            bail!("invalid identity k={} in {}", identity_k, path.display());
        }
        let identity_seeds = read_identity_seed_map(&mut r, n_segments, path)?;
        let terminal_k = read_u32(&mut r)? as usize;
        if !(5..=31).contains(&terminal_k) {
            bail!("invalid terminal k={} in {}", terminal_k, path.display());
        }
        let v_terminal_seeds = read_seed_occurrence_map(&mut r, n_segments, path)?;
        let j_terminal_seeds = read_seed_occurrence_map(&mut r, n_segments, path)?;
        let reference = VdjReference { segments };
        let mapper = Self {
            reference,
            config,
            coverage_bins,
            coverage_seeds,
            identity_k,
            identity_seeds,
            terminal_k,
            v_terminal_seeds,
            j_terminal_seeds,
        };
        let stats = mapper.seed_index_stats();
        eprintln!(
            "Loaded VDJ index {}:\n  version = 4\n  segments = {}\n  coverage k = {}\n  coverage bins = {}\n  concrete kmers = {}\n  uniquely mapping kmers = {}\n  ambiguous kmers = {}\n  max segments/kmer = {}\n  identity k = {}\n  identity kmers = {}\n  terminal k = {}\n  V-terminal kmers = {}\n  J-terminal kmers = {}",
            path.display(),
            mapper.reference.len(),
            stats.k,
            mapper.coverage_bins,
            stats.concrete_kmers,
            stats.unique_kmers,
            stats.ambiguous_kmers,
            stats.max_segments_per_kmer,
            mapper.identity_k,
            mapper.identity_seeds.len(),
            mapper.terminal_k,
            mapper.v_terminal_seeds.len(),
            mapper.j_terminal_seeds.len(),
        );
        Ok(mapper)
    }

    pub fn seed_index_stats(&self) -> SeedIndexStats {
        coverage_seed_index_stats(&self.coverage_seeds, self.config.k)
    }

    pub fn coverage_bins(&self) -> u8 {
        self.coverage_bins
    }

    /// Return the best convincing V+J candidate in either read orientation.
    ///
    /// This is a diversion/filtering call, not final clonotype assignment.
    pub fn map(&self, read: &[u8]) -> Option<VdjCandidate> {
        let forward = self.map_oriented(read, Orientation::Forward);
        let reverse = reverse_complement(read);
        let rev = self.map_oriented(&reverse, Orientation::ReverseComplement);

        match (forward, rev) {
            (Some(a), Some(b)) => {
                let sa = a.v.score + a.j.score + a.c.as_ref().map_or(0, |x| x.score);
                let sb = b.v.score + b.j.score + b.c.as_ref().map_or(0, |x| x.score);
                Some(if sa >= sb { a } else { b })
            }
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (None, None) => None,
        }
    }

    pub fn is_vdj(&self, read: &[u8]) -> bool {
        self.map(read).is_some()
    }

    /// Streaming prefilter matching the posterior's raw seed threshold. A read
    /// that cannot give any biological segment `min_hits` seed hits cannot reach
    /// local alignment and therefore cannot affect a rearrangement call.
    pub fn has_seed_candidate(&self, read: &[u8], min_hits: u32) -> bool {
        self.segment_seed_scores_ranked(read)
            .into_iter()
            .any(|(_, hits, _)| hits >= min_hits)
    }

    /// Seed evidence ranked by discriminative power. `hits` preserves the old
    /// biological-segment hit count while `weighted` strongly favors kmers that
    /// occur in few germline segments. A unique seed contributes 1024 points; a
    /// seed shared by N segments contributes roughly 1024/N.
    pub fn segment_seed_scores_ranked(&self, read: &[u8]) -> Vec<(usize, u32, u32)> {
        let reverse = reverse_complement(read);
        self.segment_seed_scores_ranked_oriented(read, &reverse)
    }

    /// Sensitive family-level coverage evidence. Repetitive kmers are retained
    /// deliberately: a seed shared by many IGHV genes is still strong IGH-V
    /// evidence if it does not also occur in other chain/kind families. V, J,
    /// and C are accumulated independently, and each hit records which relative
    /// segment bins were observed so long V segments cannot drown short J/C
    /// segments merely by offering more kmers.
    pub fn coverage_family_scores_ranked_oriented(
        &self,
        read: &[u8],
        reverse: &[u8],
    ) -> Vec<CoverageFamilyScore> {
        let mut totals: HashMap<(Chain, SegmentKind), CoverageFamilyScore> = HashMap::new();
        let mut seen_kmers = HashSet::new();
        for oriented in [read, reverse] {
            if oriented.len() < self.config.k {
                continue;
            }
            for pos in 0..=oriented.len() - self.config.k {
                let Some(kmer) = encode_kmer(&oriented[pos..pos + self.config.k]) else {
                    continue;
                };
                if !seen_kmers.insert(kmer) {
                    continue;
                }
                let Some(entry) = self.coverage_seeds.get(&kmer) else {
                    continue;
                };
                let weight = (1024u32 / entry.family_count.max(1) as u32).max(1);
                for group in &entry.groups {
                    let key = (group.chain, group.kind);
                    let score = totals.entry(key).or_insert(CoverageFamilyScore {
                        chain: group.chain,
                        kind: group.kind,
                        distinct_kmers: 0,
                        weighted_score: 0,
                        bin_mask: 0,
                    });
                    score.distinct_kmers += 1;
                    score.weighted_score = score.weighted_score.saturating_add(weight);
                    if group.bin_mask.count_ones() == 1 {
                        score.bin_mask |= group.bin_mask;
                    }
                }
            }
        }
        let mut out: Vec<_> = totals.into_values().collect();
        out.sort_by(|a, b| {
            b.covered_bins()
                .cmp(&a.covered_bins())
                .then_with(|| b.weighted_score.cmp(&a.weighted_score))
                .then_with(|| b.distinct_kmers.cmp(&a.distinct_kmers))
                .then_with(|| a.chain.cmp(&b.chain))
                .then_with(|| a.kind.cmp(&b.kind))
        });
        out
    }

    /// Dedicated recombination-facing rescue evidence. These 13-mers come
    /// only from the V 3' terminal window or J 5' terminal window and may
    /// nominate a weak junction fragment even when distributed family
    /// coverage is insufficient. They never participate in C or D inference.
    pub fn terminal_seed_scores_ranked_oriented(
        &self,
        read: &[u8],
        reverse: &[u8],
        kind: SegmentKind,
    ) -> Vec<(usize, u32, u32)> {
        let seeds = match kind {
            SegmentKind::V => &self.v_terminal_seeds,
            SegmentKind::J => &self.j_terminal_seeds,
            _ => return Vec::new(),
        };
        ranked_seed_scores(read, reverse, self.terminal_k, seeds)
    }

    /// Discriminative long-seed evidence within one already-recognized
    /// chain/kind family. Common within-family 13-mers remain available but are
    /// downweighted by their biological segment multiplicity. Relative-position
    /// bins are retained and become alignment-window hints downstream.
    pub fn identity_seed_scores_ranked_oriented(
        &self,
        read: &[u8],
        reverse: &[u8],
        chain: Chain,
        kind: SegmentKind,
    ) -> Vec<IdentitySeedScore> {
        let mut totals: HashMap<usize, IdentitySeedScore> = HashMap::new();
        let mut seen_kmers = HashSet::new();
        for oriented in [read, reverse] {
            if oriented.len() < self.identity_k {
                continue;
            }
            for pos in 0..=oriented.len() - self.identity_k {
                let Some(kmer) = encode_kmer(&oriented[pos..pos + self.identity_k]) else {
                    continue;
                };
                if !seen_kmers.insert(kmer) {
                    continue;
                }
                let Some(occurrences) = self.identity_seeds.get(&kmer) else {
                    continue;
                };
                for occurrence in occurrences {
                    let segment = &self.reference.segments[occurrence.segment_index];
                    if segment.chain != chain || segment.kind != kind {
                        continue;
                    }
                    let weight =
                        (4096u32 / occurrence.family_multiplicity.max(1) as u32).max(1);
                    let score = totals
                        .entry(occurrence.segment_index)
                        .or_insert(IdentitySeedScore {
                            segment_index: occurrence.segment_index,
                            distinct_kmers: 0,
                            weighted_score: 0,
                            bin_mask: 0,
                        });
                    score.distinct_kmers += 1;
                    score.weighted_score = score.weighted_score.saturating_add(weight);
                    if occurrence.bin_mask.count_ones() == 1 {
                        score.bin_mask |= occurrence.bin_mask;
                    }
                }
            }
        }
        let last_bin = self.coverage_bins.saturating_sub(1);
        let mut out: Vec<_> = totals.into_values().collect();
        out.sort_by(|a, b| {
            let terminal_a = match kind {
                SegmentKind::V => u8::from(a.bin_mask & (1u8 << last_bin) != 0),
                SegmentKind::J => u8::from(a.bin_mask & 1 != 0),
                _ => 0,
            };
            let terminal_b = match kind {
                SegmentKind::V => u8::from(b.bin_mask & (1u8 << last_bin) != 0),
                SegmentKind::J => u8::from(b.bin_mask & 1 != 0),
                _ => 0,
            };
            terminal_b
                .cmp(&terminal_a)
                .then_with(|| b.covered_bins().cmp(&a.covered_bins()))
                .then_with(|| b.weighted_score.cmp(&a.weighted_score))
                .then_with(|| b.distinct_kmers.cmp(&a.distinct_kmers))
                .then_with(|| a.segment_index.cmp(&b.segment_index))
        });
        out
    }

    pub fn segment_seed_scores_ranked_oriented(
        &self,
        read: &[u8],
        reverse: &[u8],
    ) -> Vec<(usize, u32, u32)> {
        ranked_coverage_segment_scores(read, reverse, self.config.k, &self.coverage_seeds)
    }

    pub fn segment_seed_scores(&self, read: &[u8]) -> Vec<(usize, u32)> {
        let mut totals: HashMap<usize, u32> = HashMap::new();
        for oriented in [read.to_vec(), reverse_complement(read)] {
            if oriented.len() < self.config.k {
                continue;
            }
            let mut seen: HashSet<(usize, u64)> = HashSet::new();
            for pos in 0..=oriented.len() - self.config.k {
                let Some(kmer) = encode_kmer(&oriented[pos..pos + self.config.k]) else {
                    continue;
                };
                if let Some(occ) = self.coverage_seeds.get(&kmer).map(|entry| &entry.segments) {
                    for o in occ {
                        if seen.insert((o.segment_index, kmer)) {
                            *totals.entry(o.segment_index).or_insert(0) += 1;
                        }
                    }
                }
            }
        }
        let mut v: Vec<_> = totals.into_iter().collect();
        v.sort_by_key(|x| std::cmp::Reverse(x.1));
        v
    }

    fn map_oriented(&self, read: &[u8], orientation: Orientation) -> Option<VdjCandidate> {
        if read.len() < self.config.k {
            return None;
        }

        let mut hits: HashMap<usize, HitAccumulator> = HashMap::new();
        // Repetitive query k-mers should contribute once per segment. This keeps
        // low-complexity sequence from manufacturing huge scores.
        let mut seen: HashSet<(usize, u64)> = HashSet::new();

        for query_pos in 0..=read.len() - self.config.k {
            let Some(kmer) = encode_kmer(&read[query_pos..query_pos + self.config.k]) else {
                continue;
            };
            let Some(occurrences) = self.coverage_seeds.get(&kmer).map(|entry| &entry.segments) else {
                continue;
            };

            for occ in occurrences {
                if !seen.insert((occ.segment_index, kmer)) {
                    continue;
                }
                let hit = hits.entry(occ.segment_index).or_insert(HitAccumulator {
                    score: 0,
                    first_query_pos: query_pos,
                    last_query_pos: query_pos,
                });
                hit.score += 1;
                hit.first_query_pos = hit.first_query_pos.min(query_pos);
                hit.last_query_pos = hit.last_query_pos.max(query_pos + self.config.k);
            }
        }

        let mut best: Option<VdjCandidate> = None;
        for chain in [
            Chain::Igh,
            Chain::Igk,
            Chain::Igl,
            Chain::Tra,
            Chain::Trb,
            Chain::Trg,
            Chain::Trd,
        ] {
            let v = best_hit(
                &self.reference,
                &hits,
                chain,
                SegmentKind::V,
                self.config.min_v_hits,
            );
            let j = best_hit(
                &self.reference,
                &hits,
                chain,
                SegmentKind::J,
                self.config.min_j_hits,
            );
            let (Some(v), Some(j)) = (v, j) else { continue };

            // In transcript orientation, useful V evidence should start before J evidence.
            // Allow overlap because short reads and shared seeds can blur boundaries.
            if v.first_query_pos > j.last_query_pos {
                continue;
            }

            let c = best_hit(
                &self.reference,
                &hits,
                chain,
                SegmentKind::C,
                self.config.min_c_hits,
            );

            let candidate = VdjCandidate {
                orientation,
                chain,
                v,
                j,
                c,
            };

            let score =
                candidate.v.score + candidate.j.score + candidate.c.as_ref().map_or(0, |x| x.score);
            let replace = best.as_ref().map_or(true, |old| {
                score > old.v.score + old.j.score + old.c.as_ref().map_or(0, |x| x.score)
            });
            if replace {
                best = Some(candidate);
            }
        }

        best
    }
}

fn relative_bin(pos: usize, sequence_len: usize, k: usize, bins: u8) -> u8 {
    if bins <= 1 || sequence_len <= k {
        return 0;
    }
    let starts = sequence_len - k + 1;
    ((pos.saturating_mul(bins as usize) / starts).min(bins as usize - 1)) as u8
}

fn ranked_coverage_segment_scores(
    read: &[u8],
    reverse: &[u8],
    k: usize,
    seeds: &HashMap<u64, CoverageSeedEntry>,
) -> Vec<(usize, u32, u32)> {
    let mut totals: HashMap<usize, (u32, u32)> = HashMap::new();
    for oriented in [read, reverse] {
        if oriented.len() < k {
            continue;
        }
        let mut seen: HashSet<(usize, u64)> = HashSet::new();
        for pos in 0..=oriented.len() - k {
            let Some(kmer) = encode_kmer(&oriented[pos..pos + k]) else {
                continue;
            };
            let Some(entry) = seeds.get(&kmer) else {
                continue;
            };
            let weight = (1024u32 / entry.unique_segment_count.max(1) as u32).max(1);
            for occurrence in &entry.segments {
                if seen.insert((occurrence.segment_index, kmer)) {
                    let total = totals.entry(occurrence.segment_index).or_insert((0, 0));
                    total.0 += 1;
                    total.1 = total.1.saturating_add(weight);
                }
            }
        }
    }
    let mut out: Vec<_> = totals
        .into_iter()
        .map(|(idx, (hits, weighted))| (idx, hits, weighted))
        .collect();
    out.sort_by(|a, b| {
        b.2.cmp(&a.2)
            .then_with(|| b.1.cmp(&a.1))
            .then_with(|| a.0.cmp(&b.0))
    });
    out
}

fn coverage_seed_index_stats(
    seeds: &HashMap<u64, CoverageSeedEntry>,
    k: usize,
) -> SeedIndexStats {
    let mut stats = SeedIndexStats {
        k,
        concrete_kmers: seeds.len(),
        ..SeedIndexStats::default()
    };
    for entry in seeds.values() {
        let n = entry.unique_segment_count as usize;
        if n == 1 {
            stats.unique_kmers += 1;
        } else if n > 1 {
            stats.ambiguous_kmers += 1;
        }
        stats.max_segments_per_kmer = stats.max_segments_per_kmer.max(n);
    }
    stats
}

fn build_coverage_seed_index(
    reference: &VdjReference,
    k: usize,
    bins: u8,
) -> HashMap<u64, CoverageSeedEntry> {
    let mut segment_sets: HashMap<u64, HashSet<usize>> = HashMap::new();
    let mut family_bins: HashMap<u64, HashMap<(Chain, SegmentKind), u8>> = HashMap::new();

    for (segment_index, segment) in reference.segments.iter().enumerate() {
        if segment.sequence.len() < k {
            continue;
        }
        for pos in 0..=segment.sequence.len() - k {
            let bin = relative_bin(pos, segment.sequence.len(), k, bins);
            for kmer in encode_reference_kmers(&segment.sequence[pos..pos + k], 8) {
                segment_sets.entry(kmer).or_default().insert(segment_index);
                family_bins
                    .entry(kmer)
                    .or_default()
                    .entry((segment.chain, segment.kind))
                    .and_modify(|mask| *mask |= 1u8 << bin)
                    .or_insert(1u8 << bin);
            }
        }
        // `segment_sets` enforces one biological occurrence per concrete kmer
        // regardless of repeats within this germline sequence.
    }

    segment_sets
        .into_iter()
        .map(|(kmer, segment_indices)| {
            let mut segments: Vec<_> = segment_indices
                .into_iter()
                .map(|segment_index| SeedOccurrence { segment_index })
                .collect();
            segments.sort_by_key(|occurrence| occurrence.segment_index);
            let mut groups: Vec<_> = family_bins
                .remove(&kmer)
                .unwrap_or_default()
                .into_iter()
                .map(|((chain, kind), bin_mask)| CoverageGroup {
                    chain,
                    kind,
                    bin_mask,
                })
                .collect();
            groups.sort_by_key(|group| (group.chain, group.kind));
            let family_count = groups.len().min(u8::MAX as usize) as u8;
            let unique_segment_count = segments.len().min(u16::MAX as usize) as u16;
            (
                kmer,
                CoverageSeedEntry {
                    segments,
                    groups,
                    family_count,
                    unique_segment_count,
                },
            )
        })
        .collect()
}

fn build_identity_seed_index(
    reference: &VdjReference,
    k: usize,
    bins: u8,
) -> HashMap<u64, Vec<IdentityOccurrence>> {
    let mut per_kmer: HashMap<u64, HashMap<usize, u8>> = HashMap::new();
    for (segment_index, segment) in reference.segments.iter().enumerate() {
        if segment.sequence.len() < k {
            continue;
        }
        for pos in 0..=segment.sequence.len() - k {
            let bin = relative_bin(pos, segment.sequence.len(), k, bins);
            for kmer in encode_reference_kmers(&segment.sequence[pos..pos + k], 8) {
                per_kmer
                    .entry(kmer)
                    .or_default()
                    .entry(segment_index)
                    .and_modify(|mask| *mask |= 1u8 << bin)
                    .or_insert(1u8 << bin);
            }
        }
    }

    per_kmer
        .into_iter()
        .map(|(kmer, positions)| {
            let mut family_counts: HashMap<(Chain, SegmentKind), u16> = HashMap::new();
            for &segment_index in positions.keys() {
                let segment = &reference.segments[segment_index];
                *family_counts.entry((segment.chain, segment.kind)).or_insert(0) += 1;
            }
            let mut occurrences: Vec<_> = positions
                .into_iter()
                .map(|(segment_index, bin_mask)| {
                    let segment = &reference.segments[segment_index];
                    IdentityOccurrence {
                        segment_index,
                        bin_mask,
                        family_multiplicity: family_counts[&(segment.chain, segment.kind)],
                    }
                })
                .collect();
            occurrences.sort_by_key(|occurrence| occurrence.segment_index);
            (kmer, occurrences)
        })
        .collect()
}

fn build_terminal_seed_index(
    reference: &VdjReference,
    kind: SegmentKind,
    k: usize,
    window: usize,
) -> HashMap<u64, Vec<SeedOccurrence>> {
    let mut seeds: HashMap<u64, Vec<SeedOccurrence>> = HashMap::new();
    for (segment_index, segment) in reference.segments.iter().enumerate() {
        if segment.kind != kind || segment.sequence.len() < k {
            continue;
        }
        let (start, end) = match kind {
            SegmentKind::V => (segment.sequence.len().saturating_sub(window), segment.sequence.len()),
            SegmentKind::J => (0, segment.sequence.len().min(window)),
            _ => continue,
        };
        if end - start < k {
            continue;
        }
        let mut seen = HashSet::new();
        for pos in start..=end - k {
            for kmer in encode_reference_kmers(&segment.sequence[pos..pos + k], 8) {
                if seen.insert(kmer) {
                    seeds.entry(kmer).or_default().push(SeedOccurrence { segment_index });
                }
            }
        }
    }
    seeds
}

fn ranked_seed_scores(
    read: &[u8],
    reverse: &[u8],
    k: usize,
    seeds: &HashMap<u64, Vec<SeedOccurrence>>,
) -> Vec<(usize, u32, u32)> {
    let mut totals: HashMap<usize, (u32, u32)> = HashMap::new();
    let mut seen: HashSet<(usize, u64)> = HashSet::new();
    for oriented in [read, reverse] {
        if oriented.len() < k {
            continue;
        }
        for pos in 0..=oriented.len() - k {
            let Some(kmer) = encode_kmer(&oriented[pos..pos + k]) else { continue; };
            let Some(occurrences) = seeds.get(&kmer) else { continue; };
            let weight = (4096u32 / occurrences.len().max(1) as u32).max(1);
            for occurrence in occurrences {
                if seen.insert((occurrence.segment_index, kmer)) {
                    let total = totals.entry(occurrence.segment_index).or_insert((0, 0));
                    total.0 += 1;
                    total.1 = total.1.saturating_add(weight);
                }
            }
        }
    }
    let mut out: Vec<_> = totals.into_iter().map(|(idx, (hits, score))| (idx, hits, score)).collect();
    out.sort_by(|a, b| b.2.cmp(&a.2).then_with(|| b.1.cmp(&a.1)).then_with(|| a.0.cmp(&b.0)));
    out
}

fn write_seed_occurrence_map<W: Write>(
    w: &mut W,
    seeds: &HashMap<u64, Vec<SeedOccurrence>>,
) -> Result<()> {
    write_u64(w, seeds.len() as u64)?;
    let mut keys: Vec<_> = seeds.keys().copied().collect();
    keys.sort_unstable();
    for key in keys {
        write_u64(w, key)?;
        let occurrences = &seeds[&key];
        write_u32(w, occurrences.len() as u32)?;
        for occurrence in occurrences {
            write_u32(w, occurrence.segment_index as u32)?;
        }
    }
    Ok(())
}

fn read_seed_occurrence_map<R: Read>(
    r: &mut R,
    n_segments: usize,
    path: &Path,
) -> Result<HashMap<u64, Vec<SeedOccurrence>>> {
    let n_seeds = read_u64(r)? as usize;
    let mut seeds = HashMap::with_capacity(n_seeds);
    for _ in 0..n_seeds {
        let key = read_u64(r)?;
        let n = read_u32(r)? as usize;
        let mut occurrences = Vec::with_capacity(n);
        for _ in 0..n {
            let segment_index = read_u32(r)? as usize;
            if segment_index >= n_segments {
                bail!("corrupt VDJ index {}: terminal seed references segment {} of {}", path.display(), segment_index, n_segments);
            }
            occurrences.push(SeedOccurrence { segment_index });
        }
        seeds.insert(key, occurrences);
    }
    Ok(seeds)
}

fn write_coverage_seed_map<W: Write>(
    w: &mut W,
    seeds: &HashMap<u64, CoverageSeedEntry>,
) -> Result<()> {
    write_u64(w, seeds.len() as u64)?;
    let mut keys: Vec<_> = seeds.keys().copied().collect();
    keys.sort_unstable();
    for key in keys {
        write_u64(w, key)?;
        let entry = &seeds[&key];
        write_u16(w, entry.unique_segment_count)?;
        write_u8(w, entry.family_count)?;
        write_u32(w, entry.segments.len() as u32)?;
        for occurrence in &entry.segments {
            write_u32(w, occurrence.segment_index as u32)?;
        }
        write_u32(w, entry.groups.len() as u32)?;
        for group in &entry.groups {
            write_u8(w, chain_code(group.chain))?;
            write_u8(w, kind_code(group.kind))?;
            write_u8(w, group.bin_mask)?;
        }
    }
    Ok(())
}

fn read_coverage_seed_map<R: Read>(
    r: &mut R,
    n_segments: usize,
    path: &Path,
) -> Result<HashMap<u64, CoverageSeedEntry>> {
    let n_seeds = read_u64(r)? as usize;
    let mut seeds = HashMap::with_capacity(n_seeds);
    for _ in 0..n_seeds {
        let key = read_u64(r)?;
        let unique_segment_count = read_u16(r)?;
        let family_count = read_u8(r)?;
        let n_segments_for_seed = read_u32(r)? as usize;
        let mut segments = Vec::with_capacity(n_segments_for_seed);
        for _ in 0..n_segments_for_seed {
            let segment_index = read_u32(r)? as usize;
            if segment_index >= n_segments {
                bail!(
                    "corrupt VDJ index {}: coverage seed references segment {} of {}",
                    path.display(),
                    segment_index,
                    n_segments
                );
            }
            segments.push(SeedOccurrence { segment_index });
        }
        let n_groups = read_u32(r)? as usize;
        let mut groups = Vec::with_capacity(n_groups);
        for _ in 0..n_groups {
            groups.push(CoverageGroup {
                chain: chain_from_code(read_u8(r)?)?,
                kind: kind_from_code(read_u8(r)?)?,
                bin_mask: read_u8(r)?,
            });
        }
        if unique_segment_count as usize != segments.len() {
            bail!(
                "corrupt VDJ index {}: coverage seed segment count mismatch",
                path.display()
            );
        }
        if family_count as usize != groups.len() {
            bail!(
                "corrupt VDJ index {}: coverage seed family count mismatch",
                path.display()
            );
        }
        seeds.insert(
            key,
            CoverageSeedEntry {
                segments,
                groups,
                family_count,
                unique_segment_count,
            },
        );
    }
    Ok(seeds)
}

fn write_identity_seed_map<W: Write>(
    w: &mut W,
    seeds: &HashMap<u64, Vec<IdentityOccurrence>>,
) -> Result<()> {
    write_u64(w, seeds.len() as u64)?;
    let mut keys: Vec<_> = seeds.keys().copied().collect();
    keys.sort_unstable();
    for key in keys {
        write_u64(w, key)?;
        let occurrences = &seeds[&key];
        write_u32(w, occurrences.len() as u32)?;
        for occurrence in occurrences {
            write_u32(w, occurrence.segment_index as u32)?;
            write_u8(w, occurrence.bin_mask)?;
            write_u16(w, occurrence.family_multiplicity)?;
        }
    }
    Ok(())
}

fn read_identity_seed_map<R: Read>(
    r: &mut R,
    n_segments: usize,
    path: &Path,
) -> Result<HashMap<u64, Vec<IdentityOccurrence>>> {
    let n_seeds = read_u64(r)? as usize;
    let mut seeds = HashMap::with_capacity(n_seeds);
    for _ in 0..n_seeds {
        let key = read_u64(r)?;
        let n = read_u32(r)? as usize;
        let mut occurrences = Vec::with_capacity(n);
        for _ in 0..n {
            let segment_index = read_u32(r)? as usize;
            if segment_index >= n_segments {
                bail!(
                    "corrupt VDJ index {}: identity seed references segment {} of {}",
                    path.display(),
                    segment_index,
                    n_segments
                );
            }
            occurrences.push(IdentityOccurrence {
                segment_index,
                bin_mask: read_u8(r)?,
                family_multiplicity: read_u16(r)?,
            });
        }
        seeds.insert(key, occurrences);
    }
    Ok(seeds)
}

fn best_hit(
    reference: &VdjReference,
    hits: &HashMap<usize, HitAccumulator>,
    chain: Chain,
    kind: SegmentKind,
    min_score: u32,
) -> Option<SegmentHit> {
    hits.iter()
        .filter_map(|(&segment_index, hit)| {
            let segment = reference.segments.get(segment_index)?;
            if segment.chain != chain || segment.kind != kind || hit.score < min_score {
                return None;
            }
            Some(SegmentHit {
                segment_index,
                score: hit.score,
                first_query_pos: hit.first_query_pos,
                last_query_pos: hit.last_query_pos,
            })
        })
        .max_by_key(|hit| hit.score)
}

fn write_u8<W: Write>(w: &mut W, v: u8) -> Result<()> { w.write_all(&[v])?; Ok(()) }
fn write_u16<W: Write>(w: &mut W, v: u16) -> Result<()> { w.write_all(&v.to_le_bytes())?; Ok(()) }
fn write_u32<W: Write>(w: &mut W, v: u32) -> Result<()> { w.write_all(&v.to_le_bytes())?; Ok(()) }
fn write_u64<W: Write>(w: &mut W, v: u64) -> Result<()> { w.write_all(&v.to_le_bytes())?; Ok(()) }
fn write_f64<W: Write>(w: &mut W, v: f64) -> Result<()> { w.write_all(&v.to_le_bytes())?; Ok(()) }
fn read_u8<R: Read>(r: &mut R) -> Result<u8> { let mut b=[0;1]; r.read_exact(&mut b)?; Ok(b[0]) }
fn read_u16<R: Read>(r: &mut R) -> Result<u16> { let mut b=[0;2]; r.read_exact(&mut b)?; Ok(u16::from_le_bytes(b)) }
fn read_u32<R: Read>(r: &mut R) -> Result<u32> { let mut b=[0;4]; r.read_exact(&mut b)?; Ok(u32::from_le_bytes(b)) }
fn read_u64<R: Read>(r: &mut R) -> Result<u64> { let mut b=[0;8]; r.read_exact(&mut b)?; Ok(u64::from_le_bytes(b)) }
fn read_f64<R: Read>(r: &mut R) -> Result<f64> { let mut b=[0;8]; r.read_exact(&mut b)?; Ok(f64::from_le_bytes(b)) }
fn write_bytes<W: Write>(w: &mut W, bytes: &[u8]) -> Result<()> { write_u32(w, bytes.len() as u32)?; w.write_all(bytes)?; Ok(()) }
fn read_bytes<R: Read>(r: &mut R) -> Result<Vec<u8>> { let n=read_u32(r)? as usize; if n > 100_000_000 { bail!("unreasonable VDJ index field length {n}"); } let mut out=vec![0;n]; r.read_exact(&mut out)?; Ok(out) }
fn write_string<W: Write>(w: &mut W, s: &str) -> Result<()> { write_bytes(w, s.as_bytes()) }
fn read_string<R: Read>(r: &mut R) -> Result<String> { String::from_utf8(read_bytes(r)?).context("non-UTF8 string in VDJ index") }

fn chain_code(c: Chain) -> u8 { match c { Chain::Igh=>0, Chain::Igk=>1, Chain::Igl=>2, Chain::Tra=>3, Chain::Trb=>4, Chain::Trg=>5, Chain::Trd=>6 } }
fn chain_from_code(v:u8)->Result<Chain> { Ok(match v {0=>Chain::Igh,1=>Chain::Igk,2=>Chain::Igl,3=>Chain::Tra,4=>Chain::Trb,5=>Chain::Trg,6=>Chain::Trd,_=>bail!("invalid chain code {v} in VDJ index")}) }
fn kind_code(k: SegmentKind)->u8 { match k { SegmentKind::V=>0, SegmentKind::D=>1, SegmentKind::J=>2, SegmentKind::C=>3 } }
fn kind_from_code(v:u8)->Result<SegmentKind> { Ok(match v {0=>SegmentKind::V,1=>SegmentKind::D,2=>SegmentKind::J,3=>SegmentKind::C,_=>bail!("invalid segment kind code {v} in VDJ index")}) }

fn write_segment<W: Write>(w:&mut W, s:&VdjSegment)->Result<()> {
    write_string(w,&s.name)?; write_string(w,&s.transcript_id)?; write_string(w,&s.gene_id)?;
    write_u8(w,chain_code(s.chain))?; write_u8(w,kind_code(s.kind))?; write_string(w,&s.chr)?;
    write_u32(w,s.start)?; write_u32(w,s.end)?; write_u8(w,u8::from(s.strand_minus))?;
    write_u64(w,s.locus_rank as u64)?; write_f64(w,s.locus_fraction)?; write_u64(w,s.distance_to_recombination_center)?;
    write_bytes(w,&s.sequence)?; Ok(())
}
fn read_segment<R: Read>(r:&mut R)->Result<VdjSegment> {
    Ok(VdjSegment { name:read_string(r)?, transcript_id:read_string(r)?, gene_id:read_string(r)?, chain:chain_from_code(read_u8(r)?)?, kind:kind_from_code(read_u8(r)?)?, chr:read_string(r)?, start:read_u32(r)?, end:read_u32(r)?, strand_minus:read_u8(r)? != 0, locus_rank:read_u64(r)? as usize, locus_fraction:read_f64(r)?, distance_to_recombination_center:read_u64(r)?, sequence:read_bytes(r)? })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::VdjSegment;

    fn seg(name: &str, kind: SegmentKind, sequence: &[u8]) -> VdjSegment {
        VdjSegment {
            name: name.into(),
            transcript_id: format!("{name}_tx"),
            gene_id: name.into(),
            chain: Chain::Trb,
            kind,
            chr: "chr1".into(),
            start: 0,
            end: sequence.len() as u32,
            strand_minus: false,
            locus_rank: 0,
            locus_fraction: 0.0,
            distance_to_recombination_center: 0,
            sequence: sequence.to_vec(),
        }
    }

    fn seg_chain(name: &str, chain: Chain, kind: SegmentKind, sequence: &[u8]) -> VdjSegment {
        let mut segment = seg(name, kind, sequence);
        segment.chain = chain;
        segment
    }

    #[test]
    fn requires_both_v_and_j() {
        let reference = VdjReference {
            segments: vec![
                seg("TRBV1", SegmentKind::V, b"AACCGGTTAACCGGTT"),
                seg("TRBJ1", SegmentKind::J, b"TTGGAACCTTGGAACC"),
            ],
        };
        let mapper = VdjMapper::new(
            reference,
            VdjMapperConfig {
                k: 5,
                min_v_hits: 2,
                min_j_hits: 2,
                min_c_hits: 2,
            },
        );

        assert!(mapper.is_vdj(b"AACCGGTTAACCGGTTNNNNNTTGGAACCTTGGAACC"));
        assert!(!mapper.is_vdj(b"AACCGGTTAACCGGTTNNNNNNNNNNNNNNNNNNNN"));
    }

    #[test]
    fn detects_reverse_complement_orientation() {
        let reference = VdjReference {
            segments: vec![
                seg("TRBV1", SegmentKind::V, b"AACCGGTTAACCGGTT"),
                seg("TRBJ1", SegmentKind::J, b"TTGGAACCTTGGAACC"),
            ],
        };
        let mapper = VdjMapper::new(
            reference,
            VdjMapperConfig {
                k: 5,
                min_v_hits: 2,
                min_j_hits: 2,
                min_c_hits: 2,
            },
        );

        let forward = b"AACCGGTTAACCGGTTAAAAATTGGAACCTTGGAACC";
        let rc = reverse_complement(forward);
        let hit = mapper.map(&rc).unwrap();
        assert_eq!(hit.orientation, Orientation::ReverseComplement);
    }

    #[test]
    fn compiled_index_roundtrip_preserves_reference_and_seeds() {
        let reference = VdjReference {
            segments: vec![
                seg("TRBV1", SegmentKind::V, b"AACCGGTTAACCGGTT"),
                seg("TRBJ1", SegmentKind::J, b"TTGGAACCTTGGAACC"),
            ],
        };
        let mapper = VdjMapper::new(
            reference,
            VdjMapperConfig { k: 5, min_v_hits: 2, min_j_hits: 2, min_c_hits: 2 },
        );
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.vdjidx");
        mapper.save_index(&path).unwrap();
        let loaded = VdjMapper::load_index(&path).unwrap();
        assert_eq!(mapper.reference().segments, loaded.reference().segments);
        assert_eq!(mapper.seed_index_stats(), loaded.seed_index_stats());
        assert_eq!(mapper.map(b"AACCGGTTAACCGGTTNNNNNTTGGAACCTTGGAACC"), loaded.map(b"AACCGGTTAACCGGTTNNNNNTTGGAACCTTGGAACC"));
    }

    #[test]
    fn compiled_index_rejects_bad_magic() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad.vdjidx");
        std::fs::write(&path, b"NOTINDEX1234").unwrap();
        assert!(VdjMapper::load_index(&path).is_err());
    }

    #[test]
    fn seed_prefilter_respects_posterior_threshold() {
        let reference = VdjReference {
            segments: vec![seg("TRBV1", SegmentKind::V, b"AACCGGTTAACCGGTT")],
        };
        let mapper = VdjMapper::new(
            reference,
            VdjMapperConfig {
                k: 5,
                min_v_hits: 2,
                min_j_hits: 2,
                min_c_hits: 2,
            },
        );

        assert!(mapper.has_seed_candidate(b"AACCGGTT", 2));
        assert!(!mapper.has_seed_candidate(b"AACCG", 2));
    }

    #[test]
    fn shared_seed_is_preserved_for_all_biological_segments() {
        let reference = VdjReference {
            segments: vec![
                seg("TRBV1", SegmentKind::V, b"AACCGGTTAAAA"),
                seg("TRBV2", SegmentKind::V, b"AACCGGTTCCCC"),
                seg("TRBJ1", SegmentKind::J, b"TTGGAACCTTGG"),
            ],
        };
        let mapper = VdjMapper::new(
            reference,
            VdjMapperConfig {
                k: 5,
                min_v_hits: 1,
                min_j_hits: 1,
                min_c_hits: 1,
            },
        );

        let scores = mapper.segment_seed_scores(b"AACCGGTT");
        let ids: HashSet<_> = scores.into_iter().map(|(idx, _)| idx).collect();
        assert!(ids.contains(&0));
        assert!(ids.contains(&1));
        assert!(mapper.seed_index_stats().ambiguous_kmers > 0);
    }
    #[test]
    fn coverage_scoring_separates_segment_kinds_and_rewards_family_repetition() {
        let shared_v = b"ACGTGCACTGATCGTA";
        let reference = VdjReference {
            segments: vec![
                seg_chain("IGHV1", Chain::Igh, SegmentKind::V, shared_v),
                seg_chain("IGHV2", Chain::Igh, SegmentKind::V, shared_v),
                seg_chain("IGHJ1", Chain::Igh, SegmentKind::J, b"TTGGAACCTTGGAACC"),
                seg_chain("IGKV1", Chain::Igk, SegmentKind::V, b"TTTTGGTTAACCAAAA"),
            ],
        };
        let mapper = VdjMapper::new(
            reference,
            VdjMapperConfig {
                k: 5,
                min_v_hits: 1,
                min_j_hits: 1,
                min_c_hits: 1,
            },
        );
        let reverse = reverse_complement(shared_v);
        let scores = mapper.coverage_family_scores_ranked_oriented(shared_v, &reverse);
        assert!(scores.iter().any(|score| {
            score.chain == Chain::Igh && score.kind == SegmentKind::V && score.covered_bins() > 0
        }));
        assert!(!scores.iter().any(|score| {
            score.chain == Chain::Igh && score.kind == SegmentKind::J && score.distinct_kmers > 0
        }));
    }

    #[test]
    fn identity_index_retains_position_bins() {
        let sequence = b"ACGTGCACTGATCGTACGATGCTAGCTACGTTAGCGTACG";
        let reference = VdjReference {
            segments: vec![seg_chain("IGHV1", Chain::Igh, SegmentKind::V, sequence)],
        };
        let mapper = VdjMapper::new(reference, VdjMapperConfig::default());
        let reverse = reverse_complement(sequence);
        let scores = mapper.identity_seed_scores_ranked_oriented(
            sequence,
            &reverse,
            Chain::Igh,
            SegmentKind::V,
        );
        assert_eq!(scores.first().map(|score| score.segment_index), Some(0));
        assert!(scores.first().map(|score| score.covered_bins()).unwrap_or(0) >= 2);
    }

    #[test]
    fn compiled_index_v4_preserves_coverage_identity_and_terminal_metadata() {
        let reference = VdjReference {
            segments: vec![
                seg_chain(
                    "IGHV1",
                    Chain::Igh,
                    SegmentKind::V,
                    b"ACGTGCACTGATCGTACGATGCTAGCTACGTTAGCGTACG",
                ),
                seg_chain(
                    "IGHJ1",
                    Chain::Igh,
                    SegmentKind::J,
                    b"TTGGAACCTTGGAACCTTGGAACCTTGGAACC",
                ),
            ],
        };
        let mapper = VdjMapper::new(reference, VdjMapperConfig::default());
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test-v4.vdjidx");
        mapper.save_index(&path).unwrap();
        let loaded = VdjMapper::load_index(&path).unwrap();

        let read = b"ACGTGCACTGATCGTACGATGCTAGCTACGTT";
        let reverse = reverse_complement(read);
        assert_eq!(
            mapper.coverage_family_scores_ranked_oriented(read, &reverse),
            loaded.coverage_family_scores_ranked_oriented(read, &reverse)
        );
        assert_eq!(
            mapper.identity_seed_scores_ranked_oriented(
                read,
                &reverse,
                Chain::Igh,
                SegmentKind::V,
            ),
            loaded.identity_seed_scores_ranked_oriented(
                read,
                &reverse,
                Chain::Igh,
                SegmentKind::V,
            )
        );
    }

    #[test]
    fn compiled_index_rejects_older_with_regeneration_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("old.vdjidx");
        std::fs::write(&path, b"LVDJIDX2").unwrap();
        let err = VdjMapper::load_index(&path).unwrap_err().to_string();
        assert!(err.contains("regenerate"));
    }

}