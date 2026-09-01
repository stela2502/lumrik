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
    seeds: HashMap<u64, Vec<SeedOccurrence>>,
}

impl VdjMapper {
    pub fn new(reference: VdjReference, config: VdjMapperConfig) -> Self {
        assert!((5..=31).contains(&config.k), "VDJ k must be 5..=31");
        let seeds = build_seed_index(&reference, config.k);
        let stats = seed_index_stats(&seeds, config.k);
        eprintln!(
            "Seed index:\n  k = {}\n  concrete kmers = {}\n  uniquely mapping kmers = {}\n  ambiguous kmers = {}\n  max segments/kmer = {}",
            stats.k,
            stats.concrete_kmers,
            stats.unique_kmers,
            stats.ambiguous_kmers,
            stats.max_segments_per_kmer
        );
        Self {
            reference,
            config,
            seeds,
        }
    }

    pub fn reference(&self) -> &VdjReference {
        &self.reference
    }

    /// Persist both the biological reference and the already-expanded seed index.
    pub fn save_index<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path = path.as_ref();
        let mut w = BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
        w.write_all(b"LVDJIDX1")?;
        write_u32(&mut w, 1)?;
        write_u32(&mut w, self.config.k as u32)?;
        write_u32(&mut w, self.config.min_v_hits)?;
        write_u32(&mut w, self.config.min_j_hits)?;
        write_u32(&mut w, self.config.min_c_hits)?;
        write_u32(&mut w, self.reference.segments.len() as u32)?;
        for segment in &self.reference.segments {
            write_segment(&mut w, segment)?;
        }
        write_u64(&mut w, self.seeds.len() as u64)?;
        let mut seed_keys: Vec<_> = self.seeds.keys().copied().collect();
        seed_keys.sort_unstable();
        for key in seed_keys {
            write_u64(&mut w, key)?;
            let occurrences = &self.seeds[&key];
            write_u32(&mut w, occurrences.len() as u32)?;
            for occurrence in occurrences {
                write_u32(&mut w, occurrence.segment_index as u32)?;
            }
        }
        w.flush()?;
        Ok(())
    }

    /// Load a versioned compiled VDJ index without reparsing GTF/FASTA or
    /// expanding IUPAC seed kmers.
    pub fn load_index<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let mut r = BufReader::new(File::open(path).with_context(|| format!("opening {}", path.display()))?);
        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)?;
        if &magic != b"LVDJIDX1" {
            bail!("{} is not a lumrik VDJ index (bad magic)", path.display());
        }
        let version = read_u32(&mut r)?;
        if version != 1 {
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
        let n_seeds = read_u64(&mut r)? as usize;
        let mut seeds = HashMap::with_capacity(n_seeds);
        for _ in 0..n_seeds {
            let key = read_u64(&mut r)?;
            let n = read_u32(&mut r)? as usize;
            let mut occurrences = Vec::with_capacity(n);
            for _ in 0..n {
                let segment_index = read_u32(&mut r)? as usize;
                if segment_index >= segments.len() {
                    bail!("corrupt VDJ index {}: seed references segment {} of {}", path.display(), segment_index, segments.len());
                }
                occurrences.push(SeedOccurrence { segment_index });
            }
            seeds.insert(key, occurrences);
        }
        let mapper = Self { reference: VdjReference { segments }, config, seeds };
        let stats = mapper.seed_index_stats();
        eprintln!(
            "Loaded VDJ index {}:\n  segments = {}\n  k = {}\n  concrete kmers = {}\n  uniquely mapping kmers = {}\n  ambiguous kmers = {}\n  max segments/kmer = {}",
            path.display(), mapper.reference.len(), stats.k, stats.concrete_kmers, stats.unique_kmers, stats.ambiguous_kmers, stats.max_segments_per_kmer
        );
        Ok(mapper)
    }

    pub fn seed_index_stats(&self) -> SeedIndexStats {
        seed_index_stats(&self.seeds, self.config.k)
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


    /// Seed evidence ranked by discriminative power. `hits` preserves the old
    /// biological-segment hit count while `weighted` strongly favors kmers that
    /// occur in few germline segments. A unique seed contributes 1024 points; a
    /// seed shared by N segments contributes roughly 1024/N.
    pub fn segment_seed_scores_ranked(&self, read: &[u8]) -> Vec<(usize, u32, u32)> {
        let reverse = reverse_complement(read);
        self.segment_seed_scores_ranked_oriented(read, &reverse)
    }

    pub fn segment_seed_scores_ranked_oriented(
        &self,
        read: &[u8],
        reverse: &[u8],
    ) -> Vec<(usize, u32, u32)> {
        let mut totals: HashMap<usize, (u32, u32)> = HashMap::new();
        for oriented in [read, reverse] {
            if oriented.len() < self.config.k {
                continue;
            }
            let mut seen: HashSet<(usize, u64)> = HashSet::new();
            for pos in 0..=oriented.len() - self.config.k {
                let Some(kmer) = encode_kmer(&oriented[pos..pos + self.config.k]) else {
                    continue;
                };
                let Some(occ) = self.seeds.get(&kmer) else {
                    continue;
                };
                let weight = (1024u32 / occ.len().max(1) as u32).max(1);
                for o in occ {
                    if seen.insert((o.segment_index, kmer)) {
                        let entry = totals.entry(o.segment_index).or_insert((0, 0));
                        entry.0 += 1;
                        entry.1 += weight;
                    }
                }
            }
        }
        let mut out: Vec<_> = totals
            .into_iter()
            .map(|(idx, (hits, weighted))| (idx, hits, weighted))
            .collect();
        out.sort_by(|a, b| b.2.cmp(&a.2).then_with(|| b.1.cmp(&a.1)).then_with(|| a.0.cmp(&b.0)));
        out
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
                if let Some(occ) = self.seeds.get(&kmer) {
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
            let Some(occurrences) = self.seeds.get(&kmer) else {
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

fn seed_index_stats(seeds: &HashMap<u64, Vec<SeedOccurrence>>, k: usize) -> SeedIndexStats {
    let mut stats = SeedIndexStats {
        k,
        concrete_kmers: seeds.len(),
        ..SeedIndexStats::default()
    };
    for occurrences in seeds.values() {
        let n = occurrences.len();
        if n == 1 {
            stats.unique_kmers += 1;
        } else if n > 1 {
            stats.ambiguous_kmers += 1;
        }
        stats.max_segments_per_kmer = stats.max_segments_per_kmer.max(n);
    }
    stats
}

fn build_seed_index(reference: &VdjReference, k: usize) -> HashMap<u64, Vec<SeedOccurrence>> {
    let mut seeds: HashMap<u64, Vec<SeedOccurrence>> = HashMap::new();

    for (segment_index, segment) in reference.segments.iter().enumerate() {
        if segment.sequence.len() < k {
            continue;
        }
        let mut seen = HashSet::new();
        for pos in 0..=segment.sequence.len() - k {
            for kmer in encode_reference_kmers(&segment.sequence[pos..pos + k], 8) {
                if seen.insert(kmer) {
                    seeds
                        .entry(kmer)
                        .or_default()
                        .push(SeedOccurrence { segment_index });
                }
            }
        }
    }

    seeds
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
fn write_u32<W: Write>(w: &mut W, v: u32) -> Result<()> { w.write_all(&v.to_le_bytes())?; Ok(()) }
fn write_u64<W: Write>(w: &mut W, v: u64) -> Result<()> { w.write_all(&v.to_le_bytes())?; Ok(()) }
fn write_f64<W: Write>(w: &mut W, v: f64) -> Result<()> { w.write_all(&v.to_le_bytes())?; Ok(()) }
fn read_u8<R: Read>(r: &mut R) -> Result<u8> { let mut b=[0;1]; r.read_exact(&mut b)?; Ok(b[0]) }
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
}
