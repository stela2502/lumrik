use anyhow::{anyhow, bail, Context, Result};
use flate2::read::MultiGzDecoder;
use gtf_splice_index::{IdNameKeys, SpliceIndex, Strand as GtfStrand};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

mod format;

pub type SegmentId = u16;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum SegmentKind {
    V,
    D,
    J,
    C,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Chain {
    Igh,
    Igk,
    Igl,
    Tra,
    Trb,
    Trg,
    Trd,
}

impl Chain {
    pub const ALL: [Self; 7] = [
        Self::Igh,
        Self::Igk,
        Self::Igl,
        Self::Tra,
        Self::Trb,
        Self::Trg,
        Self::Trd,
    ];
    pub fn has_d(self) -> bool {
        matches!(self, Self::Igh | Self::Trb | Self::Trd)
    }
    pub fn is_bcr(self) -> bool {
        matches!(self, Self::Igh | Self::Igk | Self::Igl)
    }
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Igh => "IGH",
            Self::Igk => "IGK",
            Self::Igl => "IGL",
            Self::Tra => "TRA",
            Self::Trb => "TRB",
            Self::Trg => "TRG",
            Self::Trd => "TRD",
        }
    }
    pub(crate) fn code(self) -> u8 {
        match self {
            Self::Igh => 0,
            Self::Igk => 1,
            Self::Igl => 2,
            Self::Tra => 3,
            Self::Trb => 4,
            Self::Trg => 5,
            Self::Trd => 6,
        }
    }
    pub(crate) fn from_code(x: u8) -> Result<Self> {
        Ok(match x {
            0 => Self::Igh,
            1 => Self::Igk,
            2 => Self::Igl,
            3 => Self::Tra,
            4 => Self::Trb,
            5 => Self::Trg,
            6 => Self::Trd,
            _ => bail!("invalid VDJ chain code {x}"),
        })
    }
}
impl fmt::Display for Chain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Plus,
    Minus,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VdjSegment {
    pub id: SegmentId,
    pub name: String,
    pub transcript_id: String,
    pub gene_id: String,
    pub chain: Chain,
    pub kind: SegmentKind,
    pub chromosome: String,
    pub start: u32,
    pub end: u32,
    pub strand: Strand,
    /// Genomic exon blocks for this transcript. `start..end` is the full
    /// transcript span and may contain large introns, especially for constant
    /// genes; only these blocks are bona-fide exonic segment sequence.
    pub exon_blocks: Vec<(u32, u32)>,
    /// Zero-based coding start in `sequence` for V segments, projected from the
    /// transcript CDS annotation. Other segment kinds do not need to carry the
    /// full CDS geometry.
    pub coding_start: Option<u32>,
    /// Mature transcript-oriented sequence.
    pub sequence: Vec<u8>,
}

impl VdjSegment {
    /// Zero-based coding start in the mature transcript-oriented segment
    /// sequence. `None` is a legitimate property of a V annotation without a
    /// usable CDS anchor; callers must not invent a frame in that case.
    pub fn coding_start(&self) -> Option<usize> {
        self.coding_start.map(|x| x as usize)
    }
}

#[derive(Debug, Clone)]
pub struct VdjIndex {
    pub segments: Vec<VdjSegment>,
    by_chromosome: HashMap<String, Vec<SegmentId>>,
    local_ordinals: HashMap<(Chain, SegmentKind), Vec<SegmentId>>,
    precise_exon_blocks: bool,
}

impl VdjIndex {
    pub fn from_segments(segments: Vec<VdjSegment>) -> Result<Self> {
        Self::from_segments_with_exon_precision(segments, true)
    }

    pub(crate) fn from_segments_with_exon_precision(
        mut segments: Vec<VdjSegment>,
        precise_exon_blocks: bool,
    ) -> Result<Self> {
        segments.sort_by(|a, b| {
            (a.chain, a.kind, &a.name, &a.transcript_id, a.start, a.end).cmp(&(
                b.chain,
                b.kind,
                &b.name,
                &b.transcript_id,
                b.start,
                b.end,
            ))
        });
        segments.dedup_by(|a, b| {
            a.chain == b.chain && a.kind == b.kind && a.name == b.name && a.sequence == b.sequence
        });
        if segments.len() > u16::MAX as usize {
            bail!(
                "VDJ index has {} segments; u16 SegmentId capacity exceeded",
                segments.len()
            );
        }
        for (i, s) in segments.iter_mut().enumerate() {
            s.id = i as SegmentId;
        }
        let mut out = Self {
            segments,
            by_chromosome: HashMap::new(),
            local_ordinals: HashMap::new(),
            precise_exon_blocks,
        };
        out.rebuild_lookups();
        Ok(out)
    }
    fn rebuild_lookups(&mut self) {
        self.by_chromosome.clear();
        self.local_ordinals.clear();
        for s in &self.segments {
            self.by_chromosome
                .entry(s.chromosome.clone())
                .or_default()
                .push(s.id);
            self.local_ordinals
                .entry((s.chain, s.kind))
                .or_default()
                .push(s.id);
        }
        for ids in self.by_chromosome.values_mut() {
            ids.sort_by_key(|id| {
                let s = &self.segments[*id as usize];
                (s.start, s.end, s.kind, s.name.clone())
            });
        }
    }
    pub fn len(&self) -> usize {
        self.segments.len()
    }
    pub fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }
    pub fn has_precise_exon_blocks(&self) -> bool {
        self.precise_exon_blocks
    }
    pub fn has_v_coding_starts(&self) -> bool {
        self.segments
            .iter()
            .filter(|s| s.kind == SegmentKind::V)
            .all(|s| s.coding_start.is_some())
    }
    pub fn segment(&self, id: SegmentId) -> Option<&VdjSegment> {
        self.segments.get(id as usize)
    }
    pub fn counts(&self) -> BTreeMap<(Chain, SegmentKind), usize> {
        let mut x = BTreeMap::new();
        for s in &self.segments {
            *x.entry((s.chain, s.kind)).or_insert(0) += 1;
        }
        x
    }
    pub fn segments_for(
        &self,
        chain: Chain,
        kind: SegmentKind,
    ) -> impl Iterator<Item = &VdjSegment> {
        self.local_ordinals
            .get(&(chain, kind))
            .into_iter()
            .flatten()
            .filter_map(|id| self.segment(*id))
    }
    pub fn overlapping(&self, chromosome: &str, blocks: &[(u32, u32)]) -> Vec<SegmentId> {
        let Some(ids) = self.by_chromosome.get(chromosome) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for id in ids {
            let s = &self.segments[*id as usize];
            if blocks.iter().any(|&(a, b)| {
                s.exon_blocks
                    .iter()
                    .any(|&(x, y)| a < y && b > x)
            }) {
                out.push(*id);
            }
        }
        out
    }

    /// Constant-gene spans overlapped by aligned reference sequence that does
    /// not overlap any annotated exon of that constant transcript. This is
    /// genuine intronic evidence, not ordinary C-exon coverage.
    pub fn intronic_constant_overlapping(
        &self,
        chromosome: &str,
        blocks: &[(u32, u32)],
        min_bases: u32,
    ) -> Vec<SegmentId> {
        let Some(ids) = self.by_chromosome.get(chromosome) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for id in ids {
            let s = &self.segments[*id as usize];
            if s.kind != SegmentKind::C {
                continue;
            }
            let mut intronic = 0u32;
            for &(a, b) in blocks {
                let lo = a.max(s.start);
                let hi = b.min(s.end);
                if lo >= hi {
                    continue;
                }
                let span_overlap = hi - lo;
                let exonic_overlap: u32 = s
                    .exon_blocks
                    .iter()
                    .map(|&(x, y)| {
                        let ex_lo = lo.max(x);
                        let ex_hi = hi.min(y);
                        ex_hi.saturating_sub(ex_lo)
                    })
                    .sum();
                intronic = intronic.saturating_add(span_overlap.saturating_sub(exonic_overlap));
            }
            if intronic >= min_bases {
                out.push(*id);
            }
        }
        out
    }
    pub fn local_ordinal(&self, id: SegmentId) -> Option<usize> {
        let s = self.segment(id)?;
        self.local_ordinals
            .get(&(s.chain, s.kind))?
            .iter()
            .position(|x| *x == id)
    }
    pub fn segment_by_local_ordinal(
        &self,
        chain: Chain,
        kind: SegmentKind,
        ordinal: usize,
    ) -> Option<&VdjSegment> {
        let id = *self.local_ordinals.get(&(chain, kind))?.get(ordinal)?;
        self.segment(id)
    }
    pub fn local_count(&self, chain: Chain, kind: SegmentKind) -> usize {
        self.local_ordinals.get(&(chain, kind)).map_or(0, Vec::len)
    }
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        format::save(self, path.as_ref())
    }
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        format::load(path.as_ref())
    }
}

#[derive(Debug, Clone)]
pub struct VdjIndexBuilder {
    max_feature_span: u32,
    gene_names: Option<HashSet<String>>,
}
impl Default for VdjIndexBuilder {
    fn default() -> Self {
        Self {
            max_feature_span: 1_000_000,
            gene_names: None,
        }
    }
}
impl VdjIndexBuilder {
    pub fn new(max_feature_span: u32) -> Self {
        Self {
            max_feature_span,
            gene_names: None,
        }
    }

    /// Restrict index construction to these gene names. This changes only
    /// which normal V/D/J/C segments are retained; retained genes follow the
    /// exact same construction path as an unrestricted VDJ index.
    pub fn with_gene_names<I, S>(mut self, gene_names: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.gene_names = Some(gene_names.into_iter().map(Into::into).collect());
        self
    }
    pub fn build<P: AsRef<Path>, Q: AsRef<Path>>(
        &self,
        gtf: P,
        genome_fasta: Q,
    ) -> Result<VdjIndex> {
        let gtf = gtf.as_ref();
        let genome_fasta = genome_fasta.as_ref();
        let splice = SpliceIndex::from_path(gtf, self.max_feature_span, IdNameKeys::default())
            .with_context(|| format!("building splice index from {}", gtf.display()))?;
        let immune = discover_immune_transcripts(gtf)?;
        let genome = read_fasta(genome_fasta)?;
        let mut segments = Vec::with_capacity(immune.len());
        for entry in immune.values() {
            if self
                .gene_names
                .as_ref()
                .is_some_and(|wanted| !wanted.contains(&entry.gene_name))
            {
                continue;
            }
            let tx = splice
                .transcript_by_name(&entry.transcript_id)
                .map_err(|e| {
                    anyhow!(
                        "immune transcript {} ({}) missing from splice index: {}",
                        entry.transcript_id,
                        entry.gene_name,
                        e
                    )
                })?;
            let chr_name = splice.chr_names.get(tx.chr_id).ok_or_else(|| {
                anyhow!("invalid chr_id {} for {}", tx.chr_id, entry.transcript_id)
            })?;
            let chr = genome.get(chr_name).ok_or_else(|| {
                anyhow!(
                    "chromosome {} used by {} missing from genome FASTA",
                    chr_name,
                    entry.transcript_id
                )
            })?;
            let mut sequence = Vec::new();
            let mut exon_blocks = Vec::new();
            let mut start = u32::MAX;
            let mut end = 0;
            for exon in tx.exons() {
                start = start.min(exon.start);
                end = end.max(exon.end);
                exon_blocks.push((exon.start, exon.end));
                let a = exon.start as usize;
                let b = exon.end as usize;
                sequence.extend_from_slice(chr.get(a..b).ok_or_else(|| {
                    anyhow!(
                        "exon [{a},{b}) of {} outside {}",
                        entry.transcript_id,
                        chr_name
                    )
                })?);
            }
            let strand = match tx.strand {
                GtfStrand::Plus => Strand::Plus,
                GtfStrand::Minus => {
                    sequence = reverse_complement(&sequence);
                    Strand::Minus
                }
                _ => bail!(
                    "immune transcript {} has unknown strand",
                    entry.transcript_id
                ),
            };
            if !sequence.is_empty() {
                segments.push(VdjSegment {
                    id: 0,
                    name: entry.gene_name.clone(),
                    transcript_id: entry.transcript_id.clone(),
                    gene_id: entry.gene_id.clone(),
                    chain: entry.chain,
                    kind: entry.kind,
                    chromosome: chr_name.clone(),
                    start,
                    end,
                    strand,
                    coding_start: if entry.kind == SegmentKind::V {
                        tx.cds_span().map(|span| project_cds_start(tx.exons(), tx.strand, span))
                    } else {
                        None
                    },
                    exon_blocks,
                    sequence: normalize_reference_dna(&sequence)?,
                });
            }
        }
        if let Some(wanted) = &self.gene_names {
            let found: HashSet<_> = segments.iter().map(|s| s.name.as_str()).collect();
            let mut missing: Vec<_> = wanted
                .iter()
                .filter(|name| !found.contains(name.as_str()))
                .cloned()
                .collect();
            missing.sort();
            if !missing.is_empty() {
                bail!("requested VDJ gene(s) not found in annotation: {}", missing.join(", "));
            }
        }
        if segments.is_empty() {
            bail!("no IG/TR V(D)J segments recovered from {}", gtf.display())
        }
        VdjIndex::from_segments(segments)
    }
}

fn project_cds_start(
    exons: &[gtf_splice_index::RefBlock],
    strand: GtfStrand,
    cds: (u32, u32),
) -> u32 {
    match strand {
        GtfStrand::Plus => exons
            .iter()
            .map(|e| e.end.min(cds.0).saturating_sub(e.start))
            .sum(),
        GtfStrand::Minus => exons
            .iter()
            .map(|e| e.end.saturating_sub(e.start.max(cds.1)))
            .sum(),
        _ => 0,
    }
}

#[derive(Debug, Clone)]
struct ImmuneTranscript {
    transcript_id: String,
    gene_id: String,
    gene_name: String,
    chain: Chain,
    kind: SegmentKind,
}

fn open_maybe_gz(path: &Path) -> Result<Box<dyn Read>> {
    let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    if path.extension().and_then(|x| x.to_str()) == Some("gz") {
        Ok(Box::new(MultiGzDecoder::new(f)))
    } else {
        Ok(Box::new(f))
    }
}
fn read_fasta(path: &Path) -> Result<HashMap<String, Vec<u8>>> {
    let mut out = HashMap::new();
    let mut name = None::<String>;
    let mut seq = Vec::new();
    for line in BufReader::new(open_maybe_gz(path)?).lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(n) = name.take() {
                out.insert(n, normalize_reference_dna(&seq)?);
                seq.clear();
            }
            name = Some(rest.split_whitespace().next().unwrap_or("").to_string());
        } else {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some(n) = name {
        out.insert(n, normalize_reference_dna(&seq)?);
    }
    Ok(out)
}

fn discover_immune_transcripts(path: &Path) -> Result<HashMap<String, ImmuneTranscript>> {
    let mut out = HashMap::new();
    for line in BufReader::new(open_maybe_gz(path)?).lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 || fields[2] != "exon" {
            continue;
        }
        let attrs = parse_attrs(fields[8]);
        let Some(tx) = attrs.get("transcript_id") else {
            continue;
        };
        let gene_id = attrs.get("gene_id").cloned().unwrap_or_default();
        let gene_name = attrs
            .get("gene_name")
            .or_else(|| attrs.get("gene"))
            .cloned()
            .unwrap_or_else(|| gene_id.clone());
        let Some((chain, kind)) = classify_gene(&gene_name) else {
            continue;
        };
        out.entry(tx.clone()).or_insert(ImmuneTranscript {
            transcript_id: tx.clone(),
            gene_id,
            gene_name,
            chain,
            kind,
        });
    }
    Ok(out)
}
fn parse_attrs(s: &str) -> HashMap<String, String> {
    let mut out = HashMap::new();
    for part in s.split(';') {
        let p = part.trim();
        if p.is_empty() {
            continue;
        }
        if let Some((k, v)) = p.split_once(' ') {
            out.insert(k.trim().to_string(), v.trim().trim_matches('"').to_string());
        } else if let Some((k, v)) = p.split_once('=') {
            out.insert(k.trim().to_string(), v.trim().trim_matches('"').to_string());
        }
    }
    out
}
pub fn classify_gene(name: &str) -> Option<(Chain, SegmentKind)> {
    let n = name.to_ascii_uppercase();
    let chain = if n.starts_with("IGH") {
        Chain::Igh
    } else if n.starts_with("IGK") {
        Chain::Igk
    } else if n.starts_with("IGL") {
        Chain::Igl
    } else if n.starts_with("TRA") {
        Chain::Tra
    } else if n.starts_with("TRB") {
        Chain::Trb
    } else if n.starts_with("TRG") {
        Chain::Trg
    } else if n.starts_with("TRD") {
        Chain::Trd
    } else {
        return None;
    };
    let kind = if n == "IGHD" {
        SegmentKind::C
    } else if variable_prefix(&n) {
        SegmentKind::V
    } else if diversity_prefix(&n) {
        SegmentKind::D
    } else if joining_prefix(&n) {
        SegmentKind::J
    } else if constant_prefix(&n) {
        SegmentKind::C
    } else {
        return None;
    };
    Some((chain, kind))
}
fn variable_prefix(n: &str) -> bool {
    ["IGHV", "IGKV", "IGLV", "TRAV", "TRBV", "TRGV", "TRDV"]
        .iter()
        .any(|p| n.starts_with(p))
}
fn diversity_prefix(n: &str) -> bool {
    ["IGHD", "TRBD", "TRDD"].iter().any(|p| n.starts_with(p))
}
fn joining_prefix(n: &str) -> bool {
    ["IGHJ", "IGKJ", "IGLJ", "TRAJ", "TRBJ", "TRGJ", "TRDJ"]
        .iter()
        .any(|p| n.starts_with(p))
}
fn constant_prefix(n: &str) -> bool {
    [
        "TRAC", "TRBC", "TRGC", "TRDC", "IGKC", "IGLC", "IGHM", "IGHG", "IGHA", "IGHE",
    ]
    .iter()
    .any(|p| exact_or_numbered_family(n, p))
}
fn exact_or_numbered_family(n: &str, prefix: &str) -> bool {
    let Some(rest) = n.strip_prefix(prefix) else {
        return false;
    };
    rest.is_empty()
        || (rest.as_bytes().first().is_some_and(|b| b.is_ascii_digit())
            && rest.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'-'))
}

pub(crate) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'D' => b'H',
            b'H' => b'D',
            b'V' => b'B',
            _ => b'N',
        })
        .collect()
}
pub(crate) fn normalize_reference_dna(seq: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(seq.len());
    for &raw in seq {
        let b = raw.to_ascii_uppercase();
        match b {
            b'A' | b'C' | b'G' | b'T' | b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B'
            | b'D' | b'H' | b'V' => out.push(b),
            _ if b.is_ascii_whitespace() => {}
            _ => bail!("unsupported reference base {:?}", b as char),
        }
    }
    Ok(out)
}
pub(crate) fn reference_base_matches(r: u8, q: u8) -> bool {
    let q = q.to_ascii_uppercase();
    match r.to_ascii_uppercase() {
        b'A' => q == b'A',
        b'C' => q == b'C',
        b'G' => q == b'G',
        b'T' => q == b'T',
        b'R' => matches!(q, b'A' | b'G'),
        b'Y' => matches!(q, b'C' | b'T'),
        b'S' => matches!(q, b'C' | b'G'),
        b'W' => matches!(q, b'A' | b'T'),
        b'K' => matches!(q, b'G' | b'T'),
        b'M' => matches!(q, b'A' | b'C'),
        b'B' => matches!(q, b'C' | b'G' | b'T'),
        b'D' => matches!(q, b'A' | b'G' | b'T'),
        b'H' => matches!(q, b'A' | b'C' | b'T'),
        b'V' => matches!(q, b'A' | b'C' | b'G'),
        b'N' => matches!(q, b'A' | b'C' | b'G' | b'T'),
        _ => false,
    }
}

#[cfg(test)]
mod index_contract_tests {
    use super::*;
    use std::fs;
    use tempfile::tempdir;

    fn write_reference() -> Result<(tempfile::TempDir, std::path::PathBuf, std::path::PathBuf)> {
        let dir = tempdir()?;
        let gtf = dir.path().join("tiny.gtf");
        let fasta = dir.path().join("tiny.fa");

        let mut genome = vec![b'A'; 240];
        // Give the V transcript an explicit coding start at genomic position 3
        // (GTF position 4), so the projected mature-transcript coding start is 3.
        genome[3..6].copy_from_slice(b"ATG");
        fs::write(&fasta, format!(">chr1\n{}\n", String::from_utf8(genome).unwrap()))?;
        fs::write(
            &gtf,
            concat!(
                "chr1\ttest\texon\t1\t90\t.\t+\t.\tgene_id \"IGHV1\"; gene_name \"Ighv1-64\"; transcript_id \"IGHV1-T1\";\n",
                "chr1\ttest\tCDS\t4\t90\t.\t+\t0\tgene_id \"IGHV1\"; gene_name \"Ighv1-64\"; transcript_id \"IGHV1-T1\";\n",
                "chr1\ttest\texon\t121\t150\t.\t+\t.\tgene_id \"IGHJ3\"; gene_name \"Ighj3\"; transcript_id \"IGHJ3-T1\";\n",
            ),
        )?;
        Ok((dir, gtf, fasta))
    }

    #[test]
    fn filtered_builder_is_the_same_index_contract_as_full_builder() -> Result<()> {
        let (dir, gtf, fasta) = write_reference()?;
        let full = VdjIndexBuilder::default().build(&gtf, &fasta)?;
        let filtered = VdjIndexBuilder::default()
            .with_gene_names(["Ighv1-64"])
            .build(&gtf, &fasta)?;

        let full_v = full.segments.iter().find(|s| s.name == "Ighv1-64").unwrap();
        let tiny_v = filtered.segments.iter().find(|s| s.name == "Ighv1-64").unwrap();
        assert_eq!(tiny_v.sequence, full_v.sequence);
        assert_eq!(tiny_v.exon_blocks, full_v.exon_blocks);
        assert_eq!(tiny_v.coding_start(), full_v.coding_start());
        assert_eq!(tiny_v.coding_start(), Some(3));

        let path = dir.path().join("filtered.vdjidx");
        filtered.save(&path)?;
        let loaded = VdjIndex::load(&path)?;
        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded.segments[0].name, "Ighv1-64");
        assert_eq!(loaded.segments[0].coding_start(), Some(3));
        assert!(loaded.has_precise_exon_blocks());
        Ok(())
    }

    #[test]
    fn vdj_index_roundtrip_allows_v_without_cds() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("no-cds.vdjidx");
        let index = VdjIndex::from_segments(vec![VdjSegment {
            id: 0,
            name: "IghvPseudo".into(),
            transcript_id: "tx".into(),
            gene_id: "gene".into(),
            chain: Chain::Igh,
            kind: SegmentKind::V,
            chromosome: "chr1".into(),
            start: 0,
            end: 9,
            strand: Strand::Plus,
            exon_blocks: vec![(0, 9)],
            coding_start: None,
            sequence: b"ACGTACGTA".to_vec(),
        }])?;
        index.save(&path)?;
        let loaded = VdjIndex::load(&path)?;
        assert_eq!(loaded.segments[0].coding_start(), None);
        Ok(())
    }
}
