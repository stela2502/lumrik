use anyhow::{anyhow, bail, Context, Result};
use flate2::read::MultiGzDecoder;
use gtf_splice_index::{IdNameKeys, SpliceIndex, Strand as GtfStrand};
use std::collections::{BTreeMap, HashMap};
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
    /// Mature transcript-oriented sequence.
    pub sequence: Vec<u8>,
}

#[derive(Debug, Clone)]
pub struct VdjIndex {
    pub segments: Vec<VdjSegment>,
    by_chromosome: HashMap<String, Vec<SegmentId>>,
    local_ordinals: HashMap<(Chain, SegmentKind), Vec<SegmentId>>,
}

impl VdjIndex {
    pub fn from_segments(mut segments: Vec<VdjSegment>) -> Result<Self> {
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
            if blocks.iter().any(|&(a, b)| a < s.end && b > s.start) {
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
}
impl Default for VdjIndexBuilder {
    fn default() -> Self {
        Self {
            max_feature_span: 1_000_000,
        }
    }
}
impl VdjIndexBuilder {
    pub fn new(max_feature_span: u32) -> Self {
        Self { max_feature_span }
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
            let mut start = u32::MAX;
            let mut end = 0;
            for exon in tx.exons() {
                start = start.min(exon.start);
                end = end.max(exon.end);
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
                    sequence: normalize_reference_dna(&sequence)?,
                });
            }
        }
        if segments.is_empty() {
            bail!("no IG/TR V(D)J segments recovered from {}", gtf.display())
        }
        VdjIndex::from_segments(segments)
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
