use std::fmt;

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
    pub fn is_bcr(self) -> bool {
        matches!(self, Self::Igh | Self::Igk | Self::Igl)
    }
    pub fn is_tcr(self) -> bool {
        !self.is_bcr()
    }
    pub fn has_d(self) -> bool {
        matches!(self, Self::Igh | Self::Trb | Self::Trd)
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
}
impl fmt::Display for Chain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Orientation {
    Forward,
    ReverseComplement,
}

#[derive(Debug, Clone, PartialEq)]
pub struct VdjSegment {
    pub name: String,
    pub transcript_id: String,
    pub gene_id: String,
    pub chain: Chain,
    pub kind: SegmentKind,
    pub chr: String,
    pub start: u32,
    pub end: u32,
    pub strand_minus: bool,
    /// 0 is nearest the recombination centre among segments of this kind.
    pub locus_rank: usize,
    /// 0.0 is nearest the recombination centre, 1.0 is most distal.
    pub locus_fraction: f64,
    /// Genomic distance in bp to the chain-specific J/DJ recombination centre.
    pub distance_to_recombination_center: u64,
    /// Mature transcript-oriented sequence assembled from annotated exons.
    pub sequence: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SegmentHit {
    pub segment_index: usize,
    pub score: u32,
    pub first_query_pos: usize,
    pub last_query_pos: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VdjCandidate {
    pub orientation: Orientation,
    pub chain: Chain,
    pub v: SegmentHit,
    pub j: SegmentHit,
    pub c: Option<SegmentHit>,
}
