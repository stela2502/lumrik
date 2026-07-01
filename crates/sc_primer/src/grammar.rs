use crate::anchor::AnchorSearch;
use crate::error::{PrimerError, PrimerResult};
use crate::single_cell_systems::*;

use int_to_str::IntToStr;

#[derive(Debug, Clone)]
pub struct Grammar {
    pub name: String,
    pub ops: Vec<GrammarOp>,
    cell_len: usize,
    umi_len: usize,
    system: Option<SingleCellSystem>,
    anchor_search: Option<AnchorSearch>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GrammarOp {
    Fixed { seq: Vec<u8>, mismatches: usize },
    Cell { len: usize },
    Umi { len: usize },
    PolyT { min: usize },
    Insert,
    Skip { len: usize },
    Search { start: usize, end: usize },

    BdCell { version: BdCellVersion },
    TenxCell { version: TenxVersion },

    Tag { len: usize },
    Feature { len: usize },
}

impl Grammar {
    pub fn new(name: impl Into<String>, ops: Vec<GrammarOp>) -> PrimerResult<Self> {
        let mut cell_len = 0usize;
        let mut umi_len = 0usize;
        let mut system: Option<SingleCellSystem> = None;
        let mut anchor_search: Option<AnchorSearch> = None;

        for op in &ops {
            match op {
                GrammarOp::Cell { len } => {
                    cell_len += *len;
                }

                GrammarOp::Umi { len } => {
                    umi_len += *len;
                }

                GrammarOp::BdCell { version } => {
                    cell_len += version.cell_len();
                    umi_len += version.umi_len();

                    system = Some(SingleCellSystem::Rhapsody(match version {
                        BdCellVersion::V1 => RhapsodyWhitelist::bd_v1(),
                        BdCellVersion::V2_96 => RhapsodyWhitelist::bd_v2_96(),
                        BdCellVersion::V2_384 => RhapsodyWhitelist::bd_v2_384(),
                    }));
                }

                GrammarOp::TenxCell { version } => {
                    cell_len += version.cell_len();

                    system = Some(SingleCellSystem::Tenx(TenxWhitelist::builtin(*version)?));
                }

                GrammarOp::Fixed { seq, mismatches } if anchor_search.is_none() => {
                    anchor_search = AnchorSearch::new(seq, *mismatches);
                }

                _ => {}
            }
        }

        Ok(Self {
            name: name.into(),
            ops,
            cell_len,
            umi_len,
            system,
            anchor_search,
        })
    }

    pub fn cell_len(&self) -> usize {
        self.cell_len
    }

    pub fn umi_len(&self) -> usize {
        self.umi_len
    }

    pub fn umi_from_u64(&self, id: u64) -> Vec<u8> {
        IntToStr::from_u64(id).to_string(self.umi_len).into_bytes()
    }

    pub fn anchor_search(&self) -> Option<&AnchorSearch> {
        self.anchor_search.as_ref()
    }

    pub fn system(&self) -> Option<SingleCellSystem> {
        self.system.clone()
    }

    pub fn parse(name: impl Into<String>, structure: &str) -> PrimerResult<Self> {
        let mut ops = Vec::new();

        for raw in structure.split('+') {
            let token = raw.trim();
            if token.is_empty() {
                continue;
            }
            ops.push(GrammarOp::parse_token(token)?);
        }

        if ops.is_empty() {
            return Err(PrimerError::invalid_grammar("empty grammar"));
        }

        Self::new(name, ops)
    }

    /// Create a full primer prefix from the grammar using an externally supplied
    /// cell sequence/cassette and UMI sequence.
    pub fn synthesize(&self, cell_seq: &[u8], umi_seq: &[u8]) -> PrimerResult<Vec<u8>> {
        let mut seq = Vec::new();

        for op in &self.ops {
            match op {
                GrammarOp::Fixed { seq: fixed, .. } => {
                    seq.extend_from_slice(fixed);
                }

                GrammarOp::Cell { len } => {
                    if cell_seq.len() != *len {
                        return Err(PrimerError::invalid_coordinates(
                            "cell sequence length does not match grammar",
                        ));
                    }
                    seq.extend_from_slice(cell_seq);
                }

                GrammarOp::TenxCell { version } => {
                    if cell_seq.len() != version.cell_len() {
                        return Err(PrimerError::invalid_coordinates(
                            "10x cell sequence length does not match grammar",
                        ));
                    }
                    seq.extend_from_slice(cell_seq);
                }

                GrammarOp::BdCell { version } => {
                    if cell_seq.len() != version.cell_len() {
                        return Err(PrimerError::invalid_coordinates(
                            "BD cell cassette length does not match grammar",
                        ));
                    }
                    if umi_seq.len() != version.umi_len() {
                        return Err(PrimerError::invalid_coordinates(
                            "BD UMI length does not match grammar",
                        ));
                    }

                    seq.extend_from_slice(cell_seq);
                    seq.extend_from_slice(umi_seq);
                }

                GrammarOp::Umi { len } => {
                    if umi_seq.len() != *len {
                        return Err(PrimerError::invalid_coordinates(
                            "UMI length does not match grammar",
                        ));
                    }
                    seq.extend_from_slice(umi_seq);
                }

                GrammarOp::PolyT { min } => {
                    seq.extend(std::iter::repeat_n(b'T', *min));
                }

                GrammarOp::Skip { len } => {
                    seq.extend(std::iter::repeat_n(b'A', *len));
                }

                GrammarOp::Search { .. } => {
                    // Deterministic synthetic placeholder for a flexible search window.
                    // This keeps generated primers reproducible while still reserving
                    // space for grammars that contain SEARCH.
                    seq.extend(std::iter::repeat_n(b'A', 4));
                }

                GrammarOp::Insert => break,

                GrammarOp::Tag { len } | GrammarOp::Feature { len } => {
                    seq.extend(std::iter::repeat_n(b'A', *len));
                }
            }
        }

        Ok(seq)
    }
}

impl GrammarOp {
    pub fn parse_token(token: &str) -> PrimerResult<Self> {
        if let Some(rest) = token.strip_prefix("FIXED:") {
            return Self::parse_fixed(rest);
        }
        if let Some(rest) = token.strip_prefix("CELL:") {
            return Self::parse_len(rest).map(|len| Self::Cell { len });
        }
        if let Some(rest) = token.strip_prefix("TENX_CELL:") {
            return TenxVersion::parse(rest).map(|version| Self::TenxCell { version });
        }
        if let Some(rest) = token.strip_prefix("UMI:") {
            return Self::parse_len(rest).map(|len| Self::Umi { len });
        }
        if let Some(rest) = token.strip_prefix("POLYT:") {
            return Self::parse_poly_t(rest);
        }
        if token == "INSERT" {
            return Ok(Self::Insert);
        }
        if let Some(rest) = token.strip_prefix("SKIP:") {
            return Self::parse_len(rest).map(|len| Self::Skip { len });
        }
        if let Some(rest) = token.strip_prefix("SEARCH:") {
            return Self::parse_search(rest);
        }
        if let Some(rest) = token.strip_prefix("BD_CELL:") {
            return BdCellVersion::parse(rest).map(|version| Self::BdCell { version });
        }
        if let Some(rest) = token.strip_prefix("TAG:") {
            return Self::parse_len(rest).map(|len| Self::Tag { len });
        }
        if let Some(rest) = token.strip_prefix("FEATURE:") {
            return Self::parse_len(rest).map(|len| Self::Feature { len });
        }

        Err(PrimerError::invalid_grammar(format!(
            "unknown token '{token}'"
        )))
    }

    pub fn parse_fixed(rest: &str) -> PrimerResult<Self> {
        let mut fields = rest.split(':');
        let seq = fields.next().unwrap_or_default().as_bytes().to_vec();
        if seq.is_empty() {
            return Err(PrimerError::invalid_grammar("FIXED sequence is empty"));
        }

        let mut mismatches = 0usize;
        for field in fields {
            if let Some(mm) = field.strip_prefix("mm=") {
                mismatches = Self::parse_len(mm)?;
            } else if let Some(mm) = field.strip_prefix("mm") {
                mismatches = Self::parse_len(mm)?;
            } else {
                return Err(PrimerError::invalid_grammar(format!(
                    "unknown FIXED option '{field}'"
                )));
            }
        }

        Ok(Self::Fixed { seq, mismatches })
    }

    pub fn parse_poly_t(rest: &str) -> PrimerResult<Self> {
        if let Some(min) = rest.strip_prefix("min=") {
            return Self::parse_len(min).map(|min| Self::PolyT { min });
        }

        Err(PrimerError::invalid_grammar(format!(
            "bad POLYT token '{rest}'"
        )))
    }

    pub fn parse_search(rest: &str) -> PrimerResult<Self> {
        let Some((start, end)) = rest.split_once("..") else {
            return Err(PrimerError::invalid_grammar(format!(
                "bad SEARCH range '{rest}'"
            )));
        };

        let start = Self::parse_len(start)?;
        let end = Self::parse_len(end)?;
        if end < start {
            return Err(PrimerError::invalid_grammar(format!(
                "bad SEARCH range '{rest}'"
            )));
        }

        Ok(Self::Search { start, end })
    }

    pub fn parse_len(raw: &str) -> PrimerResult<usize> {
        raw.parse::<usize>()
            .map_err(|_| PrimerError::invalid_grammar(format!("bad integer '{raw}'")))
    }
}
