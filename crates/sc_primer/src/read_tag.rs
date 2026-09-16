use anyhow::Result;
use std::fmt;
use serde::{Deserialize, Serialize};

use crate::GrammarType;

/// Cell/UMI metadata carried with a mapper-facing molecule.
///
/// sc_primer owns this record because primer detection is also the point where
/// the assay provenance (`GrammarType`) becomes known.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ReadTagRecord {
    pub read_id: String,
    pub original_read_id: Option<String>,
    pub cell_seq: Vec<u8>,
    pub cell_qual: Vec<u8>,
    pub umi_seq: Vec<u8>,
    pub umi_qual: Vec<u8>,
    #[serde(default)]
    pub grammar_type: GrammarType,
}

impl ReadTagRecord {
    /// Backwards-compatible constructor. External/legacy tags are ordinary GEX
    /// unless the caller explicitly supplies provenance.
    pub fn new(
        read_id: String,
        original_read_id: Option<String>,
        cell_seq: impl AsRef<[u8]>,
        cell_qual: impl AsRef<[u8]>,
        umi_seq: impl AsRef<[u8]>,
        umi_qual: impl AsRef<[u8]>,
    ) -> Self {
        Self::new_with_grammar_type(
            read_id, original_read_id, cell_seq, cell_qual, umi_seq, umi_qual,
            GrammarType::Gex,
        )
    }

    pub fn new_with_grammar_type(
        read_id: String,
        original_read_id: Option<String>,
        cell_seq: impl AsRef<[u8]>,
        cell_qual: impl AsRef<[u8]>,
        umi_seq: impl AsRef<[u8]>,
        umi_qual: impl AsRef<[u8]>,
        grammar_type: GrammarType,
    ) -> Self {
        Self {
            read_id,
            original_read_id,
            cell_seq: cell_seq.as_ref().to_vec(),
            cell_qual: cell_qual.as_ref().to_vec(),
            umi_seq: umi_seq.as_ref().to_vec(),
            umi_qual: umi_qual.as_ref().to_vec(),
            grammar_type,
        }
    }

    /// Mapper QNAME format v2:
    /// `read|cell|cellqual|umi|umiqual|grammar|`
    ///
    /// The v1 five-field form is still accepted and is interpreted as GEX.
    pub fn extend_qname(&self, qname: &str) -> String {
        format!(
            "{}|{}|{}|{}|{}|{}|",
            qname,
            hex_encode(&self.cell_seq),
            hex_encode(&self.cell_qual),
            hex_encode(&self.umi_seq),
            hex_encode(&self.umi_qual),
            self.grammar_type.code(),
        )
    }

    pub fn extend_fastq_qnames(&self, q1: &str, q2: &str) -> (String, String) {
        (self.extend_qname(q1), self.extend_qname(q2))
    }

    pub fn from_qname(qname: &str) -> Result<Self> {
        let mut fields = qname.split('|');
        let read_id = fields.next().filter(|id| !id.is_empty())
            .ok_or_else(|| anyhow::anyhow!("QNAME contains an empty read ID: '{qname}'"))?;
        let cell_seq = fields.next().ok_or_else(|| anyhow::anyhow!("QNAME is missing cell sequence: '{qname}'"))?;
        let cell_qual = fields.next().ok_or_else(|| anyhow::anyhow!("QNAME is missing cell quality: '{qname}'"))?;
        let umi_seq = fields.next().ok_or_else(|| anyhow::anyhow!("QNAME is missing UMI sequence: '{qname}'"))?;
        let umi_qual = fields.next().ok_or_else(|| anyhow::anyhow!("QNAME is missing UMI quality: '{qname}'"))?;
        let next = fields.next().ok_or_else(|| anyhow::anyhow!("QNAME metadata is not terminated: '{qname}'"))?;
        let grammar_type = if next.is_empty() {
            GrammarType::Gex
        } else {
            let gt = GrammarType::from_code(next)?;
            match (fields.next(), fields.next()) {
                (Some(""), None) => {}
                _ => anyhow::bail!("QNAME contains malformed ReadTagRecord metadata: '{qname}'"),
            }
            gt
        };
        if next.is_empty() && fields.next().is_some() {
            anyhow::bail!("QNAME contains malformed ReadTagRecord metadata: '{qname}'");
        }
        Ok(Self {
            read_id: read_id.to_string(), original_read_id: None,
            cell_seq: hex_decode(cell_seq)?, cell_qual: hex_decode(cell_qual)?,
            umi_seq: hex_decode(umi_seq)?, umi_qual: hex_decode(umi_qual)?, grammar_type,
        })
    }

    pub fn from_slices(
        read_id: impl Into<String>, original_read_id: Option<String>,
        cell_seq: &[u8], cell_qual: &[u8], umi_seq: &[u8], umi_qual: &[u8],
    ) -> Self {
        Self::new(read_id.into(), original_read_id, cell_seq, cell_qual, umi_seq, umi_qual)
    }

    pub fn from_tsv_fields(
        read_id: impl Into<String>, original_read_id: Option<String>, cell: &str,
        cell_qual: Option<&str>, umi: &str, umi_qual: Option<&str>,
    ) -> Self {
        Self::new(
            read_id.into(), original_read_id, cell.as_bytes(),
            cell_qual.map(phred_from_ascii).unwrap_or_default(), umi.as_bytes(),
            umi_qual.map(phred_from_ascii).unwrap_or_default(),
        )
    }

    pub fn cell_string(&self) -> String { String::from_utf8_lossy(&self.cell_seq).into_owned() }
    pub fn umi_string(&self) -> String { String::from_utf8_lossy(&self.umi_seq).into_owned() }
    pub fn cell_qual_string(&self) -> String { phred_to_ascii(&self.cell_qual) }
    pub fn umi_qual_string(&self) -> String { phred_to_ascii(&self.umi_qual) }
}

fn hex_encode(data: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut out = String::with_capacity(data.len() * 2);
    for &byte in data { out.push(HEX[(byte >> 4) as usize] as char); out.push(HEX[(byte & 0x0f) as usize] as char); }
    out
}
fn hex_decode(s: &str) -> Result<Vec<u8>> {
    if s.len() % 2 != 0 { anyhow::bail!("invalid hex string with odd length"); }
    let mut out = Vec::with_capacity(s.len()/2);
    for pair in s.as_bytes().chunks_exact(2) { out.push(u8::from_str_radix(std::str::from_utf8(pair)?, 16)?); }
    Ok(out)
}

fn phred_from_ascii(text: &str) -> Vec<u8> { text.as_bytes().iter().map(|q| q.saturating_sub(33)).collect() }
fn phred_to_ascii(qual: &[u8]) -> String { qual.iter().map(|q| q.saturating_add(33) as char).collect() }

impl fmt::Display for ReadTagRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "read_id={}, original_read_id={}, cell={}, cell_qual={}, umi={}, umi_qual={}, grammar_type={:?}",
            self.read_id, self.original_read_id.as_deref().unwrap_or("-"), self.cell_string(),
            self.cell_qual_string(), self.umi_string(), self.umi_qual_string(), self.grammar_type)
    }
}
