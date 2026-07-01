use std::collections::HashMap;
use std::fmt;
use std::io::Read;

use flate2::read::GzDecoder;

use crate::error::{PrimerError, PrimerResult};
use crate::single_cell_systems::models::Range;
use crate::single_cell_systems::CellIdGenerator;
use onehot_dna::{OneHot, OneHotSet};

static TENX_3M_FEBRUARY_2018: &[u8] =
    include_bytes!("whitelists/3M-february-2018.txt.gz");

static TENX_3M_FEBRUARY_2018_TRU: &[u8] =
    include_bytes!("whitelists/3M-february-2018_TRU.txt.gz");

static TENX_3M_3PGEX_MAY_2023_TRU: &[u8] =
    include_bytes!("whitelists/3M-3pgex-may-2023_TRU.txt.gz");

static TENX_3M_5PGEX_JAN_2023: &[u8] =
    include_bytes!("whitelists/3M-5pgex-jan-2023.txt.gz");


static TENX_737K_APRIL_2014_RC: &[u8] =
    include_bytes!("whitelists/737k-april-2014_rc.txt.gz");

static TENX_737K_ARC_V1: &[u8] =
    include_bytes!("whitelists/737K-arc-v1.txt.gz");


#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TenxVersion {
    ThreePrimeV1,
    ThreePrimeV2,
    ThreePrimeV3,
    ThreePrimeV4,
    FivePrime,
    MultiomeArcV1,
}

pub struct TenxCoords {
    pub cell: Range,
    pub umi: Range,
    pub consumed: usize,
}

#[derive(Debug, Clone)]
pub struct TenxCellCall {
    pub version: TenxVersion,
    pub cell_id: u64,
    pub cell_seq: Vec<u8>,
    pub cell_qual: Vec<u8>,
    pub umi_seq: Vec<u8>,
    pub umi_qual: Vec<u8>,
    pub consumed: usize,
    pub cell: (usize, usize),
    pub umi: (usize, usize),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TenxWhitelist {
    version: TenxVersion,
    cells: Vec<Vec<u8>>,
    exact: HashMap<Vec<u8>, u64>,
    fuzzy: OneHotSet<16>,
}

impl TenxVersion {
    pub fn parse(raw: &str) -> PrimerResult<Self> {
        match raw.to_ascii_lowercase().as_str() {
            "3pv1" | "3p-v1" | "chromium-single-cell-3-prime-v1" => Ok(Self::ThreePrimeV1),

            "3pv2" | "3p-v2" | "chromium-single-cell-3-prime-v2" => Ok(Self::ThreePrimeV2),

            "3pv3" | "3p-v3" | "chromium-single-cell-3-prime-v3" => Ok(Self::ThreePrimeV3),

            "3pv4" | "3p-v4" | "chromium-single-cell-3-prime-v4" => Ok(Self::ThreePrimeV4),

            "5p" |  "chromium-single-cell-5-prime" => Ok(Self::FivePrime),

            "arc" | "arc-v1" | "chromium-single-cell-multiome-atac-gene-expression" => {
                Ok(Self::MultiomeArcV1)
            }

            other => Err(PrimerError::invalid_grammar(format!(
                "unknown 10x chemistry '{other}'"
            ))),
        }
    }

    pub fn cell_len(self) -> usize {
        16
    }

    pub fn umi_len(self) -> usize {
        match self {
            Self::ThreePrimeV1 => 10,
            Self::ThreePrimeV2 => 10,
            Self::ThreePrimeV3 => 12,
            Self::ThreePrimeV4 => 12,

            Self::FivePrime => 10,

            Self::MultiomeArcV1 => 12,
        }
    }

    pub fn unshifted_consumed_len(self) -> usize {
        self.cell_len() + self.umi_len()
    }

    fn whitelist_gz(self) -> &'static [u8] {
        match self {
            Self::ThreePrimeV1 => TENX_737K_APRIL_2014_RC,
            Self::ThreePrimeV2 => TENX_737K_APRIL_2014_RC,

            Self::ThreePrimeV3 => TENX_3M_FEBRUARY_2018,
            Self::ThreePrimeV4 => TENX_3M_3PGEX_MAY_2023_TRU,

            Self::FivePrime => TENX_3M_5PGEX_JAN_2023,

            Self::MultiomeArcV1 => TENX_737K_ARC_V1,
        }
    }
}

impl fmt::Display for TenxCellCall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "10x {:?} cell_id={} consumed={} cell={}..{} umi={}..{} cell={} umi={}",
            self.version,
            self.cell_id,
            self.consumed,
            self.cell.0,
            self.cell.1,
            self.umi.0,
            self.umi.1,
            String::from_utf8_lossy(&self.cell_seq),
            String::from_utf8_lossy(&self.umi_seq),
        )
    }
}

impl TenxWhitelist {
    pub fn builtin(version: TenxVersion) -> PrimerResult<Self> {
        Self::from_gz_bytes(version, version.whitelist_gz())
    }

    pub fn from_gz_bytes(version: TenxVersion, bytes: &'static [u8]) -> PrimerResult<Self> {
        let mut decoder = GzDecoder::new(bytes);
        let mut text = String::new();

        decoder.read_to_string(&mut text).map_err(|e| {
            PrimerError::invalid_grammar(format!("failed to read 10x whitelist: {e}"))
        })?;

        Ok(Self::from_text(version, &text))
    }

    pub fn from_text(version: TenxVersion, text: &str) -> Self {
        let cells = text
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .map(|line| line.as_bytes().to_vec())
            .collect::<Vec<_>>();

        Self::new(version, cells)
    }

    pub fn new(version: TenxVersion, cells: Vec<Vec<u8>>) -> Self {
        let exact = cells
            .iter()
            .enumerate()
            .map(|(idx, seq)| (seq.clone(), idx as u64))
            .collect();

        let refs = cells
            .iter()
            .map(|v| {
                <&[u8; 16]>::try_from(v.as_slice())
                    .expect("10x whitelist cell must be 16 bp")
            })
            .collect::<Vec<_>>();

        let fuzzy =
            OneHotSet::<16>::from_sequences(&refs)
                .expect("10x whitelist must encode");

        Self {
            version,
            cells,
            exact,
            fuzzy,
        }
    }

    pub fn version(&self) -> TenxVersion {
        self.version
    }

    pub fn len(&self) -> usize {
        self.cells.len()
    }

    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    pub fn index_cell(&self, seq: &[u8]) -> Option<u64> {
        if let Some(idx) = self.exact.get(seq) {
            return Some(*idx);
        }

        let obs = OneHot::<16>::from_bytes(seq).ok()?;
        let (idx, _dist) = self.fuzzy.best_match(&obs, 1)?;

        Some(idx as u64)
    }

    pub fn cell_id_to_seq(&self, cell_id: u64) -> Option<Vec<u8>> {
        if cell_id == 0 {
            return None;
        }

        let idx = (cell_id - 1) as usize;
        self.cells.get(idx).cloned()
    }

    pub fn coords(&self, base: usize) -> Option<TenxCoords> {
        let cell = (base, base + self.version.cell_len());
        let umi = (cell.1, cell.1 + self.version.umi_len());

        Some(TenxCoords {
            cell,
            umi,
            consumed: umi.1,
        })
    }

    pub fn call(&self, seq: &[u8], qual: &[u8], offset: usize) -> Option<TenxCellCall> {
        let coords = self.coords(offset)?;

        let cell = coords.cell;
        let umi = coords.umi;
        let consumed = coords.consumed;

        if seq.len() < umi.1 || qual.len() < umi.1 {
            return None;
        }

        let observed_cell = &seq[cell.0..cell.1];
        let cell_idx = self.index_cell(observed_cell)?;
        let cell_id = cell_idx + 1;

        Some(TenxCellCall {
            version: self.version,
            cell_id,
            cell_seq: self.cells.get(cell_idx as usize)?.clone(),
            cell_qual: qual[cell.0..cell.1].to_vec(),
            umi_seq: seq[umi.0..umi.1].to_vec(),
            umi_qual: qual[umi.0..umi.1].to_vec(),
            consumed,
            cell,
            umi,
        })
    }
}

impl CellIdGenerator for TenxWhitelist {
    fn cell_seq_for_index(&self, allocation_index: u64) -> Option<Vec<u8>> {
        self.cells.get(allocation_index as usize).cloned()
    }

    fn cell_index_for_seq(&self, cell_seq: &[u8]) -> Option<u64> {
        self.exact.get(cell_seq).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn qual(len: usize) -> Vec<u8> {
        vec![b'I'; len]
    }

    #[test]
    fn tenx_from_text_detects_exact_cell() {
        let wl = TenxWhitelist::from_text(
            TenxVersion::ThreePrimeV3,
            "AAACCCAAGAAACACT\nAAACCCAAGAAACCAT\n",
        );

        let mut seq = Vec::new();
        seq.extend_from_slice(b"AAACCCAAGAAACACT");
        seq.extend_from_slice(b"ACGTACGTACGT");
        seq.extend_from_slice(b"TTTTTTTTTTTT");

        let qual = qual(seq.len());

        let call = wl.call(&seq, &qual, 0).expect("10x call failed");

        assert_eq!(call.version, TenxVersion::ThreePrimeV3);
        assert_eq!(call.cell_id, 1);
        assert_eq!(call.cell_seq, b"AAACCCAAGAAACACT".to_vec());
        assert_eq!(call.umi_seq, b"ACGTACGTACGT".to_vec());
        assert_eq!(call.consumed, 28);
    }

    #[test]
    fn tenx_cell_generator_is_zero_based() {
        let wl = TenxWhitelist::from_text(
            TenxVersion::ThreePrimeV3,
            "AAACCCAAGAAACACT\nAAACCCAAGAAACCAT\n",
        );

        assert_eq!(
            wl.cell_seq_for_index(0).unwrap(),
            b"AAACCCAAGAAACACT".to_vec()
        );

        assert_eq!(
            wl.cell_seq_for_index(1).unwrap(),
            b"AAACCCAAGAAACCAT".to_vec()
        );
    }

    #[test]
    fn tenx_cell_id_is_one_based() {
        let wl = TenxWhitelist::from_text(TenxVersion::ThreePrimeV3, "AAACCCAAGAAACACT\n");

        assert!(wl.cell_id_to_seq(0).is_none());

        assert_eq!(wl.cell_id_to_seq(1).unwrap(), b"AAACCCAAGAAACACT".to_vec());
    }
}
