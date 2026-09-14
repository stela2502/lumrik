use std::fmt;
use std::io::Read;

use flate2::read::GzDecoder;

use crate::error::{PrimerError, PrimerResult};
use crate::single_cell_systems::models::Range;
use crate::single_cell_systems::CellIdGenerator;
use crate::whitelist_hash::WhitelistHash;

static TENX_3M_FEBRUARY_2018: &[u8] = include_bytes!("whitelists/3M-february-2018.txt.gz");

#[allow(unused)]
static TENX_3M_FEBRUARY_2018_TRU: &[u8] = include_bytes!("whitelists/3M-february-2018_TRU.txt.gz");

static TENX_3M_3PGEX_MAY_2023_TRU: &[u8] =
    include_bytes!("whitelists/3M-3pgex-may-2023_TRU.txt.gz");

static TENX_3M_5PGEX_JAN_2023: &[u8] = include_bytes!("whitelists/3M-5pgex-jan-2023.txt.gz");

static TENX_737K_APRIL_2014_RC: &[u8] = include_bytes!("whitelists/737k-april-2014_rc.txt.gz");

static TENX_737K_ARC_V1: &[u8] = include_bytes!("whitelists/737K-arc-v1.txt.gz");

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
    hash: WhitelistHash<16>,
}

impl TenxVersion {
    pub fn parse(raw: &str) -> PrimerResult<Self> {
        match raw.to_ascii_lowercase().as_str() {
            "3pv1" | "3p-v1" | "chromium-single-cell-3-prime-v1" => Ok(Self::ThreePrimeV1),

            "3pv2" | "3p-v2" | "chromium-single-cell-3-prime-v2" => Ok(Self::ThreePrimeV2),

            "3pv3" | "3p-v3" | "chromium-single-cell-3-prime-v3" => Ok(Self::ThreePrimeV3),

            "3pv4" | "3p-v4" | "chromium-single-cell-3-prime-v4" => Ok(Self::ThreePrimeV4),

            "5p" | "chromium-single-cell-5-prime" => Ok(Self::FivePrime),

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
        Self::from_text_with_mismatches(version, text, 1)
            .expect("10x whitelist must contain unique 16 bp A/C/G/T sequences")
    }

    pub fn from_text_with_mismatches(
        version: TenxVersion,
        text: &str,
        max_mismatches: u32,
    ) -> Result<Self, String> {
        let hash = WhitelistHash::<16>::from_text(text, max_mismatches)?;
        Ok(Self { version, hash })
    }

    pub fn new(version: TenxVersion, cells: Vec<Vec<u8>>) -> Self {
        let hash = WhitelistHash::<16>::from_sequences(&cells, 1)
            .expect("10x whitelist cell must be a unique 16 bp A/C/G/T sequence");
        Self { version, hash }
    }

    pub fn version(&self) -> TenxVersion {
        self.version
    }

    pub fn len(&self) -> usize {
        self.hash.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hash.is_empty()
    }

    pub fn index_cell(&self, seq: &[u8]) -> Option<u64> {
        self.hash
            .best_match_default(seq)
            .map(|(index, _distance)| index as u64)
    }

    pub fn cell_id_to_seq(&self, cell_id: u64) -> Option<Vec<u8>> {
        if cell_id == 0 {
            return None;
        }
        self.hash.sequence((cell_id - 1) as usize)
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
            cell_seq: self.hash.sequence(cell_idx as usize)?,
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
        self.hash.sequence(allocation_index as usize)
    }

    fn cell_index_for_seq(&self, cell_seq: &[u8]) -> Option<u64> {
        self.hash
            .best_match_default(cell_seq)
            .map(|(index, _distance)| index as u64)
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
    fn tenx_one_mismatch_uses_same_exact_map() {
        let wl = TenxWhitelist::from_text(
            TenxVersion::ThreePrimeV3,
            "AAACCCAAGAAACACT\nAAACCCAAGAAACCAT\n",
        );

        // One mismatch from the first whitelist barcode.
        assert_eq!(wl.index_cell(b"TAACCCAAGAAACACT"), Some(0));
    }

    #[test]
    fn tenx_one_mismatch_tie_is_rejected() {
        let wl = TenxWhitelist::from_text(
            TenxVersion::ThreePrimeV3,
            "ACGTTGCAACGTTGCA\nACGATGCAACGTTGCA\n",
        );

        assert_eq!(wl.index_cell(b"ACGCTGCAACGTTGCA"), None);
    }

    #[test]
    fn tenx_single_n_is_one_mismatch_and_must_be_unique() {
        let unique = TenxWhitelist::from_text(TenxVersion::ThreePrimeV3, "ACGTTGCAACGTTGCA\n");
        assert_eq!(unique.index_cell(b"NCGTTGCAACGTTGCA"), Some(0));

        let tied = TenxWhitelist::from_text(
            TenxVersion::ThreePrimeV3,
            "ACGTTGCAACGTTGCA\nCCGTTGCAACGTTGCA\n",
        );
        assert_eq!(tied.index_cell(b"NCGTTGCAACGTTGCA"), None);
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
