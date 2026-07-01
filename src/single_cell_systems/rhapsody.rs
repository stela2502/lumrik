use std::collections::HashMap;
use std::fmt;

use crate::error::{PrimerError, PrimerResult};
use crate::single_cell_systems::models::Range;
use crate::single_cell_systems::traits::CellIdGenerator;

use crate::single_cell_systems::whitelists::bd_const_blocks::{
    BD_V2_384_C1, BD_V2_384_C2, BD_V2_384_C3, BD_V2_96_C1, BD_V2_96_C2, BD_V2_96_C3,
};
use onehot_dna::{OneHot, OneHotSet};

pub struct BdCoords {
    pub c1: Range,
    pub c2: Range,
    pub c3: Range,
    pub umi: Range,
    pub consumed: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BdCellVersion {
    V2_384,
    V2_96,
    V1,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RhapsodyCellCall {
    pub version: BdCellVersion,
    pub cell_id: u64,
    pub cell_seq: Vec<u8>,
    pub cell_qual: Vec<u8>,
    pub umi_seq: Vec<u8>,
    pub umi_qual: Vec<u8>,
    pub shift: usize,
    pub consumed: usize,
    pub c1: (usize, usize),
    pub c2: (usize, usize),
    pub c3: (usize, usize),
    pub umi: (usize, usize),
}

impl fmt::Display for RhapsodyCellCall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "BD {:?} cell_id={} shift={} consumed={} \
             c1={}..{} c2={}..{} c3={}..{} umi={}..{} \
             cell={} umi={}",
            self.version,
            self.cell_id,
            self.shift,
            self.consumed,
            self.c1.0,
            self.c1.1,
            self.c2.0,
            self.c2.1,
            self.c3.0,
            self.c3.1,
            self.umi.0,
            self.umi.1,
            String::from_utf8_lossy(&self.cell_seq),
            String::from_utf8_lossy(&self.umi_seq),
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RhapsodyWhitelist {
    version: BdCellVersion,
    block_size: u64,

    c1: &'static [&'static [u8; 9]],
    c2: &'static [&'static [u8; 9]],
    c3: &'static [&'static [u8; 9]],

    c1_exact: HashMap<Vec<u8>, u64>,
    c2_exact: HashMap<Vec<u8>, u64>,
    c3_exact: HashMap<Vec<u8>, u64>,

    c1_fuzzy: OneHotSet<9>,
    c2_fuzzy: OneHotSet<9>,
    c3_fuzzy: OneHotSet<9>,
}

impl BdCellVersion {
    pub fn parse(raw: &str) -> PrimerResult<Self> {
        match raw {
            "v1" => Ok(Self::V1),
            "v2.96" => Ok(Self::V2_96),
            "v2.384" => Ok(Self::V2_384),
            other => Err(PrimerError::rhapsody(format!(
                "unknown BD cell version '{other}'"
            ))),
        }
    }

    pub fn cell_len(self) -> usize {
        match self {
            Self::V1 => 52,
            Self::V2_96 | Self::V2_384 => 36,
        }
    }

    pub fn block_size(self) -> u64 {
        match self {
            Self::V1 => 96,
            Self::V2_96 => 96,
            Self::V2_384 => 384,
        }
    }

    pub fn umi_len(self) -> usize {
        match self {
            Self::V1 => 8,
            Self::V2_96 | Self::V2_384 => 6,
        }
    }

    pub fn unshifted_consumed_len(self) -> usize {
        match self {
            Self::V1 => 60,
            Self::V2_96 | Self::V2_384 => 42,
        }
    }
}

impl RhapsodyCellCall {
    pub fn empty(version: BdCellVersion) -> Self {
        Self {
            version,
            cell_id: 0,
            cell_seq: Vec::new(),
            cell_qual: Vec::new(),
            umi_seq: Vec::new(),
            umi_qual: Vec::new(),
            shift: 0,
            consumed: 0,
            c1: (0, 0),
            c2: (0, 0),
            c3: (0, 0),
            umi: (0, 0),
        }
    }
}

impl CellIdGenerator for RhapsodyWhitelist {
    fn cell_seq_for_index(&self, allocation_index: u64) -> Option<Vec<u8>> {
        let cell_id = allocation_index.checked_add(1)?;
        self.cell_id_to_cassette(cell_id)
    }

    fn cell_index_for_seq(&self, seq: &[u8]) -> Option<u64> {
        let qual = vec![b'I'; seq.len()];

        let call = self.call(seq, &qual, 0, 0, 0)?.cell_id as u64;

        Some(call - 1)
    }
}

impl RhapsodyWhitelist {
    pub fn new(
        version: BdCellVersion,
        c1s: &'static [&'static [u8; 9]],
        c2s: &'static [&'static [u8; 9]],
        c3s: &'static [&'static [u8; 9]],
    ) -> Self {
        Self {
            version,
            block_size: version.block_size(),

            c1: c1s,
            c2: c2s,
            c3: c3s,

            c1_exact: Self::make_map(c1s),
            c2_exact: Self::make_map(c2s),
            c3_exact: Self::make_map(c3s),

            c1_fuzzy: OneHotSet::<9>::from_sequences(c1s)
                .expect("builtin C1 whitelist must encode"),
            c2_fuzzy: OneHotSet::<9>::from_sequences(c2s)
                .expect("builtin C2 whitelist must encode"),
            c3_fuzzy: OneHotSet::<9>::from_sequences(c3s)
                .expect("builtin C3 whitelist must encode"),
        }
    }

    pub fn cell_len(&self) -> usize {
        self.version.cell_len()
    }

    fn make_map(entries: &[&[u8; 9]]) -> HashMap<Vec<u8>, u64> {
        entries
            .iter()
            .enumerate()
            .map(|(idx, seq)| (seq.to_vec(), idx as u64))
            .collect()
    }

    pub fn builtin(version: BdCellVersion) -> Self {
        match version {
            BdCellVersion::V1 => Self::bd_v1(),
            BdCellVersion::V2_96 => Self::bd_v2_96(),
            BdCellVersion::V2_384 => Self::bd_v2_384(),
        }
    }

    pub fn bd_v1() -> Self {
        Self::new(BdCellVersion::V1, BD_V2_96_C1, BD_V2_96_C2, BD_V2_96_C3)
    }

    pub fn bd_v2_96() -> Self {
        Self::new(BdCellVersion::V2_96, BD_V2_96_C1, BD_V2_96_C2, BD_V2_96_C3)
    }

    pub fn bd_v2_384() -> Self {
        Self::new(
            BdCellVersion::V2_384,
            BD_V2_384_C1,
            BD_V2_384_C2,
            BD_V2_384_C3,
        )
    }

    pub fn version(&self) -> BdCellVersion {
        self.version
    }

    pub fn call(
        &self,
        seq: &[u8],
        qual: &[u8],
        offset: usize,
        shift_start: usize,
        shift_end: usize,
    ) -> Option<RhapsodyCellCall> {
        for shift in shift_start..=shift_end {
            if let Some(call) = self.call_exact_shift(seq, qual, offset, shift) {
                return Some(call);
            }
        }
        None
    }

    pub fn cell_id_to_parts_ids(&self, cell_id: u64) -> Option<(usize, usize, usize)> {
        if cell_id == 0 {
            return None;
        }

        let id = cell_id - 1;
        let bs = self.block_size;

        let c1_idx = (id / (bs * bs)) as usize;

        let rem = id % (bs * bs);

        let c2_idx = (rem / bs) as usize;
        let c3_idx = (rem % bs) as usize;

        Some((c1_idx, c2_idx, c3_idx))
    }

    pub fn cell_id_to_seq(&self, cell_id: u64) -> Option<Vec<u8>> {
        let (c1_idx, c2_idx, c3_idx) = self.cell_id_to_parts_ids(cell_id)?;

        let c1 = self.c1.get(c1_idx)?;
        let c2 = self.c2.get(c2_idx)?;
        let c3 = self.c3.get(c3_idx)?;

        let mut seq = Vec::with_capacity(27);
        seq.extend_from_slice(*c1);
        seq.extend_from_slice(*c2);
        seq.extend_from_slice(*c3);

        Some(seq)
    }

    pub fn cell_id_to_cassette(&self, cell_id: u64) -> Option<Vec<u8>> {
        let (c1_idx, c2_idx, c3_idx) = self.cell_id_to_parts_ids(cell_id)?;

        self.c1.get(c1_idx)?;
        self.c2.get(c2_idx)?;
        self.c3.get(c3_idx)?;

        Some(self.create_cell_cassette(c1_idx, c2_idx, c3_idx))
    }

    pub fn call_exact_shift(
        &self,
        seq: &[u8],
        qual: &[u8],
        offset: usize,
        shift: usize,
    ) -> Option<RhapsodyCellCall> {
        let base = offset.checked_add(shift)?;
        let coords = self.coords(base)?;

        let c1 = coords.c1;
        let c2 = coords.c2;
        let c3 = coords.c3;
        let umi = coords.umi;
        let consumed = coords.consumed;

        if seq.len() < umi.1 || qual.len() < umi.1 {
            return None;
        }

        let c1_exact = Self::index_block_fast(&seq[c1.0..c1.1], &self.c1_exact);

        let c2_exact = Self::index_block_fast(&seq[c2.0..c2.1], &self.c2_exact);

        let c3_exact = Self::index_block_fast(&seq[c3.0..c3.1], &self.c3_exact);

        let missing =
            c1_exact.is_none() as u8 + c2_exact.is_none() as u8 + c3_exact.is_none() as u8;

        let (c1_idx, c2_idx, c3_idx) = match missing {
            0 => (c1_exact.unwrap(), c2_exact.unwrap(), c3_exact.unwrap()),

            1 => (
                match c1_exact {
                    Some(v) => v,
                    None => self.index_c1(&seq[c1.0..c1.1])?,
                },
                match c2_exact {
                    Some(v) => v,
                    None => self.index_c2(&seq[c2.0..c2.1])?,
                },
                match c3_exact {
                    Some(v) => v,
                    None => self.index_c3(&seq[c3.0..c3.1])?,
                },
            ),

            _ => return None,
        };

        let cell_id =
            c1_idx * self.block_size * self.block_size + c2_idx * self.block_size + c3_idx + 1;

        let mut cell_seq = Vec::with_capacity(27);
        let mut cell_qual = Vec::with_capacity(27);

        self.extend_part(
            &mut cell_seq,
            &mut cell_qual,
            seq,
            qual,
            c1,
            Some(self.c1[c1_idx as usize]),
        );
        self.extend_part(
            &mut cell_seq,
            &mut cell_qual,
            seq,
            qual,
            c2,
            Some(self.c2[c2_idx as usize]),
        );
        self.extend_part(
            &mut cell_seq,
            &mut cell_qual,
            seq,
            qual,
            c3,
            Some(self.c3[c3_idx as usize]),
        );

        Some(RhapsodyCellCall {
            version: self.version,
            cell_id,
            cell_seq,
            cell_qual,
            umi_seq: seq[umi.0..umi.1].to_vec(),
            umi_qual: qual[umi.0..umi.1].to_vec(),
            shift,
            consumed,
            c1,
            c2,
            c3,
            umi,
        })
    }

    pub fn expected_id(&self, c1: u64, c2: u64, c3: u64) -> u64 {
        c1 * self.block_size * self.block_size + c2 * self.block_size + c3 + 1
    }

    #[inline]
    fn index_block_fast(seq: &[u8], exact: &HashMap<Vec<u8>, u64>) -> Option<u64> {
        if let Some(idx) = exact.get(seq) {
            return Some(*idx);
        }
        None
    }

    #[inline]
    fn index_block_slow(seq: &[u8], fuzzy: &OneHotSet<9>, max_mismatches: u32) -> Option<u64> {
        let obs = OneHot::<9>::from_bytes(seq).ok()?;
        let (idx, _dist) = fuzzy.best_match(&obs, max_mismatches)?;

        Some(idx as u64)
    }

    pub fn index_c1(&self, seq: &[u8]) -> Option<u64> {
        Self::index_block_slow(seq, &self.c1_fuzzy, 1)
    }

    pub fn index_c2(&self, seq: &[u8]) -> Option<u64> {
        Self::index_block_slow(seq, &self.c2_fuzzy, 1)
    }

    pub fn index_c3(&self, seq: &[u8]) -> Option<u64> {
        Self::index_block_slow(seq, &self.c3_fuzzy, 1)
    }

    pub fn create_cell_cassette(&self, c1_idx: usize, c2_idx: usize, c3_idx: usize) -> Vec<u8> {
        let (c1s, c2s, c3s) = match self.version {
            BdCellVersion::V1 => (BD_V2_96_C1, BD_V2_96_C2, BD_V2_96_C3),
            BdCellVersion::V2_96 => (BD_V2_96_C1, BD_V2_96_C2, BD_V2_96_C3),
            BdCellVersion::V2_384 => (BD_V2_384_C1, BD_V2_384_C2, BD_V2_384_C3),
        };

        let mut seq = Vec::new();

        match self.version {
            BdCellVersion::V1 => {
                seq.extend_from_slice(c1s[c1_idx]);
                seq.extend_from_slice(b"AAAAAAAAAAAA");
                seq.extend_from_slice(c2s[c2_idx]);
                seq.extend_from_slice(b"AAAAAAAAAAAAA");
                seq.extend_from_slice(c3s[c3_idx]);
                seq.extend_from_slice(b"A");
            }

            BdCellVersion::V2_96 | BdCellVersion::V2_384 => {
                seq.extend_from_slice(c1s[c1_idx]);
                seq.extend_from_slice(b"AAAA");
                seq.extend_from_slice(c2s[c2_idx]);
                seq.extend_from_slice(b"AAAA");
                seq.extend_from_slice(c3s[c3_idx]);
                seq.extend_from_slice(b"A");
            }
        }

        seq
    }

    pub fn coords(&self, base: usize) -> Option<BdCoords> {
        match self.version {
            BdCellVersion::V1 => Some(BdCoords {
                c1: (base, base + 9),
                c2: (base + 21, base + 30),
                c3: (base + 43, base + 52),
                umi: (base + 52, base + 60),
                consumed: base + 60,
            }),
            BdCellVersion::V2_96 | BdCellVersion::V2_384 => Some(BdCoords {
                c1: (base, base + 9),
                c2: (base + 13, base + 22),
                c3: (base + 26, base + 35),
                umi: (base + 36, base + 42),
                consumed: base + 42,
            }),
        }
    }

    pub fn extend_part(
        &self,
        cell_seq: &mut Vec<u8>,
        cell_qual: &mut Vec<u8>,
        seq: &[u8],
        qual: &[u8],
        range: (usize, usize),
        corrected: Option<&[u8; 9]>,
    ) {
        match corrected {
            Some(block) => cell_seq.extend_from_slice(block),
            None => cell_seq.extend_from_slice(&seq[range.0..range.1]),
        }

        cell_qual.extend_from_slice(&qual[range.0..range.1]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Chemistry, PrimerDetector};

    fn qual(len: usize) -> Vec<u8> {
        vec![40; len]
    }

    #[test]
    fn bd_v2_384_detects_real_r1_like_read_at_shift_1() {
        let wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);

        let seq = b"TNACGGAGAGATGTGAGCGCCATATGACAGCGGAGCATTGAACCTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let qual = qual(seq.len());

        let call = wl
            .call(seq, &qual, 0, 0, 4)
            .expect("BD v2.384 call should be detected");

        assert_eq!(call.version, BdCellVersion::V2_384);
        assert_eq!(call.shift, 3);
        assert_eq!(
            &seq[call.c1.0..call.c1.1],
            b"CGGAGAGAT",
            "expected CGGAGAGAT from {} to {}",
            call.c1.0,
            call.c1.1
        );
        assert_eq!(
            &seq[call.c2.0..call.c2.1],
            b"GCGCCATAT",
            "expected GCGCCATAT from {} to {}",
            call.c2.0,
            call.c2.1
        );
        assert_eq!(
            &seq[call.c3.0..call.c3.1],
            b"GCGGAGCAT",
            "expected GCGGAGCAT from {} to {}",
            call.c3.0,
            call.c3.1
        );

        assert_eq!(call.cell_id, 45928512,);

        assert_eq!(call.cell_seq, b"CGGAGAGATGCGCCATATGCGGAGCAT".to_vec(),);
    }

    #[test]
    fn bd_v2_384_fails_without_required_shift() {
        let wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);

        let seq = b"TNACGGAGAGATGTGAGCGCCATATGACAGCGGAGCATTGAACCTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let qual = qual(seq.len());

        assert!(
            wl.call(seq, &qual, 0, 0, 0).is_none(),
            "shift 0 should not match this read"
        );
    }

    #[test]
    fn bd_v2_384_detects_real_r1_against_builtin_whitelist() {
        let wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);

        let seq = b"TNACGGAGAGATGTGAGCGCCATATGACAGCGGAGCATTGAACCTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let qual = vec![40; seq.len()];

        let call = wl
            .call(seq, &qual, 0, 0, 4)
            .expect("BD v2.384 builtin whitelist should detect this read");

        eprintln!("shift: {}", call.shift);
        eprintln!("cell_id: {}", call.cell_id);
        eprintln!("cell_seq: {}", String::from_utf8_lossy(&call.cell_seq));
        eprintln!("umi: {}", String::from_utf8_lossy(&call.umi_seq));

        assert_eq!(call.version, BdCellVersion::V2_384);
        assert_eq!(call.umi_seq.len(), 6);
        assert_eq!(call.cell_seq.len(), 27);
    }

    #[test]
    fn bd_v2_384_false_positive_stress_test_detect_all() {
        let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();

        let mut seq = Vec::new();

        // Build a worst-case read consisting entirely of valid
        // whitelist entries but never an intentionally constructed
        // BD primer.

        for i in 0..2000 {
            seq.extend_from_slice(BD_V2_384_C1[i % BD_V2_384_C1.len()]);
            seq.extend_from_slice(BD_V2_384_C2[(i * 7) % BD_V2_384_C2.len()]);
            seq.extend_from_slice(BD_V2_384_C3[(i * 13) % BD_V2_384_C3.len()]);
        }

        let qual = vec![b'I'; seq.len()];

        let hits = detector
            .detect_all(&seq, &qual)
            .expect("detect_all should not fail");

        let top10 = hits
            .iter()
            .take(10)
            .enumerate()
            .map(|(i, h)| {
                let start = h.primer_start.saturating_sub(20);
                let end = (h.primer_end + 20).min(seq.len());

                format!(
                    "{i}: {h}\n    seq={}",
                    String::from_utf8_lossy(&seq[start..end])
                )
            })
            .collect::<Vec<_>>()
            .join("\n");

        assert!(
            hits.is_empty(),
            "false positives detected: n={} top 10:\n{}",
            hits.len(),
            top10,
        );
    }

    #[test]
    fn bd_v2_384_detect_all_multimer_with_mutations() {
        let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();

        let mut seq = Vec::new();
        let mut qual = Vec::new();

        let mut expected_ids = Vec::new();

        for i in 0..100 {
            let cell_id = (i + 1) as u64;
            let umi = format!("{:06}", i)
                .replace('0', "A")
                .replace('1', "C")
                .replace('2', "G")
                .replace('3', "T")
                .into_bytes();

            let wl = RhapsodyWhitelist::bd_v2_384();

            let (c1, c2, c3) = wl
                .cell_id_to_parts_ids(cell_id)
                .expect("Used a wrong cell id - lib error!");

            let cell_seq_primer = wl.create_cell_cassette(c1, c2, c3);

            let mut primer = detector
                .grammar()
                .synthesize(&cell_seq_primer, &umi)
                .expect("Primer creation failed!");

            // Add some insert sequence.
            primer.extend_from_slice(b"GATCGATCGATCGATCGATCGATCGATCG");

            seq.extend_from_slice(&primer);
            qual.extend(std::iter::repeat_n(b'I', primer.len()));

            expected_ids.push(cell_id);
        }

        let hits = detector.detect_all(&seq, &qual).unwrap();

        assert_eq!(
            hits.len(),
            expected_ids.len(),
            "expected {} hits, got {}",
            expected_ids.len(),
            hits.len()
        );

        for (i, hit) in hits.iter().enumerate() {
            assert_eq!(
                hit.bd_cell_id,
                Some(expected_ids[i]),
                "wrong cell id for hit {i}"
            );
        }
    }
}
