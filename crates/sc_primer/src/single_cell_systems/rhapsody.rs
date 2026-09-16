use std::collections::HashMap;
use std::fmt;

use onehot_dna::{OneHot, OneHotSequence};

use crate::error::{PrimerError, PrimerResult};
use crate::single_cell_systems::models::Range;
use crate::single_cell_systems::traits::CellIdGenerator;

use crate::single_cell_systems::whitelists::bd_const_blocks::{
    BD_V2_384_C1, BD_V2_384_C2, BD_V2_384_C3, BD_V2_96_C1, BD_V2_96_C2, BD_V2_96_C3,
};
use crate::whitelist_hash::WhitelistHash;

const BD_V2_LINKER_1: &[u8; 4] = b"GTGA";
const BD_V2_LINKER_2: &[u8; 4] = b"GACA";
const BD_V2_LINKERS: &[u8; 8] = b"GTGAGACA";
const BD_V2_VDJ_LINKER_1: &[u8; 4] = b"AATG";
const BD_V2_VDJ_LINKER_2: &[u8; 4] = b"CCAC";
const BD_V2_VDJ_LINKERS: &[u8; 8] = b"AATGCCAC";
const BD_V2_MAX_LINKER_MISMATCHES: u32 = 2;

fn nearest_hamming(query: &[u8], whitelist: &[&'static [u8; 9]]) -> u32 {
    whitelist
        .iter()
        .map(|candidate| {
            query
                .iter()
                .zip(candidate.iter())
                .filter(|(left, right)| left != right)
                .count() as u32
        })
        .min()
        .unwrap_or(query.len() as u32)
}

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
    V2_384Vdj,
    V2_96,
    V1,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RhapsodyCallDiagnostics {
    pub linker_signature: [u8; 8],
    pub linker_mismatches: u8,
    pub c1_mismatches: u8,
    pub c2_mismatches: u8,
    pub c3_mismatches: u8,
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
    pub(crate) diagnostics: Option<RhapsodyCallDiagnostics>,
}

impl RhapsodyCellCall {
    pub fn diagnostics(&self) -> Option<RhapsodyCallDiagnostics> {
        self.diagnostics
    }
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BdMismatchProfile {
    pub shift: usize,
    pub c1_mismatches: u32,
    pub c2_mismatches: u32,
    pub c3_mismatches: u32,
}

impl BdMismatchProfile {
    pub fn total_mismatches(self) -> u32 {
        self.c1_mismatches + self.c2_mismatches + self.c3_mismatches
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

    c1_fuzzy: WhitelistHash<9>,
    c2_fuzzy: WhitelistHash<9>,
    c3_fuzzy: WhitelistHash<9>,
}

impl BdCellVersion {
    pub fn parse(raw: &str) -> PrimerResult<Self> {
        match raw {
            "v1" => Ok(Self::V1),
            "v2.96" => Ok(Self::V2_96),
            "v2.384" => Ok(Self::V2_384),
            "v2.384-vdj" => Ok(Self::V2_384Vdj),
            other => Err(PrimerError::rhapsody(format!(
                "unknown BD cell version '{other}'"
            ))),
        }
    }

    pub fn cell_len(self) -> usize {
        match self {
            Self::V1 => 52,
            Self::V2_96 | Self::V2_384 | Self::V2_384Vdj => 36,
        }
    }

    pub fn block_size(self) -> u64 {
        match self {
            Self::V1 => 96,
            Self::V2_96 => 96,
            Self::V2_384 | Self::V2_384Vdj => 384,
        }
    }

    pub fn umi_len(self) -> usize {
        match self {
            Self::V1 => 8,
            Self::V2_96 | Self::V2_384 | Self::V2_384Vdj => 6,
        }
    }

    pub fn unshifted_consumed_len(self) -> usize {
        match self {
            Self::V1 => 60,
            Self::V2_96 | Self::V2_384 | Self::V2_384Vdj => 42,
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
            diagnostics: None,
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

            c1_fuzzy: WhitelistHash::<9>::from_sequences(c1s, 1)
                .expect("builtin C1 whitelist must encode"),
            c2_fuzzy: WhitelistHash::<9>::from_sequences(c2s, 1)
                .expect("builtin C2 whitelist must encode"),
            c3_fuzzy: WhitelistHash::<9>::from_sequences(c3s, 1)
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
            BdCellVersion::V2_384Vdj => Self::bd_v2_384_vdj(),
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

    pub fn bd_v2_384_vdj() -> Self {
        Self::new(
            BdCellVersion::V2_384Vdj,
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
        if matches!(
            self.version,
            BdCellVersion::V2_96 | BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj
        ) {
            let packed = OneHotSequence::from_bytes(seq);
            return self.call_with_packed(seq, qual, &packed, offset, shift_start, shift_end);
        }

        for shift in shift_start..=shift_end {
            if let Some(call) = self.call_exact_shift(seq, qual, offset, shift) {
                return Some(call);
            }
        }
        None
    }

    /// BD v2 call using a read-sized OneHotSequence prepared by the caller.
    /// This avoids re-packing the same read at every candidate start.
    pub fn call_with_packed(
        &self,
        seq: &[u8],
        qual: &[u8],
        packed: &OneHotSequence,
        offset: usize,
        shift_start: usize,
        shift_end: usize,
    ) -> Option<RhapsodyCellCall> {
        for shift in shift_start..=shift_end {
            if let Some(call) = self.call_exact_shift_packed(seq, qual, packed, offset, shift) {
                return Some(call);
            }
        }
        None
    }

    /// Find the next plausible BD v2 cassette start using only the packed
    /// linker windows. This is the cheap read-wide search primitive shared by
    /// the production detector and diagnostics. Whitelist resolution is not
    /// attempted here.
    pub fn next_candidate_start_packed(
        &self,
        packed: &OneHotSequence,
        from: usize,
        shift_start: usize,
        shift_end: usize,
    ) -> Option<usize> {
        if !matches!(
            self.version,
            BdCellVersion::V2_96 | BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj
        ) {
            return None;
        }

        let last = packed
            .len()
            .checked_sub(self.version.unshifted_consumed_len())?;
        for outer in from..=last {
            for shift in shift_start..=shift_end {
                let Some(base) = outer.checked_add(shift) else {
                    continue;
                };
                let Some(distance) = self.v2_linker_mismatches_at(packed, base) else {
                    continue;
                };
                if distance <= BD_V2_MAX_LINKER_MISMATCHES {
                    return Some(outer);
                }
            }
        }
        None
    }

    /// Scan a complete BD v2 read with the cheap linker gate and only resolve
    /// C1/C2/C3 for offsets that pass it. The read is packed once.
    pub fn scan_first(
        &self,
        seq: &[u8],
        qual: &[u8],
        shift_start: usize,
        shift_end: usize,
    ) -> Option<RhapsodyCellCall> {
        if seq.len() != qual.len() {
            return None;
        }
        if !matches!(
            self.version,
            BdCellVersion::V2_96 | BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj
        ) {
            return self.call(seq, qual, 0, shift_start, shift_end);
        }

        let packed = OneHotSequence::from_bytes(seq);
        let mut cursor = 0usize;
        while let Some(offset) =
            self.next_candidate_start_packed(&packed, cursor, shift_start, shift_end)
        {
            if let Some(call) =
                self.call_with_packed(seq, qual, &packed, offset, shift_start, shift_end)
            {
                return Some(call);
            }
            cursor = offset.saturating_add(1);
        }
        None
    }

    /// Return the next outer grammar start for BD v2 without using linker or
    /// exact-whitelist gates. The grammar itself will try the configured
    /// SEARCH shifts and the fuzzy whitelist matcher decides whether a cell is
    /// valid.
    pub fn next_candidate_start(
        &self,
        seq: &[u8],
        from: usize,
        _shift_start: usize,
        _shift_end: usize,
    ) -> Option<usize> {
        if !matches!(
            self.version,
            BdCellVersion::V2_96 | BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj
        ) {
            return None;
        }

        (from < seq.len()).then_some(from)
    }

    pub fn explain_call_failure(
        &self,
        seq: &[u8],
        offset: usize,
        shift_start: usize,
        shift_end: usize,
    ) -> String {
        if !matches!(
            self.version,
            BdCellVersion::V2_96 | BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj
        ) {
            return "BD_CELL: no complete whitelist match".to_string();
        }

        let mut saw_full_length = false;
        let mut whitelist_reason = None;

        for shift in shift_start..=shift_end {
            let Some(base) = offset.checked_add(shift) else {
                continue;
            };
            let Some(coords) = self.coords(base) else {
                continue;
            };

            if seq.len() < coords.umi.1 {
                continue;
            }

            saw_full_length = true;

            let c1 = &seq[coords.c1.0..coords.c1.1];
            let c2 = &seq[coords.c2.0..coords.c2.1];
            let c3 = &seq[coords.c3.0..coords.c3.1];

            if !self.v2_linkers_pass(seq, &coords) {
                whitelist_reason =
                    Some("BD_CELL: combined 8 bp linker has more than two mismatches".to_string());
                continue;
            }

            if self.index_c1_unique(c1).is_none() {
                whitelist_reason = Some(
                    "BD_CELL: C1 unique-nearest whitelist assignment was ambiguous".to_string(),
                );
                continue;
            }
            if self.index_c2_unique(c2).is_none() {
                whitelist_reason = Some(
                    "BD_CELL: C2 unique-nearest whitelist assignment was ambiguous".to_string(),
                );
                continue;
            }
            if self.index_c3_unique(c3).is_none() {
                whitelist_reason = Some(
                    "BD_CELL: C3 unique-nearest whitelist assignment was ambiguous".to_string(),
                );
                continue;
            }

            return "BD_CELL: cassette passed fuzzy whitelist checks but the full grammar failed later"
                .to_string();
        }

        if !saw_full_length {
            return "BD_CELL: sequence is too short for the cassette and UMI in the SEARCH window"
                .to_string();
        }

        whitelist_reason.unwrap_or_else(|| "BD_CELL: barcode whitelist match failed".to_string())
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

    /// Convert a canonical 27-base (C1+C2+C3) barcode sequence to the official
    /// one-based positional BD/Rustody cell id. This is intentionally exact:
    /// upstream primer calling is responsible for barcode correction.
    pub fn cell_id_for_seq(&self, seq: &[u8]) -> Option<u64> {
        if seq.len() != 27 {
            return None;
        }
        let c1 = *self.c1_exact.get(&seq[0..9])?;
        let c2 = *self.c2_exact.get(&seq[9..18])?;
        let c3 = *self.c3_exact.get(&seq[18..27])?;
        Some(c1 * self.block_size * self.block_size + c2 * self.block_size + c3 + 1)
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

    /// For diagnostics only: find the closest possible BD v2 barcode blocks
    /// at any allowed SEARCH shift, without applying the production mismatch
    /// threshold or requiring a unique whitelist hit.
    pub fn best_mismatch_profile(
        &self,
        seq: &[u8],
        offset: usize,
        shift_start: usize,
        shift_end: usize,
    ) -> Option<BdMismatchProfile> {
        if !matches!(
            self.version,
            BdCellVersion::V2_96 | BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj
        ) {
            return None;
        }

        let mut best = None;
        for shift in shift_start..=shift_end {
            let base = offset.checked_add(shift)?;
            let Some(coords) = self.coords(base) else {
                continue;
            };
            if seq.len() < coords.c3.1 {
                continue;
            }

            let profile = BdMismatchProfile {
                shift,
                c1_mismatches: nearest_hamming(&seq[coords.c1.0..coords.c1.1], &self.c1),
                c2_mismatches: nearest_hamming(&seq[coords.c2.0..coords.c2.1], &self.c2),
                c3_mismatches: nearest_hamming(&seq[coords.c3.0..coords.c3.1], &self.c3),
            };

            if best.is_none_or(|current: BdMismatchProfile| {
                profile.total_mismatches() < current.total_mismatches()
                    || (profile.total_mismatches() == current.total_mismatches()
                        && profile.shift < current.shift)
            }) {
                best = Some(profile);
            }
        }
        best
    }

    fn v2_call_diagnostics(
        &self,
        seq: &[u8],
        coords: &BdCoords,
        c1_idx: u64,
        c2_idx: u64,
        c3_idx: u64,
    ) -> Option<RhapsodyCallDiagnostics> {
        let linker1 = seq.get(coords.c1.1..coords.c2.0)?;
        let linker2 = seq.get(coords.c2.1..coords.c3.0)?;
        if linker1.len() != 4 || linker2.len() != 4 {
            return None;
        }

        let mut linker_signature = [0u8; 8];
        linker_signature[..4].copy_from_slice(linker1);
        linker_signature[4..].copy_from_slice(linker2);

        let observed_linker = OneHot::<8>::from_bytes(&linker_signature).ok()?;
        let expected_linker = OneHot::<8>::from_bytes(self.v2_linkers()).ok()?;
        let c1 = OneHot::<9>::from_bytes(seq.get(coords.c1.0..coords.c1.1)?).ok()?;
        let c2 = OneHot::<9>::from_bytes(seq.get(coords.c2.0..coords.c2.1)?).ok()?;
        let c3 = OneHot::<9>::from_bytes(seq.get(coords.c3.0..coords.c3.1)?).ok()?;
        let c1_expected = OneHot::<9>::from_bytes(*self.c1.get(c1_idx as usize)?).ok()?;
        let c2_expected = OneHot::<9>::from_bytes(*self.c2.get(c2_idx as usize)?).ok()?;
        let c3_expected = OneHot::<9>::from_bytes(*self.c3.get(c3_idx as usize)?).ok()?;

        Some(RhapsodyCallDiagnostics {
            linker_signature,
            linker_mismatches: observed_linker.mismatches(expected_linker) as u8,
            c1_mismatches: c1.mismatches(c1_expected) as u8,
            c2_mismatches: c2.mismatches(c2_expected) as u8,
            c3_mismatches: c3.mismatches(c3_expected) as u8,
        })
    }

    #[inline]
    fn v2_linker_mismatches_at(&self, packed: &OneHotSequence, base: usize) -> Option<u32> {
        let linker1 = packed.window::<4>(base.checked_add(9)?).ok()?;
        let linker2 = packed.window::<4>(base.checked_add(22)?).ok()?;
        let expected1 = OneHot::<4>::from_bytes(self.v2_linker_1()).ok()?;
        let expected2 = OneHot::<4>::from_bytes(self.v2_linker_2()).ok()?;
        Some(linker1.mismatches(expected1) + linker2.mismatches(expected2))
    }

    #[inline]
    fn index_block_unique_onehot(
        query: OneHot<9>,
        whitelist: &WhitelistHash<9>,
    ) -> Option<(u64, u32)> {
        let (idx, distance) = whitelist.unique_nearest_routed_onehot(query)?;
        Some((idx as u64, distance))
    }

    fn call_exact_shift_packed(
        &self,
        seq: &[u8],
        qual: &[u8],
        packed: &OneHotSequence,
        offset: usize,
        shift: usize,
    ) -> Option<RhapsodyCellCall> {
        let base = offset.checked_add(shift)?;
        let coords = self.coords(base)?;
        let consumed = coords.consumed.checked_sub(offset)?;
        if seq.len() < coords.umi.1 || qual.len() < coords.umi.1 {
            return None;
        }

        if !matches!(
            self.version,
            BdCellVersion::V2_96 | BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj
        ) {
            return self.call_exact_shift(seq, qual, offset, shift);
        }

        let linker_mismatches = self.v2_linker_mismatches_at(packed, base)?;
        if linker_mismatches > BD_V2_MAX_LINKER_MISMATCHES {
            return None;
        }

        let (c1_idx, c1_mismatches) =
            Self::index_block_unique_onehot(packed.window::<9>(coords.c1.0).ok()?, &self.c1_fuzzy)?;
        let (c2_idx, c2_mismatches) =
            Self::index_block_unique_onehot(packed.window::<9>(coords.c2.0).ok()?, &self.c2_fuzzy)?;
        let (c3_idx, c3_mismatches) =
            Self::index_block_unique_onehot(packed.window::<9>(coords.c3.0).ok()?, &self.c3_fuzzy)?;

        let mut linker_signature = [0u8; 8];
        linker_signature[..4].copy_from_slice(seq.get(coords.c1.1..coords.c2.0)?);
        linker_signature[4..].copy_from_slice(seq.get(coords.c2.1..coords.c3.0)?);
        let diagnostics = Some(RhapsodyCallDiagnostics {
            linker_signature,
            linker_mismatches: linker_mismatches as u8,
            c1_mismatches: c1_mismatches as u8,
            c2_mismatches: c2_mismatches as u8,
            c3_mismatches: c3_mismatches as u8,
        });

        let cell_id =
            c1_idx * self.block_size * self.block_size + c2_idx * self.block_size + c3_idx + 1;

        let mut cell_seq = Vec::with_capacity(27);
        let mut cell_qual = Vec::with_capacity(27);
        self.extend_part(
            &mut cell_seq,
            &mut cell_qual,
            seq,
            qual,
            coords.c1,
            Some(self.c1[c1_idx as usize]),
        );
        self.extend_part(
            &mut cell_seq,
            &mut cell_qual,
            seq,
            qual,
            coords.c2,
            Some(self.c2[c2_idx as usize]),
        );
        self.extend_part(
            &mut cell_seq,
            &mut cell_qual,
            seq,
            qual,
            coords.c3,
            Some(self.c3[c3_idx as usize]),
        );

        Some(RhapsodyCellCall {
            version: self.version,
            cell_id,
            cell_seq,
            cell_qual,
            umi_seq: seq[coords.umi.0..coords.umi.1].to_vec(),
            umi_qual: qual[coords.umi.0..coords.umi.1].to_vec(),
            shift,
            consumed,
            c1: coords.c1,
            c2: coords.c2,
            c3: coords.c3,
            umi: coords.umi,
            diagnostics,
        })
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

        // consumed must be relative to the grammar start, not absolute
        let consumed = coords.consumed.checked_sub(offset)?;

        if seq.len() < coords.umi.1 || qual.len() < coords.umi.1 {
            return None;
        }

        // The two 4 bp BD v2 linkers are our cheap structural gate. Reads with
        // more than two mismatches across the combined 8 bp linker are not
        // worth attempting to rescue. Once that gate passes, the cassette is
        // considered bona fide and each 9 bp cell block is assigned to its
        // unique nearest whitelist entry. There is deliberately no additional
        // barcode mismatch cutoff: only a tied nearest neighbour is rejected.
        let (c1_idx, c2_idx, c3_idx) = match self.version {
            BdCellVersion::V2_96 | BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj => {
                if !self.v2_linkers_pass(seq, &coords) {
                    return None;
                }
                (
                    self.index_c1_unique(&seq[coords.c1.0..coords.c1.1])?,
                    self.index_c2_unique(&seq[coords.c2.0..coords.c2.1])?,
                    self.index_c3_unique(&seq[coords.c3.0..coords.c3.1])?,
                )
            }
            BdCellVersion::V1 => (
                self.index_c1(&seq[coords.c1.0..coords.c1.1])?,
                self.index_c2(&seq[coords.c2.0..coords.c2.1])?,
                self.index_c3(&seq[coords.c3.0..coords.c3.1])?,
            ),
        };

        let diagnostics = match self.version {
            BdCellVersion::V2_96 | BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj => {
                self.v2_call_diagnostics(seq, &coords, c1_idx, c2_idx, c3_idx)
            }
            BdCellVersion::V1 => None,
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
            coords.c1,
            Some(self.c1[c1_idx as usize]),
        );

        self.extend_part(
            &mut cell_seq,
            &mut cell_qual,
            seq,
            qual,
            coords.c2,
            Some(self.c2[c2_idx as usize]),
        );

        self.extend_part(
            &mut cell_seq,
            &mut cell_qual,
            seq,
            qual,
            coords.c3,
            Some(self.c3[c3_idx as usize]),
        );

        Some(RhapsodyCellCall {
            version: self.version,
            cell_id,
            cell_seq,
            cell_qual,
            umi_seq: seq[coords.umi.0..coords.umi.1].to_vec(),
            umi_qual: qual[coords.umi.0..coords.umi.1].to_vec(),
            shift,
            consumed,
            c1: coords.c1,
            c2: coords.c2,
            c3: coords.c3,
            umi: coords.umi,
            diagnostics,
        })
    }

    pub fn expected_id(&self, c1: u64, c2: u64, c3: u64) -> u64 {
        c1 * self.block_size * self.block_size + c2 * self.block_size + c3 + 1
    }

    #[inline]
    fn index_block_slow(seq: &[u8], fuzzy: &WhitelistHash<9>, max_mismatches: u32) -> Option<u64> {
        // The terminal 4-base keys only select a small whitelist candidate set.
        // The complete 9-base OneHot comparison still decides the unique best
        // hit, so this preserves correction semantics without scanning all 384
        // entries after every non-exact barcode block.
        let (idx, _dist) = fuzzy.best_match(seq, max_mismatches)?;
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

    #[inline]
    fn index_block_unique(seq: &[u8], whitelist: &WhitelistHash<9>) -> Option<u64> {
        let (idx, _dist) = whitelist.unique_nearest(seq)?;
        Some(idx as u64)
    }

    #[inline]
    fn index_c1_unique(&self, seq: &[u8]) -> Option<u64> {
        Self::index_block_unique(seq, &self.c1_fuzzy)
    }

    #[inline]
    fn index_c2_unique(&self, seq: &[u8]) -> Option<u64> {
        Self::index_block_unique(seq, &self.c2_fuzzy)
    }

    #[inline]
    fn index_c3_unique(&self, seq: &[u8]) -> Option<u64> {
        Self::index_block_unique(seq, &self.c3_fuzzy)
    }

    #[inline]
    fn v2_linkers_pass(&self, seq: &[u8], coords: &BdCoords) -> bool {
        let linker1_start = coords.c1.1;
        let linker1_end = coords.c2.0;
        let linker2_start = coords.c2.1;
        let linker2_end = coords.c3.0;
        if linker1_end - linker1_start != 4 || linker2_end - linker2_start != 4 {
            return false;
        }

        let mut observed = [0u8; 8];
        observed[..4].copy_from_slice(&seq[linker1_start..linker1_end]);
        observed[4..].copy_from_slice(&seq[linker2_start..linker2_end]);

        let Ok(observed) = OneHot::<8>::from_bytes(&observed) else {
            return false;
        };
        let expected =
            OneHot::<8>::from_bytes(self.v2_linkers()).expect("builtin BD v2 linker must encode");
        observed.mismatches(expected) <= BD_V2_MAX_LINKER_MISMATCHES
    }

    #[inline]
    fn v2_linker_1(&self) -> &'static [u8; 4] {
        match self.version {
            BdCellVersion::V2_384Vdj => BD_V2_VDJ_LINKER_1,
            _ => BD_V2_LINKER_1,
        }
    }

    #[inline]
    fn v2_linker_2(&self) -> &'static [u8; 4] {
        match self.version {
            BdCellVersion::V2_384Vdj => BD_V2_VDJ_LINKER_2,
            _ => BD_V2_LINKER_2,
        }
    }

    #[inline]
    fn v2_linkers(&self) -> &'static [u8; 8] {
        match self.version {
            BdCellVersion::V2_384Vdj => BD_V2_VDJ_LINKERS,
            _ => BD_V2_LINKERS,
        }
    }

    pub fn create_cell_cassette(&self, c1_idx: usize, c2_idx: usize, c3_idx: usize) -> Vec<u8> {
        let (c1s, c2s, c3s) = match self.version {
            BdCellVersion::V1 => (BD_V2_96_C1, BD_V2_96_C2, BD_V2_96_C3),
            BdCellVersion::V2_96 => (BD_V2_96_C1, BD_V2_96_C2, BD_V2_96_C3),
            BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj => {
                (BD_V2_384_C1, BD_V2_384_C2, BD_V2_384_C3)
            }
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

            BdCellVersion::V2_96 | BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj => {
                seq.extend_from_slice(c1s[c1_idx]);
                seq.extend_from_slice(self.v2_linker_1());
                seq.extend_from_slice(c2s[c2_idx]);
                seq.extend_from_slice(self.v2_linker_2());
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
            BdCellVersion::V2_96 | BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj => {
                Some(BdCoords {
                    c1: (base, base + 9),
                    c2: (base + 13, base + 22),
                    c3: (base + 26, base + 35),
                    umi: (base + 36, base + 42),
                    consumed: base + 42,
                })
            }
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
    use crate::Orientation;
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

    fn bd_v2_384_real_r1_with_prefix(prefix_len: usize, corrupt_linkers: bool) -> Vec<u8> {
        assert!(prefix_len <= 4);

        // Real whitelist blocks from the R1 fixture used below:
        //   C1 = CGGAGAGAT
        //   C2 = GCGCCATAT
        //   C3 = GCGGAGCAT
        // Keep the cassette otherwise realistic, including a 6 bp UMI and
        // poly-T tail, while varying only the number of leading bases.
        let mut seq = vec![b'N'; prefix_len];
        seq.extend_from_slice(b"CGGAGAGAT");
        seq.extend_from_slice(if corrupt_linkers { b"GTAA" } else { b"GTGA" });
        seq.extend_from_slice(b"GCGCCATAT");
        seq.extend_from_slice(if corrupt_linkers { b"GATA" } else { b"GACA" });
        seq.extend_from_slice(b"GCGGAGCAT");
        seq.extend_from_slice(b"TGAACC");
        seq.extend_from_slice(b"TTTTTTTTTTTTTTTTTTTT");
        seq
    }

    #[test]
    fn bd_v2_384_builtin_chemistry_accepts_real_variable_starts_0_through_4() {
        let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();

        for prefix_len in 0..=4 {
            let seq = bd_v2_384_real_r1_with_prefix(prefix_len, false);
            let qual = qual(seq.len());
            let hit = detector
                .detect_first(&seq, &qual)
                .unwrap()
                .unwrap_or_else(|| panic!("BD v2.384 failed at real R1 start offset {prefix_len}"));

            assert_eq!(hit.bd_cell_id, Some(45928512), "offset {prefix_len}");
        }
    }

    #[test]
    fn bd_v2_384_builtin_chemistry_tolerates_crappy_linkers_at_variable_starts() {
        let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();

        for prefix_len in 0..=4 {
            let seq = bd_v2_384_real_r1_with_prefix(prefix_len, true);
            let qual = qual(seq.len());
            let hit = detector
                .detect_first(&seq, &qual)
                .unwrap()
                .unwrap_or_else(|| {
                    panic!(
                        "BD v2.384 rejected a whitelist-valid cell at start offset {prefix_len} \
                         solely because both linkers contain one sequencing error"
                    )
                });

            assert_eq!(hit.bd_cell_id, Some(45928512), "offset {prefix_len}");
        }
    }

    #[test]
    fn bd_v2_384_rejects_when_combined_linkers_have_three_mismatches() {
        let wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);
        let mut seq = bd_v2_384_real_r1_with_prefix(0, false);
        // Linker1 starts at 9, linker2 at 22. Introduce three total errors.
        seq[9] = b'A';
        seq[10] = b'A';
        seq[22] = b'A';
        let qual = qual(seq.len());
        assert!(wl.call(&seq, &qual, 0, 0, 0).is_none());
    }

    #[test]
    fn bd_v2_384_recovers_unique_cell_block_beyond_one_mismatch() {
        let wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);
        let mut seq = bd_v2_384_real_r1_with_prefix(0, false);
        let expected = wl.call(&seq, &qual(seq.len()), 0, 0, 0).unwrap();

        // Damage two positions in C1. The linker gate remains perfect and the
        // unique nearest whitelist entry must still recover the same cell.
        seq[0] = if seq[0] == b'A' { b'C' } else { b'A' };
        seq[1] = if seq[1] == b'A' { b'C' } else { b'A' };
        let rescued = wl
            .call(&seq, &qual(seq.len()), 0, 0, 0)
            .expect("unique nearest C1 should be rescued beyond one mismatch");
        assert_eq!(rescued.cell_id, expected.cell_id);
        assert_eq!(rescued.cell_seq, expected.cell_seq);
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
    fn bd_v2_384_detect_all_is_positional() {
        let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();
        let wl = RhapsodyWhitelist::bd_v2_384();

        let mut seq = Vec::new();
        let mut qual = Vec::new();

        for cell_id in [1u64, 2u64] {
            let (c1, c2, c3) = wl
                .cell_id_to_parts_ids(cell_id)
                .expect("Used a wrong cell id - lib error!");
            let cell_seq_primer = wl.create_cell_cassette(c1, c2, c3);
            let mut primer = detector
                .grammar()
                .synthesize(&cell_seq_primer, b"ACGTAC")
                .expect("Primer creation failed!");
            primer.extend_from_slice(b"GATCGATCGATCGATCGATCGATCGATCG");
            seq.extend_from_slice(&primer);
            qual.extend(std::iter::repeat_n(b'I', primer.len()));
        }

        let hits = detector.detect_all(&seq, &qual).unwrap();
        assert_eq!(hits.len(), 1, "BD detection must not wander downstream");
        assert_eq!(hits[0].bd_cell_id, Some(1));
    }

    #[test]
    fn bd_v2_384_accepts_shift_four_in_both_orientations() {
        let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();
        let wl = RhapsodyWhitelist::bd_v2_384();
        let (c1, c2, c3) = wl.cell_id_to_parts_ids(1).unwrap();
        let cell = wl.create_cell_cassette(c1, c2, c3);
        let primer = detector.grammar().synthesize(&cell, b"ACGTAC").unwrap();

        // Grammar::synthesize() materializes SEARCH as four deterministic
        // placeholder bases, so `primer` already represents the maximum
        // allowed shift (4). Do not prepend another four bases here.
        let mut shifted = primer;
        shifted.extend_from_slice(b"GATCGATCGATC");
        let q = qual(shifted.len());
        let forward = detector.detect_first(&shifted, &q).unwrap().unwrap();
        assert_eq!(forward.orientation, Orientation::Forward);
        assert_eq!(forward.bd_cell_id, Some(1));

        let rc = PrimerDetector::reverse_complement(&shifted);
        let rcq = qual(rc.len());
        let reverse = detector.detect_first(&rc, &rcq).unwrap().unwrap();
        assert_eq!(reverse.orientation, Orientation::ReverseComplement);
        assert_eq!(reverse.bd_cell_id, Some(1));
    }

    #[test]
    fn mouse_igk_transcript_must_not_match_bd_v2_384() {
        let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();
        let seq = b"ANAGGAAACTCATGGTGCGTGGATCTGGCAATGAGCCTGCCGCCACTATCAGTCGTGGCATATGTGAGTCGTGATTATAGAGAGAGAGACCAAAATTCAAAGAGAAAATGGATTTTCAGGTGCAGATTTTCAGCTTCCTGCTAATCAGTGC";
        let q = qual(seq.len());

        let hit = detector.detect_first(seq, &q).unwrap();
        assert!(
            hit.is_none(),
            "BLAST-confirmed mouse Igk transcript must never be accepted as a BD cell cassette: {hit:?}"
        );
    }

    #[test]
    fn canonical_sequence_roundtrips_to_positional_cell_id() {
        let wl = RhapsodyWhitelist::bd_v2_384();
        for id in [1u64, 384, 385, 147_457] {
            let seq = wl.cell_id_to_seq(id).unwrap();
            assert_eq!(wl.cell_id_for_seq(&seq), Some(id));
        }
    }
}
