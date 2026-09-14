use std::collections::HashMap;
use std::fs;
use std::path::Path;

use onehot_dna::OneHot;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HashPart {
    pub start: usize,
    pub len: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WhitelistLayout {
    parts: Vec<HashPart>,
}

impl WhitelistLayout {
    /// Materialize non-overlapping positional hash parts.
    ///
    /// Four bases fill one u8. Any remainder (0..=3 bases) is deliberately
    /// left as gap sequence between full 4 bp parts; the complete OneHot
    /// comparison remains authoritative for those bases.
    pub fn automatic(len: usize) -> Result<Self, String> {
        if len == 0 || len > OneHot::<1>::MAX_LEN {
            return Err(format!(
                "whitelist barcode length must be in 1..={}, observed {len}",
                OneHot::<1>::MAX_LEN
            ));
        }

        if len <= 4 {
            return Ok(Self {
                parts: vec![HashPart { start: 0, len }],
            });
        }

        let part_count = len / 4;
        let gap_bases = len - part_count * 4;
        let gap_slots = part_count.saturating_sub(1);

        let mut parts = Vec::with_capacity(part_count);
        let mut start = 0usize;
        let base_gap = if gap_slots == 0 {
            0
        } else {
            gap_bases / gap_slots
        };
        let extra = if gap_slots == 0 {
            0
        } else {
            gap_bases % gap_slots
        };

        for part in 0..part_count {
            parts.push(HashPart { start, len: 4 });
            start += 4;
            if part < gap_slots {
                start += base_gap + usize::from(part < extra);
            }
        }

        debug_assert!(parts
            .last()
            .is_some_and(|part| part.start + part.len <= len));

        Ok(Self { parts })
    }

    pub fn parts(&self) -> &[HashPart] {
        &self.parts
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct WhitelistBucket {
    /// One vector per positional hash part. Entries are canonical whitelist IDs.
    positions: Vec<Vec<u32>>,
}

impl WhitelistBucket {
    fn new(parts: usize) -> Self {
        Self {
            positions: (0..parts).map(|_| Vec::new()).collect(),
        }
    }
}

/// Error-correcting fixed-length DNA whitelist.
///
/// Each whitelist sequence is stored exactly once as `OneHot<N>`. Candidate
/// routing uses a single 256-way u8 hash; every bucket contains one ID vector
/// per positional hash part. A sequence found at the wrong position therefore
/// contributes no evidence. The complete OneHot sequence decides the unique
/// best match after routing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WhitelistHash<const N: usize> {
    entries: Vec<OneHot<N>>,
    exact: HashMap<u64, u32>,
    exact_onehot: HashMap<u128, u32>,
    data: [WhitelistBucket; 256],
    layout: WhitelistLayout,
    max_mismatches: u32,
}

impl<const N: usize> WhitelistHash<N> {
    fn with_capacity(capacity: usize, max_mismatches: u32) -> Result<Self, String> {
        let layout = WhitelistLayout::automatic(N)?;
        Ok(Self {
            entries: Vec::with_capacity(capacity),
            exact: HashMap::with_capacity(capacity),
            exact_onehot: HashMap::with_capacity(capacity),
            data: std::array::from_fn(|_| WhitelistBucket::new(layout.parts.len())),
            layout,
            max_mismatches,
        })
    }

    pub fn from_sequences<S: AsRef<[u8]>>(seqs: &[S], max_mismatches: u32) -> Result<Self, String> {
        let mut out = Self::with_capacity(seqs.len(), max_mismatches)?;
        for seq in seqs {
            out.insert(seq.as_ref())?;
        }
        Ok(out)
    }

    pub fn from_text(text: &str, max_mismatches: u32) -> Result<Self, String> {
        // Do not materialize millions of temporary Vec<u8>s for large 10x
        // whitelists. Count once for capacity, then insert directly from lines.
        let capacity = text
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .count();
        let mut out = Self::with_capacity(capacity, max_mismatches)?;
        for line in text
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
        {
            out.insert(line.as_bytes())?;
        }
        Ok(out)
    }

    pub fn from_path(path: impl AsRef<Path>, max_mismatches: u32) -> Result<Self, String> {
        let path = path.as_ref();
        let text = fs::read_to_string(path)
            .map_err(|e| format!("failed to read whitelist '{}': {e}", path.display()))?;
        Self::from_text(&text, max_mismatches)
    }

    fn insert(&mut self, seq: &[u8]) -> Result<(), String> {
        if seq.len() != N {
            return Err(format!(
                "whitelist sequence has length {}, expected {N}: {}",
                seq.len(),
                String::from_utf8_lossy(seq)
            ));
        }
        if !seq.iter().all(|base| Self::base_bits(*base).is_some()) {
            return Err(format!(
                "whitelist sequence contains non-ACGT base: {}",
                String::from_utf8_lossy(seq)
            ));
        }

        let one_hot = OneHot::<N>::from_bytes(seq).map_err(|e| e.to_string())?;
        let packed = Self::pack_clean(seq).expect("validated ACGT whitelist sequence");
        if self.exact.contains_key(&packed) {
            return Err(format!(
                "duplicate whitelist sequence: {}",
                String::from_utf8_lossy(seq)
            ));
        }

        let id = u32::try_from(self.entries.len())
            .map_err(|_| "whitelist has more than u32::MAX entries".to_string())?;
        self.entries.push(one_hot);
        self.exact.insert(packed, id);
        self.exact_onehot.insert(one_hot.bits(), id);

        for (position, part) in self.layout.parts.iter().copied().enumerate() {
            let key =
                Self::encode_part(seq, part).expect("validated ACGT whitelist part must encode");
            self.data[key as usize].positions[position].push(id);
        }

        Ok(())
    }

    #[inline]
    fn base_bits(base: u8) -> Option<u8> {
        match base {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    }

    #[inline]
    fn pack_clean(seq: &[u8]) -> Option<u64> {
        if seq.len() > 32 {
            return None;
        }
        let mut packed = 0u64;
        for &base in seq {
            packed = (packed << 2) | u64::from(Self::base_bits(base)?);
        }
        Some(packed)
    }

    #[inline]
    fn encode_part(seq: &[u8], part: HashPart) -> Option<u8> {
        let end = part.start.checked_add(part.len)?;
        if end > seq.len() || part.len == 0 || part.len > 4 {
            return None;
        }
        let mut packed = 0u8;
        for &base in &seq[part.start..end] {
            packed = (packed << 2) | Self::base_bits(base)?;
        }
        Some(packed)
    }

    #[inline]
    fn encode_onehot_part(query: OneHot<N>, part: HashPart) -> Option<u8> {
        let end = part.start.checked_add(part.len)?;
        if end > N || part.len == 0 || part.len > 4 {
            return None;
        }

        let mut packed = 0u8;
        for position in part.start..end {
            let nibble = ((query.bits() >> (position * 4)) & 0b1111) as u8;
            let base = match nibble {
                0b0001 => 0,
                0b0010 => 1,
                0b0100 => 2,
                0b1000 => 3,
                _ => return None,
            };
            packed = (packed << 2) | base;
        }
        Some(packed)
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn max_mismatches(&self) -> u32 {
        self.max_mismatches
    }

    pub fn layout(&self) -> &WhitelistLayout {
        &self.layout
    }

    pub fn exact_index(&self, seq: &[u8]) -> Option<usize> {
        let packed = Self::pack_clean(seq)?;
        self.exact.get(&packed).copied().map(|id| id as usize)
    }

    pub fn sequence(&self, index: usize) -> Option<Vec<u8>> {
        self.entries
            .get(index)
            .copied()
            .map(|entry| entry.to_dna_string().into_bytes())
    }

    pub fn best_match_default(&self, seq: &[u8]) -> Option<(usize, u32)> {
        self.best_match(seq, self.max_mismatches)
    }

    /// Find the unique nearest whitelist entry without imposing a mismatch
    /// radius. Exact A/C/G/T queries retain the O(1) packed lookup; otherwise
    /// every whitelist entry is compared with the complete OneHot sequence.
    ///
    /// This is intended for small, strongly pre-gated whitelists such as the
    /// BD Rhapsody C1/C2/C3 blocks. A tied nearest neighbour is rejected.
    pub fn unique_nearest(&self, seq: &[u8]) -> Option<(usize, u32)> {
        if seq.len() != N {
            return None;
        }

        let query = OneHot::<N>::from_bytes(seq).ok()?;
        self.unique_nearest_onehot(query)
    }

    /// Find the unique nearest whitelist entry from an already packed query.
    /// This is the zero-reencoding path used by read-sized OneHotSequence
    /// scanners. Exact hits remain O(1); only non-exact queries scan the
    /// complete small whitelist.
    #[inline]
    pub fn unique_nearest_onehot(&self, query: OneHot<N>) -> Option<(usize, u32)> {
        if let Some(index) = self.exact_onehot.get(&query.bits()).copied() {
            return Some((index as usize, 0));
        }
        self.scan_all(query, N as u32)
    }

    /// Find the unique nearest whitelist entry using the positional hash as a
    /// candidate router, without imposing a final mismatch radius.
    ///
    /// At least one positional hash part must match exactly. This is the
    /// CellHash-style rescue path used by BD after its linker gate: a damaged
    /// barcode may still be rescued from candidates routed by an intact 4 bp
    /// part, but a query that damages every routing part is rejected rather
    /// than triggering a full whitelist scan.
    #[inline]
    pub fn unique_nearest_routed_onehot(&self, query: OneHot<N>) -> Option<(usize, u32)> {
        if let Some(index) = self.exact_onehot.get(&query.bits()).copied() {
            return Some((index as usize, 0));
        }

        const MAX_PARTS: usize = 8;
        let mut lists: [Option<&[u32]>; MAX_PARTS] = [None; MAX_PARTS];
        let mut active = 0usize;

        for (position, part) in self.layout.parts.iter().copied().enumerate() {
            let Some(key) = Self::encode_onehot_part(query, part) else {
                continue;
            };
            let list = &self.data[key as usize].positions[position];
            if list.is_empty() {
                continue;
            }
            lists[active] = Some(list);
            active += 1;
        }

        // No intact positional hash part means no cheap candidate evidence.
        // Deliberately do not fall back to an all-entry nearest-neighbour scan.
        if active == 0 {
            return None;
        }

        let mut offsets = [0usize; MAX_PARTS];
        let mut best_index = None;
        let mut best_dist = N as u32 + 1;
        let mut tied = false;

        loop {
            let mut next_id = None;
            for slot in 0..active {
                let list = lists[slot].expect("active list");
                if let Some(&id) = list.get(offsets[slot]) {
                    next_id = Some(next_id.map_or(id, |current: u32| current.min(id)));
                }
            }
            let Some(id) = next_id else {
                break;
            };

            // Advance every routing list containing this id so each candidate
            // is evaluated exactly once even when both 4 bp parts vote for it.
            for slot in 0..active {
                let list = lists[slot].expect("active list");
                if list.get(offsets[slot]).copied() == Some(id) {
                    offsets[slot] += 1;
                }
            }

            let candidate = *self.entries.get(id as usize)?;
            let dist = query.mismatches(candidate);
            if dist < best_dist {
                best_index = Some(id as usize);
                best_dist = dist;
                tied = false;
            } else if dist == best_dist && best_index != Some(id as usize) {
                tied = true;
            }
        }

        if tied {
            None
        } else {
            best_index.map(|index| (index, best_dist))
        }
    }

    /// Find the unique best whitelist entry within `max_mismatches`.
    ///
    /// Exact A/C/G/T queries take the O(1) packed path. Fuzzy queries consult
    /// only the position-correct u8 buckets. IDs that cannot collect enough
    /// positional votes are never subjected to a full OneHot comparison.
    pub fn best_match(&self, seq: &[u8], max_mismatches: u32) -> Option<(usize, u32)> {
        if seq.len() != N || max_mismatches > N as u32 {
            return None;
        }

        if let Some(index) = self.exact_index(seq) {
            return Some((index, 0));
        }

        // For very large whitelists and the common radius-one case, probing
        // the 3*N packed neighbours is cheaper than proving an intersection of
        // several crowded positional buckets at every failed read offset.
        // Small whitelists (BD) stay on the positional CellHash-style path.
        if self.entries.len() > 4096 && max_mismatches == 1 {
            return self.best_match_large_radius_one(seq);
        }

        let query = OneHot::<N>::from_bytes(seq).ok()?;
        let required_votes = self
            .layout
            .parts
            .len()
            .saturating_sub(max_mismatches as usize);
        if required_votes == 0 {
            return self.scan_all(query, max_mismatches);
        }

        // N <= 32 and parts are non-overlapping, so there can be at most 8
        // four-base routing parts. Unknown bases simply disable the affected
        // routing part; the final OneHot distance still counts them.
        const MAX_PARTS: usize = 8;
        let mut lists: [Option<&[u32]>; MAX_PARTS] = [None; MAX_PARTS];
        let mut active = 0usize;

        for (position, part) in self.layout.parts.iter().copied().enumerate() {
            let Some(key) = Self::encode_part(seq, part) else {
                continue;
            };
            lists[active] = Some(&self.data[key as usize].positions[position]);
            active += 1;
        }

        if active < required_votes {
            return None;
        }

        let mut offsets = [0usize; MAX_PARTS];
        let mut best_index = None;
        let mut best_dist = max_mismatches + 1;
        let mut tied = false;

        loop {
            let mut next_id = None;
            for slot in 0..active {
                let list = lists[slot].expect("active list");
                if let Some(&id) = list.get(offsets[slot]) {
                    next_id = Some(next_id.map_or(id, |current: u32| current.min(id)));
                }
            }
            let Some(id) = next_id else {
                break;
            };

            let mut votes = 0usize;
            for slot in 0..active {
                let list = lists[slot].expect("active list");
                if list.get(offsets[slot]).copied() == Some(id) {
                    votes += 1;
                    offsets[slot] += 1;
                }
            }

            if votes < required_votes {
                continue;
            }

            let candidate = *self.entries.get(id as usize)?;
            let dist = query.mismatches(candidate);
            if dist > max_mismatches {
                continue;
            }

            if dist < best_dist {
                best_index = Some(id as usize);
                best_dist = dist;
                tied = false;
            } else if dist == best_dist && best_index != Some(id as usize) {
                tied = true;
            }
        }

        if tied {
            None
        } else {
            best_index.map(|index| (index, best_dist))
        }
    }

    fn best_match_large_radius_one(&self, seq: &[u8]) -> Option<(usize, u32)> {
        let mut packed = 0u64;
        let mut unknown = None;

        for (position, &base) in seq.iter().enumerate() {
            packed <<= 2;
            match Self::base_bits(base) {
                Some(bits) => packed |= u64::from(bits),
                None if unknown.is_none() => unknown = Some(position),
                None => return None,
            }
        }

        let mut hit = None;
        let mut consider = |candidate: u64| -> bool {
            let Some(&id) = self.exact.get(&candidate) else {
                return true;
            };
            match hit {
                None => {
                    hit = Some(id);
                    true
                }
                Some(current) if current == id => true,
                Some(_) => false,
            }
        };

        if let Some(position) = unknown {
            let shift = 2 * (N - 1 - position);
            let mask = 0b11u64 << shift;
            for replacement in 0..4u64 {
                if !consider((packed & !mask) | (replacement << shift)) {
                    return None;
                }
            }
            return hit.map(|id| (id as usize, 1));
        }

        for position in 0..N {
            let shift = 2 * (N - 1 - position);
            let mask = 0b11u64 << shift;
            let observed = (packed & mask) >> shift;
            for replacement in 0..4u64 {
                if replacement == observed {
                    continue;
                }
                if !consider((packed & !mask) | (replacement << shift)) {
                    return None;
                }
            }
        }

        hit.map(|id| (id as usize, 1))
    }

    fn scan_all(&self, query: OneHot<N>, max_mismatches: u32) -> Option<(usize, u32)> {
        let mut best = None;
        let mut best_dist = max_mismatches + 1;
        let mut tied = false;
        for (index, candidate) in self.entries.iter().copied().enumerate() {
            let dist = query.mismatches(candidate);
            if dist > max_mismatches {
                continue;
            }
            if dist < best_dist {
                best = Some(index);
                best_dist = dist;
                tied = false;
            } else if dist == best_dist {
                tied = true;
            }
        }
        if tied {
            None
        } else {
            best.map(|index| (index, best_dist))
        }
    }
}

/// Runtime-length facade used by the CLI. The hot matcher remains the same
/// const-generic `WhitelistHash<N>` used by built-in chemistries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RuntimeWhitelistHash {
    N1(WhitelistHash<1>),
    N2(WhitelistHash<2>),
    N3(WhitelistHash<3>),
    N4(WhitelistHash<4>),
    N5(WhitelistHash<5>),
    N6(WhitelistHash<6>),
    N7(WhitelistHash<7>),
    N8(WhitelistHash<8>),
    N9(WhitelistHash<9>),
    N10(WhitelistHash<10>),
    N11(WhitelistHash<11>),
    N12(WhitelistHash<12>),
    N13(WhitelistHash<13>),
    N14(WhitelistHash<14>),
    N15(WhitelistHash<15>),
    N16(WhitelistHash<16>),
    N17(WhitelistHash<17>),
    N18(WhitelistHash<18>),
    N19(WhitelistHash<19>),
    N20(WhitelistHash<20>),
    N21(WhitelistHash<21>),
    N22(WhitelistHash<22>),
    N23(WhitelistHash<23>),
    N24(WhitelistHash<24>),
    N25(WhitelistHash<25>),
    N26(WhitelistHash<26>),
    N27(WhitelistHash<27>),
    N28(WhitelistHash<28>),
    N29(WhitelistHash<29>),
    N30(WhitelistHash<30>),
    N31(WhitelistHash<31>),
    N32(WhitelistHash<32>),
}

macro_rules! runtime_dispatch {
    ($len:expr, $text:expr, $mm:expr; $(($n:literal, $variant:ident)),+ $(,)?) => {
        match $len {
            $($n => WhitelistHash::<$n>::from_text($text, $mm).map(RuntimeWhitelistHash::$variant),)+
            other => Err(format!("runtime whitelist length must be in 1..=32, observed {other}")),
        }
    };
}

impl RuntimeWhitelistHash {
    pub fn from_text(len: usize, text: &str, max_mismatches: u32) -> Result<Self, String> {
        runtime_dispatch!(len, text, max_mismatches;
            (1,N1),(2,N2),(3,N3),(4,N4),(5,N5),(6,N6),(7,N7),(8,N8),
            (9,N9),(10,N10),(11,N11),(12,N12),(13,N13),(14,N14),(15,N15),(16,N16),
            (17,N17),(18,N18),(19,N19),(20,N20),(21,N21),(22,N22),(23,N23),(24,N24),
            (25,N25),(26,N26),(27,N27),(28,N28),(29,N29),(30,N30),(31,N31),(32,N32)
        )
    }

    pub fn from_path(
        len: usize,
        path: impl AsRef<Path>,
        max_mismatches: u32,
    ) -> Result<Self, String> {
        let path = path.as_ref();
        let text = fs::read_to_string(path)
            .map_err(|e| format!("failed to read whitelist '{}': {e}", path.display()))?;
        Self::from_text(len, &text, max_mismatches)
    }

    pub fn best_match(&self, seq: &[u8]) -> Option<(usize, u32)> {
        macro_rules! call { ($($variant:ident),+ $(,)?) => { match self { $(Self::$variant(hash) => hash.best_match_default(seq),)+ } }; }
        call!(
            N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19,
            N20, N21, N22, N23, N24, N25, N26, N27, N28, N29, N30, N31, N32
        )
    }

    pub fn sequence(&self, index: usize) -> Option<Vec<u8>> {
        macro_rules! call { ($($variant:ident),+ $(,)?) => { match self { $(Self::$variant(hash) => hash.sequence(index),)+ } }; }
        call!(
            N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19,
            N20, N21, N22, N23, N24, N25, N26, N27, N28, N29, N30, N31, N32
        )
    }

    pub fn len(&self) -> usize {
        macro_rules! call { ($($variant:ident),+ $(,)?) => { match self { $(Self::$variant(hash) => hash.len(),)+ } }; }
        call!(
            N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19,
            N20, N21, N22, N23, N24, N25, N26, N27, N28, N29, N30, N31, N32
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nine_base_layout_leaves_one_base_gap() {
        let layout = WhitelistLayout::automatic(9).unwrap();
        assert_eq!(
            layout.parts(),
            &[HashPart { start: 0, len: 4 }, HashPart { start: 5, len: 4 }]
        );
    }

    #[test]
    fn sixteen_base_layout_has_four_positional_chunks() {
        let layout = WhitelistLayout::automatic(16).unwrap();
        assert_eq!(
            layout.parts(),
            &[
                HashPart { start: 0, len: 4 },
                HashPart { start: 4, len: 4 },
                HashPart { start: 8, len: 4 },
                HashPart { start: 12, len: 4 },
            ]
        );
    }

    #[test]
    fn one_mutated_chunk_is_corrected_by_full_onehot() {
        let hash = WhitelistHash::<16>::from_sequences(&[b"ACGTTGCAGGCCAATT"], 1).unwrap();
        assert_eq!(hash.best_match_default(b"ACGTTGCAGGCTAATT"), Some((0, 1)));
    }

    #[test]
    fn mutation_in_gap_is_still_seen_by_onehot() {
        let hash = WhitelistHash::<9>::from_sequences(&[b"ACGTATGCA"], 1).unwrap();
        assert_eq!(hash.best_match_default(b"ACGTCTGCA"), Some((0, 1)));
    }

    #[test]
    fn wrong_position_does_not_vote() {
        let hash =
            WhitelistHash::<16>::from_sequences(&[b"ACGTTGCAGGCCAATT", b"TGCAACGTGGCCAATT"], 1)
                .unwrap();
        assert_eq!(hash.best_match_default(b"ACGTTGCAGGCTAATT"), Some((0, 1)));
    }

    #[test]
    fn equal_distance_tie_is_rejected() {
        let hash = WhitelistHash::<8>::from_sequences(&[b"ACGTTGCA", b"ACGATGCA"], 1).unwrap();
        assert_eq!(hash.best_match_default(b"ACGCTGCA"), None);
    }

    #[test]
    fn unique_nearest_can_rescue_beyond_configured_radius() {
        let hash = WhitelistHash::<8>::from_sequences(&[b"ACGTTGCA", b"TTTTTTTT"], 1).unwrap();
        assert_eq!(hash.best_match_default(b"ACGATGTA"), None);
        assert_eq!(hash.unique_nearest(b"ACGATGTA"), Some((0, 2)));
    }

    #[test]
    fn unique_nearest_rejects_ties_without_a_radius() {
        let hash = WhitelistHash::<8>::from_sequences(&[b"ACGTTGCA", b"ACGATGCA"], 1).unwrap();
        assert_eq!(hash.unique_nearest(b"ACGCTGCA"), None);
    }

    #[test]
    fn unique_nearest_onehot_matches_byte_api() {
        let hash = WhitelistHash::<9>::from_sequences(&[b"ACGTTGCAA", b"TTTTTTTTT"], 1).unwrap();
        let query = OneHot::<9>::from_bytes(b"ACGATGTAA").unwrap();
        assert_eq!(
            hash.unique_nearest_onehot(query),
            hash.unique_nearest(b"ACGATGTAA")
        );
    }

    #[test]
    fn routed_unique_nearest_rescues_from_one_intact_hash_part() {
        let hash = WhitelistHash::<9>::from_sequences(&[b"ACGTATGCA", b"TTTTATTTT"], 1).unwrap();
        let query = OneHot::<9>::from_bytes(b"ACGTCTGTA").unwrap();

        // Two mismatches are allowed here because the intact leading ACGT
        // routing part restricts the candidate set before full-distance
        // comparison.
        assert_eq!(hash.unique_nearest_routed_onehot(query), Some((0, 2)));
    }

    #[test]
    fn routed_unique_nearest_rejects_when_all_hash_parts_are_damaged() {
        let hash = WhitelistHash::<9>::from_sequences(&[b"ACGTATGCA"], 1).unwrap();
        let query = OneHot::<9>::from_bytes(b"ACGAATGTA").unwrap();

        // N=9 routes on bases 0..4 and 5..9. Both four-base routing parts
        // contain a mismatch here, so there is deliberately no expensive
        // all-whitelist fallback.
        assert_eq!(hash.unique_nearest_routed_onehot(query), None);
    }

    #[test]
    fn one_n_is_preserved_as_onehot_mismatch() {
        let hash = WhitelistHash::<16>::from_sequences(&[b"ACGTTGCAGGCCAATT"], 1).unwrap();
        assert_eq!(hash.best_match_default(b"NCGTTGCAGGCCAATT"), Some((0, 1)));
    }
}
