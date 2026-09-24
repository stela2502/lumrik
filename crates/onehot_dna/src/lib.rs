//! Fixed-length one-hot DNA encoding for small barcode/primer matching.
//!
//! Encoding:
//!
//! - `A = 0001`
//! - `C = 0010`
//! - `G = 0100`
//! - `T = 1000`
//! Fixed-size barcode encodings keep their historical strict A/C/G/T behaviour:
//! non-ACGT input becomes `0000` and therefore counts as a mismatch.
//!
//! [`OneHotSequence::from_iupac_bytes`] is the biological-sequence representation.
//! It preserves the complete IUPAC DNA alphabet as four-bit possibility masks, so
//! compatibility is simply `a & b != 0`.
//!
//! `OneHot<N>` stores `N` bases in a `u128`, using four bits per base. Therefore
//! `N <= 32`.

use core::fmt;
use core::str::FromStr;

/// Error returned when constructing a [`OneHot`] sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OneHotError {
    /// The input length did not match the const-generic sequence length.
    WrongLength { expected: usize, observed: usize },

    /// `OneHot<N>` stores four bits per base in a `u128`, so `N` must be <= 32.
    TooLong { max: usize, observed: usize },
}

impl fmt::Display for OneHotError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongLength { expected, observed } => {
                write!(
                    f,
                    "wrong sequence length: expected {expected}, observed {observed}"
                )
            }
            Self::TooLong { max, observed } => {
                write!(
                    f,
                    "sequence too long for u128 one-hot encoding: max {max}, observed {observed}"
                )
            }
        }
    }
}

impl std::error::Error for OneHotError {}

/// Fixed-length one-hot DNA sequence.
///
/// Unknown/non-ACGT bases are encoded as zero and therefore match nothing.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct OneHot<const N: usize> {
    bits: u128,
}

/// Convenience alias for BD Rhapsody C1/C2/C3 9 bp blocks.
pub type OneHot9 = OneHot<9>;

impl<const N: usize> OneHot<N> {
    /// Maximum supported length for this implementation.
    pub const MAX_LEN: usize = 32;

    /// Create from raw one-hot bits.
    ///
    /// This does not validate that every nibble contains only one bit.
    #[inline]
    pub const fn from_bits(bits: u128) -> Self {
        Self { bits }
    }

    /// Return raw packed one-hot bits.
    #[inline]
    pub const fn bits(self) -> u128 {
        self.bits
    }

    /// Encode a sequence of exactly `N` bases.
    ///
    /// A/C/G/T are one-hot encoded. Everything else becomes zero.
    pub fn from_bytes(seq: &[u8]) -> Result<Self, OneHotError> {
        if N > Self::MAX_LEN {
            return Err(OneHotError::TooLong {
                max: Self::MAX_LEN,
                observed: N,
            });
        }

        if seq.len() != N {
            return Err(OneHotError::WrongLength {
                expected: N,
                observed: seq.len(),
            });
        }

        let mut bits = 0u128;

        for (i, base) in seq.iter().copied().enumerate() {
            bits |= encode_base(base) << (i * 4);
        }

        Ok(Self { bits })
    }

    /// Encode a sequence of at least `start + N` bases from a larger read.
    pub fn from_window(seq: &[u8], start: usize) -> Result<Self, OneHotError> {
        let end = start.checked_add(N).ok_or(OneHotError::WrongLength {
            expected: N,
            observed: 0,
        })?;

        if end > seq.len() {
            return Err(OneHotError::WrongLength {
                expected: N,
                observed: seq.len().saturating_sub(start),
            });
        }

        Self::from_bytes(&seq[start..end])
    }

    /// Convert back to ASCII DNA.
    ///
    /// Unknown/zero nibbles become `N`.
    pub fn to_dna_string(self) -> String {
        let mut out = String::with_capacity(N);
        for i in 0..N {
            out.push(decode_nibble(((self.bits >> (i * 4)) & 0b1111) as u8) as char);
        }
        out
    }

    /// Count barcode-style mismatches.
    ///
    /// Two bases match if their one-hot masks overlap. Because `N`/unknown is
    /// encoded as zero, unknown bases match nothing and count as mismatches.
    #[inline]
    pub fn mismatches(self, other: Self) -> u32 {
        (self.bits ^ other.bits).count_ones().div_ceil(2)
    }

    /// Exact bit equality.
    #[inline]
    pub fn exact_match(self, other: Self) -> bool {
        self.bits == other.bits
    }

    /// Returns true if the barcode-style mismatch count is at most `max_mismatches`.
    #[inline]
    pub fn within(self, other: Self, max_mismatches: u32) -> bool {
        self.mismatches(other) <= max_mismatches
    }

    /// Reverse-complement this packed sequence without decoding it.
    ///
    /// The one-hot nibble encoding is deliberately symmetric:
    /// `A=0001`, `C=0010`, `G=0100`, `T=1000`. Reversing all meaningful
    /// bits therefore both reverses the order of the base nibbles and
    /// complements each base in one operation. Unknown/zero nibbles remain
    /// zero.
    #[inline]
    pub fn reverse_complement(self) -> Self {
        if N == 0 {
            return self;
        }

        let meaningful_bits = N * 4;
        Self {
            bits: self.bits.reverse_bits() >> (128 - meaningful_bits),
        }
    }
}

/// Packed one-hot representation of an arbitrarily long biological DNA sequence.
///
/// Two four-bit IUPAC possibility masks are stored per byte. Even positions use
/// the low nibble and odd positions use the high nibble. This keeps arbitrary
/// offsets cheap while retaining the full IUPAC alphabet.
#[derive(Clone, PartialEq, Eq)]
pub struct OneHotSequence {
    packed: Vec<u8>,
    len: usize,
}

/// One biological DNA state represented as an IUPAC possibility mask.
pub type OneHotNibble = u8;

impl OneHotSequence {
    pub const BASES_PER_BYTE: usize = 2;

    /// Pack DNA directly into the two-bases-per-byte representation.
    ///
    /// Every byte is decoded through a 256-entry IUPAC lookup table. Invalid
    /// sequence symbols are rejected instead of being silently converted to N.
    #[inline]
    pub fn from_bytes(seq: &[u8]) -> Self {
        Self::try_from_iupac_bytes(seq).expect("invalid IUPAC DNA sequence")
    }

    /// Fallible constructor for external/input boundaries. Invalid sequence
    /// symbols are returned as `anyhow::Error`; no private error hierarchy is
    /// introduced for arbitrary biological sequences.
    #[inline]
    pub fn try_from_bytes(seq: &[u8]) -> anyhow::Result<Self> {
        Self::try_from_iupac_bytes(seq)
    }

    /// Pack a biological DNA sequence while preserving the complete IUPAC alphabet.
    /// This is a single-pass, single-allocation conversion. Invalid input panics;
    /// use `try_from_iupac_bytes` at untrusted input boundaries.
    #[inline]
    pub fn from_iupac_bytes(seq: &[u8]) -> Self {
        Self::try_from_iupac_bytes(seq).expect("invalid IUPAC DNA sequence")
    }

    pub fn try_from_iupac_bytes(seq: &[u8]) -> anyhow::Result<Self> {
        let mut packed = Vec::with_capacity(seq.len().div_ceil(Self::BASES_PER_BYTE));
        let mut i = 0usize;
        while i + 1 < seq.len() {
            let lo = IUPAC_LUT[seq[i] as usize];
            let hi = IUPAC_LUT[seq[i + 1] as usize];
            if lo == INVALID_IUPAC {
                anyhow::bail!(
                    "invalid IUPAC base {:?} at sequence position {}",
                    seq[i] as char,
                    i
                );
            }
            if hi == INVALID_IUPAC {
                anyhow::bail!(
                    "invalid IUPAC base {:?} at sequence position {}",
                    seq[i + 1] as char,
                    i + 1
                );
            }
            packed.push(lo | (hi << 4));
            i += 2;
        }
        if i < seq.len() {
            let lo = IUPAC_LUT[seq[i] as usize];
            if lo == INVALID_IUPAC {
                anyhow::bail!(
                    "invalid IUPAC base {:?} at sequence position {}",
                    seq[i] as char,
                    i
                );
            }
            packed.push(lo);
        }
        Ok(Self {
            packed,
            len: seq.len(),
        })
    }

    /// Build directly from four-bit biological possibility masks.
    pub fn from_masks(masks: &[u8]) -> Self {
        let mut packed = vec![0u8; masks.len().div_ceil(Self::BASES_PER_BYTE)];
        for (i, &mask) in masks.iter().enumerate() {
            let shift = (i & 1) * 4;
            packed[i >> 1] |= (mask & 0x0f) << shift;
        }
        Self {
            packed,
            len: masks.len(),
        }
    }

    /// Expand an existing four-bases-per-byte 2-bit sequence directly into the
    /// nibble-packed OneHot layout. One input byte becomes exactly two output
    /// bytes through a 256-entry lookup table; no DNA text is materialized.
    pub fn from_2bit_bytes(encoded: &[u8], len: usize) -> Self {
        let mut packed = Vec::with_capacity(len.div_ceil(Self::BASES_PER_BYTE));
        for &byte in encoded.iter().take(len.div_ceil(4)) {
            let expanded = TWO_BIT_TO_ONEHOT[byte as usize];
            packed.push(expanded as u8);
            packed.push((expanded >> 8) as u8);
        }
        packed.truncate(len.div_ceil(Self::BASES_PER_BYTE));
        if len & 1 != 0 {
            if let Some(last) = packed.last_mut() {
                *last &= 0x0f;
            }
        }
        Self { packed, len }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Expose the two-bases-per-byte packed storage.
    #[inline]
    pub fn packed_bytes(&self) -> &[u8] {
        &self.packed
    }

    /// Return the four-bit IUPAC possibility mask at one base position.
    #[inline(always)]
    pub fn mask_at(&self, pos: usize) -> Option<OneHotNibble> {
        if pos >= self.len {
            return None;
        }
        let byte = self.packed[pos >> 1];
        Some(if pos & 1 == 0 { byte & 0x0f } else { byte >> 4 })
    }

    /// Return two consecutive OneHot nibbles packed into one byte. The first
    /// base is in the low nibble. At the final odd base the high nibble is zero.
    #[inline(always)]
    pub fn packed_pair_at(&self, pos: usize) -> Option<u8> {
        if pos >= self.len {
            return None;
        }
        if pos & 1 == 0 {
            return self.packed.get(pos >> 1).copied();
        }
        let lo = self.packed[pos >> 1] >> 4;
        let hi = self.packed.get((pos >> 1) + 1).copied().unwrap_or(0) & 0x0f;
        Some(lo | (hi << 4))
    }

    #[inline(always)]
    pub fn compatible_at(&self, pos: usize, other: &Self, other_pos: usize) -> bool {
        match (self.mask_at(pos), other.mask_at(other_pos)) {
            (Some(a), Some(b)) => compatible_masks(a, b),
            _ => false,
        }
    }

    /// Find the first position at or after `from` whose lookup nibble and every
    /// following evidence nibble are IUPAC-compatible.  At least one external
    /// nibble is mandatory: a one-base hit is deliberately never sufficient.
    ///
    /// The sequence stays nibble-packed; candidates are rejected with a single
    /// mask AND per base and the scan stops at the first incompatible nibble.
    pub fn find_next_with_external_after(
        &self,
        lookup: OneHotNibble,
        external: &[OneHotNibble],
        from: usize,
    ) -> Option<usize> {
        if lookup == 0 || external.is_empty() {
            return None;
        }
        let width = 1usize.checked_add(external.len())?;
        let last = self.len.checked_sub(width)?;
        if from > last {
            return None;
        }
        'candidate: for pos in from..=last {
            if !compatible_masks(self.mask_at(pos)?, lookup) {
                continue;
            }
            for (i, &expected) in external.iter().enumerate() {
                if expected == 0 || !compatible_masks(self.mask_at(pos + 1 + i)?, expected) {
                    continue 'candidate;
                }
            }
            return Some(pos);
        }
        None
    }

    #[inline]
    pub fn find_with_external_after(
        &self,
        lookup: OneHotNibble,
        external: &[OneHotNibble],
    ) -> Option<usize> {
        self.find_next_with_external_after(lookup, external, 0)
    }

    /// Mirror of `find_with_external_after`: `external` occurs immediately
    /// before the lookup nibble.  External evidence is ordered left-to-right.
    pub fn find_next_with_external_before(
        &self,
        lookup: OneHotNibble,
        external: &[OneHotNibble],
        from: usize,
    ) -> Option<usize> {
        if lookup == 0 || external.is_empty() {
            return None;
        }
        let first = from.max(external.len());
        if first >= self.len {
            return None;
        }
        'candidate: for pos in first..self.len {
            if !compatible_masks(self.mask_at(pos)?, lookup) {
                continue;
            }
            let start = pos - external.len();
            for (i, &expected) in external.iter().enumerate() {
                if expected == 0 || !compatible_masks(self.mask_at(start + i)?, expected) {
                    continue 'candidate;
                }
            }
            return Some(pos);
        }
        None
    }

    #[inline]
    pub fn find_with_external_before(
        &self,
        lookup: OneHotNibble,
        external: &[OneHotNibble],
    ) -> Option<usize> {
        self.find_next_with_external_before(lookup, external, external.len())
    }

    /// Find a compatible pattern using a strong packed seed before checking the
    /// complete pattern. `seed_len` is normally eight bases (four packed bytes).
    /// Shorter seeds are supported for callers operating at physical boundaries.
    pub fn find_next_compatible_seed(
        &self,
        pattern: &Self,
        from: usize,
        seed_len: usize,
    ) -> Option<usize> {
        self.find_next_compatible_seed_with_mismatches(pattern, from, seed_len, 0)
    }

    /// Find a pattern candidate using a packed seed while allowing a bounded
    /// number of incompatible seed positions. The complete caller-specific
    /// mismatch policy must still be checked after this fast candidate gate.
    pub fn find_next_compatible_seed_with_mismatches(
        &self,
        pattern: &Self,
        from: usize,
        seed_len: usize,
        max_mismatches: usize,
    ) -> Option<usize> {
        if pattern.is_empty() || seed_len == 0 || seed_len > pattern.len() {
            return None;
        }
        let last = self.len.checked_sub(pattern.len())?;
        if from > last {
            return None;
        }

        'candidate: for pos in from..=last {
            let mut done = 0usize;
            let mut mismatches = 0usize;
            while done < seed_len {
                let observed = self.packed_pair_at(pos + done)?;
                let expected = pattern.packed_pair_at(done)?;
                let take = (seed_len - done).min(2);
                for shift in [0, 4].into_iter().take(take) {
                    if !compatible_masks((observed >> shift) & 0x0f, (expected >> shift) & 0x0f) {
                        mismatches += 1;
                        if mismatches > max_mismatches {
                            continue 'candidate;
                        }
                    }
                }
                done += take;
            }

            return Some(pos);
        }
        None
    }

    /// Count informative and compatible positions across an aligned overlap.
    /// Two bases are consumed per iteration using a single packed byte from
    /// each sequence, including arbitrary odd starts.
    pub fn compatibility_counts(
        &self,
        start: usize,
        other: &Self,
        other_start: usize,
        len: usize,
    ) -> Option<(usize, usize)> {
        if start.checked_add(len)? > self.len || other_start.checked_add(len)? > other.len {
            return None;
        }
        let mut informative = 0usize;
        let mut compatible = 0usize;
        let mut done = 0usize;
        while done < len {
            let a = self.packed_pair_at(start + done)?;
            let b = other.packed_pair_at(other_start + done)?;
            let take = (len - done).min(2);
            for shift in [0, 4].into_iter().take(take) {
                let am = (a >> shift) & 0x0f;
                let bm = (b >> shift) & 0x0f;
                if am != 0 && bm != 0 {
                    informative += 1;
                    if compatible_masks(am, bm) {
                        compatible += 1;
                    }
                }
            }
            done += take;
        }
        Some((informative, compatible))
    }

    pub fn to_dna_vec(&self) -> Vec<u8> {
        (0..self.len)
            .map(|pos| decode_nibble(self.mask_at(pos).unwrap_or(0)))
            .collect()
    }

    /// Extract a fixed-size window from the packed sequence.
    #[inline]
    pub fn window<const N: usize>(&self, start: usize) -> Result<OneHot<N>, OneHotError> {
        if N > OneHot::<N>::MAX_LEN {
            return Err(OneHotError::TooLong {
                max: OneHot::<N>::MAX_LEN,
                observed: N,
            });
        }
        let end = start.checked_add(N).ok_or(OneHotError::WrongLength {
            expected: N,
            observed: 0,
        })?;
        if end > self.len {
            return Err(OneHotError::WrongLength {
                expected: N,
                observed: self.len.saturating_sub(start),
            });
        }
        let mut bits = 0u128;
        for i in 0..N {
            bits |= (self.mask_at(start + i).unwrap_or(0) as u128) << (i * 4);
        }
        Ok(OneHot::from_bits(bits))
    }

    #[inline]
    pub fn reverse_complement_window<const N: usize>(
        &self,
        start: usize,
    ) -> Result<OneHot<N>, OneHotError> {
        let end = start.checked_add(N).ok_or(OneHotError::WrongLength {
            expected: N,
            observed: 0,
        })?;
        if end > self.len {
            return Err(OneHotError::WrongLength {
                expected: N,
                observed: self.len.saturating_sub(start),
            });
        }
        let forward_start = self.len - end;
        Ok(self.window::<N>(forward_start)?.reverse_complement())
    }

    #[inline]
    pub fn mismatches_at<const N: usize>(
        &self,
        start: usize,
        target: OneHot<N>,
    ) -> Result<u32, OneHotError> {
        Ok(self.window::<N>(start)?.mismatches(target))
    }

    #[inline]
    pub fn reverse_complement_mismatches_at<const N: usize>(
        &self,
        start: usize,
        target: OneHot<N>,
    ) -> Result<u32, OneHotError> {
        Ok(self
            .reverse_complement_window::<N>(start)?
            .mismatches(target))
    }

    pub fn to_dna_string(&self) -> String {
        String::from_utf8(self.to_dna_vec()).expect("IUPAC DNA is ASCII")
    }
}

impl fmt::Debug for OneHotSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("OneHotSequence")
            .field("len", &self.len)
            .field("packed_bytes", &self.packed.len())
            .field("seq", &self.to_dna_string())
            .finish()
    }
}

impl<const N: usize> fmt::Debug for OneHot<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("OneHot")
            .field("N", &N)
            .field("seq", &self.to_dna_string())
            .field("bits", &format_args!("0x{:x}", self.bits))
            .finish()
    }
}

impl<const N: usize> fmt::Display for OneHot<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_dna_string())
    }
}

impl<const N: usize> FromStr for OneHot<N> {
    type Err = OneHotError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_bytes(s.as_bytes())
    }
}

/// Encode one DNA base into a four-bit one-hot mask.
///
/// Unknown/non-ACGT bases are encoded as zero.
#[inline]
const fn encode_base(base: u8) -> u128 {
    match base {
        b'A' | b'a' => 0b0001,
        b'C' | b'c' => 0b0010,
        b'G' | b'g' => 0b0100,
        b'T' | b't' => 0b1000,
        _ => 0b0000,
    }
}

const INVALID_IUPAC: u8 = 0xff;

const fn build_iupac_lut() -> [u8; 256] {
    let mut lut = [INVALID_IUPAC; 256];
    lut[b'A' as usize] = 0b0001;
    lut[b'a' as usize] = 0b0001;
    lut[b'C' as usize] = 0b0010;
    lut[b'c' as usize] = 0b0010;
    lut[b'G' as usize] = 0b0100;
    lut[b'g' as usize] = 0b0100;
    lut[b'T' as usize] = 0b1000;
    lut[b't' as usize] = 0b1000;
    lut[b'U' as usize] = 0b1000;
    lut[b'u' as usize] = 0b1000;
    lut[b'R' as usize] = 0b0101;
    lut[b'r' as usize] = 0b0101;
    lut[b'Y' as usize] = 0b1010;
    lut[b'y' as usize] = 0b1010;
    lut[b'S' as usize] = 0b0110;
    lut[b's' as usize] = 0b0110;
    lut[b'W' as usize] = 0b1001;
    lut[b'w' as usize] = 0b1001;
    lut[b'K' as usize] = 0b1100;
    lut[b'k' as usize] = 0b1100;
    lut[b'M' as usize] = 0b0011;
    lut[b'm' as usize] = 0b0011;
    lut[b'B' as usize] = 0b1110;
    lut[b'b' as usize] = 0b1110;
    lut[b'D' as usize] = 0b1101;
    lut[b'd' as usize] = 0b1101;
    lut[b'H' as usize] = 0b1011;
    lut[b'h' as usize] = 0b1011;
    lut[b'V' as usize] = 0b0111;
    lut[b'v' as usize] = 0b0111;
    lut[b'N' as usize] = 0b1111;
    lut[b'n' as usize] = 0b1111;
    lut
}

const IUPAC_LUT: [u8; 256] = build_iupac_lut();

/// Four-bit IUPAC possibility mask. Each bit denotes one canonical base.
#[inline(always)]
pub const fn iupac_mask(base: u8) -> u8 {
    match base {
        b'A' | b'a' => 0b0001,
        b'C' | b'c' => 0b0010,
        b'G' | b'g' => 0b0100,
        b'T' | b't' | b'U' | b'u' => 0b1000,
        b'R' | b'r' => 0b0101,
        b'Y' | b'y' => 0b1010,
        b'S' | b's' => 0b0110,
        b'W' | b'w' => 0b1001,
        b'K' | b'k' => 0b1100,
        b'M' | b'm' => 0b0011,
        b'B' | b'b' => 0b1110,
        b'D' | b'd' => 0b1101,
        b'H' | b'h' => 0b1011,
        b'V' | b'v' => 0b0111,
        b'N' | b'n' => 0b1111,
        _ => 0,
    }
}

#[inline(always)]
pub const fn compatible_masks(a: u8, b: u8) -> bool {
    (a & b) != 0
}

#[inline]
const fn decode_nibble(nibble: u8) -> u8 {
    match nibble & 0x0f {
        0b0001 => b'A',
        0b0010 => b'C',
        0b0100 => b'G',
        0b1000 => b'T',
        0b0101 => b'R',
        0b1010 => b'Y',
        0b0110 => b'S',
        0b1001 => b'W',
        0b1100 => b'K',
        0b0011 => b'M',
        0b1110 => b'B',
        0b1101 => b'D',
        0b1011 => b'H',
        0b0111 => b'V',
        0b1111 => b'N',
        _ => b'N',
    }
}

const fn expand_2bit_byte(packed: u8) -> u16 {
    let mut out = 0u16;
    let mut i = 0usize;
    while i < 4 {
        let base = (packed >> (i * 2)) & 0b11;
        out |= (1u16 << base) << (i * 4);
        i += 1;
    }
    out
}

const fn build_2bit_to_onehot() -> [u16; 256] {
    let mut table = [0u16; 256];
    let mut i = 0usize;
    while i < 256 {
        table[i] = expand_2bit_byte(i as u8);
        i += 1;
    }
    table
}

const TWO_BIT_TO_ONEHOT: [u16; 256] = build_2bit_to_onehot();

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct OneHotSet<const N: usize> {
    data: Vec<OneHot<N>>,
}

impl<const N: usize> OneHotSet<N> {
    pub fn from_sequences<S: AsRef<[u8]>>(seqs: &[S]) -> Result<Self, OneHotError> {
        let data = seqs
            .iter()
            .map(|s| OneHot::<N>::from_bytes(s.as_ref()))
            .collect::<Result<Vec<_>, _>>()?;

        Ok(Self { data })
    }

    pub fn min_pairwise_mismatches(&self) -> Option<usize> {
        if self.data.len() < 2 {
            return None;
        }

        let mut best = usize::MAX;

        for i in 0..self.data.len() {
            for j in (i + 1)..self.data.len() {
                best = best.min(self.data[i].mismatches(self.data[j]) as usize);
            }
        }

        Some(best)
    }

    pub fn correction_radius(&self) -> Option<usize> {
        self.min_pairwise_mismatches()
            .map(|d| d.saturating_sub(1) / 2)
    }

    pub fn as_slice(&self) -> &[OneHot<N>] {
        &self.data
    }

    pub fn into_vec(self) -> Vec<OneHot<N>> {
        self.data
    }

    /// Find the unique best match against this candidate table.
    ///
    /// Returns `(index, distance)` if the best match is unique and within
    /// `max_mismatches`. Returns `None` for no hit or ties.
    pub fn best_match(&self, query: &OneHot<N>, max_mismatches: u32) -> Option<(usize, u32)> {
        let mut best_index = None;
        let mut best_dist = max_mismatches + 1;
        let mut ties = 0u32;

        for (i, candidate) in self.data.iter().copied().enumerate() {
            let d = query.mismatches(candidate);

            if d < best_dist {
                best_index = Some(i);
                best_dist = d;
                ties = 1;
            } else if d == best_dist {
                ties += 1;
            }
        }

        match (best_index, best_dist <= max_mismatches, ties == 1) {
            (Some(i), true, true) => Some((i, best_dist)),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodes_and_decodes() {
        let x = OneHot::<9>::from_bytes(b"ACGTACGTN").unwrap();
        assert_eq!(x.to_dna_string(), "ACGTACGTN");
    }

    #[test]
    fn counts_mismatches() {
        let a = OneHot::<9>::from_bytes(b"AAAAAAAAA").unwrap();
        let b = OneHot::<9>::from_bytes(b"AAAAAAAAC").unwrap();
        assert_eq!(a.mismatches(b), 1);
    }

    #[test]
    fn n_counts_as_mismatch() {
        let a = OneHot::<9>::from_bytes(b"AAAAAAAAA").unwrap();
        let n = OneHot::<9>::from_bytes(b"AAAAAAAAN").unwrap();
        assert_eq!(a.mismatches(n), 1);
        assert!(!a.within(n, 0));
        assert!(a.within(n, 1));
    }

    #[test]
    fn best_match_requires_unique_hit() {
        let candidates = OneHotSet::<9>::from_sequences(&[
            b"AAAAAAAAA".as_slice(),
            b"CCCCCCCCC".as_slice(),
            b"GGGGGGGGG".as_slice(),
        ])
        .unwrap();

        let obs = OneHot9::from_bytes(b"AAAAAAAAC").unwrap();
        assert_eq!(candidates.best_match(&obs, 1), Some((0, 1)));
    }

    #[test]
    fn best_match_rejects_ties() {
        let candidates =
            OneHotSet::<9>::from_sequences(&[b"AAAAAAAAA".as_slice(), b"AAAAAAAAC".as_slice()])
                .unwrap();

        let obs = OneHot9::from_bytes(b"AAAAAAAAN").unwrap();
        assert_eq!(candidates.best_match(&obs, 1), None);
    }

    #[test]
    fn window_extracts_from_read() {
        let read = b"XXACGTACGTNYY";
        let x = OneHot::<9>::from_window(read, 2).unwrap();
        assert_eq!(x.to_dna_string(), "ACGTACGTN");
    }

    #[test]
    fn reverse_complement_is_pure_bit_operation() {
        let x = OneHot::<8>::from_bytes(b"ACGTNAGT").unwrap();
        assert_eq!(x.reverse_complement().to_dna_string(), "ACTNACGT");
        assert_eq!(x.reverse_complement().reverse_complement(), x);
    }

    #[test]
    fn packed_sequence_extracts_windows_across_word_boundaries() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTTGCATGCA";
        let packed = OneHotSequence::from_bytes(seq);

        assert_eq!(packed.len(), seq.len());
        assert_eq!(packed.to_dna_string().as_bytes(), seq);
        assert_eq!(
            packed.window::<9>(28).unwrap().to_dna_string().as_bytes(),
            &seq[28..37]
        );
    }

    #[test]
    fn packed_sequence_reverse_complement_windows_use_mirrored_coordinates() {
        let seq = b"AAAACCCCGGGGTTTTACGTN";
        let packed = OneHotSequence::from_bytes(seq);
        let rc = b"NACGTAAAACCCCGGGGTTTT";

        for start in 0..=(rc.len() - 8) {
            assert_eq!(
                packed
                    .reverse_complement_window::<8>(start)
                    .unwrap()
                    .to_dna_string()
                    .as_bytes(),
                &rc[start..start + 8]
            );
        }
    }

    #[test]
    fn packed_mismatch_scans_reuse_one_encoding() {
        let packed = OneHotSequence::from_bytes(b"TTTTGTGAGACAAAAA");
        let linker = OneHot::<8>::from_bytes(b"GTGAGACA").unwrap();

        assert_eq!(packed.mismatches_at(4, linker).unwrap(), 0);
        assert!(packed.mismatches_at(3, linker).unwrap() > 0);
    }
    #[test]
    fn iupac_sequence_roundtrips_complete_alphabet() {
        let seq = b"ACGTRYSWKMBDHVN";
        let packed = OneHotSequence::from_iupac_bytes(seq);
        assert_eq!(packed.to_dna_vec(), seq);
    }

    #[test]
    fn iupac_compatibility_is_mask_intersection() {
        let ambiguity = OneHotSequence::from_iupac_bytes(b"RYN");
        let agt = OneHotSequence::from_iupac_bytes(b"AGT");
        assert!(ambiguity.compatible_at(0, &agt, 0)); // R with A
        assert!(!ambiguity.compatible_at(1, &agt, 1)); // Y with G
        assert!(ambiguity.compatible_at(2, &agt, 2)); // N with T
    }

    #[test]
    fn two_bit_expansion_matches_strict_onehot() {
        // IntToDna byte layout: low two bits are the first base.
        let encoded = [0b11_10_01_00u8]; // A C G T
        let packed = OneHotSequence::from_2bit_bytes(&encoded, 4);
        assert_eq!(packed.to_dna_string(), "ACGT");
    }

    #[test]
    fn evidence_masks_preserve_ambiguity_and_zero_means_unobserved() {
        let packed = OneHotSequence::from_masks(&[0b0001, 0b0101, 0, 0b1111]);
        assert_eq!(packed.mask_at(0), Some(0b0001));
        assert_eq!(packed.mask_at(1), Some(0b0101));
        assert_eq!(packed.mask_at(2), Some(0));
        assert_eq!(packed.mask_at(3), Some(0b1111));
    }

    #[test]
    fn packed_compatibility_counts_ambiguity_in_32_base_chunks() {
        let a = OneHotSequence::from_masks(&[0b0001, 0b0101, 0, 0b1010, 0b1111]);
        let b = OneHotSequence::from_iupac_bytes(b"AGGCTA");
        // A/A yes, R/G yes, zero/G uninformative, Y/C yes, N/T yes.
        assert_eq!(a.compatibility_counts(0, &b, 0, 5), Some((4, 4)));
    }

    #[test]
    fn packed_pair_is_two_nibbles_in_one_byte() {
        let seq = OneHotSequence::from_iupac_bytes(b"ACR");
        assert_eq!(seq.packed_pair_at(0), Some(0x21)); // A, C
        assert_eq!(seq.packed_pair_at(1), Some(0x52)); // C, R
        assert_eq!(seq.packed_pair_at(2), Some(0x05)); // R, end
    }

    #[test]
    fn packed_compatibility_handles_different_odd_even_offsets() {
        let a = OneHotSequence::from_iupac_bytes(b"TACGTRYSWKMBDHVNACGTACGTACGTACGTACGT");
        let b = OneHotSequence::from_iupac_bytes(b"GGACGTRYSWKMBDHVNACGTACGTACGTACGTACGT");
        assert_eq!(a.compatibility_counts(1, &b, 2, 35), Some((35, 35)));
    }

    #[test]
    fn packed_compatibility_handles_unaligned_word_boundaries() {
        let a = OneHotSequence::from_iupac_bytes(
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAARYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY",
        );
        let b = OneHotSequence::from_iupac_bytes(
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        );
        assert_eq!(a.compatibility_counts(31, &b, 31, 33), Some((33, 33)));
    }

    #[test]
    fn lut_packing_rejects_non_iupac_input() {
        let err = OneHotSequence::try_from_iupac_bytes(b"ACGT|N").unwrap_err();
        assert!(err.to_string().contains("position 4"));
    }

    #[test]
    fn eight_base_seed_finds_even_and_odd_packed_offsets() {
        let pattern = OneHotSequence::from_iupac_bytes(b"ACGTRYSWKM");
        let even = OneHotSequence::from_iupac_bytes(b"TTACGTRYSWKMCC");
        let odd = OneHotSequence::from_iupac_bytes(b"TTTACGTRYSWKMCC");
        assert_eq!(even.find_next_compatible_seed(&pattern, 0, 8), Some(2));
        assert_eq!(odd.find_next_compatible_seed(&pattern, 0, 8), Some(3));
    }

    #[test]
    fn eight_base_seed_strict_rejects_one_mismatch() {
        let observed = OneHotSequence::from_iupac_bytes(b"TTACGTACGTTT");
        let pattern = OneHotSequence::from_iupac_bytes(b"ACGTTCGT");
        assert_eq!(observed.find_next_compatible_seed(&pattern, 0, 8), None);
    }

    #[test]
    fn eight_base_seed_permissive_accepts_one_but_not_two_mismatches() {
        let observed = OneHotSequence::from_iupac_bytes(b"TTACGTACGTTT");
        let one_error = OneHotSequence::from_iupac_bytes(b"ACGTTCGT");
        let two_errors = OneHotSequence::from_iupac_bytes(b"ACGTTGGT");
        assert_eq!(
            observed.find_next_compatible_seed_with_mismatches(&one_error, 0, 8, 1),
            Some(2)
        );
        assert_eq!(
            observed.find_next_compatible_seed_with_mismatches(&two_errors, 0, 8, 1),
            None
        );
    }

    #[test]
    fn onehot_external_search_requires_and_matches_all_evidence() {
        let seq = OneHotSequence::from_iupac_bytes(b"TTACGTRNAA");
        let pattern = OneHotSequence::from_iupac_bytes(b"ACGTR");
        let lookup = pattern.mask_at(0).unwrap();
        let external: Vec<_> = (1..pattern.len())
            .map(|i| pattern.mask_at(i).unwrap())
            .collect();
        assert_eq!(seq.find_with_external_after(lookup, &external), Some(2));
        assert_eq!(seq.find_with_external_after(lookup, &[]), None);
        let wrong = OneHotSequence::from_iupac_bytes(b"ACGTC");
        let wrong_external: Vec<_> = (1..wrong.len())
            .map(|i| wrong.mask_at(i).unwrap())
            .collect();
        assert_eq!(
            seq.find_with_external_after(wrong.mask_at(0).unwrap(), &wrong_external),
            None
        );
    }

    #[test]
    fn onehot_external_search_honours_iupac_and_before_direction() {
        let seq = OneHotSequence::from_iupac_bytes(b"CCAGTCC");
        let expected = OneHotSequence::from_iupac_bytes(b"RGT");
        let after = [expected.mask_at(1).unwrap(), expected.mask_at(2).unwrap()];
        assert_eq!(
            seq.find_with_external_after(expected.mask_at(0).unwrap(), &after),
            Some(2)
        );
        let before = [expected.mask_at(0).unwrap(), expected.mask_at(1).unwrap()];
        assert_eq!(
            seq.find_with_external_before(expected.mask_at(2).unwrap(), &before),
            Some(4)
        );
    }
}
