//! Fixed-length one-hot DNA encoding for small barcode/primer matching.
//!
//! Encoding:
//!
//! - `A = 0001`
//! - `C = 0010`
//! - `G = 0100`
//! - `T = 1000`
//! - anything else, including `N`, becomes `0000`
//!
//! This means unknown bases match nothing and count as mismatches. That is usually
//! what you want for barcode correction.
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

#[inline]
const fn decode_nibble(nibble: u8) -> u8 {
    match nibble {
        0b0001 => b'A',
        0b0010 => b'C',
        0b0100 => b'G',
        0b1000 => b'T',
        _ => b'N',
    }
}

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
}
