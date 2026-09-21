use std::fmt;

/// Five-bit amino-acid identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AminoAcid(u8);

/// Intrinsic amino-acid chemistry flags. These describe residue chemistry, not
/// biological context or model-derived properties.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Chemistry(u16);

impl Chemistry {
    pub const HYDROPHOBIC: Self = Self(1 << 0);
    pub const POLAR: Self = Self(1 << 1);
    pub const POSITIVE: Self = Self(1 << 2);
    pub const NEGATIVE: Self = Self(1 << 3);
    pub const AROMATIC: Self = Self(1 << 4);
    pub const ALIPHATIC: Self = Self(1 << 5);
    pub const SMALL: Self = Self(1 << 6);
    pub const TINY: Self = Self(1 << 7);
    pub const SULFUR: Self = Self(1 << 8);
    pub const HYDROXYL: Self = Self(1 << 9);
    pub const AMIDE: Self = Self(1 << 10);
    pub const PROLINE: Self = Self(1 << 11);
    pub const GLYCINE: Self = Self(1 << 12);

    #[inline]
    pub const fn bits(self) -> u16 { self.0 }
}

impl AminoAcid {
    #[inline]
    pub const fn code(self) -> u8 { self.0 }

    #[inline]
    pub fn residue(self) -> Option<u8> { IntToProt::decode_binary(self) }

    #[inline]
    pub const fn chemistry(self) -> u16 { chemistry_bits(self.0) }

    #[inline]
    pub const fn is(self, property: Chemistry) -> bool {
        self.chemistry() & property.bits() != 0
    }

    /// Intrinsic chemistry distance in [0, 1]. This is deliberately independent
    /// of biological context. Equal residues have distance 0; otherwise the
    /// distance is the normalized Hamming distance between chemistry flags.
    /// Unknown/non-residue symbols are maximally distant from ordinary residues.
    pub fn dist(self, other: Self) -> f32 {
        if self == other { return 0.0; }
        let a = self.chemistry();
        let b = other.chemistry();
        if a == 0 || b == 0 { return 1.0; }
        let union = (a | b).count_ones();
        let different = (a ^ b).count_ones();
        different as f32 / union as f32
    }
}

const fn chemistry_bits(code: u8) -> u16 {
    let h = Chemistry::HYDROPHOBIC.0; let p = Chemistry::POLAR.0;
    let pos = Chemistry::POSITIVE.0; let neg = Chemistry::NEGATIVE.0;
    let ar = Chemistry::AROMATIC.0; let al = Chemistry::ALIPHATIC.0;
    let sm = Chemistry::SMALL.0; let ti = Chemistry::TINY.0;
    let su = Chemistry::SULFUR.0; let oh = Chemistry::HYDROXYL.0;
    let am = Chemistry::AMIDE.0; let pro = Chemistry::PROLINE.0; let gly = Chemistry::GLYCINE.0;
    match code {
        0 => sm | ti,                         // A
        1 => p | sm | su,                     // C
        2 => p | neg | sm,                    // D
        3 => p | neg,                         // E
        4 => h | ar,                          // F
        5 => sm | ti | gly,                   // G
        6 => p | pos | ar,                    // H
        7 => h | al,                          // I
        8 => p | pos,                         // K
        9 => h | al,                          // L
        10 => h | su,                         // M
        11 => p | sm | am,                    // N
        12 => h | sm | pro,                   // P
        13 => p | am,                         // Q
        14 => p | pos,                        // R
        15 => p | sm | ti | oh,               // S
        16 => p | sm | oh,                    // T
        17 => h | al | sm,                    // V
        18 => h | ar,                         // W
        19 => p | ar | oh,                    // Y
        22 => p | neg | sm | am,              // B: D/N ambiguity
        23 => p | neg | am,                   // Z: E/Q ambiguity
        24 => h | al,                         // J: I/L ambiguity
        25 => p | sm | su,                    // U: selenocysteine, C-like
        26 => p | pos,                        // O: pyrrolysine, K-like
        _ => 0,                               // X, stop, gap, reserved
    }
}

pub const A: AminoAcid = AminoAcid(0);
pub const C: AminoAcid = AminoAcid(1);
pub const D: AminoAcid = AminoAcid(2);
pub const E: AminoAcid = AminoAcid(3);
pub const F: AminoAcid = AminoAcid(4);
pub const G: AminoAcid = AminoAcid(5);
pub const H: AminoAcid = AminoAcid(6);
pub const I: AminoAcid = AminoAcid(7);
pub const K: AminoAcid = AminoAcid(8);
pub const L: AminoAcid = AminoAcid(9);
pub const M: AminoAcid = AminoAcid(10);
pub const N: AminoAcid = AminoAcid(11);
pub const P: AminoAcid = AminoAcid(12);
pub const Q: AminoAcid = AminoAcid(13);
pub const R: AminoAcid = AminoAcid(14);
pub const S: AminoAcid = AminoAcid(15);
pub const T: AminoAcid = AminoAcid(16);
pub const V: AminoAcid = AminoAcid(17);
pub const W: AminoAcid = AminoAcid(18);
pub const Y: AminoAcid = AminoAcid(19);
pub const X: AminoAcid = AminoAcid(20);
pub const STOP: AminoAcid = AminoAcid(21);
pub const B: AminoAcid = AminoAcid(22);
pub const Z: AminoAcid = AminoAcid(23);
pub const J: AminoAcid = AminoAcid(24);
pub const U: AminoAcid = AminoAcid(25);
pub const O: AminoAcid = AminoAcid(26);
pub const GAP: AminoAcid = AminoAcid(27);

const BITS_PER_RESIDUE: usize = 5;
const CODE_MASK: u64 = 0b1_1111;

/// Compact exact protein sequence representation.
///
/// Residues are stored as consecutive 5-bit codes. The representation is
/// deliberately concerned only with residue identity; biochemical feature
/// models belong in separate layers.
#[derive(Default, Debug, PartialEq, Eq, Clone, Hash)]
pub struct IntToProt {
    packed: Vec<u64>,
    size: usize,
}

impl IntToProt {
    pub fn new<T: AsRef<[u8]>>(input: T) -> Self {
        Self::try_new(input).unwrap_or_else(|e| panic!("failed to encode protein sequence: {e}"))
    }

    pub fn try_new<T: AsRef<[u8]>>(input: T) -> Result<Self, String> {
        let seq = input.as_ref();
        let mut ret = Self {
            packed: vec![0; (seq.len() * BITS_PER_RESIDUE).div_ceil(64)],
            size: seq.len(),
        };

        for (pos, &residue) in seq.iter().enumerate() {
            ret.set_code(pos, Self::encode_binary(residue)?);
        }

        Ok(ret)
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.size
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    #[inline]
    pub fn get(&self, pos: usize) -> Option<AminoAcid> {
        if pos >= self.size {
            return None;
        }

        let bit = pos * BITS_PER_RESIDUE;
        let word = bit / 64;
        let shift = bit % 64;
        let mut value = self.packed[word] >> shift;

        if shift > 64 - BITS_PER_RESIDUE {
            value |= self.packed.get(word + 1).copied().unwrap_or(0) << (64 - shift);
        }

        Some(AminoAcid((value & CODE_MASK) as u8))
    }

    pub fn encode_binary(residue: u8) -> Result<AminoAcid, String> {
        match residue.to_ascii_uppercase() {
            b'A' => Ok(A),
            b'C' => Ok(C),
            b'D' => Ok(D),
            b'E' => Ok(E),
            b'F' => Ok(F),
            b'G' => Ok(G),
            b'H' => Ok(H),
            b'I' => Ok(I),
            b'K' => Ok(K),
            b'L' => Ok(L),
            b'M' => Ok(M),
            b'N' => Ok(N),
            b'P' => Ok(P),
            b'Q' => Ok(Q),
            b'R' => Ok(R),
            b'S' => Ok(S),
            b'T' => Ok(T),
            b'V' => Ok(V),
            b'W' => Ok(W),
            b'Y' => Ok(Y),
            b'X' => Ok(X),
            b'*' => Ok(STOP),
            b'B' => Ok(B),
            b'Z' => Ok(Z),
            b'J' => Ok(J),
            b'U' => Ok(U),
            b'O' => Ok(O),
            b'-' => Ok(GAP),
            other => Err(format!("cannot encode '{}' as an amino acid", other as char)),
        }
    }

    #[inline]
    pub fn decode_binary(code: AminoAcid) -> Option<u8> {
        const ALPHABET: &[u8; 28] = b"ACDEFGHIKLMNPQRSTVWYX*BZJUO-";
        ALPHABET.get(code.0 as usize).copied()
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        (0..self.size)
            .map(|pos| Self::decode_binary(self.get(pos).unwrap()).unwrap())
            .collect()
    }

    fn set_code(&mut self, pos: usize, code: AminoAcid) {
        let bit = pos * BITS_PER_RESIDUE;
        let word = bit / 64;
        let shift = bit % 64;
        self.packed[word] |= (code.0 as u64) << shift;

        if shift > 64 - BITS_PER_RESIDUE {
            self.packed[word + 1] |= (code.0 as u64) >> (64 - shift);
        }
    }
}

impl fmt::Display for IntToProt {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for residue in self.to_bytes() {
            write!(f, "{}", residue as char)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_canonical_protein() {
        let seq = b"ACDEFGHIKLMNPQRSTVWY";
        let protein = IntToProt::new(seq);
        assert_eq!(protein.len(), seq.len());
        assert_eq!(protein.to_bytes(), seq);
        assert_eq!(protein.to_string(), "ACDEFGHIKLMNPQRSTVWY");
    }

    #[test]
    fn round_trip_extended_alphabet_across_word_boundaries() {
        let seq = b"MPEPTIDEACDEFGHIKLMNPQRSTVWYX*BZJUO-";
        let protein = IntToProt::new(seq);
        assert_eq!(protein.to_bytes(), seq);
    }

    #[test]
    fn accepts_lowercase() {
        assert_eq!(IntToProt::new(b"mkwv").to_string(), "MKWV");
    }

    #[test]
    fn rejects_non_protein_symbols() {
        assert!(IntToProt::try_new(b"M?K").is_err());
    }
}

/// A protease supported by the built-in virtual digestor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Protease {
    Trypsin,
    LysC,
    ArgC,
    GluC,
    AspN,
    Chymotrypsin,
}

/// A peptide produced by virtual proteolysis. Coordinates are 0-based,
/// half-open positions in the parent protein.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Peptide {
    pub sequence: IntToProt,
    pub start: usize,
    pub end: usize,
    pub missed_cleavages: usize,
}

impl IntToProt {
    /// Return an exact packed subsequence using 0-based, half-open coordinates.
    pub fn slice(&self, start: usize, end: usize) -> Self {
        assert!(start <= end, "protein slice start must not exceed end");
        assert!(end <= self.size, "protein slice end exceeds sequence length");
        let bytes: Vec<u8> = (start..end)
            .map(|pos| Self::decode_binary(self.get(pos).unwrap()).unwrap())
            .collect();
        Self::new(bytes)
    }

    /// Perform an in-silico digest and retain parent-protein coordinates.
    /// `max_missed_cleavages = 0` returns only fully cleaved peptides.
    pub fn digest(&self, protease: Protease, max_missed_cleavages: usize) -> Vec<Peptide> {
        if self.is_empty() {
            return Vec::new();
        }

        let seq = self.to_bytes();
        let mut cuts = vec![0usize];

        match protease {
            Protease::AspN => {
                for (i, &aa) in seq.iter().enumerate() {
                    if aa == b'D' && i != 0 {
                        cuts.push(i);
                    }
                }
            }
            _ => {
                for i in 0..seq.len() {
                    let aa = seq[i];
                    let next = seq.get(i + 1).copied();
                    let cut = match protease {
                        Protease::Trypsin => (aa == b'K' || aa == b'R') && next != Some(b'P'),
                        Protease::LysC => aa == b'K',
                        Protease::ArgC => aa == b'R',
                        Protease::GluC => aa == b'E',
                        Protease::Chymotrypsin => {
                            matches!(aa, b'F' | b'W' | b'Y') && next != Some(b'P')
                        }
                        Protease::AspN => unreachable!(),
                    };
                    if cut && i + 1 < seq.len() {
                        cuts.push(i + 1);
                    }
                }
            }
        }

        if cuts.last().copied() != Some(seq.len()) {
            cuts.push(seq.len());
        }
        cuts.sort_unstable();
        cuts.dedup();

        let mut peptides = Vec::new();
        for start_cut in 0..cuts.len() - 1 {
            let max_end_cut = (start_cut + max_missed_cleavages + 1).min(cuts.len() - 1);
            for end_cut in (start_cut + 1)..=max_end_cut {
                let start = cuts[start_cut];
                let end = cuts[end_cut];
                peptides.push(Peptide {
                    sequence: self.slice(start, end),
                    start,
                    end,
                    missed_cleavages: end_cut - start_cut - 1,
                });
            }
        }
        peptides
    }
}

#[cfg(test)]
mod chemistry_tests {
    use super::*;

    #[test]
    fn chemistry_flags_are_intrinsic_and_overlapping() {
        assert!(L.is(Chemistry::HYDROPHOBIC));
        assert!(L.is(Chemistry::ALIPHATIC));
        assert!(!L.is(Chemistry::NEGATIVE));
        assert!(D.is(Chemistry::POLAR));
        assert!(D.is(Chemistry::NEGATIVE));
        assert!(F.is(Chemistry::AROMATIC));
    }

    #[test]
    fn chemistry_distance_respects_similarity() {
        assert_eq!(L.dist(L), 0.0);
        assert_eq!(L.dist(I), 0.0);
        assert!(L.dist(V) < L.dist(D));
        assert!(K.dist(R) < K.dist(E));
        assert_eq!(X.dist(L), 1.0);
    }

    #[test]
    fn protein_exposes_amino_acids_with_chemistry() {
        let protein = IntToProt::new(b"LID");
        let leucine = protein.get(0).unwrap();
        let isoleucine = protein.get(1).unwrap();
        let aspartate = protein.get(2).unwrap();
        assert_eq!(leucine.residue(), Some(b'L'));
        assert!(leucine.is(Chemistry::HYDROPHOBIC));
        assert!(leucine.dist(isoleucine) < leucine.dist(aspartate));
    }
}

#[cfg(test)]
mod sequence_operation_tests {
    use super::*;

    #[test]
    fn protein_slice_preserves_exact_sequence() {
        let protein = IntToProt::new(b"MPEPTIDEKTAIL");
        assert_eq!(protein.slice(1, 8).to_string(), "PEPTIDE");
    }

    #[test]
    fn trypsin_digest_retains_parent_coordinates() {
        let protein = IntToProt::new(b"MPEPTIDEKTAIL");
        let peptides = protein.digest(Protease::Trypsin, 0);
        let observed: Vec<(usize, usize, String)> = peptides
            .iter()
            .map(|p| (p.start, p.end, p.sequence.to_string()))
            .collect();
        assert_eq!(
            observed,
            vec![
                (0, 9, "MPEPTIDEK".to_string()),
                (9, 13, "TAIL".to_string()),
            ]
        );
    }

    #[test]
    fn trypsin_respects_proline_exception_and_missed_cleavages() {
        let protein = IntToProt::new(b"AKPQRTAIL");
        let peptides = protein.digest(Protease::Trypsin, 1);
        assert!(peptides.iter().any(|p| p.sequence.to_string() == "AKPQR"));
        assert!(peptides.iter().any(|p| p.sequence.to_string() == "AKPQRTAIL"));
    }
}
