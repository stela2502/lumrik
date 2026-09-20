use std::fmt;

/// Five-bit amino-acid code.
pub type AminoAcid = u8;

pub const A: AminoAcid = 0;
pub const C: AminoAcid = 1;
pub const D: AminoAcid = 2;
pub const E: AminoAcid = 3;
pub const F: AminoAcid = 4;
pub const G: AminoAcid = 5;
pub const H: AminoAcid = 6;
pub const I: AminoAcid = 7;
pub const K: AminoAcid = 8;
pub const L: AminoAcid = 9;
pub const M: AminoAcid = 10;
pub const N: AminoAcid = 11;
pub const P: AminoAcid = 12;
pub const Q: AminoAcid = 13;
pub const R: AminoAcid = 14;
pub const S: AminoAcid = 15;
pub const T: AminoAcid = 16;
pub const V: AminoAcid = 17;
pub const W: AminoAcid = 18;
pub const Y: AminoAcid = 19;
pub const X: AminoAcid = 20;
pub const STOP: AminoAcid = 21;
pub const B: AminoAcid = 22;
pub const Z: AminoAcid = 23;
pub const J: AminoAcid = 24;
pub const U: AminoAcid = 25;
pub const O: AminoAcid = 26;
pub const GAP: AminoAcid = 27;

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

        Some((value & CODE_MASK) as AminoAcid)
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
        ALPHABET.get(code as usize).copied()
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
        self.packed[word] |= (code as u64) << shift;

        if shift > 64 - BITS_PER_RESIDUE {
            self.packed[word + 1] |= (code as u64) >> (64 - shift);
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
