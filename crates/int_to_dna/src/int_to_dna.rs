use onehot_dna::{OneHot, OneHotError, OneHotSequence};
use std::collections::BTreeMap;
use std::fmt;
use int_to_prot::IntToProt;
//use crate::errors::SeqError;
//use crate::traits::BinaryMatcher;

// logics copied from https://github.com/COMBINE-lab/kmers/
pub type Base = u8;
pub const A: Base = 0;
pub const C: Base = 1;
pub const G: Base = 2;
pub const T: Base = 3;

#[derive(Default, Debug, PartialEq, Clone)]
pub struct IntToDna {
    pub u8_encoded: Vec<u8>, // the 2bit encoded array (4 times compressed)
    pub lost: usize,         // how many times have I lost 4bp?
    pub size: usize,         // how long was the original sequence
    pub kmer_size: usize,
    step_size: usize,
    checker: BTreeMap<u8, usize>,
    mask: u64, //a mask to fill matching sequences to match the index's kmer_len
    current_position: usize,
}

trait ToLeBytes {
    fn to_le_byte_vec(&self) -> Vec<u8>;
}

impl ToLeBytes for u8 {
    fn to_le_byte_vec(&self) -> Vec<u8> {
        vec![*self]
    }
}

impl ToLeBytes for u16 {
    fn to_le_byte_vec(&self) -> Vec<u8> {
        self.to_le_bytes().to_vec()
    }
}

impl ToLeBytes for u32 {
    fn to_le_byte_vec(&self) -> Vec<u8> {
        self.to_le_bytes().to_vec()
    }
}

impl ToLeBytes for u64 {
    fn to_le_byte_vec(&self) -> Vec<u8> {
        self.to_le_bytes().to_vec()
    }
}

impl ToLeBytes for u128 {
    fn to_le_byte_vec(&self) -> Vec<u8> {
        self.to_le_bytes().to_vec()
    }
}

impl fmt::Display for IntToDna {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Helper: format u64 or usize into grouped binary string
        fn to_bin_grouped<T: Into<u64>>(num: T, bits: usize) -> String {
            let raw = format!("{:0width$b}", num.into(), width = bits);
            raw.chars()
                .collect::<Vec<_>>()
                .chunks(4)
                .map(|chunk| chunk.iter().collect::<String>())
                .collect::<Vec<_>>()
                .join("_")
        }

        // Helper: format u8 vector as binary
        fn vec_to_bin(vec: &[u8]) -> Vec<String> {
            vec.iter()
                .map(|b| format!("0b{}", to_bin_grouped(*b, 8)))
                .collect()
        }

        let dna_str = self.to_string(self.size);

        writeln!(f, "IntToDna {{")?;
        writeln!(f, "  stored human readable: \"{}\",", dna_str)?;
        writeln!(
            f,
            "  u8_encoded: {},",
            vec_to_bin(&self.u8_encoded).join(", ")
        )?;
        writeln!(f, "  lost: {},", self.lost)?;
        writeln!(f, "  size: {},", self.size)?;
        writeln!(f, "  checker: {:?},", self.checker)?;
        writeln!(f, "  mask: 0x{},", to_bin_grouped(self.mask, 64))?; // hex for clarity
        writeln!(f, "  current_position: {}", self.current_position)?;
        write!(f, "}}")
    }
}

// Implement the Index trait for MyClass
use std::ops::Index;

impl Index<usize> for IntToDna {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        self.u8_encoded.get(index).unwrap_or(&0_u8)
    }
}

impl PartialEq<Vec<u8>> for IntToDna {
    fn eq(&self, other: &Vec<u8>) -> bool {
        &self.u8_encoded == other
    }
}

/// This iterates over the bytes of the u8_encoded and therfore returns a u16, u64 pair at every 4 bp.
/// If you wnat more you need to create slices of this object first.
impl Iterator for IntToDna {
    type Item = (u16, u64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.size < (self.current_position + 4) * 4 {
            None
        } else {
            let start = self.current_position;
            let max = self.current_position + 10;
            self.current_position += 1;
            Some((
                Self::into_uint::<u16, 2>(&self.u8_encoded[start..max], false),
                Self::into_uint::<u64, 2>(&self.u8_encoded[(start + 2)..max], false),
            ))
        }
    }
}

/// Here I have my accessorie functions that more or less any of the classes would have use of.
/// I possibly learn a better way to have them...
impl IntToDna {
    pub fn new<T: AsRef<[u8]>>(input: T) -> Self {
        let seq = input.as_ref();
        // 4 of the array u8 fit into one result u8
        //eprintln!("Somtimes I die?! -> processed seq: {:?}", seq);

        let u8_encoded =
            Self::enc_bytes(seq).unwrap_or_else(|e| panic!("failed to encode sequence: {e}"));

        Self {
            u8_encoded,
            lost: 0,
            size: seq.len(),
            kmer_size: 16, //deprecated
            checker: BTreeMap::<u8, usize>::new(),
            mask: 0, // useless
            current_position: 0,
            step_size: 1,
        }
    }

    /// Construct directly from this crate's native packed 2-bit bytes.
    ///
    /// Base 0 occupies the least-significant two bits of byte 0. The final
    /// byte may contain unused high pairs when `size` is not divisible by 4.
    /// This is primarily useful for binary genome-format adapters that have
    /// already converted their alphabet/layout into IntToDna representation.
    pub fn from_packed_2bit(u8_encoded: Vec<u8>, size: usize) -> Self {
        assert_eq!(u8_encoded.len(), size.div_ceil(4), "packed byte count does not match sequence size");
        Self {
            u8_encoded,
            lost: 0,
            size,
            kmer_size: 16,
            checker: BTreeMap::<u8, usize>::new(),
            mask: 0,
            current_position: 0,
            step_size: 1,
        }
    }

    pub fn try_new<T: AsRef<[u8]>>(input: T) -> Result<Self, String> {
        let seq = input.as_ref();

        let u8_encoded = Self::enc_bytes(seq)?;

        Ok(Self {
            u8_encoded,
            lost: 0,
            size: seq.len(),
            kmer_size: 16,
            checker: BTreeMap::<u8, usize>::new(),
            mask: 0,
            current_position: 0,
            step_size: 1,
        })
    }

    pub fn enc_bytes(seq: &[u8]) -> Result<Vec<u8>, String> {
        let num_bytes = seq.len().div_ceil(4);
        let mut u8_encoded = Vec::<u8>::with_capacity(num_bytes);

        for id in 0..num_bytes {
            let mut packed = 0_u8;

            for add in (0..4).rev() {
                let current_utf8_id = id * 4 + add;
                packed <<= 2;

                if seq.len() > current_utf8_id {
                    let bits = Self::encode_binary(seq[current_utf8_id])?;
                    packed |= bits;
                }
            }

            u8_encoded.push(packed);
        }

        Ok(u8_encoded)
    }

    pub fn encode_binary(c: u8) -> Result<Base, String> {
        // might have to play some tricks for lookup in a const
        // array at some point
        match c {
            b'A' | b'a' => Ok(A),
            b'C' | b'c' => Ok(C),
            b'G' | b'g' => Ok(G),
            b'T' | b't' => Ok(T),
            b'N' | b'n' => Ok(A), // this is necessary as we can not even load a N containing sequence
            _ => Err("cannot encode {c} into 2 bit encoding".to_string()),
        }
    }

    /// Encode a DNA string (ACGTN) into a compact u64 using 2-bit encoding.
    ///
    /// A=0, C=1, G=2, T=3, N→A (0)
    ///
    /// # Limits
    /// Maximum length is 32 bases (64 bits).
    pub fn str_to_u64(seq: &str) -> Result<u64, String> {
        if seq.len() > 32 {
            return Err(format!(
                "sequence too long for u64 encoding: {} > 32",
                seq.len()
            ));
        }

        let mut value = 0u64;

        for &b in seq.as_bytes() {
            let base = Self::encode_binary(b)?;
            value = (value << 2) | (base as u64);
        }

        Ok(value)
    }

    pub fn is_empty(&self) -> bool {
        self.u8_encoded.is_empty()
    }

    pub fn len(&self) -> usize {
        self.u8_encoded.len()
    }

    /// Return the first packed 4-base byte.
    ///
    /// Empty sequences return zero, matching the existing `Index` behaviour.
    #[inline]
    pub fn first_u8(&self) -> u8 {
        self.u8_encoded.first().copied().unwrap_or(0)
    }

    /// Pack four bases starting at an arbitrary base offset.
    ///
    /// This reads directly from the existing 2-bit representation, so callers
    /// can scan a packed sequence without allocating 4-base DNA slices or
    /// re-encoding them. Returns `None` when fewer than four bases remain.
    #[inline]
    pub fn packed_u8_at(&self, base_offset: usize) -> Option<u8> {
        if base_offset.checked_add(4)? > self.size {
            return None;
        }

        let byte = base_offset / 4;
        let shift = (base_offset % 4) * 2;
        if shift == 0 {
            return self.u8_encoded.get(byte).copied();
        }

        let low = self.u8_encoded[byte] >> shift;
        let high = self.u8_encoded.get(byte + 1).copied().unwrap_or(0) << (8 - shift);
        Some(low | high)
    }

    /// Find the first packed 4-base anchor whose immediately following packed
    /// evidence matches `external`. Returns the anchor's exact base position.
    ///
    /// The scan advances one base at a time, so all four packed reading frames
    /// are covered without constructing shifted `IntToDna` objects. Empty
    /// external evidence is allowed and reduces this to an anchor search.
    #[inline]
    pub fn find_with_external_after(&self, lookup_id: u8, external: &[u8]) -> Option<usize> {
        let needed = 4usize.checked_mul(external.len().checked_add(1)?)?;
        if self.size < needed {
            return None;
        }

        for pos in 0..=self.size - needed {
            if self.packed_u8_at(pos) != Some(lookup_id) {
                continue;
            }
            if external
                .iter()
                .enumerate()
                .all(|(i, expected)| self.packed_u8_at(pos + (i + 1) * 4) == Some(*expected))
            {
                return Some(pos);
            }
        }
        None
    }

    /// Find the first packed 4-base anchor whose immediately preceding packed
    /// evidence matches `external`. Returns the anchor's exact base position.
    ///
    /// `external` is supplied in sequence order, i.e. `[a, b]` matches
    /// `[a][b][lookup_id]`. Any evidence length is accepted.
    #[inline]
    pub fn find_with_external_before(&self, lookup_id: u8, external: &[u8]) -> Option<usize> {
        let prefix_bases = 4usize.checked_mul(external.len())?;
        if self.size < prefix_bases + 4 {
            return None;
        }

        for anchor_pos in prefix_bases..=self.size - 4 {
            if self.packed_u8_at(anchor_pos) != Some(lookup_id) {
                continue;
            }
            let start = anchor_pos - prefix_bases;
            if external
                .iter()
                .enumerate()
                .all(|(i, expected)| self.packed_u8_at(start + i * 4) == Some(*expected))
            {
                return Some(anchor_pos);
            }
        }
        None
    }

    /// Return the maximally informative packed byte for the final four bases.
    ///
    /// `u8_encoded.last()` is only maximally informative when the sequence
    /// length is divisible by four. For e.g. a 9 bp barcode its final packed
    /// byte contains one real base plus three padded A bases. This method
    /// instead repacks the last four *real* bases, without converting back to
    /// a DNA string.
    ///
    /// For sequences shorter than four bases this returns their existing
    /// packed byte. Empty sequences return zero.
    #[inline]
    pub fn last_informative_u8(&self) -> u8 {
        if self.size == 0 {
            return 0;
        }

        if self.size <= 4 {
            return self.first_u8();
        }

        let start = self.size - 4;
        let mut packed = 0u8;

        for target_pos in 0..4 {
            let source_pos = start + target_pos;
            let source_byte = self.u8_encoded[source_pos / 4];
            let base = (source_byte >> ((source_pos % 4) * 2)) & 0b11;
            packed |= base << (target_pos * 2);
        }

        packed
    }

    /// Convert this already 2-bit-packed sequence directly to `OneHot<N>`.
    ///
    /// This avoids `IntToDna -> DNA bytes -> OneHot`. The sequence must have
    /// exactly `N` meaningful bases. Note that `IntToDna` historically maps
    /// `N` to `A` during 2-bit encoding; callers that need unknown bases to
    /// remain mismatches must reject them before constructing `IntToDna`.
    /// Expand this 2-bit sequence directly into a read-sized OneHotSequence.
    /// No intermediate A/C/G/T byte string is allocated.
    #[inline]
    pub fn as_one_hot_sequence(&self) -> OneHotSequence {
        OneHotSequence::from_2bit_bytes(&self.u8_encoded, self.size)
    }

    pub fn as_one_hot<const N: usize>(&self) -> Result<OneHot<N>, OneHotError> {
        if N > OneHot::<N>::MAX_LEN {
            return Err(OneHotError::TooLong {
                max: OneHot::<N>::MAX_LEN,
                observed: N,
            });
        }

        if self.size != N {
            return Err(OneHotError::WrongLength {
                expected: N,
                observed: self.size,
            });
        }

        let mut bits = 0u128;

        for pos in 0..N {
            let packed = self.u8_encoded[pos / 4];
            let base = (packed >> ((pos % 4) * 2)) & 0b11;
            bits |= (1u128 << base) << (pos * 4);
        }

        Ok(OneHot::<N>::from_bits(bits))
    }

    fn reverse_bits_in_byte(b: u8) -> u8 {
        let mut b = b;
        b = b.rotate_left(4);
        b = ((b & 0b1100_1100) >> 2) | ((b & 0b0011_0011) << 2);
        b = ((b & 0b1010_1010) >> 1) | ((b & 0b0101_0101) << 1);
        b
    }

    ///helper function for crating the unsigned intergers from the internal u8 array.
    fn into_uint<T, const N: usize>(encoded: &[u8], reverse_bits: bool) -> T
    where
        T: Default + From<u8> + std::ops::Shl<u8, Output = T> + std::ops::BitOr<Output = T>,
    {
        let mut ret = T::default();

        for i in (0..N).rev() {
            let mut byte = encoded.get(i).copied().unwrap_or(0);
            if reverse_bits {
                byte = Self::reverse_bits_in_byte(byte);
            }

            ret = (ret << 8) | T::from(byte);
        }

        ret
    }
    /// this is needed for the initial mappers
    /// takes the self encoded seqences and converts the last entries into one u16
    pub fn into_u16(&self) -> u16 {
        Self::into_uint::<u16, 2>(&self.u8_encoded, false)
    }
    /// this is needed for completing the set
    pub fn into_u32(&self) -> u32 {
        Self::into_uint::<u32, 4>(&self.u8_encoded, false)
    }
    /// needed for the secondary mappings
    /// takes the UTF8 encoded sequence and encodes the first 32 into a u64
    pub fn into_u64(&self) -> u64 {
        Self::into_uint::<u64, 8>(&self.u8_encoded, false)
    }
    /// this is needed for completing the set
    pub fn into_u128(&self) -> u128 {
        Self::into_uint::<u128, 16>(&self.u8_encoded, false)
    }

    /// reverse into_uint - create the internal structure from the numbers
    fn from_uint<T, const N: usize>(value: T, reverse_bits: bool) -> Self
    where
        T: Copy
            + ToLeBytes
            + std::ops::Shr<u8, Output = T>
            + std::ops::BitAnd<Output = T>
            + From<u8>
            + Into<u128>, // For size handling
    {
        let bytes = value.to_le_byte_vec();

        // Ensure we have exactly N bytes, padding with 0s if needed
        let mut buf = vec![0u8; N];
        for (i, b) in bytes.iter().enumerate().take(N) {
            buf[i] = *b;
        }

        if reverse_bits {
            for b in &mut buf {
                *b = Self::reverse_bits_in_byte(*b);
            }
        }

        IntToDna {
            u8_encoded: buf,
            lost: 0,
            kmer_size: 16,
            size: bytes.len() * 4,
            step_size: 3,
            checker: BTreeMap::new(),
            mask: 0_u64,
            current_position: 0,
        }
    }

    pub fn from_u16(val: u16) -> Self {
        Self::from_uint::<u16, 2>(val, false)
    }

    pub fn from_u32(val: u32) -> Self {
        Self::from_uint::<u32, 4>(val, false)
    }

    pub fn from_u64(val: u64) -> Self {
        Self::from_uint::<u64, 8>(val, false)
    }

    pub fn from_u128(val: u128) -> Self {
        Self::from_uint::<u128, 16>(val, false)
    }

    pub fn ints_at_position(&self, position: usize) -> Result<(u16, u64), String> {
        if !position.is_multiple_of(4) {
            return Err("get yourself a sliced IntToSeq!".to_string());
        }
        let start = position / 4;
        if start + 3 < self.len() {
            return Err("I do not have enough data to do this!".to_string());
        }
        let max = (start + 10).min(self.len());
        Ok((
            Self::into_uint::<u16, 2>(&self.u8_encoded[start..max], false),
            Self::into_uint::<u64, 2>(&self.u8_encoded[(start + 2)..max], false),
        ))
    }

    /// Translate this DNA sequence with the standard genetic code in frame 0.
    pub fn translate(&self) -> IntToProt {
        self.translate_frame(0)
    }

    /// Translate this DNA sequence with the standard genetic code in frame 0, 1 or 2.
    /// Trailing bases that do not form a complete codon are ignored.
    pub fn translate_frame(&self, frame: usize) -> IntToProt {
        assert!(frame < 3, "translation frame must be 0, 1 or 2");
        let dna = self.to_string(self.size);
        let bytes = dna.as_bytes();
        let mut protein = Vec::with_capacity(bytes.len().saturating_sub(frame) / 3);
        let mut pos = frame;
        while pos + 3 <= bytes.len() {
            protein.push(translate_codon(&bytes[pos..pos + 3]));
            pos += 3;
        }
        IntToProt::new(protein)
    }

    /// get a new IntToSeq objects from start..end of the DNA sequence.
    pub fn slice(&self, start: Option<usize>, end: Option<usize>) -> Self {
        let start: usize = start.unwrap_or(0);
        let end: usize = end.unwrap_or(self.size);

        let u8_enc = &self.u8_encoded[(start / 4)..end];
        let seq = Self::u8_array_to_str(u8_enc);
        //let shift = start % 4;

        Self::new(&seq.as_bytes()[start..end])
    }

    /// Convert any integer value obtained from any of the self.into_u16.to_le_bytes() like functions
    /// and adds them at the end of the mutable data string.
    pub fn to_string(&self, bases: usize) -> String {
        let num_bytes = bases.div_ceil(4).min(self.u8_encoded.len());
        let mut data = Self::u8_array_to_str(&self.u8_encoded[0..num_bytes]);
        data.truncate(bases);
        data
    }

    pub fn u8_array_to_str(u8_encoded: &[u8]) -> String {
        let mut data = String::default();
        for u8_4bp in u8_encoded.iter() {
            Self::u8_to_str(u8_4bp, &mut data);
        }
        data
    }

    pub fn u8_to_str(u8_rep: &u8, data: &mut String) {
        let mut loc: u8 = *u8_rep;
        //println!("converting u8 {loc:b} to string with {kmer_size} bp.");
        for _i in (0..4).rev() {
            let ch = match loc & 0b11 {
                0 => "A", //0b00
                1 => "C", // 0b01
                2 => "G", // 0b10
                3 => "T", // 0b11
                _ => "N",
            };
            *data += ch;
            loc >>= 2;
            //println!("{ch} and loc {loc:b}");
        }

        //println!("\nMakes sense? {:?}", data);
    }

    pub fn print(&self) {
        let max_nuc = self.len() * 4;
        println!(
            ">seq (n={} [{} meaningful bp])\n{}",
            max_nuc,
            self.size,
            self.to_string(self.size)
        );
        print!("[ ");
        for en in &self.u8_encoded {
            print!(", {en:b}");
        }
        println!("]");
    }
}

#[cfg(test)]
mod tests {
    use crate::IntToDna;

    fn format_bytes_binary(bytes: &[u8]) -> String {
        bytes
            .iter()
            .map(|b| format!("{:08b}", b))
            .collect::<Vec<_>>()
            .join(", ")
    }
    fn roundtrip_test(input: &str) {
        let used_str = if input.len() < 64 {
            let padding = "A".repeat(64 - input.len());
            format!("{input}{padding}")
        } else {
            input.to_string()
        };
        println!("==> Testing input: {used_str}");
        let encoder = IntToDna::new(used_str.as_bytes());
        let original = encoder.u8_encoded.clone();
        println!("Original encoded: {:?}", original);

        // Encode
        let val16 = encoder.into_u16();
        let val32 = encoder.into_u32();
        let val64 = encoder.into_u64();
        let val128 = encoder.into_u128();

        println!("Encoded u16:   {:b}", val16);
        println!("Encoded u32:   {:b}", val32);
        println!("Encoded u64:   {:b}", val64);
        println!("Encoded u128:  {:b}", val128);

        //panic!("Will die anyhow");
        // Decode
        let decoded16 = IntToDna::from_u16(val16);
        let decoded32 = IntToDna::from_u32(val32);
        let decoded64 = IntToDna::from_u64(val64);
        let decoded128 = IntToDna::from_u128(val128);

        println!(
            "Decoded u16:   {}",
            format_bytes_binary(&decoded16.u8_encoded)
        );
        println!(
            "Decoded u32:   {}",
            format_bytes_binary(&decoded32.u8_encoded)
        );
        println!(
            "Decoded u64:   {}",
            format_bytes_binary(&decoded64.u8_encoded)
        );
        println!(
            "Decoded u128:  {}",
            format_bytes_binary(&decoded128.u8_encoded)
        );

        // Assertions
        assert_eq!(&original[..2], &decoded16.u8_encoded[..2], "u16 mismatch");
        assert_eq!(&original[..4], &decoded32.u8_encoded[..4], "u32 mismatch");
        assert_eq!(&original[..8], &decoded64.u8_encoded[..8], "u64 mismatch");
        assert_eq!(
            &original[..16],
            &decoded128.u8_encoded[..16],
            "u128 mismatch"
        );

        //assert_eq!(&encoder.storage[..8], &decoded16.storage[..8], "u16 mismatch");
        //assert_eq!(&encoder.storage[..16], &decoded32.storage[..16], "u32 mismatch");
        //assert_eq!(&encoder.storage[..32], &decoded64.storage[..32], "u64 mismatch");
        assert_eq!(&decoded128.to_string(input.len()), input);
    }
    #[test]
    fn test_encoding() {
        let encoder = IntToDna::new(b"ACGTTC");
        assert_eq!(
            encoder.into_u16(),
            0b011111100100,
            "encoded {:b} is not the expected {:b}",
            encoder.into_u16(),
            0b011111100100
        );
        assert_eq!(&encoder.to_string(5), "ACGTT");
        //this will only overreach by max 3 bp of A's (0)
        assert_eq!(&encoder.to_string(15), "ACGTTCAA");
    }

    #[test]
    fn test_from_uint() {
        // the from_uint function is the 'mother' of all fromn16 to from_u128.
        // So testing one should also test all others.
        let encoder = IntToDna::from_u16(0b1111100100_u16);
        let exp = IntToDna::new(b"ACGTT");
        assert_eq!(encoder.u8_encoded, exp.u8_encoded);
    }

    #[test]
    fn test_uxx_roundtrip_aaaaaaacaaagaaat() {
        roundtrip_test("AAAAAAACAAAGAAAT");
    }

    #[test]
    fn test_uxx_roundtrip_aaaaaacaagaataaa() {
        roundtrip_test("AAAAAACAAGAATAAA");
    }

    #[test]
    fn packed_u8_at_reads_unaligned_four_base_windows() {
        let encoded = IntToDna::new(b"AACGTTGCA");

        assert_eq!(encoded.packed_u8_at(0), Some(IntToDna::new(b"AACG").first_u8()));
        assert_eq!(encoded.packed_u8_at(1), Some(IntToDna::new(b"ACGT").first_u8()));
        assert_eq!(encoded.packed_u8_at(4), Some(IntToDna::new(b"TTGC").first_u8()));
        assert_eq!(encoded.packed_u8_at(5), Some(IntToDna::new(b"TGCA").first_u8()));
        assert_eq!(encoded.packed_u8_at(6), None);
    }

    #[test]
    fn packed_external_search_covers_all_frames_and_returns_anchor_position() {
        let encoded = IntToDna::new(b"TTACGTTGCAGGAACC");
        let pattern = IntToDna::new(b"ACGTTGCAGGAA");
        let external = [pattern.packed_u8_at(0).unwrap(), pattern.packed_u8_at(4).unwrap()];
        let lookup = pattern.packed_u8_at(8).unwrap();

        assert_eq!(encoded.find_with_external_before(lookup, &external), Some(10));

        let after_external = [pattern.packed_u8_at(4).unwrap(), pattern.packed_u8_at(8).unwrap()];
        let first = pattern.packed_u8_at(0).unwrap();
        assert_eq!(encoded.find_with_external_after(first, &after_external), Some(2));
    }

    #[test]
    fn packed_external_search_accepts_only_the_evidence_supplied() {
        let encoded = IntToDna::new(b"GGGGACGTCCCCACGT");
        let lookup = IntToDna::new(b"ACGT").first_u8();
        let gggg = IntToDna::new(b"GGGG").first_u8();
        let cccc = IntToDna::new(b"CCCC").first_u8();

        assert_eq!(encoded.find_with_external_before(lookup, &[]), Some(4));
        assert_eq!(encoded.find_with_external_before(lookup, &[gggg]), Some(4));
        assert_eq!(encoded.find_with_external_before(lookup, &[cccc]), Some(12));
        assert_eq!(encoded.find_with_external_after(lookup, &[cccc]), Some(4));
    }

    #[test]
    fn last_informative_u8_uses_last_four_real_bases() {
        let encoded = IntToDna::new(b"GTCAGCTAC");
        let expected = IntToDna::new(b"CTAC");

        assert_eq!(encoded.last_informative_u8(), expected.first_u8());
        assert_ne!(
            encoded.last_informative_u8(),
            *encoded.u8_encoded.last().unwrap()
        );
    }

    #[test]
    fn first_and_last_u8_can_be_identical_without_padding_artifacts() {
        let encoded = IntToDna::new(b"AGCTAGCT");
        let agct = IntToDna::new(b"AGCT").first_u8();

        assert_eq!(encoded.first_u8(), agct);
        assert_eq!(encoded.last_informative_u8(), agct);
    }

    #[test]
    fn as_one_hot_matches_direct_encoding_without_string_roundtrip() {
        let seq = b"GTCAGCTAC";
        let encoded = IntToDna::new(seq);
        let direct = onehot_dna::OneHot::<9>::from_bytes(seq).unwrap();

        assert_eq!(encoded.as_one_hot::<9>().unwrap(), direct);
    }

    #[test]
    fn as_one_hot_rejects_wrong_const_length() {
        let encoded = IntToDna::new(b"GTCAGCTAC");
        assert!(matches!(
            encoded.as_one_hot::<8>(),
            Err(onehot_dna::OneHotError::WrongLength {
                expected: 8,
                observed: 9
            })
        ));
    }
}

#[inline]
fn translate_codon(codon: &[u8]) -> u8 {
    match codon {
        b"TTT" | b"TTC" => b'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"ATG" => b'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"TAT" | b"TAC" => b'Y',
        b"TAA" | b"TAG" | b"TGA" => b'*',
        b"CAT" | b"CAC" => b'H',
        b"CAA" | b"CAG" => b'Q',
        b"AAT" | b"AAC" => b'N',
        b"AAA" | b"AAG" => b'K',
        b"GAT" | b"GAC" => b'D',
        b"GAA" | b"GAG" => b'E',
        b"TGT" | b"TGC" => b'C',
        b"TGG" => b'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
        _ => b'X',
    }
}

#[cfg(test)]
mod translation_tests {
    use super::*;

    #[test]
    fn translates_standard_code_in_frame_zero() {
        let dna = IntToDna::new(b"ATGGCTGAATTTTAA");
        assert_eq!(dna.translate().to_string(), "MAEF*");
    }

    #[test]
    fn translates_requested_frame() {
        let dna = IntToDna::new(b"AATGGCTGAATTT");
        assert_eq!(dna.translate_frame(1).to_string(), "MAEF");
    }
}
