use std::collections::BTreeMap;
use std::fmt;
//use crate::errors::SeqError;
//use crate::traits::BinaryMatcher;

// logics copied from https://github.com/COMBINE-lab/kmers/
pub type Base = u8;
pub const A: Base = 0;
pub const C: Base = 1;
pub const G: Base = 2;
pub const T: Base = 3;

#[derive(Default, Debug, PartialEq, Clone)]
pub struct IntToStr {
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

impl fmt::Display for IntToStr {
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

        writeln!(f, "IntToStr {{")?;
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

impl Index<usize> for IntToStr {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        self.u8_encoded.get(index).unwrap_or(&0_u8)
    }
}

impl PartialEq<Vec<u8>> for IntToStr {
    fn eq(&self, other: &Vec<u8>) -> bool {
        &self.u8_encoded == other
    }
}

/// This iterates over the bytes of the u8_encoded and therfore returns a u16, u64 pair at every 4 bp.
/// If you wnat more you need to create slices of this object first.
impl Iterator for IntToStr {
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
impl IntToStr {
    pub fn new<T: AsRef<[u8]>>(input: T) -> Self {
        let seq = input.as_ref();
        // 4 of the array u8 fit into one result u8
        //eprintln!("Somtimes I die?! -> processed seq: {:?}", seq);
        
        let u8_encoded = Self::enc_bytes(seq)
        .unwrap_or_else(|e| panic!("failed to encode sequence: {e}"));

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

        IntToStr {
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
    use crate::IntToStr;

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
        let encoder = IntToStr::new(used_str.as_bytes());
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
        let decoded16 = IntToStr::from_u16(val16);
        let decoded32 = IntToStr::from_u32(val32);
        let decoded64 = IntToStr::from_u64(val64);
        let decoded128 = IntToStr::from_u128(val128);

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
        let encoder = IntToStr::new(b"ACGTTC");
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
        let encoder = IntToStr::from_u16(0b1111100100_u16);
        let exp = IntToStr::new(b"ACGTT");
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
}
