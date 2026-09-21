use crate::IntToDna;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

const TWO_BIT_SIGNATURE: u32 = 0x1A41_2743;

#[derive(Debug, Clone, Copy)]
enum ByteOrder {
    Big,
    Little,
}

impl ByteOrder {
    fn u32(self, bytes: [u8; 4]) -> u32 {
        match self {
            Self::Big => u32::from_be_bytes(bytes),
            Self::Little => u32::from_le_bytes(bytes),
        }
    }
}

#[derive(Debug, Clone)]
struct SequenceRecord {
    dna_size: u32,
    n_blocks: Vec<(u32, u32)>,
    packed_offset: u64,
}

/// Random-access reader for the UCSC `.2bit` genome format.
///
/// UCSC stores T=00, C=01, A=10, G=11 with the first base in the most
/// significant pair of each byte. `IntToDna` stores A=00, C=01, G=10, T=11
/// with the first base in the least significant pair. Conversion therefore
/// changes both the alphabet and the pair order; it is performed a byte at a
/// time without creating an intermediate DNA string.
pub struct TwoBitReader {
    reader: BufReader<File>,
    order: ByteOrder,
    index: HashMap<String, u64>,
    records: HashMap<String, SequenceRecord>,
}

impl TwoBitReader {
    pub fn open(path: impl AsRef<Path>) -> io::Result<Self> {
        let mut reader = BufReader::new(File::open(path)?);
        let mut signature = [0u8; 4];
        reader.read_exact(&mut signature)?;
        let order = if u32::from_be_bytes(signature) == TWO_BIT_SIGNATURE {
            ByteOrder::Big
        } else if u32::from_le_bytes(signature) == TWO_BIT_SIGNATURE {
            ByteOrder::Little
        } else {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "not a UCSC twoBit file"));
        };

        let version = read_u32(&mut reader, order)?;
        if version != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("unsupported UCSC twoBit version {version}"),
            ));
        }
        let sequence_count = read_u32(&mut reader, order)?;
        let reserved = read_u32(&mut reader, order)?;
        if reserved != 0 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "twoBit header reserved field is not zero"));
        }

        let mut index = HashMap::with_capacity(sequence_count as usize);
        for _ in 0..sequence_count {
            let mut name_len = [0u8; 1];
            reader.read_exact(&mut name_len)?;
            let mut name = vec![0u8; name_len[0] as usize];
            reader.read_exact(&mut name)?;
            let name = String::from_utf8(name)
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "twoBit sequence name is not UTF-8"))?;
            let offset = read_u32(&mut reader, order)? as u64;
            index.insert(name, offset);
        }

        Ok(Self { reader, order, index, records: HashMap::new() })
    }

    pub fn sequence_names(&self) -> impl Iterator<Item = &str> {
        self.index.keys().map(String::as_str)
    }

    pub fn sequence_len(&mut self, name: &str) -> io::Result<u32> {
        Ok(self.record(name)?.dna_size)
    }

    /// Read a zero-based, half-open genomic interval directly into `IntToDna`.
    ///
    /// UCSC N-blocks are represented as `A`, matching the existing `IntToDna`
    /// policy for ambiguous N bases. Soft-mask blocks deliberately do not alter
    /// the result because `IntToDna` has no case/masking state.
    pub fn sequence(&mut self, name: &str, start: u32, end: u32) -> io::Result<IntToDna> {
        let record = self.record(name)?.clone();
        if start > end || end > record.dna_size {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("invalid interval {name}:{start}-{end}; sequence length is {}", record.dna_size),
            ));
        }
        if start == end {
            return Ok(IntToDna::new(b""));
        }

        let first_byte = start / 4;
        let last_byte = end.div_ceil(4);
        let mut packed = vec![0u8; (last_byte - first_byte) as usize];
        self.reader.seek(SeekFrom::Start(record.packed_offset + first_byte as u64))?;
        self.reader.read_exact(&mut packed)?;

        // Convert complete UCSC bytes to IntToDna bytes first. Boundary bases
        // are then selected from this packed representation without strings.
        for byte in &mut packed {
            *byte = ucsc_byte_to_int_to_dna(*byte);
        }

        let wanted = (end - start) as usize;
        let offset = (start % 4) as usize;
        let mut out = Vec::with_capacity(wanted.div_ceil(4));
        for i in 0..wanted {
            let absolute = offset + i;
            let src = packed[absolute / 4];
            let base = (src >> (2 * (absolute % 4))) & 0b11;
            set_packed_base(&mut out, i, base);
        }

        // N has no native representation in IntToDna. Preserve its established
        // N->A policy, but use the authoritative UCSC N-block table rather than
        // whatever placeholder bits happen to be present in the packed DNA.
        for &(n_start, n_size) in &record.n_blocks {
            let n_end = n_start.saturating_add(n_size);
            let overlap_start = start.max(n_start);
            let overlap_end = end.min(n_end);
            for pos in overlap_start..overlap_end {
                set_existing_packed_base(&mut out, (pos - start) as usize, 0);
            }
        }

        Ok(IntToDna::from_packed_2bit(out, wanted))
    }

    fn record(&mut self, name: &str) -> io::Result<&SequenceRecord> {
        if !self.records.contains_key(name) {
            let offset = *self.index.get(name).ok_or_else(|| {
                io::Error::new(io::ErrorKind::NotFound, format!("sequence {name:?} not present in twoBit file"))
            })?;
            self.reader.seek(SeekFrom::Start(offset))?;
            let dna_size = read_u32(&mut self.reader, self.order)?;
            let n_count = read_u32(&mut self.reader, self.order)? as usize;
            let mut n_starts = Vec::with_capacity(n_count);
            for _ in 0..n_count { n_starts.push(read_u32(&mut self.reader, self.order)?); }
            let mut n_sizes = Vec::with_capacity(n_count);
            for _ in 0..n_count { n_sizes.push(read_u32(&mut self.reader, self.order)?); }
            let n_blocks = n_starts.into_iter().zip(n_sizes).collect();

            let mask_count = read_u32(&mut self.reader, self.order)? as u64;
            self.reader.seek(SeekFrom::Current((mask_count * 4) as i64))?; // mask starts
            self.reader.seek(SeekFrom::Current((mask_count * 4) as i64))?; // mask sizes
            let reserved = read_u32(&mut self.reader, self.order)?;
            if reserved != 0 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "twoBit sequence reserved field is not zero"));
            }
            let packed_offset = self.reader.stream_position()?;
            self.records.insert(name.to_owned(), SequenceRecord { dna_size, n_blocks, packed_offset });
        }
        Ok(self.records.get(name).expect("record was inserted"))
    }
}

fn read_u32(reader: &mut impl Read, order: ByteOrder) -> io::Result<u32> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(order.u32(bytes))
}

const fn build_ucsc_byte_lut() -> [u8; 256] {
    let mut lut = [0u8; 256];
    let mut byte = 0usize;
    while byte < 256 {
        let mut out = 0u8;
        let mut i = 0usize;
        while i < 4 {
            let ucsc = ((byte as u8) >> (6 - 2 * i)) & 0b11;
            let base = match ucsc {
                0 => 3, // T
                1 => 1, // C
                2 => 0, // A
                3 => 2, // G
                _ => 0,
            };
            out |= base << (2 * i);
            i += 1;
        }
        lut[byte] = out;
        byte += 1;
    }
    lut
}

const UCSC_BYTE_TO_INT_TO_DNA: [u8; 256] = build_ucsc_byte_lut();

#[inline]
fn ucsc_byte_to_int_to_dna(byte: u8) -> u8 {
    UCSC_BYTE_TO_INT_TO_DNA[byte as usize]
}

#[inline]
fn set_packed_base(out: &mut Vec<u8>, position: usize, base: u8) {
    let byte = position / 4;
    if byte == out.len() { out.push(0); }
    out[byte] |= base << (2 * (position % 4));
}

#[inline]
fn set_existing_packed_base(out: &mut [u8], position: usize, base: u8) {
    let byte = position / 4;
    let shift = 2 * (position % 4);
    out[byte] = (out[byte] & !(0b11 << shift)) | (base << shift);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;

    fn write_u32(out: &mut Vec<u8>, value: u32, little: bool) {
        let bytes = if little { value.to_le_bytes() } else { value.to_be_bytes() };
        out.extend_from_slice(&bytes);
    }

    fn tiny_two_bit(path: &Path, little: bool) {
        // chrTest = ACGTNNGTA. Packed placeholder bases for Ns are T; the N
        // block must override those bits when read.
        let name = b"chrTest";
        let header_and_index = 16 + 1 + name.len() + 4;
        let record_offset = header_and_index as u32;
        let mut data = Vec::new();
        write_u32(&mut data, TWO_BIT_SIGNATURE, little);
        write_u32(&mut data, 0, little);
        write_u32(&mut data, 1, little);
        write_u32(&mut data, 0, little);
        data.push(name.len() as u8);
        data.extend_from_slice(name);
        write_u32(&mut data, record_offset, little);

        write_u32(&mut data, 9, little); // dna size
        write_u32(&mut data, 1, little); // N blocks
        write_u32(&mut data, 4, little); // N start
        write_u32(&mut data, 2, little); // N size
        write_u32(&mut data, 0, little); // mask blocks
        write_u32(&mut data, 0, little); // reserved
        // A C G T | T T G T | A in UCSC encoding, high pair first.
        data.extend_from_slice(&[0b10_01_11_00, 0b00_00_11_00, 0b10_00_00_00]);

        let mut file = File::create(path).unwrap();
        file.write_all(&data).unwrap();
    }

    #[test]
    fn ucsc_byte_conversion_matches_int_to_dna_layout() {
        let converted = ucsc_byte_to_int_to_dna(0b10_01_11_00); // ACGT
        assert_eq!(IntToDna::from_packed_2bit(vec![converted], 4).to_string(4), "ACGT");
    }

    #[test]
    fn reads_full_and_unaligned_ranges_and_applies_n_blocks() {
        let path = std::env::temp_dir().join(format!("int-to-dna-{}.2bit", std::process::id()));
        tiny_two_bit(&path, false);
        let mut reader = TwoBitReader::open(&path).unwrap();
        assert_eq!(reader.sequence_len("chrTest").unwrap(), 9);
        assert_eq!(reader.sequence("chrTest", 0, 9).unwrap().to_string(9), "ACGTAAGTA");
        assert_eq!(reader.sequence("chrTest", 1, 8).unwrap().to_string(7), "CGTAAGT");
        assert!(reader.sequence("chrTest", 0, 10).is_err());
        fs::remove_file(path).ok();
    }

    #[test]
    fn accepts_byte_swapped_two_bit_files() {
        let path = std::env::temp_dir().join(format!("int-to-dna-le-{}.2bit", std::process::id()));
        tiny_two_bit(&path, true);
        let mut reader = TwoBitReader::open(&path).unwrap();
        assert_eq!(reader.sequence("chrTest", 0, 4).unwrap().to_string(4), "ACGT");
        fs::remove_file(path).ok();
    }
}
