use anyhow::{Result, bail};
use rust_htslib::bam::Record;

/// All consecutive BAM records belonging to one query name.
///
/// This is deliberately alignment-preserving: primary, secondary and
/// supplementary records remain in the group so consumers such as a future
/// fusion detector can inspect the complete fragment evidence.
#[derive(Debug)]
pub struct ReadGroup {
    qname: Vec<u8>,
    records: Vec<Record>,
}

impl ReadGroup {
    pub fn new(first: Record) -> Self {
        Self {
            qname: first.qname().to_vec(),
            records: vec![first],
        }
    }

    pub fn qname(&self) -> &[u8] {
        &self.qname
    }

    pub fn records(&self) -> &[Record] {
        &self.records
    }

    pub fn records_mut(&mut self) -> &mut [Record] {
        &mut self.records
    }

    pub fn push(&mut self, record: Record) -> Result<()> {
        if record.qname() != self.qname.as_slice() {
            bail!(
                "attempted to add QNAME '{}' to read group '{}'",
                String::from_utf8_lossy(record.qname()),
                String::from_utf8_lossy(&self.qname),
            );
        }
        self.records.push(record);
        Ok(())
    }

    /// Recover the original sequenced R1/R2 sequences from the primary BAM
    /// records. BAM stores SEQ reverse-complemented for reverse-strand
    /// alignments, therefore we normalize it back before molecule identity is
    /// calculated. This keeps FASTQ and BAM identity generation equivalent.
    pub fn primary_pair_sequences(&self) -> Result<(Vec<u8>, Vec<u8>)> {
        let mut r1: Option<Vec<u8>> = None;
        let mut r2: Option<Vec<u8>> = None;

        for record in &self.records {
            if record.is_secondary() || record.is_supplementary() {
                continue;
            }

            let target = if record.is_first_in_template() {
                &mut r1
            } else if record.is_last_in_template() {
                &mut r2
            } else {
                continue;
            };

            if target.is_some() {
                bail!(
                    "QNAME '{}' contains more than one primary alignment for the same mate",
                    String::from_utf8_lossy(&self.qname),
                );
            }

            *target = Some(original_sequence(record));
        }

        let r1 = r1.ok_or_else(|| {
            anyhow::anyhow!(
                "QNAME '{}' has no primary R1 record; unbarcoded BAM quantification requires paired reads",
                String::from_utf8_lossy(&self.qname),
            )
        })?;
        let r2 = r2.ok_or_else(|| {
            anyhow::anyhow!(
                "QNAME '{}' has no primary R2 record; unbarcoded BAM quantification requires paired reads",
                String::from_utf8_lossy(&self.qname),
            )
        })?;

        Ok((r1, r2))
    }
}

pub fn original_sequence(record: &Record) -> Vec<u8> {
    let mut seq = record.seq().as_bytes();
    if record.is_reverse() {
        reverse_complement_in_place(&mut seq);
    }
    seq
}

fn reverse_complement_in_place(seq: &mut [u8]) {
    seq.reverse();
    for base in seq {
        *base = match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            other => other,
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    fn record(qname: &[u8], seq: &[u8], flags: u16) -> Record {
        let mut record = Record::new();
        let cigar = CigarString(vec![]);
        let qual = vec![b'I' - 33; seq.len()];
        record.set(qname, Some(&cigar), seq, &qual);
        record.set_flags(flags);
        record
    }

    #[test]
    fn read_group_keeps_all_alignments_but_pair_identity_uses_primaries() {
        // paired + first
        let mut group = ReadGroup::new(record(b"frag1", b"ACGTACGTACGTACGT", 0x1 | 0x40));
        // supplementary R1 must remain in the group but not replace primary R1
        group
            .push(record(b"frag1", b"TTTTTTTTTTTTTTTT", 0x1 | 0x40 | 0x800))
            .unwrap();
        // paired + second
        group
            .push(record(b"frag1", b"CCCCCCCCCCCCCCCC", 0x1 | 0x80))
            .unwrap();

        assert_eq!(group.records().len(), 3);
        let (r1, r2) = group.primary_pair_sequences().unwrap();
        assert_eq!(r1, b"ACGTACGTACGTACGT");
        assert_eq!(r2, b"CCCCCCCCCCCCCCCC");
    }

    #[test]
    fn reverse_alignment_is_normalized_to_original_sequence() {
        // BAM reverse-strand SEQ is the reverse complement of the original
        // sequenced string. TGCA -> original TGCA here is palindromic-ish, so
        // use an asymmetric example.
        let record = record(b"frag", b"CGTT", 0x10);
        assert_eq!(original_sequence(&record), b"AACG");
    }

    #[test]
    fn missing_mate_is_an_error_for_unbarcoded_pair_identity() {
        let group = ReadGroup::new(record(b"frag1", b"ACGTACGTACGTACGT", 0x1 | 0x40));
        let err = group.primary_pair_sequences().unwrap_err().to_string();
        assert!(err.contains("no primary R2"));
    }
    #[test]
    fn bam_pair_identity_matches_direct_fastq_identity_for_none_grammar() {
        let grammar = sc_primer::Grammar::parse("none", "NONE").unwrap();
        let r1_original = b"AAAACCCCGGGGTTTT";
        let r2_original = b"TTTTGGGGCCCCAAAA";

        let mut group = ReadGroup::new(record(b"frag1", r1_original, 0x1 | 0x40));
        group
            .push(record(b"frag1", r2_original, 0x1 | 0x80))
            .unwrap();

        let (bam_r1, bam_r2) = group.primary_pair_sequences().unwrap();
        let from_bam = grammar
            .molecule_identity(None, None, &bam_r1, &bam_r2)
            .unwrap();
        let from_fastq = grammar
            .molecule_identity(None, None, r1_original, r2_original)
            .unwrap();

        assert_eq!(from_bam, from_fastq);
    }

    #[test]
    fn different_qnames_with_same_pair_sequence_have_same_none_identity() {
        let grammar = sc_primer::Grammar::parse("none", "NONE").unwrap();
        let r1 = b"AAAACCCCGGGGTTTT";
        let r2 = b"TTTTGGGGCCCCAAAA";

        let mut first = ReadGroup::new(record(b"frag1", r1, 0x1 | 0x40));
        first.push(record(b"frag1", r2, 0x1 | 0x80)).unwrap();
        let mut second = ReadGroup::new(record(b"frag2", r1, 0x1 | 0x40));
        second.push(record(b"frag2", r2, 0x1 | 0x80)).unwrap();

        let (a1, a2) = first.primary_pair_sequences().unwrap();
        let (b1, b2) = second.primary_pair_sequences().unwrap();

        assert_eq!(
            grammar.molecule_identity(None, None, &a1, &a2).unwrap(),
            grammar.molecule_identity(None, None, &b1, &b2).unwrap(),
        );
    }
}
