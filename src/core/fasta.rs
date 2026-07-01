use anyhow::{anyhow, Result, Context};
use std::io::Write;

use gtf_splice_index::{Transcript, Strand};
use snp_index::Genome;


#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    pub id: String,
    pub seq: Vec<u8>,
}

impl FastaRecord {
    pub fn new(id: impl Into<String>, seq: Vec<u8>) -> Self {
        Self {
            id: id.into(),
            seq,
        }
    }

    pub fn revcomp(mut self) -> Self {
        self.seq.reverse();

        for base in self.seq.iter_mut() {
            *base = complement(*base);
        }

        self
    }

    pub fn from_transcript(
        genome: &Genome,
        tx: &Transcript,
        chr_name: &str,
    ) -> Result<Self> {
        let tx_name = tx
            .primary_name()
            .ok_or_else(|| anyhow!("transcript has no primary name"))?;

        let genome_chr_id = genome
            .chr_id(chr_name)
            .with_context(|| format!("chromosome missing from genome FASTA: {chr_name}"))?;

        let mut seq = Vec::<u8>::new();

        for exon in tx.exons() {
            let slice = genome
                .slice(genome_chr_id, exon.start, exon.end)
                .with_context(|| {
                    format!(
                        "failed to fetch {}:{}-{} for transcript {}",
                        chr_name,
                        exon.start + 1,
                        exon.end,
                        tx_name
                    )
                })?;

            seq.extend_from_slice(slice);
        }

        let record = Self::new(tx_name, seq);

        Ok(if matches!(tx.strand, Strand::Minus) {
            record.revcomp()
        } else {
            record
        })
    }

    pub fn write<W: Write>(&self, writer: &mut W, line_width: usize) -> Result<()> {
        if line_width == 0 {
            return Err(anyhow!("line_width must be > 0"));
        }

        writeln!(writer, ">{}", self.id)?;

        for chunk in self.seq.chunks(line_width) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }

        Ok(())
    }
}

#[inline(always)]
pub fn complement(base: u8) -> u8 {
    DNA_COMPLEMENT[base as usize]
}

const DNA_COMPLEMENT: [u8; 256] = {
    let mut table = [0u8; 256];

    let mut i = 0;
    while i < 256 {
        table[i] = i as u8;
        i += 1;
    }

    table[b'A' as usize] = b'T';
    table[b'a' as usize] = b't';
    table[b'T' as usize] = b'A';
    table[b't' as usize] = b'a';
    table[b'G' as usize] = b'C';
    table[b'g' as usize] = b'c';
    table[b'C' as usize] = b'G';
    table[b'c' as usize] = b'g';
    table[b'N' as usize] = b'N';
    table[b'n' as usize] = b'n';

    table
};

#[test]
fn fasta_record_from_transcript_builds_strand_oriented_mrna() {
    use crate::core::fasta::FastaRecord;
    use snp_index::Genome;
    use gtf_splice_index::types::{RefBlock, Strand};
    use gtf_splice_index::Transcript;

    let genome = Genome::new(vec![(
        "chr1".to_string(),
        b"AAAACCGTGGCAtttt".to_vec(),
    )])
    .unwrap();

    let mut plus = Transcript::new(0, 0, "tx_plus", 0, Strand::Plus);
    plus.add_exon(RefBlock::new(4, 8));   // CCGT
    plus.add_exon(RefBlock::new(8, 12));  // GGCA
    plus.finalize();

    let plus_record = FastaRecord::from_transcript(&genome, &plus, "chr1").unwrap();

    assert_eq!(plus_record.id, "tx_plus");
    assert_eq!(plus_record.seq, b"CCGTGGCA");

    let mut minus = Transcript::new(1, 0, "tx_minus", 0, Strand::Minus);
    minus.add_exon(RefBlock::new(4, 8));   // CCGT
    minus.add_exon(RefBlock::new(8, 12));  // GGCA
    minus.finalize();

    let minus_record = FastaRecord::from_transcript(&genome, &minus, "chr1").unwrap();

    // genomic concatenation: CCGTGGCA
    // reverse-complement: TGCCACGG
    assert_eq!(minus_record.id, "tx_minus");
    assert_eq!(minus_record.seq, b"TGCCACGG");
}
