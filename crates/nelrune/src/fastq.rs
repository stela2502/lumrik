use anyhow::{anyhow, bail, Result};

use bam_tide::fastq::{
    FastqPairReader,
    FastqRecord,
};

use int_to_str::IntToStr;
use sc_primer::PrimerDetector;

pub struct ParsedPair {
    pub read_id: String,

    pub cell_seq: Vec<u8>,
    pub umi_seq: Vec<u8>,

    pub cell_id: u64,
    pub umi_id: u64,

    pub r2: FastqRecord,
}

pub fn next_parsed_pair(
    reader: &mut FastqPairReader,
    primer: &PrimerDetector,
) -> Result<Option<ParsedPair>> {
    loop {
        let Some((r1, r2)) =
            reader.next_pair()?
        else {
            return Ok(None);
        };

        let r1_id =
            r1.clean_id();

        let r2_id =
            r2.clean_id();

        if r1_id != r2_id {
            bail!(
                "FASTQ pair ID mismatch: R1='{r1_id}', R2='{r2_id}'"
            );
        }

        let Some(hit) =
            primer
                .detect_first(
                    &r1.seq,
                    &r1.qual,
                )
                .map_err(|e| {
                    anyhow!(
                        "primer detection failed for read '{}': {e}",
                        r1_id
                    )
                })?
        else {
            continue;
        };

        let cell =
            hit.get_cell(
                &r1.seq,
                &r1.qual,
            )
            .map_err(|e| {
                anyhow!(
                    "failed to extract cell barcode from '{}': {e}",
                    r1_id
                )
            })?;

        let umi =
            hit.get_umi(
                &r1.seq,
                &r1.qual,
            )
            .map_err(|e| {
                anyhow!(
                    "failed to extract UMI from '{}': {e}",
                    r1_id
                )
            })?;

        let cell_id =
            dna_to_u64(
                &cell.seq,
                "cell barcode",
            )?;

        let umi_id =
            dna_to_u64(
                &umi.seq,
                "UMI",
            )?;

        return Ok(
            Some(
                ParsedPair {
                    read_id: r2_id,

                    cell_seq:
                        cell.seq,

                    umi_seq:
                        umi.seq,

                    cell_id,
                    umi_id,

                    r2,
                }
            )
        );
    }
}

fn dna_to_u64(
    seq: &[u8],
    label: &str,
) -> Result<u64> {
    IntToStr::try_new(seq)
        .map(|encoded| encoded.into_u64())
        .map_err(|e| {
            anyhow!(
                "failed to encode {label}: {e}"
            )
        })
}