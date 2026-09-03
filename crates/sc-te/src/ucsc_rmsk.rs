use std::io::{BufRead, Write};

use anyhow::{bail, Context, Result};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RmskConversionSummary {
    pub input_rows: u64,
    pub output_rows: u64,
}

/// Convert a UCSC `rmsk.txt` table into a TE-oriented GTF.
///
/// UCSC `rmsk` rows contain 17 tab-separated columns. The fields used here are:
/// chromosome (6), zero-based half-open start/end (7/8), strand (10), repeat
/// name (11), class (12), family (13). Each RepeatMasker row is emitted as one
/// GTF exon. The repeat name is the GTF gene identity and each input row gets a
/// deterministic unique transcript identity.
pub fn convert_ucsc_rmsk_to_gtf<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
) -> Result<RmskConversionSummary> {
    let mut summary = RmskConversionSummary::default();

    for (line_idx, line) in reader.lines().enumerate() {
        let line_no = line_idx + 1;
        let line = line.with_context(|| format!("failed reading UCSC rmsk line {line_no}"))?;
        if line.trim().is_empty() {
            continue;
        }

        summary.input_rows += 1;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 17 {
            bail!(
                "UCSC rmsk line {line_no}: expected at least 17 tab-separated columns, found {}",
                fields.len()
            );
        }

        let chrom = fields[5];
        let start0: u64 = fields[6]
            .parse()
            .with_context(|| format!("UCSC rmsk line {line_no}: invalid genoStart {:?}", fields[6]))?;
        let end0: u64 = fields[7]
            .parse()
            .with_context(|| format!("UCSC rmsk line {line_no}: invalid genoEnd {:?}", fields[7]))?;
        if end0 <= start0 {
            bail!(
                "UCSC rmsk line {line_no}: invalid half-open interval {chrom}:{start0}-{end0}"
            );
        }

        let strand = match fields[9] {
            "+" => "+",
            "-" => "-",
            other => bail!("UCSC rmsk line {line_no}: invalid strand {other:?}"),
        };
        let rep_name = fields[10];
        let rep_class = fields[11];
        let rep_family = fields[12];

        ensure_gtf_attribute_value(rep_name, "repName", line_no)?;
        ensure_gtf_attribute_value(rep_class, "repClass", line_no)?;
        ensure_gtf_attribute_value(rep_family, "repFamily", line_no)?;

        // UCSC coordinates are 0-based half-open; GTF is 1-based inclusive.
        let gtf_start = start0 + 1;
        let gtf_end = end0;
        let transcript_id = format!("{rep_name}_ucsc_rmsk_{line_no}");

        writeln!(
            writer,
            "{chrom}\tRepeatMasker\texon\t{gtf_start}\t{gtf_end}\t.\t{strand}\t.\tgene_id \"{rep_name}\"; transcript_id \"{transcript_id}\"; gene_name \"{rep_name}\"; family_id \"{rep_family}\"; class_id \"{rep_class}\";"
        )
        .with_context(|| format!("failed writing GTF row for UCSC rmsk line {line_no}"))?;

        summary.output_rows += 1;
    }

    Ok(summary)
}

fn ensure_gtf_attribute_value(value: &str, field: &str, line_no: usize) -> Result<()> {
    if value.contains(['\t', '\n', '\r', '"']) {
        bail!(
            "UCSC rmsk line {line_no}: {field} contains a character unsafe for a GTF attribute: {value:?}"
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn converts_ucsc_coordinates_and_te_metadata() {
        let input = b"585\t4637\t195\t23\t0\tchr1\t3000000\t3000123\t-190000000\t-\tL1Md_A\tLINE\tL1\t1\t123\t0\t42\n";
        let mut output = Vec::new();

        let summary = convert_ucsc_rmsk_to_gtf(Cursor::new(input), &mut output).unwrap();
        let output = String::from_utf8(output).unwrap();

        assert_eq!(summary.input_rows, 1);
        assert_eq!(summary.output_rows, 1);
        assert_eq!(
            output,
            "chr1\tRepeatMasker\texon\t3000001\t3000123\t.\t-\t.\tgene_id \"L1Md_A\"; transcript_id \"L1Md_A_ucsc_rmsk_1\"; gene_name \"L1Md_A\"; family_id \"L1\"; class_id \"LINE\";\n"
        );
    }

    #[test]
    fn rejects_wrong_table_shape() {
        let mut output = Vec::new();
        let err = convert_ucsc_rmsk_to_gtf(Cursor::new(b"too\tfew\tcolumns\n"), &mut output)
            .unwrap_err();
        assert!(err.to_string().contains("expected at least 17"));
    }

    #[test]
    fn rejects_invalid_interval() {
        let input = b"0\t0\t0\t0\t0\tchr1\t100\t100\t0\t+\tB1\tSINE\tB4\t0\t0\t0\t1\n";
        let mut output = Vec::new();
        let err = convert_ucsc_rmsk_to_gtf(Cursor::new(input), &mut output).unwrap_err();
        assert!(err.to_string().contains("invalid half-open interval"));
    }
}
