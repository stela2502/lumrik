use anyhow::{bail, Context, Result};
use lumrik_clap::{Parser, ValueEnum, ReadableFile, OutputFile};
use flate2::read::MultiGzDecoder;
use sc_primer::{BdCellVersion, RhapsodyWhitelist};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Chemistry {
    BdV1,
    BdV2_96,
    BdV2_384,
    BdV2_384vdj,
}

impl Chemistry {
    fn version(self) -> BdCellVersion {
        match self {
            Self::BdV1 => BdCellVersion::V1,
            Self::BdV2_96 => BdCellVersion::V2_96,
            Self::BdV2_384 => BdCellVersion::V2_384,
            Self::BdV2_384vdj => BdCellVersion::V2_384Vdj,
        }
    }
}

#[derive(Parser, Debug)]
#[command(version = env!("CARGO_PKG_VERSION"))]
#[command(about = "Convert a nucleotide barcode column in a TSV to numeric BD/Rustody cell IDs")]
struct Args {
    #[arg(long)]
    input: ReadableFile,
    #[arg(long)]
    output: OutputFile,
    #[arg(long, value_enum, default_value_t = Chemistry::BdV2_384)]
    chemistry: Chemistry,
    #[arg(long, default_value = "cell")]
    column: String,
    #[arg(long, default_value = "rustody_cell_id")]
    output_column: String,
    /// Allow replacing nonempty existing IDs; conflicts are counted and reported.
    #[arg(long)]
    overwrite: bool,
}

#[derive(Debug)]
struct ConvertOptions {
    chemistry: Chemistry,
    column: String,
    output_column: String,
    overwrite: bool,
}

impl From<&Args> for ConvertOptions {
    fn from(args: &Args) -> Self {
        Self {
            chemistry: args.chemistry,
            column: args.column.clone(),
            output_column: args.output_column.clone(),
            overwrite: args.overwrite,
        }
    }
}

fn convert<R: BufRead, W: Write>(mut reader: R, mut writer: W, options: &ConvertOptions) -> Result<()> {
    let whitelist = RhapsodyWhitelist::builtin(options.chemistry.version());
    let mut line = String::new();
    if reader.read_line(&mut line)? == 0 {
        bail!("empty input TSV");
    }

    let header: Vec<String> = line
        .trim_end_matches(['\r', '\n'])
        .split('\t')
        .map(str::to_owned)
        .collect();

    let columns = &header;
    let source = columns
        .iter()
        .position(|c| *c == options.column)
        .with_context(|| format!("missing barcode column {:?}", options.column))?;
    let target = columns.iter().position(|c| *c == options.output_column);
    writeln!(
        writer,
        "{}{}",
        header.join("\t"),
        if target.is_none() {
            format!("\t{}", options.output_column)
        } else {
            String::new()
        }
    )?;
    let (mut rows, mut converted, mut failed, mut preserved, mut conflicts) =
        (0u64, 0u64, 0u64, 0u64, 0u64);
    line.clear();
    while reader.read_line(&mut line)? != 0 {
        let raw = line.trim_end_matches(['\r', '\n']);
        let mut fields: Vec<&str> = raw.split('\t').collect();
        if fields.len() != columns.len() {
            bail!(
                "row {} has {} columns; expected {}",
                rows + 2,
                fields.len(),
                columns.len()
            );
        }
        rows += 1;
        let barcode = fields[source].trim().to_ascii_uppercase();
        let numeric = whitelist.cell_id_for_seq(barcode.as_bytes());
        let existing = target.map(|i| fields[i].trim()).unwrap_or("");
        let has_existing = !existing.is_empty() && !existing.eq_ignore_ascii_case("NA");
        let calculated = numeric.map(|n| n.to_string());
        if numeric.is_some() {
            converted += 1;
        } else {
            failed += 1;
        }
        if has_existing && calculated.as_deref() != Some(existing) {
            conflicts += 1;
        }
        let value = if has_existing && !options.overwrite {
            preserved += 1;
            existing.to_owned()
        } else {
            calculated.unwrap_or_default()
        };
        if let Some(i) = target {
            fields[i] = &value;
        }
        writeln!(
            writer,
            "{}{}",
            fields.join("\t"),
            if target.is_none() {
                format!("\t{value}")
            } else {
                String::new()
            }
        )?;
        line.clear();
    }
    writer.flush()?;
    eprintln!("Input rows: {rows}\nConverted barcodes: {converted}\nUnconvertible barcodes: {failed}\nPreserved existing IDs: {preserved}\nExisting ID conflicts: {conflicts}\nOutput rows: {rows}");
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    if args.input.as_path() == args.output.as_path() {
        bail!("input and output paths must differ");
    }
    args.output.prepare().with_context(|| format!("preparing --output {}", args.output.display()))?;
    let input =
        File::open(&args.input).with_context(|| format!("opening {}", args.input.display()))?;
    let reader: Box<dyn Read> = if args.input.extension().is_some_and(|e| e == "gz") {
        Box::new(MultiGzDecoder::new(input))
    } else {
        Box::new(input)
    };
    let output = File::create(&args.output)
        .with_context(|| format!("creating {}", args.output.display()))?;
    let options = ConvertOptions::from(&args);
    convert(BufReader::new(reader), BufWriter::new(output), &options)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn converts_and_preserves_rows() {
        let wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);
        let barcode = String::from_utf8(wl.cell_id_to_seq(1234).unwrap()).unwrap();
        let input = format!("cell\trustody_cell_id\tchain\n{barcode}\tNA\tIGH\nBAD\t\tIGK\n");
        let options = ConvertOptions {
            chemistry: Chemistry::BdV2_384,
            column: "cell".into(),
            output_column: "rustody_cell_id".into(),
            overwrite: false,
        };
        let mut out = Vec::new();
        convert(input.as_bytes(), &mut out, &options).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains(&format!("{barcode}\t1234\tIGH")));
        assert!(text.contains("BAD\t\tIGK"));
        assert_eq!(text.lines().count(), 3);
    }
}
