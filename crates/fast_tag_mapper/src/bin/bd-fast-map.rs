use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::PathBuf,
};

use clap::Parser;
use fast_tag_mapper::FastMapperCli;
use mapping_info::MappingInfo;

#[derive(Debug, Parser)]
#[command(author, version, about = "Map sequences to fast-tag feature ids")]
struct Cli {
    #[command(flatten)]
    mapper: FastMapperCli,

    /// Input file with one sequence per line. Reads stdin if omitted.
    #[arg(short, long)]
    input: Option<PathBuf>,
}

fn main() {
    let cli = Cli::parse();

    let mapper = match cli.mapper.mapper() {
        Ok(mapper) => mapper,
        Err(err) => {
            panic!("Failed to build FastTagMapper: {err}");
        }
    };

    let mut info = MappingInfo::new(None, 0.0, 0);

    match cli.input {
        Some(path) => {
            let reader = match File::open(&path) {
                Ok(file) => BufReader::new(file),
                Err(err) => {
                    panic!(
                        "Failed to open input file '{}': {err}",
                        path.display()
                    );
                }
            };

            if let Err(err) = run(reader, &mapper, &mut info) {
                panic!("Failed while processing '{}': {err}", path.display());
            }
        }
        None => {
            let stdin = io::stdin();

            if let Err(err) = run(stdin.lock(), &mapper, &mut info) {
                panic!("Failed while reading stdin: {err}");
            }
        }
    }

    eprintln!("{info}");
}

fn run<R: BufRead>(
    reader: R,
    mapper: &fast_tag_mapper::FastTagMapper,
    info: &mut MappingInfo,
) -> std::io::Result<()> {
    for line in reader.lines() {
        let seq = line?;

        if seq.starts_with('@') || seq.starts_with('+') || seq.is_empty() {
            continue;
        }

        match mapper.map_feature_id(seq.as_bytes(), info) {
            Some(feature_id) => println!("{feature_id}"),
            None => println!("NA"),
        }
    }

    Ok(())
}