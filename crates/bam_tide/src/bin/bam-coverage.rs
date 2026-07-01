// src/bin/bam_coverage.rs
use bam_tide::bed_data::BedData;
use bam_tide::cli::CoverageCli;
use clap::Parser;

fn main() {
    let opts = CoverageCli::parse();

    let bed = match BedData::from_bam_with_policy(&opts) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("Error: {e}");
            std::process::exit(1);
        }
    };

    if let Err(e) = write_output(&bed, &opts.outfile) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

pub fn write_output(bed: &BedData, outfile: &str) -> Result<(), String> {
    if outfile.ends_with(".bw") || outfile.ends_with(".bigwig") {
        bed.write_bigwig(outfile)
    } else {
        bed.write_bedgraph(outfile).map_err(|e| e.to_string())
    }
}
