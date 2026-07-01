use bam_tide::bed_data::{ BedData, Normalize };
use clap::Parser;
use bam_tide::gtf_logics::{ AnalysisType };

/// A command-line tool for converting BAM coverage data to BigWig format.
#[derive(Parser, Debug)]
#[clap(version = "0.0.3", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Args {
    /// Input BAM file (sorted by chromosome position).
    #[arg(short, long)]
    bam: String,

    /// Output BigWig file.
    #[arg(short, long)]
    outfile: String,

    /// tag name for the CELL information (default CB for velocity default - change to CR for CellRanger)
    #[clap(short, long)]
    cell_tag:Option<String>,

    /// tag name for the UMI information (default UB for velocity default - change to UR for CellRanger)
    #[clap(short, long)]
    umi_tag:Option<String>,

    /// Collect single cell info or bulk
    #[clap(short, long, value_enum, default_value = "bulk")]
    analysis_type: AnalysisType,

    /// Normalize the data somehow
    #[clap(short, long, value_enum, default_value = "not")]
    normalize: Normalize,

    /// Bin width for coverage calculation.
    #[arg(short, long, default_value_t = 50)]
    width: usize,

    /// Collect only R1 areas
    #[arg( long )]
    only_r1: bool,

    /// Minimum mapping quality to include a read
    #[arg(long, default_value_t = 0)]
    min_mapping_quality: u8,

}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    println!("BamTide CLI");
    println!("Processing BAM: {}", args.bam);
    println!("Output BigWig: {}", args.outfile);
    println!("Bin Width: {} bp", args.width);

    let umi_tag: [u8; 2] = args.umi_tag.unwrap_or_else(|| "UB".to_string()).into_bytes().try_into().expect("umi-tag must be exactly 2 chars long");
    let cell_tag: [u8; 2] = args.cell_tag.unwrap_or_else(|| "CB".to_string()).into_bytes().try_into().expect("umi-tag must be exactly 2 chars long");

    let add_introns = false;
    let mut bed_data = BedData::new( &args.bam, args.width, 2, 
        &args.analysis_type, &cell_tag, &umi_tag, add_introns, 
        args.only_r1, args.min_mapping_quality );

    bed_data.normalize( &args.normalize );

    bed_data.write_bigwig( &args.outfile )?;

    println!("Conversion completed successfully.");
    Ok(())
}

