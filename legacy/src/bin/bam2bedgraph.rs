use bam_tide::bed_data::{ BedData, Normalize };
use clap::Parser;
use bam_tide::gtf_logics::{ AnalysisType };



/// A command-line tool for converting BAM coverage data to BedGraph format.
#[derive(Parser, Debug)]
#[clap(version = "0.0.1", author = "Stefan L. <stefan.lang@med.lu.se>")]
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

    /* Taken from the bamCoverage help:
      --normalizeUsing {RPKM,CPM,BPM,RPGC,None}
    Use one of the entered methods to normalize the number of reads per bin. By default, no normalization is performed. RPKM = Reads Per Kilobase per Million mapped reads; CPM = Counts Per Million mapped
    reads, same as CPM in RNA-seq; BPM = Bins Per Million mapped reads, same as TPM in RNA-seq; RPGC = reads per genomic content (1x normalization); Mapped reads are considered after blacklist filtering (if
    applied). RPKM (per bin) = number of reads per bin / (number of mapped reads (in millions) * bin length (kb)). CPM (per bin) = number of reads per bin / number of mapped reads (in millions). BPM (per bin)
    = number of reads per bin / sum of all reads per bin (in millions). RPGC (per bin) = number of reads per bin / scaling factor for 1x average coverage. None = the default and equivalent to not setting this
    option at all. This scaling factor, in turn, is determined from the sequencing depth: (total number of mapped reads * fragment length) / effective genome size. The scaling factor used is the inverse of
    the sequencing depth computed for the sample to match the 1x coverage. This option requires --effectiveGenomeSize. Each read is considered independently, if you want to only count one mate from a pair in
    paired-end data, then use the --samFlagInclude/--samFlagExclude options. (Default: None)
    */
    /// Normalize the data somehow
    #[clap(short, long, value_enum, default_value = "not")]
    normalize: Normalize,

    /// Bin width for coverage calculation (default: 50bp).
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

    match bed_data.write_bedgraph(&args.outfile, ){
        Ok(_) => {
            println!("Conversion completed successfully.");
        },
        Err(e) => {
            println!("Some error occured - check if that is a problem: {e:?}")
        }
    }

    
    Ok(())
}
