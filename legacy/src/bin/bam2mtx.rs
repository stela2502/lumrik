
//use rustody::cellids::CellIds;
//use rustody::cellids10x::CellIds10x;
//use rustody::traits::CellIndex;
use mapping_info::MappingInfo;

use bam_tide::gtf::GTF;
use bam_tide::mutation_processor::MutationProcessor;
use bam_tide::gtf_logics::{process_data, PROGRAM_NAME, AnalysisType, MatchType, GtfType };

//use rustody::ofiles::{Ofiles, Fspot};

use std::path::PathBuf;
use std::fs;
use std::fs::File;

//use std::thread;
//use rayon::prelude::*;

use std::time::SystemTime;

use clap::{Parser};


/// This tool can quantify a bam file from both single cell data as well as bulk data.
/// The default is to quantify single cell data from the BD Rhapsody tool.
/// with these default settings: --cell_tag CB --umi_tag UB --analysis_type single-cell --match_type exact --gtf_type genes.
/// Switching to --analysis_type "bulk" will create the mtx out files for one single cell id: "1"
/// Switching to --match_type "overlap" will quantify reads in a sticky way - any overlap will be called a match - even paired ones!!
/// Switching to --gtf_type "exons" will quantify exons instead of genes - make sure our exons are uniquely named!
#[derive(Parser)]
#[clap(version = "0.4.3", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the bam file to quantify
    #[clap(short, long)]
    bam: String,

    /// the gtf file fitting to the Bam file (text or gzipped)
    #[clap(short, long)]
    gtf: String,

    /// the outpath
    #[clap(short, long)]
    outpath: String,

    /// the minimum (UMI) reads per cell (sample + genes + antibody combined)
    #[clap(short, long)]
    min_umi: usize,

    /// used processor cores (default all)
    #[clap(short, long)]
    num_proc: Option<usize>,

    /// tag name for the CELL information (default CB for velocity default - change to CR for CellRanger)
    #[clap(short, long)]
    cell_tag:Option<String>,

    /// tag name for the UMI information (default UB for velocity default - change to UR for CellRanger)
    #[clap(short, long)]
    umi_tag:Option<String>,

    /// For mutation collection please give me a quality cutoff for accepting a nucl as a valid mutation (20? 30?).
    #[clap(short, long)]
    qual:Option<usize>,

    /// Collect single cell info or bulk
    #[clap(short, long, value_enum, default_value = "single-cell")]
    analysis_type: AnalysisType,

    /// Match only inside exons or overlapping?
    #[clap( long, value_enum, default_value = "exact")]
    match_type: MatchType,

    /// Group gtf info into genes or e.g. singel cell quantifications (genes) or use each exon on it's own for e.g. TE analyses (exon)
    #[clap( long, value_enum, default_value = "genes")]
    gtf_type: GtfType,

}


// Main function
fn main() {
    let now = SystemTime::now();
    let opts: Opts = Opts::parse();

    // Create output directory if needed
    if fs::metadata(&opts.outpath).is_err() {
        if let Err(err) = fs::create_dir_all(&opts.outpath) {
            eprintln!("Error creating directory {}: {}", &opts.outpath, err);
        } else {
            println!("New output directory created successfully!");
        }
    }

    // Set up logging
    let log_file_str = PathBuf::from(&opts.outpath).join("Mapping_log.txt");
    let log_file = File::create(log_file_str).expect("Failed to create log file");
    let umi_tag: [u8; 2] = opts.umi_tag.unwrap_or_else(|| "UB".to_string()).into_bytes().try_into().expect("umi-tag must be exactly 2 chars long");
    let cell_tag: [u8; 2] = opts.cell_tag.unwrap_or_else(|| "CB".to_string()).into_bytes().try_into().expect("umi-tag must be exactly 2 chars long");

    let num_threads = opts.num_proc.unwrap_or_else(rayon::current_num_threads);

    let mut mapping_info = MappingInfo::new(Some(log_file), 3.0, 0);
    mapping_info.start_counter();

    // Parse BAM and GTF
    println!("reading GTF file");
    

    let mut gtf = GTF::new();
    match &opts.gtf_type {
        GtfType::Genes => gtf.parse_gtf(&opts.gtf).unwrap(),
        GtfType::Exons => gtf.parse_gtf_only_exons(&opts.gtf, "exon" ).unwrap(),
    };

    /*let mutations: Option<MutationProcessor> = match opts.qual {
        Some(quality) => {
            Some(MutationProcessor { quality_cutoff:quality })
        },
        None => None
    };*/
    let mutations: Option<MutationProcessor> = None;

    // Process data
    let ( mut expr_results, mut mut_results ) =  match process_data(
        &opts.bam,
        &mut mapping_info,
        &gtf,
        cell_tag,
        umi_tag,
        num_threads,
        &mutations,
        &opts.analysis_type,
        &opts.match_type,
    ){
        Ok(ret) => ret,
        Err(e) => {
            panic!("{e:?}");
        }
    };

    // Final reporting and cleanup

    let file_path_sp = PathBuf::from(&opts.outpath).join( &*PROGRAM_NAME );
    println!("Writing data to path {:?}", file_path_sp);

    let info = expr_results.0.write_sparse_sub(file_path_sp, &expr_results.1, &expr_results.1.get_all_gene_names(), opts.min_umi).unwrap();
    mapping_info.write_to_log( info );

    if mut_results.0.len() > 0 {
        let file_path_sp_mut = PathBuf::from(&opts.outpath).join( &*PROGRAM_NAME ).join("_mutations");
        println!("Writing mutation data to path {:?}", file_path_sp_mut);
        let i2 = mut_results.0.write_sparse_sub(file_path_sp_mut, &mut_results.1, &mut_results.1.get_all_gene_names(), opts.min_umi).unwrap();
        mapping_info.write_to_log( i2 );
    }

    
    mapping_info.log_report( );

    println!("The total issues report:\n{}", mapping_info.report_to_string());
    println!("Runtime assessment:\n{}", mapping_info.program_states_string());

    match now.elapsed() {
        Ok(elapsed) => {
            let mut milli = elapsed.as_millis();

            let mil = milli % 1000;
            milli= (milli - mil) /1000;

            let sec = milli % 60;
            milli= (milli -sec) /60;

            let min = milli % 60;
            milli= (milli -min) /60;

            mapping_info.write_to_log( format!("\nfinished in {milli}h {min}min {sec} sec {mil}milli sec\n") );
            println!("\nI have analyzed the bam file in {milli}h {min}min {sec} sec {mil}milli sec\n");
        },
        Err(e) => {println!("Error: {e:?}");}
    }

    mapping_info.report_to_string();
}

