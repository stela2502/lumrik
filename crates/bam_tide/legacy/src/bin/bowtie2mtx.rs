
//use rustody::cellids::CellIds;
//use rustody::cellids10x::CellIds10x;
//use rustody::traits::CellIndex;
use mapping_info::MappingInfo;

use bam_tide::bed_data::BedData;
use bam_tide::mutation_processor::MutationProcessor;
use bam_tide::gtf_logics::{process_data_bowtie2, PROGRAM_NAME, AnalysisType, MatchType };
//use rust_htslib::bam::Reader;
//use rust_htslib::bam::Read;

//use rustody::ofiles::{Ofiles, Fspot};

use std::path::PathBuf;
use std::fs;
use std::fs::File;
use std::io::Write;

//use std::thread;
//use rayon::prelude::*;

use std::time::SystemTime;

use clap::{Parser};


/// Under devlopment: This tool is used to collect the coverage and mutation loads from a BAM file.
/// At the moment the results from this program are buggy!
#[derive(Parser)]
#[clap(version = "0.4.3", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the bam file to quantify
    #[clap(short, long)]
    bam: String,

    /// the fasta file needed for mutation anaylsis
    #[clap(short, long)]
    fasta: Option<String>,

    
    /// the bin width for the genome coverage info
    #[clap( long, default_value_t = 1000 )]
    bin_width: usize,

    /// the outpath
    #[clap(short, long)]
    outpath: String,

    /// the minimum (UMI) reads per cell (sample + genes + antibody combined)
    #[clap(short, long)]
    min_umi: usize,

    /// used processor cores (default all)
    #[clap(short, long)]
    num_proc: Option<usize>,

    /// For mutation collection please give me a quality cutoff for accepting a nucl as a valid mutation (20? 30?).
    #[clap(short, long)]
    qual:Option<usize>,

    /// Collect single cell info or bulk
    #[clap(short, long, value_enum, default_value = "single-cell")]
    analysis_type: AnalysisType,

    /// Match only inside exons or overlapping?
    #[clap( long, value_enum, default_value = "exact")]
    match_type: MatchType,

    /// Ignore mutations in the first and last mut_clip[%] area of the read
    #[clap( long, default_value_t = 0.0)]
    mut_clip: f32,

    /// Collect only from these chromosomes
    #[clap(long, value_parser, num_args = 1.., )]
    chromosomes: Option<Vec<String>>,


    /*
    /// Group gtf info into genes or e.g. singel cell quantifications (genes) or use each exon on it's own for e.g. TE analyses (exon)
    #[clap( long, value_enum, default_value = "genes")]
    gtf_type: GtfType,
    */
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

    let num_threads = opts.num_proc.unwrap_or_else(rayon::current_num_threads);

    let mut mapping_info = MappingInfo::new(Some(log_file), 3.0, 0);
    mapping_info.start_counter();

    // Parse BAM and GTF
    println!("creating Bed coverage info");
    
    let bed = BedData::init( &opts.bam, opts.bin_width, num_threads, opts.chromosomes );

    println!("Created {} bed areas:_\n{bed}", bed.coverage_data.len() );

    let mut mutations: Option<MutationProcessor> = match opts.fasta {
        Some(fasta_path) => {
            match MutationProcessor::new(&opts.bam, &fasta_path, opts.mut_clip) {
                Ok(processor) => Some(processor),
                Err(e) => {
                    eprintln!("Failed to initialize MutationProcessor: {}", e);
                    None
                }
            }
        },
        None => None
    };

    // Process data
    let ( mut expr_results, mut mut_results ) =  match process_data_bowtie2(
        &opts.bam,
        &mut mapping_info,
        &bed,
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
    if let Some(mproc) = mutations.as_mut() {
        mproc.hist = mapping_info.hist.clone();
    }

    let file_path_sp = PathBuf::from(&opts.outpath).join( &*PROGRAM_NAME );
    println!("Writing data to path {:?}", file_path_sp);

    let info = expr_results.0.write_sparse_sub(file_path_sp, &expr_results.1, &expr_results.1.get_all_gene_names(), opts.min_umi).unwrap();
    mapping_info.write_to_log( info );

    if mut_results.0.len() > 0 {
        let file_path_sp_mut = PathBuf::from(&opts.outpath).join( &*PROGRAM_NAME ).join("mutations");
        println!("Writing mutation data to path {:?}", file_path_sp_mut);
        let i2 = mut_results.0.write_sparse_sub(file_path_sp_mut, &mut_results.1, &mut_results.1.get_all_gene_names(), opts.min_umi).unwrap();
        mapping_info.write_to_log( i2 );
    }

    mapping_info.stop_file_io_time();
    
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

    match mutations {
        Some( mutation ) => {
            println!("{}", mutation );
            let path = PathBuf::from(&opts.outpath).join("mutations_report.txt");
            let mut file = File::create(path).expect("Could not create the mutations report file");
            writeln!(file, "{}", mutation).expect("Could not write the mutations report");
        },
        None => {
            println!("Mutaion analysis needs a matching fasta formated reference genome");
        }
    }
    mapping_info.report_to_string();
}
