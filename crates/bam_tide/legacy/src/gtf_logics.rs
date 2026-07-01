
use rayon::slice::ParallelSlice;
use rayon::iter::ParallelIterator;

use indicatif::MultiProgress;
use indicatif::ProgressStyle;
use indicatif::ProgressBar;

use crate::gtf::{ ExonIterator };
//use crate::gtf::exon_iterator::ReadResult;
use crate::feature_matcher::FeatureMatcher;
use crate::mutation_processor::MutationProcessor;
use crate::read_data::ReadData;

use scdata::{Scdata, IndexedGenes };
use scdata::scdata::MatrixValueType;

use mapping_info::MappingInfo;
//use rustody::int_to_str::IntToStr;
//use rustody::analysis::bam_flag::BamFlag;
//use rustody::genes_mapper::cigar::CigarEnum;

use std::path::Path;
use std::collections::HashMap;

use rust_htslib::bam::{Reader,Read};
use rust_htslib::bam::Header;

//use rust_htslib::header::record::LinearMap;



const BUFFER_SIZE: u64 = 1_000_000;
const ATAC_BUFFER_SIZE: u64 = 100_000;


use std::env;
use lazy_static::lazy_static;

use clap::ValueEnum;

use std::str::FromStr;

lazy_static! {
    pub static  ref PROGRAM_NAME: String = {
        if let Some(program_path) = env::args().next() {
            if let Some(program_name) = Path::new(&program_path).file_name() {
                program_name.to_string_lossy().to_string()
            } else {
                String::from("Unknown")
            }
        } else {
            String::from("Unknown")
        }
    };
}

#[derive(ValueEnum, Clone, Debug)]
pub enum GtfType {
    Genes,
    Exons,
}

impl FromStr for GtfType {
    type Err = String;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        match input.to_lowercase().as_str() {
            "genes" => Ok(GtfType::Genes),
            "exons" => Ok(GtfType::Exons),
            _ => Err(format!("Invalid value for GtfType: {}", input)),
        }
    }
}

#[derive(ValueEnum, Clone, Debug)]
pub enum AnalysisType {
    Bulk,
    SingleCell,
}

impl FromStr for AnalysisType {
    type Err = String;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        match input.to_lowercase().as_str() {
            "bulk" => Ok(AnalysisType::Bulk),
            "single-cell" => Ok(AnalysisType::SingleCell),
            _ => Err(format!("Invalid value for AnalysisType: {}", input)),
        }
    }
}

#[derive(ValueEnum, Clone, Debug)]
pub enum MatchType {
    Overlap,
    Exact,
}

impl FromStr for MatchType {
    type Err = String;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        match input.to_lowercase().as_str() {
            "overlap" => Ok(MatchType::Overlap),
            "exact" => Ok(MatchType::Exact),
            _ => Err(format!("Invalid value for MatchType: {}", input)),
        }
    }
}


// Function to create a mapping of reference ID to name
pub fn create_ref_id_to_name_hashmap(header_view: &Header ) -> HashMap<i32, String> {
    // Convert header to a hashmap
    // let header_map: HashMap<String, Vec<LinearMap<String, String>>> = header_view.to_hashmap();
    let header_map = header_view.to_hashmap();
    
    /*let mut chromosomes = Vec::new();

    header_view// Iterate over the reference sequence data
    if let Some(reference_info) = header_map.get("SQ") {
        for record in reference_info {
            if let Some(length) = record.get("LN") {
                // The "LN" tag holds the length of the chromosome
                if let Ok(length_int) = length.parse::<usize>() {
                    if let Some(name) = record.get("SN") {
                        chromosomes.push((name.to_string(), length_int));
                    }
                }
            }
        }
    }*/

    let mut ref_id_to_name = HashMap::new();

    // Get reference sequences ("SQ") from the header
    if let Some(reference_info) = header_map.get("SQ") {
        for (i, record) in reference_info.iter().enumerate() {
            // Get the chromosome name and ID
            if let Some(name) = record.get("SN") {
                ref_id_to_name.insert(i as i32, name.to_string());
            }
        }
    }

    ref_id_to_name
}

// Function to process data
pub fn process_data<T: FeatureMatcher + std::fmt::Debug>(
    bam_file: &str,
    mapping_info: &mut MappingInfo,
    gtf: &T,
    cell_tag: [u8; 2],
    umi_tag: [u8; 2],
    num_threads: usize,
    mutations: &Option<MutationProcessor>,
    analysis_type: &AnalysisType,
    match_type: &MatchType,
    ) -> Result<( (Scdata, IndexedGenes),(Scdata, IndexedGenes) ), String> {

    // Open BAM file with rust_htslib
    //let mut reader = Reader::from_path( bam_file  ).unwrap();
    //let header_view = Header::from_template(reader.header() );

    let bam_path = Path::new(bam_file);
    let mut reader = Reader::from_path(bam_path).map_err(|e| format!("Error opening BAM file: {}", e))?;
    reader.set_threads(4).expect("Failed to set threads");
    let header = Header::from_template(reader.header());
    let ref_id_to_name = create_ref_id_to_name_hashmap( &header );

    let m = MultiProgress::new();
    let pb = m.add(ProgressBar::new(5000));
    pb.set_style(ProgressStyle::default_bar().template("{prefix:.bold.dim} {spinner} {wide_msg}").unwrap());
    pb.set_message("");

    let mut lines = 0_u64;
    let split = BUFFER_SIZE as usize * num_threads;

    println!("Using {} processors and processing {} reads per batch", num_threads, split);

    let mut buffer = Vec::<(ReadData, Option<ReadData>)>::with_capacity(split);
    let mut expr_gex = Scdata::new(1, MatrixValueType::Integer);
    let mut expr_idx = IndexedGenes::empty(Some(0));

    let mut mut_gex = Scdata::new(1, MatrixValueType::Integer);
    let mut mut_idx = IndexedGenes::empty(Some(0));

    let mut singlets = HashMap::<String, ReadData>::new();

    for r in reader.records() {
        // Read a record from BAM file
        let record = match r {
            Ok(r) => r,
            Err(e) => panic!("I could not collect a read: {e:?}"),
        }; 

        if lines % BUFFER_SIZE == 0 {
            pb.set_message(format!("{} million reads processed", lines / BUFFER_SIZE));
            pb.inc(1);
        }
        lines += 1;


        // Choose the correct function to extract the data based on the AnalysisType
        let data_tuple = match *analysis_type {
            AnalysisType::SingleCell => ReadData::from_single_cell(&record, &ref_id_to_name, &cell_tag, &umi_tag),
            AnalysisType::Bulk => ReadData::from_bulk(&record, &ref_id_to_name, &umi_tag, lines, "1"),
        };

        match data_tuple {
            Ok(res) => {
                let qname = &res.0; // Cell ID or read name as key
                let read_data = res.1.clone(); // Clone to avoid ownership issues

                if read_data.is("paired") {
                    match singlets.remove(qname) {
                        Some(first_read) => {
                            // Mate found! Process the pair
                            //println!("Found paired reads: {:?} <-> {:?}", first_read, read_data);
                            buffer.push((first_read, Some(read_data)));
                        }
                        None => {
                            // Handle orphaned reads (mate unmapped)
                            if read_data.is("mate_unmapped") {
                                //println!("Orphaned read (mate unmapped): {:?}", read_data);
                                buffer.push((read_data, None)); // Process it as a single read
                            } else {
                                // No mate found, store this read for future pairing
                                singlets.insert(qname.to_string(), read_data.clone());
                                //println!("Storing read for future pairing: {:?}", read_data);
                            }
                        }
                    }
                } else {
                    // Unpaired read (not part of a pair at all)
                    //println!("Unpaired read: {:?}", read_data);
                    buffer.push((read_data, None));
                }
            }
            Err(err) => {
                mapping_info.report(err);
            }
        }
        

        // Process buffer when it reaches the split size
        if buffer.len() >= split {
            pb.set_message(format!("{} million reads - processing", lines / BUFFER_SIZE));
            process_buffer(
                &buffer,
                num_threads,
                &mut expr_gex,
                &mut expr_idx,
                &mut mut_gex,
                &mut mut_idx,
                mapping_info,
                gtf,
                mutations,
                match_type,
            );
            pb.set_message(format!("{} million reads - processing finished", lines / BUFFER_SIZE));
            buffer.clear();
        }
    }

    // Process remaining buffer
    if !buffer.is_empty() {
        pb.set_message(format!("{} million reads - processing", lines / BUFFER_SIZE));
        process_buffer(
            &buffer,
            num_threads,
            &mut expr_gex,
            &mut expr_idx,
            &mut mut_gex,
            &mut mut_idx,
            mapping_info,
            gtf,
            mutations,
            match_type,
        );
        pb.set_message(format!("{} million reads - processing finished", lines / BUFFER_SIZE));
    }

    return Ok(((expr_gex, expr_idx), (mut_gex, mut_idx)));
}


// Function to process data
pub fn process_data_bowtie2<T: FeatureMatcher + std::fmt::Debug>(
    bam_file: &str,
    mapping_info: &mut MappingInfo,
    gtf: &T,
    num_threads: usize,
    mutations: &Option<MutationProcessor>,
    analysis_type: &AnalysisType,
    match_type: &MatchType,
    ) -> Result<( (Scdata, IndexedGenes),(Scdata, IndexedGenes) ), String> {

    // Open BAM file with rust_htslib
    //let mut reader = Reader::from_path( bam_file  ).unwrap();
    //let header_view = Header::from_template(reader.header() );

    let bam_path = Path::new(bam_file);
    let mut reader = Reader::from_path(bam_path).map_err(|e| format!("Error opening BAM file: {}", e))?;
    reader.set_threads(4).expect("Failed to set threads");
    let header = Header::from_template(reader.header());
    let ref_id_to_name = create_ref_id_to_name_hashmap( &header );

    let m = MultiProgress::new();
    let pb = m.add(ProgressBar::new(5000));
    pb.set_style(ProgressStyle::default_bar().template("{prefix:.bold.dim} {spinner} {wide_msg}").unwrap());
    pb.set_message("");

    let mut lines = 0_u64;
    let split = ATAC_BUFFER_SIZE as usize * num_threads;

    println!("Using {} processors and processing {} reads per batch", num_threads, split);

    let mut buffer = Vec::<(ReadData, Option<ReadData>)>::with_capacity(split);
    let mut expr_gex = Scdata::new(1, MatrixValueType::Integer);
    let mut expr_idx = IndexedGenes::empty(Some(0));

    let mut mut_gex = Scdata::new(1, MatrixValueType::Integer);
    let mut mut_idx = IndexedGenes::empty(Some(0));

    let mut singlets = HashMap::<String, ReadData>::new();
    let mut pseuso_umi = 0_u64;

    for r in reader.records() {
        // Read a record from BAM file
        let record = match r {
            Ok(r) => r,
            Err(e) => panic!("I could not collect a read: {e:?}"),
        }; 

        if lines % BUFFER_SIZE == 0 {
            pb.set_message(format!("{} million reads processed", lines / BUFFER_SIZE));
            pb.inc(1);
        }
        lines += 1;

        // Choose the correct function to extract the data based on the AnalysisType
        let data_tuple = match analysis_type{
            AnalysisType::SingleCell => 
                ReadData::from_singlecell_bowtie2(&record, &ref_id_to_name, pseuso_umi),
            AnalysisType::Bulk => 
                ReadData::from_bulk_bowtie2(&record, &ref_id_to_name, "1", lines),

        };
        pseuso_umi +=1;
        
   
        match data_tuple {
            Ok( (qname, read_data) ) => {
                #[cfg(debug_assertions)]
                println!("I obtained read {} with data {}",qname, read_data);
                if read_data.is("paired") {
                    match singlets.remove(&qname) {
                        Some(first_read) => {
                            // Mate found! Process the pair
                            #[cfg(debug_assertions)]
                            println!("Found paired reads: {:?} <-> {:?}", first_read, read_data);
                            buffer.push((first_read, Some(read_data)));
                        }
                        None => {
                            // Handle orphaned reads (mate unmapped)
                            if read_data.is("mate_unmapped") {
                                #[cfg(debug_assertions)]
                                println!("Orphaned read (mate unmapped): {:?}", read_data);
                                buffer.push((read_data, None)); // Process it as a single read
                            } else {
                                // No mate found, store this read for future pairing
                                singlets.insert(qname.to_string(), read_data.clone());
                                #[cfg(debug_assertions)]
                                println!("Storing read for future pairing: {:?}", read_data);
                            }
                        }
                    }
                } else {
                    // Unpaired read (not part of a pair at all)
                    #[cfg(debug_assertions)]
                    println!("Unpaired read: {:?}", read_data);
                    buffer.push((read_data, None));
                }
            }
            Err(err) => {
                mapping_info.report(err);
            }
        }
        

        // Process buffer when it reaches the split size
        if buffer.len() >= split {
            pb.set_message(format!("{} million reads - processing", lines / BUFFER_SIZE));
            mapping_info.write_to_log(format!("{} million reads - processing", lines / BUFFER_SIZE));
            process_buffer(
                &buffer,
                num_threads,
                &mut expr_gex,
                &mut expr_idx,
                &mut mut_gex,
                &mut mut_idx,
                mapping_info,
                gtf,
                mutations,
                match_type,
            );
            pb.set_message(format!("{} million reads - processing finished", lines / BUFFER_SIZE));
            mapping_info.write_to_log(format!("{} million reads - processing finished", lines / BUFFER_SIZE));
            buffer.clear();
        }
    }

    // Process remaining buffer
    if !buffer.is_empty() {
        pb.set_message(format!("{} million reads - processing", lines / BUFFER_SIZE));
        process_buffer(
            &buffer,
            num_threads,
            &mut expr_gex,
            &mut expr_idx,
            &mut mut_gex,
            &mut mut_idx,
            mapping_info,
            gtf,
            mutations,
            match_type,
        );
        pb.set_message(format!("{} million reads - processing finished", lines / BUFFER_SIZE));
    }

    return Ok(((expr_gex, expr_idx), (mut_gex, mut_idx)));
}



// Function to process a buffer
fn process_buffer<T: FeatureMatcher + Sync + std::fmt::Debug>(
    buffer: &[(ReadData, Option<ReadData>)],
    num_threads: usize,
    expr_gex: &mut Scdata,
    expt_idx: &mut IndexedGenes,
    mut_gex: &mut Scdata,
    mut_idx: &mut IndexedGenes,
    mapping_info: &mut MappingInfo,
    gtf: &T,
    mutations: &Option<MutationProcessor>,
    match_type: &MatchType,
) {
    // everything outside this function is file io
    mapping_info.stop_file_io_time();
    let chunk_size = (buffer.len() / num_threads).max(1);

    let results: Vec<_> = buffer
        .par_chunks(chunk_size)
        .map(| chunk| process_chunk(chunk, gtf, mutations, match_type ))
        .collect();

    // only the above is the multiprocessor step
    mapping_info.stop_multi_processor_time();

    for (gex_results, mut_results, local_report) in results {

        // collect expr data
        let expr_trans = expt_idx.merge( &gex_results.1 );
        expr_gex.merge_re_id_genes( gex_results.0, &expr_trans );
        
        // collect mutations
        let mut_trans = mut_idx.merge( &mut_results.1 );
        mut_gex.merge_re_id_genes( mut_results.0, &mut_trans);
        
        // fill in the report
        mapping_info.merge(&local_report);
    }
    
    // and that is the main single processor part
    mapping_info.stop_single_processor_time();

}

// Function to process a chunk
fn process_chunk<T: FeatureMatcher + std::fmt::Debug>(
    chunk: &[(ReadData, Option<ReadData>)], 
    gtf: &T, 
    mutations: &Option<MutationProcessor>, 
    match_type: &MatchType, ) -> (
        ( Scdata, IndexedGenes),  
        ( Scdata, IndexedGenes), 
        MappingInfo 
    ) {
    let mut local_iterator = ExonIterator::new("part");
    // for the expression data
    let mut expr_gex = Scdata::new(1, MatrixValueType::Integer);
    let mut expr_idx = IndexedGenes::empty(Some(0)) ;
    // for the mutation counts
    let mut mut_gex = Scdata::new(1, MatrixValueType::Integer);
    let mut mut_idx = IndexedGenes::empty(Some(0)) ;

    let mut local_report = MappingInfo::new(None, 3.0, 0 );
    let mut last_chr = "unset";
    //     0        1    2      3      4    5                  6    7 
    //for (data.0, umi, start, cigar, chr, is_reverse_strand, seq, qual ) in chunk {
    for data in chunk{
        if last_chr != data.0.chromosome {
            match gtf.init_search( &data.0.chromosome, (data.0.start).try_into().unwrap(), &mut local_iterator){
                Ok(_) => {},
                Err(e) => {
                    //panic!("Does the GTF match to the bam file?! {}\n{gtf} ({e:?})",data.0.chromosome, );
                    local_report.report( &format!("{:?}", e) );
                    continue;
                    //
                },
            };
            last_chr = &data.0.chromosome ;
        }
        // I do actually not care about the error - it should have already been reported anyhow.
        let _ = gtf.process_feature(
            data,
            mutations,
            &mut local_iterator,
            &mut expr_gex,
            &mut expr_idx,
            &mut mut_gex,
            &mut mut_idx,
            &mut local_report,
            match_type,

        );
    }
    ( (expr_gex, expr_idx), (mut_gex, mut_idx), local_report )

}
