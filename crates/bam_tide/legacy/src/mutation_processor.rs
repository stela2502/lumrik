use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use rust_htslib::bam::{ Read, Reader, Header };

use crate::read_data::ReadData;
use crate::bed_data::ChrArea;


use scdata::IndexedGenes;

use scdata::scdata::MatrixValueType;
use scdata::Scdata;
#[allow(unused_variables)]
use mapping_info::MappingInfo;
use scdata::cell_data::GeneUmiHash;

use std::collections::HashSet;
use std::collections::HashMap;

use std::fmt;


fn find_area_for_pos(pos: usize, areas: &Option<HashSet<ChrArea>>) -> Option<&ChrArea> {
    if let Some(areas) = &*areas {
        areas.iter().find(|area| area.match_pos(pos))
    }else {
        None
    }
}

pub struct MutationProcessor {
    /// the offsets for each chr
    genome_offsets: Vec<usize>,
    /// the chr names (for matching to the fasta names)
    chr_names: HashMap<String, usize>,
    /// the real genome vector
    genome: Vec<u8>,
    /// The end of reads have overproportional high mutation counts - clip in percent
    mut_clip: f32,
    /// 1.0 - clip in percent
    rev_mut_clip: f32,
    /// a vec of mutation positions
    pub hist: Vec<usize>, // this is a stub - will be replaced by whatever Rustody::mapping_info provides!
} 


impl fmt::Display for MutationProcessor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "MutationProcessor {{")?;
        writeln!(f, "  chr_names: [{} entries]", self.chr_names.len())?;
        writeln!(f, "  genome: [{} bases]", self.genome.len())?;
        writeln!(f, "  mut_clip: {:.3}", self.mut_clip)?;
        writeln!(f, "  hist: [{} entries]", self.hist.len())?;
        
        // Format each entry in hist with an index: "0: 123"
        let hist_str = self.hist
            .iter()
            .enumerate()
            .map(|(i, val)| format!("    {}: {}", i, val))
            .collect::<Vec<_>>()
            .join("\n");

        writeln!(f, "  hist:\n{}", hist_str)?;
        writeln!(f, "}}")
    }
}

impl MutationProcessor {


    fn open_fasta_reader<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn BufRead>> {
        let file = File::open(&path)?;
        let reader: Box<dyn std::io::Read> = if path.as_ref().extension().map_or(false, |e| e == "gz") {
            Box::new(MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };
        Ok(Box::new(BufReader::new(reader)))
    }
    /// Initializes the MutationProcessor from a BAM header and a FASTA file
    pub fn new( bam_file: &str, fasta_path: &str, mut_clip: f32 ) -> Result<Self, String> {
        
        let reader = match Reader::from_path(bam_file) {
            Ok(r) => r,
            Err(e) => panic!("Error opening BAM file: {}", e),
        };


        let header = Header::from_template(reader.header());

        let header_map = header.to_hashmap();
        let mut total_length = 0;

        let mut chr_names = HashMap::new();
        let mut genome_offsets = Vec::new();
        let mut i = 0;

        // Parse the BAM header and extract chromosome names and lengths
        if let Some(reference_info) = header_map.get("SQ") {
            for record in reference_info {
                if let Some(length) = record.get("LN") {
                    // The "LN" tag holds the length of the chromosome
                    if let Ok(length_int) = length.parse::<usize>() {
                        if let Some(name) = record.get("SN") {
                            chr_names.insert(name.to_string(), i);
                            genome_offsets.push(total_length);
                            total_length += length_int as usize;
                            i +=1;
                        }
                    }
                }
            }
        }

        // Initialize the genome vector
        let mut genome = vec![0; total_length]; // Initialize with zeros

        // Load the FASTA file and populate the genome sequence
        let reader = Self::open_fasta_reader(&fasta_path).map_err(|e| format!("Failed to open: {e}"))?;
        let mut lines = reader.lines();

        #[allow(unused_variables)]
        let mut found_chrs = 0;
        let mut current_chr_id: Option<usize> = None;
        let mut write_pos = 0;
        let mut expected_len = 0;

        while let Some(line) = lines.next() {
            let line = line.map_err(|_| "Error reading FASTA".to_string())?;
            if line.starts_with('>') {
                // Header line
                if let Some(chr_id) = current_chr_id {
                    // Check previous chr finished
                    if write_pos != expected_len {
                        return Err(format!(
                            "Chromosome {chr_id} - sequence too short: expected {expected_len}, got {write_pos}"
                        ));
                    }
                }

                let chr_name = line[1..].split_whitespace().next().unwrap_or("");
                if let Some(&chr_id) = chr_names.get(chr_name) {
                    current_chr_id = Some(chr_id);
                    write_pos = 0;
                    let offset = genome_offsets[chr_id];
                    expected_len = if chr_id + 1 < genome_offsets.len() {
                        genome_offsets[chr_id + 1] - offset
                    } else {
                        total_length - offset
                    };
                    found_chrs += 1;
                } else {
                    current_chr_id = None; // skip unknown chromosomes
                }
            } else if let Some(chr_id) = current_chr_id {
                // Sequence line
                let offset = genome_offsets[chr_id];
                let bytes = line.as_bytes();
                let end = write_pos + bytes.len();
                if end > expected_len {
                    return Err(format!(
                        "Chromosome {} longer than expected",
                        chr_id
                    ));
                }
                genome[offset + write_pos..offset + end].copy_from_slice(bytes);
                write_pos = end;
            }
        }

        // Ensure all chromosomes from BAM header are found in the FASTA file
        // Or actually - do not do that here....
        /*
        if found_chrs != chr_names.len() {
            return Err(format!(
                "Missing {} chromosomes expected from the bam file header: {:?}",
                chr_names.len() - found_chrs,
                chr_names
            ));
        }
        */

        // Return the initialized MutationProcessor
        Ok(MutationProcessor {
            genome_offsets,
            chr_names,
            genome,
            mut_clip,
            rev_mut_clip: 1.0 - mut_clip,
            hist : vec![0, 20 ],
        })
    }

    pub fn handle_mutations(
        &self,
        data: &ReadData,
        gene_id: &str,
        mut_idx: &mut IndexedGenes,
        mut_gex: &mut Scdata,
        mapping_info: &mut MappingInfo,
        cell_id: &u64,
        read_len:usize,
        guhs: Option<HashSet::<ChrArea>>
    ) {
        
        for mutation_name in self.get_all_mutations(
            &data.chromosome,
            data.start.try_into().unwrap(),
            &data.cigar,
            read_len,
            guhs,
            mapping_info
        ) {

            let snip = format!("{gene_id}/{}", mutation_name);
            let mut_id = mut_idx.get_gene_id(&snip);
            //println!("I have added the mutation {snip} as id {mut_id}");
            let ghum = GeneUmiHash(mut_id, data.umi);

            let _ = mut_gex.try_insert(cell_id, ghum, 1_f32, mapping_info);
        }

    }


    pub fn get_all_mutations(&self, chr_name: &str, start: usize, md_tag: &str, read_len:usize, guhs: Option<HashSet::<ChrArea>>, mapping_info: &mut MappingInfo, ) -> Vec<String> {
        let mut mutations = Vec::new();
        
        // Retrieve chromosome name
        let chr_id = self.chr_names.get(chr_name).unwrap_or_else(|| panic!("I copuld not find chr {} in my data!", chr_name));
        
        let mut current_pos = start -1;
        let mut i = 0;
        
        // Iterate over the MD tag, which contains the match/mismatch information
        //let segments = md_tag.split(|c| c == 'A' || c == 'T' || c == 'C' || c == 'G');
        let mut match_len = 0;

        let rl= read_len as f32;

        let advance_if_needed = |current_pos: &mut usize, match_len: &mut usize| {
            if *match_len > 0 {
                *current_pos += *match_len;
                *match_len = 0;
            }
        };

        let ignore = |pos: usize, rl: f32| {
            let rel = pos as f32 / rl;
            /*println!( "For pos {pos} and read_length {rl}");
            println!( "I'll check this: rel < self.mut_clip || rel > self.rev_mut_clip");
            println!( "{rel} < {} || {rel} > {}", self.mut_clip, self.rev_mut_clip);*/

            rel < self.mut_clip || rel > self.rev_mut_clip
        };

        let add_to_hist = |pos: usize, rl: f32, mapping_info: &mut MappingInfo |{
            let rel = pos as f32 / rl;
            let bin = (rel * mapping_info.hist.len() as f32 ).floor() as usize;
            mapping_info.iterate_hist( bin );
        };

        for c in md_tag.chars() {
            match c {
                // Match or mismatch (number indicates match length)
                '0'..='9' => {
                    match_len = match_len * 10 + c.to_digit(10).unwrap() as usize;
                }
                // Mismatch - the character in the MD string
                'A' | 'T' | 'C' | 'G' => {
                    // First, move the current position based on the match length
                    advance_if_needed(&mut current_pos, &mut match_len);

                    if ignore(match_len, rl ) {
                        continue;
                    }
                    // When we encounter a mismatch, we construct the mutation string
                    //println!("I am getting the nucleotide at position {}: {}", self.genome_offsets[*chr_id] + current_pos, char::from(self.genome[self.genome_offsets[*chr_id] + current_pos]));
                    let ref_base = self.genome[self.genome_offsets[*chr_id] + current_pos -1 ];// Get the reference base at current position
                    let alt_base = c;

                    let addon = match find_area_for_pos( current_pos, &guhs){
                        Some(d) => format!("{}", d),
                        None => "unknown".to_string(),
                    };
                    //println!("I have obtained the addon {addon} and will return {addon}/{chr_name}:g.{current_pos}{}>{alt_base}",ref_base as char);
                    mutations.push(format!(
                        "{addon}/{chr_name}:g.{current_pos}{}>{alt_base}",
                        ref_base as char
                    ));

                    add_to_hist( match_len, rl, mapping_info );

                    // Update the current position
                    current_pos += 1;
                }
                // Handle deletions, where ^ indicates deleted bases in the reference
                '^' => {
                    advance_if_needed(&mut current_pos, &mut match_len);

                    if ignore(match_len, rl ) {
                        continue;
                    }
 
                    i += 1; // Skip the '^'
                    let mut deleted_bases = String::new();
                    while i < md_tag.len() && md_tag.chars().nth(i).unwrap().is_alphabetic() {
                        deleted_bases.push(md_tag.chars().nth(i).unwrap());
                        i += 1;
                    }

                    // We are not changing the current position since the bases were deleted in the reference
                    let addon = match find_area_for_pos( current_pos, &guhs){
                        Some(d) => format!("{}", d),
                        None => "unknown".to_string(),
                    };
                    mutations.push(format!(
                        "{}/{}:g.{}del{}",
                        addon,
                        chr_name,
                        current_pos, // Position in 1-based coordinate
                        deleted_bases
                    ));
                }
                // Skip anything else - should not be anything in fact
                _ => {}
            }
        }

        mutations
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_mutation_processor_new() {
        let fasta_path = "testData/mutation_test.fa.gz";
        let bam_path = "testData/mutation_test.bam"; // Or use a dummy path if unused

        assert!(Path::new(fasta_path).exists(), "FASTA test file not found.");
        assert!(Path::new(bam_path).exists(), "BAM test file not found.");

        match MutationProcessor::new(bam_path, fasta_path, 0.0) {
            Ok(processor) => {
                assert!(!processor.genome.is_empty(), "Genome should not be empty.");
                assert!(
                    !processor.chr_names.is_empty(),
                    "Chromosome names should be loaded."
                );
                println!("Loaded genome with {} bases.", processor.genome.len());
            }
            Err(e) => panic!("MutationProcessor failed to load: {}", e),
        }
    }
    
    use std::path::PathBuf;
    //use flate2::read::GzDecoder;
    //use std::fs::File;
    //use std::io::{BufRead, BufReader};
    use std::fs;


    #[test]
    fn test_handle_mutations() {
        // Set up dummy genome
        let fasta_path = "testData/mutation_test.fa.gz";
        /*
>chr1
00000000001111111111222222222233
01234567890123456789012345678901
ACTGACTGACTGACTGACTGACTGACTGACTG


        */
        let bam_path = "testData/mutation_test.bam"; 
        /*
@HD VN:1.6  SO:coordinate
@SQ SN:chr1 LN:32
@PG ID:samtools PN:samtools VN:1.19.2   CL:samtools view -b mutation_test.sam
@PG ID:samtools.1   PN:samtools PP:samtools VN:1.19.2   CL:samtools view -h mutation_test.bam
@PG ID:samtools.2   PN:samtools PP:samtools.1   VN:1.19.2   CL:samtools view -h mutation_test.sam
@PG ID:samtools.3   PN:samtools PP:samtools.2   VN:1.19.2   CL:samtools view -h mutation_test.bam
read1:AGCTGCTGAGATCGATACAGTAATTGGTCA    0   chr1    1   60  10M *   0   0   ACTGACTTAC* MD:Z:6C3
read2:AGCTGCTGAGATCGATACAGTTTTTGGTCA    0   chr1    4   60  10M *   0   0   ACTGACTTAC* MD:Z:3T7
read3:AGCTGCTGAGATCGATACAGTTAAAGGTCA    0   chr1    4   60  10M *   0   0   ACTGACTTAC* MD:Z:4C5
        */
        let processor = MutationProcessor::new(bam_path, fasta_path, 0.0)
            .expect("Failed to load MutationProcessor");

        let mut reader = match Reader::from_path(bam_path) {
            Ok(r) => r,
            Err(e) => panic!("Error opening BAM file: {}", e),
        };

        let mut mapping_info = MappingInfo::new(None, 3.0, 0 );
        let mut gex = Scdata::new(1, MatrixValueType::Integer);
        let mut idx = IndexedGenes::empty(Some(0));
        let cell_id = 1_u64;

        let mut id = 0_u64;
        for r in reader.records() {
            // Read a record from BAM file
            let record = match r {
                Ok(r) => r,
                Err(e) => panic!("I could not collect a read: {e:?}"),
            };
            let mut ref_id_to_name = HashMap::<i32, String>::new();
            ref_id_to_name.insert( 0, "chr1".to_string() );
            id +=1;
            match ReadData::from_singlecell_bowtie2(&record, &ref_id_to_name, id) {
                Ok((_id, read)) => {
                    // Process the read if successful
                    processor.handle_mutations(&read, "unknown", &mut idx, &mut gex, 
                        &mut mapping_info, &cell_id,30, None);
                },
                Err(e) => {
                    // Handle the error, log it, and continue processing
                    panic!("Error processing BAM record: {:?}", e);

                }
            }
            
        }

        assert_eq!( idx.get_all_gene_names().len(), 3, "detected only one mutation {:?}", idx.get_all_gene_names() );
        //MD:Z:6T3 
        //   000000000011111111112222222222
        //   012345678901234567890123456789
        //r: AGCTGCTGAGATCGATACAGTAATTGGTCA
        //g: AGCTGAGCAGATCGATACAGTAATTGGTCA

        assert_eq!( idx.get_all_gene_names()[0], "unknown/unknown/chr1:g.6A>C" );
        //start 4 MD:Z:3T7
        assert_eq!( idx.get_all_gene_names()[1], "unknown/unknown/chr1:g.6A>T" );
        // read3:AGCTGCTGAGATCGATACAGTTAAAGGTCA 0   chr1    4   60  10M *   0   0   ACTGACTTAC* MD:Z:4C5
        assert_eq!( idx.get_all_gene_names()[2], "unknown/unknown/chr1:g.7G>C" );

        // check if the cell names make sense
        let file_path = PathBuf::from("test_output/mutations");
        if file_path.exists() {
            fs::remove_dir_all(&file_path).expect("Failed to remove existing directory");
        }

        gex.write_sparse( file_path.clone(), &idx, 0 ).expect("Failed to write sparse matrix");

        let expected_files = vec![
            "barcodes.tsv.gz",
            "features.tsv.gz",
            "matrix.mtx.gz",
        ];

        for file in expected_files {
            let file = file_path.join(file);
            assert!(file.exists(), "Expected file {:?} was not created", file);
        }

    }

}
