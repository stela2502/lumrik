// /src/gtf/spliced_read.rs

use regex::Regex;
use crate::gtf::gene::Exon;


#[derive(Debug)]
pub struct SplicedRead {
	pub start: usize,
	pub end: usize,
    pub exons: Vec<Exon>,
}

impl SplicedRead {

    // Create a new SplicedRead from a CIGAR string
    pub fn new(cigar: &str, initial_position: usize) -> (Self, usize) {
        let re_start = Regex::new(r"(\d+)([MIDNSHPX=Z])").unwrap();
        let mut exons = Vec::new();
        let start = initial_position;
        let mut current_position = initial_position; // Start from the provided position
        let mut current_exon_start = current_position; // Track the start of the current exon
        let mut exon_length = 0; // To accumulate length of the current exon
        let mut mutation_count = 0; // To track the number of SNPs

        for cap in re_start.captures_iter(cigar) {
            let count: usize = cap[1].parse().unwrap(); // Get the count
            let operation: char = cap[2].chars().next().unwrap(); // Get the operation type

            match operation {
                'M' | '=' => { // Match or exact match
                    exon_length += count; // Length increases with matches
                    current_position += count; // Move current position
                }
                'I' => { // Insertion
                    exon_length += count; // Insertions can contribute to the current exon length
                    current_position += count; // Move current position
                    mutation_count += count; // Count deletions as mutations
                }
                'D' => { // Deletion
                    mutation_count += count; // Count deletions as mutations
                    current_position += count; // Move current position
                    // Finalize the exon because of deletion
                    if exon_length > 0 {
                        exons.push(Exon {
                            start: current_exon_start,
                            end: current_exon_start + exon_length - 1,
                            mutations: mutation_count,
                        });
                        // Reset for the next exon
                        exon_length = 0;
                        mutation_count = 0; // Reset mutation count for the next exon
                        current_exon_start = current_position; // Set new start for the next exon
                    }
                }
                'N' => { // Skipped regions (introns)
                    // Finalize the exon because of the skipped region
                    if exon_length > 0 {
                        exons.push(Exon {
                            start: current_exon_start,
                            end: current_exon_start + exon_length - 1,
                            mutations: mutation_count,
                        });
                        // Reset for the next exon
                        exon_length = 0;
                        mutation_count = 0; // Reset mutation count for the next exon
                    }
                    current_position += count; // Move current position
                }
                'S' => { // Soft clip
                    // Soft clipped bases are not included in the exon, just track the position
                    current_position += count; // Move current position
                }
                'H' => { // Hard clip
                    // Hard clipped bases are completely ignored, simply move the position
                    current_position += count; // Move current position
                }
                'P' => { // Padding
                    // Padding is typically ignored; just move the position
                    current_position += count; // Move current position
                }
                'X' => { // Mismatch
                    mutation_count += count; // Count mismatches as mutations
                    current_position += count; // Move current position
                    // Mismatches can be part of the current exon
                }
                _ => {
                    // Handle unsupported operations if needed
                    println!("Unsupported CIGAR operation: {}", operation);
                }
            }
        }
        
        // If we ended with an exon, finalize it
        if exon_length > 0 {
            exons.push(Exon {
                start: current_exon_start,
                end: current_exon_start + exon_length - 1,
                mutations: mutation_count,
            });
        }else{
            exons.push(Exon {
                start,
                end: current_position,
                mutations: mutation_count,
            })
        }

        ( SplicedRead { 
	        start,
	        end: current_position,
	        exons 
	       }, 
        current_position ) // Return SplicedRead and current position
    }

    pub fn new_mut(md_tag: &str, start: usize) -> (Self, usize) {

        let mut current_pos = start -1;
        let mut i = 0;
        
        // Iterate over the MD tag, which contains the match/mismatch information
        //let segments = md_tag.split(|c| c == 'A' || c == 'T' || c == 'C' || c == 'G');
        let mut match_len = 0;

        for c in md_tag.chars() {
            match c {
                // Match or mismatch (number indicates match length)
                '0'..='9' => {
                    match_len = match_len * 10 + c.to_digit(10).unwrap() as usize;
                }
                // Mismatch - the character in the MD string
                'A' | 'T' | 'C' | 'G' => {
                    // First, move the current position based on the match length
                    if match_len > 0 {
                        current_pos += match_len;
                        match_len = 0;  // Reset match length after moving
                    }

                    // Update the current position
                    current_pos += 1;
                }
                // Handle deletions, where ^ indicates deleted bases in the reference
                '^' => {
                    i += 1; // Skip the '^'
                    let mut deleted_bases = String::new();
                    while i < md_tag.len() && md_tag.chars().nth(i).unwrap().is_alphabetic() {
                        deleted_bases.push(md_tag.chars().nth(i).unwrap());
                        i += 1;
                    }
                }
                // Skip anything else - should not be anything in fact
                _ => {}
            }
        }
        ( SplicedRead { 
            start,
            end: start +i,
            exons: vec![ 
            Exon{ 
                start,
                end: start +current_pos,
                mutations :0
            }
            ] 
        }, 
        start +i )
    } 
}