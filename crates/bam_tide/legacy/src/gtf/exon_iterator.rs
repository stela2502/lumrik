//ExonIterator.rs
use core::fmt;

use std::time::{SystemTime};

use crate::gtf::{RegionStatus};

use std::collections::{VecDeque};
use std::hash::{Hash};

// Define a struct to represent the hashable index
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ReadIndex {
    start: usize,
    cigar: String,
}

#[derive(Debug, Clone)]
pub struct ReadResult {
    pub gene: String,
    pub sens_orientation: bool,
    pub match_type: RegionStatus,
}


#[derive(Debug)]
pub struct ExonIterator {
    results: VecDeque< (ReadIndex, Vec<ReadResult> ) >, // Store results in a VecDeque
    max_size: usize,
    gene_id: usize,
    exon_id: usize,
    start: SystemTime,
    name: String,
}

/// ExonInterator is the multiprocessor storage of the current search state.
/// It stores the exon_id and gene_id of the current chromosomal area we are working on.

// Implement the Display trait for ExonIterator
impl fmt::Display for ExonIterator {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {

        let (h, min, sec, mill ) = self.since_start();
        writeln!(f, "{} processing time: {}h {} min {} sec and {} millisec", 
            self.name, h, min, sec, mill )?;

        Ok(())
    }
}


impl ExonIterator {
    // Constructor to create a new ExonIterator with std 0,0 as the ids
    pub fn new( name:&str ) -> Self {
        ExonIterator {
            results: VecDeque::with_capacity(5),
            max_size: 5,
            gene_id: 0,
            exon_id: 0,
            start: SystemTime::now(),
            name: name.to_string(),

        }
    }

    pub fn last_result_matches( &mut self, cigar: &str, start_position: usize ) -> Option<Vec<ReadResult>> {

        let read_index = ReadIndex {
            start: start_position,
            cigar: cigar.to_string(),
        };
        // Check if we have already processed this read by inspecting the VecDeque
        for (index, result) in &self.results {
            if *index == read_index {
                // Return the previously computed result
                return Some( result.clone() ); // Access ReadResult
            }
        }
        return None
    }

    pub fn add_last_match (&mut self, cigar: &str, start_position: usize, result: &Vec<ReadResult> ) {
        // Check if we need to remove the oldest entry to maintain the max size
        let read_index = ReadIndex {
            start: start_position,
            cigar: cigar.to_string(),
        };
        if self.results.len() >= self.max_size {
            self.results.pop_front(); // Remove the oldest entry
        }
        // Insert the new result into the VecDeque
        self.results.push_back((read_index.clone(), result.clone()));
    }

    pub fn since_start(&self ) -> (u128, u128, u128, u128){
        let mut milli = self.start.elapsed().unwrap().as_millis();

        let mil = milli % 1000;
        milli= (milli - mil) /1000;

        let sec = milli % 60;
        milli= (milli -sec) /60;

        let min = milli % 60;
        milli= (milli -min) /60;

        (milli, min, sec, mil )
    }

    pub fn exon_id(&self) -> usize{
        self.exon_id
    }

    pub fn gene_id(&self) -> usize{
        self.gene_id
    }

    pub fn set_exon_id(&mut self, id:usize ) {
        self.exon_id = id;
    }

    pub fn set_gene_id(&mut self, id:usize ) {
        self.gene_id = id;
    }

}