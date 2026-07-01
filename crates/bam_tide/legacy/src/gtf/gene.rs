//use crate::gtf::ExonIterator;
use crate::gtf::SplicedRead;
use crate::gtf::exon_iterator::ReadResult;
use std::fmt;

#[derive(Debug, Clone, Copy)]
pub struct Exon {
    pub start: usize,
    pub end: usize,
    pub mutations: usize,
}



#[derive(Debug, Clone)]
pub struct Gene {
    pub gene_id: String,
    pub gene_name: String,
    pub start: usize,
    pub end: usize,
    pub sens_orientation: bool,
    exons: Vec<Exon>, // New field to store exons
}

impl fmt::Display for Gene{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let ex = self.exons.len();
        writeln!(f, "ID: {} name: {} {}-{} #exons: {}", self.gene_id, self.gene_name, self.start, self.end, ex )?;


        Ok(())
    }
} 

/// the order of these is from best to worst outcome
#[derive(Debug, PartialEq, Clone, PartialOrd)]
pub enum RegionStatus {
    InsideExon,
    SpanningBoundary,
    InsideIntron,
    ExtTag,
    AfterGene,
    BeforeGene,
}

impl RegionStatus {
    // Method to check if `self` is 'better' than `other` based on the enum ordering
    pub fn is_better(&self, other: &Self) -> bool {
        self < other
    }
}

impl Gene {
    // Constructor for Gene
    pub fn new(gene_id: &str, gene_name: &str, start: usize, end: usize, sens_orientation: bool) -> Self {

        Gene {
            start,
            end,
            sens_orientation,
            exons: Vec::new(), // Initialize the exons vector
            gene_id: gene_id.to_string(),
            gene_name: gene_name.to_string(),
        }
    }


    pub fn find_search_position(&self, pos: usize) -> Option<usize> {
        let mut low = 0;
        let mut high = self.exons.len();

        if pos < self.start {
            return Some(0)
        }
        
        // Perform a binary search for the first gene that starts after or overlaps with pos
        while low < high {
            let mid = (low + high) / 2;
            
            if self.exons[mid].start > pos {
                high = mid; // Narrow down to the left half
            } else {
                low = mid + 1; // Narrow down to the right half
            }
        }

        if low ==0 {
            return Some(low)
        }

        // At this point, `low` is the first index where the gene's start is >= pos or closest after `pos`
        if low > 0 && self.exons[low - 1].end >= pos {
            // If the gene just before the found position overlaps with pos, return it
            return Some(low - 1);
        }

        // Otherwise, check the gene at the found index (which should be the first gene starting after pos)
        if low < self.exons.len() && self.exons[low].end < pos {
            return Some(low);
        }

        // If no valid exon is found, return None
        None
    }

    // Method to add an exon to the gene
    pub fn add_exon(&mut self, start: usize, end: usize) {
        let mutations = 0;
        let exon = Exon { start, end, mutations };
        self.exons.push(exon); // Add the exon to the gene's exon list

        // Optionally, update the gene's start and end to reflect the new exon
        if start < self.start {
            self.start = start;
        }
        if end > self.end {
            self.end = end;
        }
    }

    /*
    InsideExon,
    SpanningBoundary,
    InsideIntron,
    AfterGene,
    BeforeGene,*/

    pub fn match_to(&self, spliced_read: &SplicedRead) -> ReadResult {
        let mut has_partial_match = false;
        #[allow(unused_variables)] //throws a warning otherwise...
        let mut fully_matched_exons = 0;
        let mut exon_matched = false;

        for read_exon in &spliced_read.exons {

            for gene_exon in &self.exons {
                if read_exon.start >= gene_exon.start && read_exon.end <= gene_exon.end {
                    // Full alignment of read exon within gene exon
                    fully_matched_exons += 1;
                    exon_matched = true;
                } else if read_exon.start < gene_exon.end && read_exon.end > gene_exon.start {
                    // Partial overlap with gene exon
                    has_partial_match = true;
                    exon_matched = true;
                }
            }
        }

        // Determine match type based on alignment results
        let match_type = if exon_matched && has_partial_match {
            RegionStatus::SpanningBoundary
        } else if exon_matched {
            RegionStatus::InsideExon
        }else if spliced_read.start >= self.start && spliced_read.end <= self.end {
            RegionStatus::InsideIntron // this should stop the tests
        } else if spliced_read.start > self.end {
            RegionStatus::AfterGene // this should remove the gene from the list of optional genes
        } else {
            RegionStatus::BeforeGene
        };

        ReadResult { gene: self.gene_name.to_string(), sens_orientation: self.sens_orientation,  match_type }
    }

    /// in comaprison with the function match_to this will check if the read overlaps with the feature with at least one base.
    pub fn match_to_overlap(&self, spliced_read: &SplicedRead) -> ReadResult {
        let mut has_partial_match = false;
        #[allow(unused_variables)] //throws a warning otherwise...
        let mut fully_matched_exons = 0;
        let mut exon_matched = false;

        for read_exon in &spliced_read.exons {

            for gene_exon in &self.exons {
                if read_exon.end >= gene_exon.start && read_exon.start <= gene_exon.end {
                    // Full alignment of read exon within gene exon
                    fully_matched_exons += 1;
                    exon_matched = true;
                } else if read_exon.start < gene_exon.end && read_exon.end > gene_exon.start {
                    // Partial overlap with gene exon
                    has_partial_match = true;
                    exon_matched = true;
                }
            }
        }

        // Determine match type based on alignment results
        let match_type = if exon_matched && has_partial_match {
            RegionStatus::SpanningBoundary
        } else if exon_matched {
            RegionStatus::InsideExon
        }else if spliced_read.start >= self.start && spliced_read.end <= self.end {
            RegionStatus::InsideIntron // this should stop the tests
        } else if spliced_read.start > self.end {
            RegionStatus::AfterGene // this should remove the gene from the list of optional genes
        } else {
            RegionStatus::BeforeGene
        };

        ReadResult { gene: self.gene_name.to_string(), sens_orientation: self.sens_orientation,  match_type }
    }
}