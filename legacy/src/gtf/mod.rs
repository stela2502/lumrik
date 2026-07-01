// src/gtf/mod.rs

pub mod exon_iterator;
pub mod gene;
pub mod gtf;
pub mod spliced_read;

pub use exon_iterator::ExonIterator; // Re-export `ExonIterator` to make it accessible from `crate::gtf`
pub use gtf::GTF;                    // Re-export `GTF` to make it accessible from `crate::gtf`
pub use gene::RegionStatus;          // Re-export `RegionStatus` enum to make it accessible from `crate::gtf`
pub use spliced_read::SplicedRead; 
