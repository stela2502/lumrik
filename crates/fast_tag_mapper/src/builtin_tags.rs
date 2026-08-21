use clap::ValueEnum;

use crate::{FastTagMapper, FeatureEntry};

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum BuiltinTagSet {
    Human,
    Mouse,
}

pub const HUMAN_SAMPLE_TAGS: [&[u8]; 12] = [
    b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG",
    b"TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG",
    b"CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT",
    b"ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT",
    b"CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG",
    b"TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC",
    b"TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT",
    b"CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT",
    b"GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG",
    b"GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC",
    b"CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC",
    b"GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG",
];

pub const MOUSE_SAMPLE_TAGS: [&[u8]; 12] = [
    b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG",
    b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC",
    b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC",
    b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT",
    b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT",
    b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC",
    b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC",
    b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG",
    b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC",
    b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC",
    b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT",
    b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC",
];

impl FastTagMapper {

    pub fn add_builtin( &mut self,  set: BuiltinTagSet ) -> usize{
        let start_id = self.feature_count() as u64;
        match set {
            BuiltinTagSet::Human => {
                for (i, seq) in HUMAN_SAMPLE_TAGS.iter().enumerate() {
                    let sample_id = i +1;
                    let ext_id = start_id + sample_id as u64;

                    self.add_feature(
                        seq,
                        FeatureEntry::bd_human(ext_id, sample_id),
                    );
                }
            }

            BuiltinTagSet::Mouse => {
                for (i, seq) in MOUSE_SAMPLE_TAGS.iter().enumerate() {
                    let sample_id = i+1;
                    let ext_id = start_id + sample_id as u64;

                    self.add_feature(
                        seq,
                        FeatureEntry::bd_mouse(ext_id, sample_id),
                    );
                }
            }
        }
        self.feature_count() - start_id as usize
    }

    pub fn human_samples() -> Self {
        let mut mapper = Self::new();
        for (i, seq) in HUMAN_SAMPLE_TAGS.iter().enumerate() {
            let id = (i + 1) as u64;
            mapper.add_feature(seq, FeatureEntry::bd_human(id, id as usize));
        }
        mapper
    }

    pub fn mouse_samples() -> Self {
        let mut mapper = Self::new();
        for (i, seq) in MOUSE_SAMPLE_TAGS.iter().enumerate() {
            let id = (i + 1) as u64;
            mapper.add_feature(seq, FeatureEntry::bd_mouse(id, id as usize));
        }
        mapper
    }

    pub fn builtin(set: BuiltinTagSet) -> Self {
        match set {
            BuiltinTagSet::Human => Self::human_samples(),
            BuiltinTagSet::Mouse => Self::mouse_samples(),
        }
    }
}
