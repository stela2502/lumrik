use rust_htslib::bam::Record;
use std::ops::Deref;

#[derive(Debug, Clone)]
pub struct MapperRecord {
    pub record: rust_htslib::bam::Record,
}

impl MapperRecord {
    pub fn new(record: Record) -> Self {
        Self { record }
    }

    pub fn clean_id(&self) -> String {
        std::str::from_utf8(self.record.qname())
            .unwrap_or("")
            .split_whitespace()
            .next()
            .unwrap_or("")
            .trim_start_matches('@')
            .to_string()
    }

    pub fn qname(&self) -> &[u8] {
        self.record.qname()
    }

    pub fn flags(&self) -> u16 {
        self.record.flags()
    }

    pub fn into_inner(self) -> Record {
        self.record
    }

}

impl Deref for MapperRecord {
    type Target = rust_htslib::bam::Record;

    fn deref(&self) -> &Self::Target {
        &self.record
    }
}


impl From<Record> for MapperRecord {
    fn from(record: Record) -> Self {
        Self::new(record)
    }
}


#[derive(Debug, Clone)]
pub struct SamReadCluster {
    pub read_id: String,
    pub records: Vec<MapperRecord>,
}
