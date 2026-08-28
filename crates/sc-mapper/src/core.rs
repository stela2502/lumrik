use anyhow::Result;
use bam_tide::fastq::FastqRecord;
use std::path::PathBuf;
use rust_htslib::bam::{Header};

use crate::process::SamReadCluster;


#[derive(Debug, Clone)]
pub struct SamRecord {
    pub qname: String,
    pub raw: String,
}


#[derive(Debug, Clone)]
pub struct MapperLaunch {
    pub mapper_bin: PathBuf,
    pub index: PathBuf,
    pub threads: usize,
    pub paired: bool,
    pub options: Vec<String>,
}


pub struct StreamingMapper {
    process: Box<dyn MapperProcessLike>,
}


impl StreamingMapper {
    pub fn new(process: Box<dyn MapperProcessLike>) -> Self {
        Self { process }
    }

    /// Submit one FASTQ read/pair to the external mapper.
    ///
    /// This does not imply that a mapping result is immediately available.
    pub fn submit(
        &mut self,
        r1: &FastqRecord,
        r2: Option<&FastqRecord>,
    ) -> Result<()> {
        self.process.write_fastq(r1, r2)
    }

    /// Return one completed mapping result if one is currently available.
    ///
    /// The returned mapping can belong to ANY previously submitted read.
    pub fn try_next(&mut self) -> Result<Option<MappingCall>> {
        let Some(cluster) = self.process.next_cluster()? else {
            return Ok(None);
        };

        Ok(Some(Self::cluster_to_call(cluster)))
    }

    /// Convenience operation:
    ///
    /// submit one read and opportunistically return one completed result.
    pub fn process(
        &mut self,
        r1: &FastqRecord,
        r2: Option<&FastqRecord>,
    ) -> Result<Option<MappingCall>> {
        self.submit(r1, r2)?;
        self.try_next()
    }

    fn cluster_to_call(cluster: SamReadCluster) -> MappingCall {
        let class = match cluster.records.len() {
            0 => MappingClass::Unmapped,
            1 => MappingClass::Unique,
            _ => MappingClass::MultiMapper,
        };

        MappingCall {
            read_id: cluster.read_id.clone(),
            class,
            gene_name: None,
            records: cluster,
        }
    }

    pub fn is_running(&mut self) -> Result<bool> {
        self.process.is_running()
    }

    pub fn header(&mut self) -> Result<&rust_htslib::bam::Header> {
        self.process.header()
    }

    pub fn header_loaded(&mut self) -> bool {
        self.process.header_loaded()
    }

    
    /// Close mapper input and return every mapping result still outstanding.
    pub fn finish(self) -> Result<Vec<MappingCall>> {
        let clusters = self.process.finish()?;

        Ok(clusters
            .into_iter()
            .map(Self::cluster_to_call)
            .collect())
    }
}

pub trait MapperProcessLike: Send {
    fn write_fastq(
        &mut self,
        r1: &FastqRecord,
        r2: Option<&FastqRecord>,
    ) -> Result<()>;

    fn next_cluster(
        &mut self,
    ) -> Result<Option<SamReadCluster>>;

    fn finish(self: Box<Self>) -> Result<Vec<SamReadCluster>>;

    fn header(&mut self) -> Result<&Header>;

    fn header_loaded(&mut self) -> bool;

    fn is_running(
        &mut self,
    ) -> Result<bool>;
}



#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MappingClass {
    Unmapped,
    Unique,
    MultiMapper,
}

pub struct MappingCall {
    pub read_id: String,
    pub class: MappingClass,
    pub gene_name: Option<String>,
    pub records: SamReadCluster,
}


