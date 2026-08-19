use crate::fastq::{FastqRecord, FastqWriter, FastqPairReader, SimpleFastqReader, FastqRead};
use crate::read_tag_table::{ReadTagRecord, ReadTagTable};
use crate::index::FastTagFeatureIndex;

use anyhow::{bail, Context, Result};
use mapping_info::MappingInfo;
use sc_primer::Orientation;
//use scdata::cell_data::GeneUmiHash;
use scdata::{MatrixValueType, Scdata};


use fast_tag_mapper::{FastTagMapper};

use std::path::{Path, PathBuf};

pub const CHUNK_SIZE: usize = 1_000_000;

#[derive(Debug, Clone)]
pub struct NormalizedMolecule {
    pub fastq: FastqRecord,
    pub original_read_id: Option<String>,
    pub orientation: Orientation,
    pub cell_seq: Vec<u8>,
    pub cell_qual: Vec<u8>,
    pub umi_seq: Vec<u8>,
    pub umi_qual: Vec<u8>,
}

impl NormalizedMolecule {
    pub fn insert_read_tag(self, read_tags: &mut ReadTagTable) -> FastqRecord {
        read_tags.insert(ReadTagRecord::new(
            self.fastq.id.clone(),
            self.original_read_id,
            &self.cell_seq,
            &self.cell_qual,
            &self.umi_seq,
            &self.umi_qual,
        ));

        self.fastq
    }

    pub fn orientation_label(&self) -> &'static str {
        NgsNormalizerSupport::orientation_label(self.orientation)
    }
}

pub struct NormalizerPartial {
    pub fastq_records: Vec<FastqRecord>,
    pub read_tags: ReadTagTable,
    pub feature_tag_table: Scdata,
    pub stats: MappingInfo,
}

impl NormalizerPartial {
    pub fn new() -> Self {
        Self {
            fastq_records: Vec::new(),
            read_tags: ReadTagTable::new(),
            feature_tag_table: Scdata::new(1, MatrixValueType::Real),
            stats: MappingInfo::new(None, 0.0, 0),
        }
    }

    pub fn push_fastq(&mut self, record: FastqRecord) {
        self.fastq_records.push(record);
    }

    pub fn push_molecule(&mut self, molecule: NormalizedMolecule) {
        let record = molecule.insert_read_tag(&mut self.read_tags);
        self.fastq_records.push(record);
    }

    pub fn merge_into(
        self,
        stats: &mut MappingInfo,
        read_tags: &mut ReadTagTable,
        feature_tag_table: &mut Scdata,
        fastq: &mut FastqWriter,
    ) -> Result<()> {
        stats.merge(&self.stats);
        read_tags.merge(self.read_tags);
        feature_tag_table.merge(&self.feature_tag_table);

        for record in self.fastq_records {
            fastq.write(&record)?;
            stats.report("fastq_reads_written");
        }

        Ok(())
    }
}

impl Default for NormalizerPartial {
    fn default() -> Self {
        Self::new()
    }
}

pub struct NgsNormalizerSupport;

impl NgsNormalizerSupport {
    pub fn new_stats() -> MappingInfo {
        let mut info = MappingInfo::new(None, 0.0, 0);
        info.start_counter();
        info
    }

    pub fn new_feature_tag_table() -> Scdata {
        Scdata::new(1, MatrixValueType::Real)
    }

    pub fn configure_rayon_threads(threads: usize) {
        if threads > 1 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .ok();
        }
    }

    pub fn orientation_label(orientation: Orientation) -> &'static str {
        match orientation {
            Orientation::Forward => "forward",
            Orientation::ReverseComplement => "reverse_complement",
        }
    }

    pub fn report_orientation(stats: &mut MappingInfo, orientation: Orientation) {
        match orientation {
            Orientation::Forward => stats.report("forward_molecules"),
            Orientation::ReverseComplement => stats.report("reverse_molecules"),
        }
    }

    pub fn normalized_molecule_id(read_id: &str, molecule_index: usize) -> String {
        format!("{read_id}/mol{molecule_index}")
    }

    pub fn encode_sequence_id(seq: &[u8]) -> u64 {
        let mut out = 0_u64;

        for &base in seq.iter().take(32) {
            out <<= 2;
            out |= match base.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 0,
            };
        }

        out
    }


    pub fn write_feature_tag_table_if_present(
        feature_tag_table: &mut Scdata,
        mapper: Option<&FastTagMapper>,
        fastq_out: &Path,
    ) -> Result<()> {
        if feature_tag_table.is_empty() {
            return Ok(());
        }

        let Some(mapper) = mapper else {
            return Ok(());
        };

        let out_dir = fastq_out
            .with_extension("")
            .join("feature_tag_table_unfiltered");

        std::fs::create_dir_all(&out_dir)
            .with_context(|| format!("failed to create {}", out_dir.display()))?;

        let feature_index = FastTagFeatureIndex::new(mapper);
        feature_tag_table
            .write_sparse(&out_dir, &feature_index)
            .map_err(|err| anyhow::anyhow!("writing feature tag table failed: {err}"))?;

        Ok(())
    }
}
