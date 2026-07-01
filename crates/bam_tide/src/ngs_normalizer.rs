use crate::fastq::{record::FastqRecord, writer::FastqWriter};
use crate::read_tag_table::{ReadTagRecord, ReadTagTable};
use crate::index::FastTagFeatureIndex;

use anyhow::{bail, Context, Result};
use mapping_info::MappingInfo;
use sc_primer::Orientation;
use scdata::cell_data::GeneUmiHash;
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

pub struct FastqPairReader {
    r1: Box<dyn FastqRead>,
    r2: Box<dyn FastqRead>,
}

impl FastqPairReader {
    pub fn from_paths(r1: &PathBuf, r2: &PathBuf) -> Result<Self> {
        Ok(Self {
            r1: Self::open_fastq_reader(r1)?,
            r2: Self::open_fastq_reader(r2)?,
        })
    }

    pub fn next_pair(&mut self) -> Result<Option<(FastqRecord, FastqRecord)>> {
        let r1 = self.r1.next_record()?;
        let r2 = self.r2.next_record()?;

        match (r1, r2) {
            (Some(r1), Some(r2)) => Ok(Some((r1, r2))),
            (None, None) => Ok(None),
            (Some(_), None) => bail!("R1 has more records than R2"),
            (None, Some(_)) => bail!("R2 has more records than R1"),
        }
    }

    fn open_fastq_reader(path: &PathBuf) -> Result<Box<dyn FastqRead>> {
        use flate2::read::MultiGzDecoder;
        use std::fs::File;
        use std::io::{BufRead, BufReader};

        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let reader: Box<dyn BufRead + Send> = if path.extension().is_some_and(|x| x == "gz") {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        Ok(Box::new(SimpleFastqReader {
            reader,
            line: String::new(),
        }))
    }
}

trait FastqRead: Send {
    fn next_record(&mut self) -> Result<Option<FastqRecord>>;
}

struct SimpleFastqReader {
    reader: Box<dyn std::io::BufRead + Send>,
    line: String,
}

impl FastqRead for SimpleFastqReader {
    fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        use std::io::BufRead;

        self.line.clear();
        if self.reader.read_line(&mut self.line)? == 0 {
            return Ok(None);
        }

        let id_line = self.line.trim_end().to_string();
        if !id_line.starts_with('@') {
            bail!("invalid FASTQ record: expected @ header, got {id_line}");
        }

        self.line.clear();
        self.reader.read_line(&mut self.line)?;
        let seq = self.line.trim_end().as_bytes().to_vec();

        self.line.clear();
        self.reader.read_line(&mut self.line)?;
        let plus = self.line.trim_end().to_string();
        if !plus.starts_with('+') {
            bail!("invalid FASTQ record: expected + line, got {plus}");
        }

        self.line.clear();
        self.reader.read_line(&mut self.line)?;
        let qual_ascii = self.line.trim_end().as_bytes().to_vec();

        if seq.len() != qual_ascii.len() {
            bail!(
                "invalid FASTQ record {}: sequence length {} != quality length {}",
                id_line,
                seq.len(),
                qual_ascii.len()
            );
        }

        let qual: Vec<u8> = qual_ascii.iter().map(|q| q.saturating_sub(33)).collect();

        Ok(Some(FastqRecord::new(
            id_line.trim_start_matches('@').to_string(),
            &seq,
            &qual,
        )))
    }
}
