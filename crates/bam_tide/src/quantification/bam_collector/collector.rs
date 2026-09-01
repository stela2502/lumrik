// crates/bam_tide/src/quantification/bam_collector/collector.rs

use std::process::ChildStdout;
use std::thread::{self, JoinHandle};

use anyhow::{Context, Result, anyhow};

use rust_htslib::bam::{self, Read, Reader, record::Aux};

use crate::quantification::{
    bam_collector::config::BamCollectorConfig,
    chunk_processor::ChunkProcessor,
    job::{Job, JobBuilder},
    processor_options::ProcessorOptions,
    snp::SnpSideChannel,
};

use crate::results::QuantData;

use read_tag_table::ReadTagRecord;

use gtf_splice_index::{MatchOptions, SpliceIndex};

use snp_index::Genome;

const CHUNK: usize = 100_000;

pub struct BamCollectorResult {
    pub data: QuantData,
    pub snp: Option<SnpSideChannel>,
}

pub struct BamCollector {
    config: BamCollectorConfig,

    index: SpliceIndex,
    genome: Option<Genome>,

    match_options: MatchOptions,
    processor_options: ProcessorOptions,
}

pub struct BamCollectorHandle {
    handle: JoinHandle<Result<BamCollectorResult>>,
}

impl BamCollector {
    /// Creates a BAM collector from the user-facing collector configuration.
    ///
    /// Resources that do not depend on the mapper/BAM header are loaded here:
    ///
    /// - splice index
    /// - optional genome FASTA
    ///
    /// Header-dependent resources such as the SNP side-channel are created
    /// later when the input SAM/BAM stream has been opened.
    pub fn from_cli(config: BamCollectorConfig) -> Result<Self> {
        let index = SpliceIndex::load(&config.index)
            .with_context(|| format!("reading splice index {}", config.index.display()))?;

        let genome = match &config.genome {
            Some(path) => Some(
                Genome::from_fasta(path)
                    .with_context(|| format!("reading genome FASTA {}", path.display()))?,
            ),

            None => None,
        };

        if config.vcf.is_some() && genome.is_none() {
            anyhow::bail!("--vcf requires --genome");
        }

        let match_options = MatchOptions {
            require_strand: config.require_strand,

            require_exact_junction_chain: config.require_exact_junction_chain,

            max_5p_overhang_bp: config.max_5p_overhang_bp,

            max_3p_overhang_bp: config.max_3p_overhang_bp,

            allowed_intronic_gap_size: config.allowed_intronic_gap_size,
        };

        let processor_options = ProcessorOptions {
            min_mapq: config.min_mapq,

            read1_only: config.read1_only,

            require_strand: config.require_strand,

            quant_mode: config.quant_mode,

            ..ProcessorOptions::default()
        };

        Ok(Self {
            config,
            index,
            genome,
            match_options,
            processor_options,
        })
    }

    /// Starts collection from a mapper stdout stream on a background thread.
    ///
    /// The returned handle can be joined with `finish()` after mapper input
    /// has been closed.
    pub fn spawn(self, stdout: ChildStdout) -> Result<BamCollectorHandle> {
        let handle = thread::spawn(move || {
            let reader = Self::reader_from_stdout(stdout)?;

            self.collect(reader)
        });

        Ok(BamCollectorHandle { handle })
    }

    pub fn run_paths(self, paths: &[std::path::PathBuf]) -> Result<BamCollectorResult> {
        let mut data = QuantData::new();
        let Some(first_path) = paths.first() else {
            anyhow::bail!("no BAM files supplied");
        };
        data.report = mapping_info::MappingInfo::new(
            None,
            self.config.min_mapq as f32,
            self.config.max_reads.unwrap_or(usize::MAX),
        );

        data.report.start_counter();

        let mut n_seen = 0usize;

        /*
         * Open the first BAM so we have a header from which the
         * SNP side-channel can be constructed.
         */
        let reader = Reader::from_path(first_path)
            .with_context(|| format!("reading BAM {}", first_path.display()))?;

        let header = reader.header().clone();

        let snp = self.load_snp_side_channel(&header)?;

        for path in paths {
            if self.config.max_reads.is_some_and(|max| n_seen >= max) {
                break;
            }

            let mut reader = Reader::from_path(path)
                .with_context(|| format!("reading BAM {}", path.display()))?;

            self.collect_reader(&mut reader, snp.as_ref(), &mut data, &mut n_seen)
                .with_context(|| format!("collecting BAM {}", path.display()))?;
        }

        Ok(BamCollectorResult { data, snp })
    }

    fn collect(&self, mut reader: Reader) -> Result<BamCollectorResult> {
        let header = reader.header().clone();

        let snp = self.load_snp_side_channel(&header)?;

        let mut data = QuantData::new();

        data.report = mapping_info::MappingInfo::new(
            None,
            self.config.min_mapq as f32,
            self.config.max_reads.unwrap_or(usize::MAX),
        );

        data.report.start_counter();

        let mut n_seen = 0usize;

        self.collect_reader(&mut reader, snp.as_ref(), &mut data, &mut n_seen)?;

        Ok(BamCollectorResult { data, snp })
    }

    fn collect_reader(
        &self,
        reader: &mut Reader,
        snp: Option<&SnpSideChannel>,
        data: &mut QuantData,
        n_seen: &mut usize,
    ) -> Result<()> {
        let header = reader.header().clone();

        let processor = ChunkProcessor::new(
            &self.index,
            snp,
            self.match_options.clone(),
            self.processor_options.clone(),
        );

        let job_builder = JobBuilder::new(&header, self.index.chr_map(), *b"CB", *b"UB")
            .with_genome(self.genome.as_ref(), !self.config.no_genome_refine)
            .with_snp_index(snp.as_ref().map(|s| &s.index))
            .with_min_mapq(self.config.min_mapq)
            .read1_only(self.config.read1_only);

        let mut jobs = Vec::<Job>::with_capacity(CHUNK);

        for record_result in reader.records() {
            let mut record = record_result.context("BAM/SAM read error")?;

            let qname =
                std::str::from_utf8(record.qname()).context("mapper returned a non-UTF8 QNAME")?;

            let read_tag = ReadTagRecord::from_qname(qname)
                .with_context(|| format!("failed to recover ReadTagRecord from QNAME '{qname}'"))?;

            record.set_qname(read_tag.read_id.as_bytes());

            let cell = std::str::from_utf8(&read_tag.cell_seq)?;

            let cell_qual = std::str::from_utf8(&read_tag.cell_qual)?;

            let umi = std::str::from_utf8(&read_tag.umi_seq)?;

            let umi_qual = std::str::from_utf8(&read_tag.umi_qual)?;

            record.push_aux(b"CB", Aux::String(cell))?;

            record.push_aux(b"CY", Aux::String(cell_qual))?;

            record.push_aux(b"UB", Aux::String(umi))?;

            record.push_aux(b"UY", Aux::String(umi_qual))?;

            if let Some(job) = job_builder.build(&record, &mut data.report)? {
                jobs.push(job);
                *n_seen += 1;
            }

            if jobs.len() >= CHUNK {
                self.flush_jobs(&processor, &mut jobs, data)?;
            }

            if self.config.max_reads.is_some_and(|max| *n_seen >= max) {
                break;
            }
        }

        self.flush_jobs(&processor, &mut jobs, data)?;

        Ok(())
    }

    fn flush_jobs(
        &self,
        processor: &ChunkProcessor<'_>,
        jobs: &mut Vec<Job>,
        data: &mut QuantData,
    ) -> Result<()> {
        if jobs.is_empty() {
            return Ok(());
        }

        data.report.stop_file_io_time();

        processor.process_into(self.config.quant_mode, jobs, data)?;

        data.report.stop_single_processor_time();

        jobs.clear();

        Ok(())
    }

    fn load_snp_side_channel(&self, header: &bam::HeaderView) -> Result<Option<SnpSideChannel>> {
        let Some(vcf) = &self.config.vcf else {
            return Ok(None);
        };

        let chr_names = (0..header.target_count())
            .map(|tid| std::str::from_utf8(header.tid2name(tid)).map(str::to_owned))
            .collect::<std::result::Result<Vec<_>, _>>()?;

        let chr_lengths = (0..header.target_count())
            .map(|tid| header.target_len(tid).unwrap_or(0) as u32)
            .collect();

        SnpSideChannel::from_vcf_path(vcf, chr_names, chr_lengths, self.config.snp_min_anchor)
            .map(Some)
    }

    #[cfg(unix)]
    fn reader_from_stdout(stdout: ChildStdout) -> Result<Reader> {
        use std::os::fd::AsRawFd;

        let path = format!("/proc/self/fd/{}", stdout.as_raw_fd());

        Reader::from_path(&path).context("opening mapper stdout as SAM/BAM")
    }
}

impl BamCollectorHandle {
    pub fn finish(self) -> Result<BamCollectorResult> {
        self.handle
            .join()
            .map_err(|_| anyhow!("BAM collector thread panicked"))?
    }
}
