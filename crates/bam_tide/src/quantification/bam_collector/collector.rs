// crates/bam_tide/src/quantification/bam_collector/collector.rs

use std::collections::HashSet;
use std::process::ChildStdout;
use std::thread::{self, JoinHandle};

use anyhow::{Context, Result, anyhow};

use rust_htslib::bam::{self, Read, Reader, record::Aux};

use crate::quantification::{
    bam_collector::{config::BamCollectorConfig, read_group::ReadGroup},
    chunk_processor::ChunkProcessor,
    job::{Job, JobBuilder},
    processor_options::ProcessorOptions,
    snp::SnpSideChannel,
};

use crate::results::QuantData;

use read_tag_table::{ReadTagRecord, ReadTagTable};
use sc_primer::Grammar;

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
    grammar: Option<Grammar>,
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
            grammar: None,
        })
    }

    /// Supply the read grammar that defined molecule identity upstream.
    ///
    /// Encoded Lumrik QNAMEs and ordinary CB/UB-tagged BAMs remain valid
    /// without a grammar. An unbarcoded grammar (NONE / INSERT-only) enables
    /// direct paired-BAM sequence deduplication when no encoded QNAME exists.
    pub fn with_grammar(mut self, grammar: Grammar) -> Self {
        self.grammar = Some(grammar);
        self
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
        let mut seen_unbarcoded = HashSet::<(u64, u64)>::new();

        /*
         * Open the first BAM so we have a header from which the
         * SNP side-channel can be constructed.
         */
        let reader = Reader::from_path(first_path)
            .with_context(|| format!("reading BAM {}", first_path.display()))?;

        let header = reader.header().clone();

        let snp = self.load_snp_side_channel(&header)?;

        if !self.config.read_tags.read_tag_table.is_empty()
            && self.config.read_tags.read_tag_table.len() != paths.len()
        {
            anyhow::bail!(
                "Number of --read-tag-table files ({}) must match number of BAM files ({})",
                self.config.read_tags.read_tag_table.len(),
                paths.len(),
            );
        }

        for (path_id, path) in paths.iter().enumerate() {
            if self.config.max_reads.is_some_and(|max| n_seen >= max) {
                break;
            }

            let read_tag_table = if self.config.read_tags.read_tag_table.is_empty() {
                None
            } else {
                Some(
                    self.config
                        .read_tags
                        .load_for_id(path_id)
                        .with_context(|| format!("reading read-tag table for {}", path.display()))?,
                )
            };

            let mut reader = Reader::from_path(path)
                .with_context(|| format!("reading BAM {}", path.display()))?;

            self.collect_reader(
                &mut reader,
                snp.as_ref(),
                &mut data,
                &mut n_seen,
                &mut seen_unbarcoded,
                read_tag_table.as_ref(),
            )
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
        let mut seen_unbarcoded = HashSet::<(u64, u64)>::new();

        self.collect_reader(
            &mut reader,
            snp.as_ref(),
            &mut data,
            &mut n_seen,
            &mut seen_unbarcoded,
            None,
        )?;

        Ok(BamCollectorResult { data, snp })
    }

    fn collect_reader(
        &self,
        reader: &mut Reader,
        snp: Option<&SnpSideChannel>,
        data: &mut QuantData,
        n_seen: &mut usize,
        seen_unbarcoded: &mut HashSet<(u64, u64)>,
        read_tag_table: Option<&ReadTagTable>,
    ) -> Result<()> {
        let header = reader.header().clone();

        let processor = ChunkProcessor::new(
            &self.index,
            snp,
            self.match_options.clone(),
            self.processor_options.clone(),
        );

        let job_builder = JobBuilder::new(
            &header,
            self.index.chr_map(),
            self.config.cell_tag.0,
            self.config.umi_tag.0,
        )
            .with_genome(self.genome.as_ref(), !self.config.no_genome_refine)
            .with_snp_index(snp.as_ref().map(|s| &s.index))
            .with_read_tag_table(read_tag_table)
            .with_min_mapq(self.config.min_mapq)
            .read1_only(self.config.read1_only);

        let mut jobs = Vec::<Job>::with_capacity(CHUNK);
        let mut current: Option<ReadGroup> = None;

        for record_result in reader.records() {
            let record = record_result.context("BAM/SAM read error")?;

            if let Some(group) = current.as_mut() {
                if group.qname() == record.qname() {
                    group.push(record)?;
                    continue;
                }
            }

            if let Some(group) = current.take() {
                if self.process_read_group(
                    group,
                    &header,
                    &job_builder,
                    &processor,
                    &mut jobs,
                    data,
                    n_seen,
                    seen_unbarcoded,
                )? {
                    break;
                }
            }

            current = Some(ReadGroup::new(record));
        }

        if !self.config.max_reads.is_some_and(|max| *n_seen >= max) {
            if let Some(group) = current.take() {
                self.process_read_group(
                    group,
                    &header,
                    &job_builder,
                    &processor,
                    &mut jobs,
                    data,
                    n_seen,
                    seen_unbarcoded,
                )?;
            }
        }

        self.flush_jobs(&processor, &mut jobs, data)?;
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn process_read_group(
        &self,
        mut group: ReadGroup,
        header: &bam::HeaderView,
        job_builder: &JobBuilder<'_>,
        processor: &ChunkProcessor<'_>,
        jobs: &mut Vec<Job>,
        data: &mut QuantData,
        n_seen: &mut usize,
        seen_unbarcoded: &mut HashSet<(u64, u64)>,
    ) -> Result<bool> {
        let qname = std::str::from_utf8(group.qname())
            .context("BAM contains a non-UTF8 QNAME")?
            .to_owned();

        // Lumrik mapper output has the complete ReadTagRecord encoded in the
        // QNAME. This takes precedence even for a NONE grammar because the
        // FASTQ normalizer has already created and deduplicated the molecule
        // identity before mapping.
        if let Ok(read_tag) = ReadTagRecord::from_qname(&qname) {
            for record in group.records_mut() {
                record.set_qname(read_tag.read_id.as_bytes());

                let cell = std::str::from_utf8(&read_tag.cell_seq)?;
                let cell_qual = std::str::from_utf8(&read_tag.cell_qual)?;
                let umi = std::str::from_utf8(&read_tag.umi_seq)?;
                let umi_qual = std::str::from_utf8(&read_tag.umi_qual)?;

                record.push_aux(b"CB", Aux::String(cell))?;
                record.push_aux(b"CY", Aux::String(cell_qual))?;
                record.push_aux(b"UB", Aux::String(umi))?;
                record.push_aux(b"UY", Aux::String(umi_qual))?;

                self.push_job(job_builder.build(record, &mut data.report)?, jobs, n_seen);
                if self.after_job(processor, jobs, data, *n_seen)? {
                    return Ok(true);
                }
            }
            return Ok(false);
        }

        if let Some(grammar) = self.grammar.as_ref().filter(|g| g.is_unbarcoded()) {
            if !header_is_queryname_sorted(header) {
                anyhow::bail!(
                    "unbarcoded BAM quantification (NONE / INSERT-only grammar) requires a query-name sorted BAM; run `samtools sort -n -o name_sorted.bam input.bam`"
                );
            }

            let (r1, r2) = group.primary_pair_sequences()?;
            let identity = grammar
                .molecule_identity(None, None, &r1, &r2)
                .map_err(anyhow::Error::msg)?;

            if !seen_unbarcoded.insert((identity.cell_id, identity.molecule_id)) {
                data.report.report("unbarcoded PCR duplicate");
                return Ok(false);
            }
            data.report.report("unbarcoded sequence identity");

            for record in group.records() {
                self.push_job(
                    job_builder.build_with_identity(
                        record,
                        &mut data.report,
                        identity.cell_id,
                        identity.molecule_id,
                    )?,
                    jobs,
                    n_seen,
                );
                if self.after_job(processor, jobs, data, *n_seen)? {
                    return Ok(true);
                }
            }
            return Ok(false);
        }

        // Generic external BAM: preserve its existing CB/UB tags. Grouping is
        // now shared by all inputs, but ordinary quantification semantics stay
        // record-for-record identical to the previous JobBuilder path.
        for record in group.records() {
            self.push_job(job_builder.build(record, &mut data.report)?, jobs, n_seen);
            if self.after_job(processor, jobs, data, *n_seen)? {
                return Ok(true);
            }
        }

        Ok(false)
    }

    fn push_job(&self, job: Option<Job>, jobs: &mut Vec<Job>, n_seen: &mut usize) {
        if let Some(job) = job {
            jobs.push(job);
            *n_seen += 1;
        }
    }

    fn after_job(
        &self,
        processor: &ChunkProcessor<'_>,
        jobs: &mut Vec<Job>,
        data: &mut QuantData,
        n_seen: usize,
    ) -> Result<bool> {
        if jobs.len() >= CHUNK {
            self.flush_jobs(processor, jobs, data)?;
        }
        Ok(self.config.max_reads.is_some_and(|max| n_seen >= max))
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

fn header_is_queryname_sorted(header: &bam::HeaderView) -> bool {
    header_text_is_queryname_sorted(&String::from_utf8_lossy(header.as_bytes()))
}

fn header_text_is_queryname_sorted(text: &str) -> bool {
    text.lines().any(|line| {
        line.starts_with("@HD")
            && line
                .split('\t')
                .any(|field| field.eq_ignore_ascii_case("SO:queryname"))
    })
}

impl BamCollectorHandle {
    pub fn finish(self) -> Result<BamCollectorResult> {
        self.handle
            .join()
            .map_err(|_| anyhow!("BAM collector thread panicked"))?
    }
}


#[cfg(test)]
mod tests {
    use super::header_text_is_queryname_sorted;

    #[test]
    fn queryname_sort_header_is_accepted() {
        assert!(header_text_is_queryname_sorted(
            "@HD\tVN:1.6\tSO:queryname\tSS:queryname:natural\n@SQ\tSN:chr1\tLN:1000\n"
        ));
    }

    #[test]
    fn coordinate_sort_header_is_rejected_for_none_path() {
        assert!(!header_text_is_queryname_sorted(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
        ));
    }
}
