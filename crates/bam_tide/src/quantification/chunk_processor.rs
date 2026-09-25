//chunk_processor.rs
use anyhow::Result;
use gtf_splice_index::{MatchClass, MatchOptions, SpliceIndex};
use rayon::prelude::*;
use scdata::cell_data::GeneUmiHash;
use mapping_info::MappingInfo;

use crate::quantification::job::Job;
use crate::quantification::processor_options::ProcessorOptions;
use crate::quantification::snp::SnpSideChannel;
use scdata::QuantData;

pub struct ChunkProcessor<'a> {
    idx: &'a SpliceIndex,
    snp: Option<&'a SnpSideChannel>,
    match_opts: MatchOptions,
    #[allow(dead_code)]
    config: ProcessorOptions,
}

impl<'a> ChunkProcessor<'a> {
    pub fn new(
        idx: &'a SpliceIndex,
        snp: Option<&'a SnpSideChannel>,
        match_opts: MatchOptions,
        config: ProcessorOptions,
    ) -> Self {
        Self {
            idx,
            snp,
            match_opts,
            config,
        }
    }

    pub fn process_into(
        &self,
        jobs: &[Job],
        merged: &mut QuantData,
        report: &mut MappingInfo,
    ) -> Result<()> {
        if jobs.is_empty() {
            return Ok(());
        }

        let threads = rayon::current_num_threads().max(1);
        let chunk_size = (jobs.len() / threads).max(10_000);

        report.start_counter();
        report.start_timer("bam_tide/multi_cpu/quantify_chunk");

        let partials: Vec<(QuantData, MappingInfo)> = jobs
            .par_chunks(chunk_size)
            .map(|chunk| self.process_partial(chunk))
            .collect();

        report.stop_timer("bam_tide/multi_cpu/quantify_chunk");
        report.stop_multi_processor_time();

        report.start_timer("bam_tide/single_cpu/merge_quantification");
        for (partial_data, partial_report) in partials {
            let merge_report = merged.merge(&partial_data);
            report.merge(&partial_report);
            report.merge(&merge_report);
        }
        report.stop_timer("bam_tide/single_cpu/merge_quantification");
        report.stop_single_processor_time();

        Ok(())
    }

    fn process_partial(&self, jobs: &[Job]) -> (QuantData, MappingInfo) {
        let mut out = QuantData::standard();
        let mut report = MappingInfo::new(None, 0.0, jobs.len());

        for job in jobs {
            self.add_hit(job, &mut out, &mut report);

            self.add_snp_hits(job, &mut out, &mut report);
        }

        (out, report)
    }

    fn record_splice_mismatch_histograms(
        &self,
        job: &Job,
        transcript: &gtf_splice_index::Transcript,
        report: &mut MappingInfo,
    ) {
        const MAX_REPORTED_OFFSET_BP: i32 = 50;
        for (donor, acceptor) in transcript
            .junction_mismatch_offsets(&job.spliced, self.match_opts.allowed_intronic_gap_size)
        {
            report.observe_signed_histogram(
                "splice donor offset [bp]",
                donor,
                -MAX_REPORTED_OFFSET_BP,
                MAX_REPORTED_OFFSET_BP,
            );
            report.observe_signed_histogram(
                "splice acceptor offset [bp]",
                acceptor,
                -MAX_REPORTED_OFFSET_BP,
                MAX_REPORTED_OFFSET_BP,
            );
            report.observe_signed_histogram(
                "splice nearest-boundary miss [bp]",
                donor.abs().max(acceptor.abs()),
                0,
                MAX_REPORTED_OFFSET_BP,
            );
        }
    }

    fn add_hit(&self, job: &Job, out: &mut QuantData, report: &mut MappingInfo) {
        let hits = self.idx.match_features(&job.spliced, self.match_opts);

        if hits.is_empty() {
            report.report("no hit");
            return;
        }

        let hit = &hits[0];
        report.report(hit.hit.class.to_string());

        if matches!(
            hit.hit.class,
            MatchClass::JunctionMismatch | MatchClass::Intronic | MatchClass::Incompatible
        ) {
            self.record_splice_mismatch_histograms(job, hit.transcript, report);
        }

        let feature_umi = GeneUmiHash(hit.feature_id, job.umi);
        if let Some(class) = hit.hit.class.quant_class() {
            out.try_insert(class.as_str(), &job.cell, feature_umi, 1.0, report);
        }
    }

    fn add_snp_hits(&self, job: &Job, out: &mut QuantData, report: &mut MappingInfo) {
        let Some(snp) = self.snp else {
            return;
        };

        snp.add_hits(
            job.aligned.as_ref(),
            job.cell,
            job.umi,
            out,
            report,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use gtf_splice_index::{IdNameKeys, RefBlock, SplicedRead, Strand};
    use std::io::Cursor;

    fn build_chr14_index() -> SpliceIndex {
        let gtf = "\
chr14\tsrc\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\";\n\
chr14\tsrc\texon\t201\t250\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\";\n";

        SpliceIndex::new(100)
            .from_reader(Cursor::new(gtf.as_bytes()), IdNameKeys::default())
            .unwrap()
    }

    fn build_14_index() -> SpliceIndex {
        let gtf = "\
14\tsrc\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\";\n\
14\tsrc\texon\t201\t250\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\";\n";

        SpliceIndex::new(100)
            .from_reader(Cursor::new(gtf.as_bytes()), IdNameKeys::default())
            .unwrap()
    }

    #[test]
    fn chunk_processor_matches_read_when_job_uses_plain_chr14_against_14_index() {
        let idx = build_14_index();

        // This is the real regression check:
        // the index was built from "chr14", but lookup by "14" must work.
        let chr_id = idx
            .chr_id("chr14")
            .expect("expected chr14 index to resolve plain chromosome alias '14'");

        let mut spliced = SplicedRead::new(
            chr_id,
            Strand::Plus,
            vec![RefBlock::new(110, 150), RefBlock::new(200, 250)],
        );
        spliced.finalize();

        let job = Job {
            cell: 1,
            umi: 1,
            spliced,
            aligned: None,
        };

        let processor = ChunkProcessor::new(
            &idx,
            None,
            MatchOptions::default(),
            ProcessorOptions::default(),
        );

        let mut merged = QuantData::standard();
        let mut report = MappingInfo::new(None, 0.0, 1);

        processor
            .process_into(&[job], &mut merged, &mut report)
            .unwrap();

        // At minimum: it must not become "no hit".
        // Depending on QuantData internals, this may need adapting to your matrix API.
        assert!(
            !merged.is_empty(QuantData::EXONIC) || !merged.is_empty(QuantData::INTRONIC),
            "expected one gene or intron count from chr14/14 alias match"
        );
    }

    #[test]
    fn chunk_processor_matches_read_when_job_uses_plain_14_against_chr14_index() {
        let idx = build_chr14_index();

        // This is the real regression check:
        // the index was built from "chr14", but lookup by "14" must work.
        let chr_id = idx
            .chr_id("14")
            .expect("expected 14 index to resolve plain chromosome alias 'chr14'");

        let mut spliced = SplicedRead::new(
            chr_id,
            Strand::Plus,
            vec![RefBlock::new(110, 150), RefBlock::new(200, 250)],
        );
        spliced.finalize();

        let job = Job {
            cell: 1,
            umi: 1,
            spliced,
            aligned: None,
        };

        let processor = ChunkProcessor::new(
            &idx,
            None,
            MatchOptions::default(),
            ProcessorOptions::default(),
        );

        let mut merged = QuantData::standard();
        let mut report = MappingInfo::new(None, 0.0, 1);

        processor
            .process_into(&[job], &mut merged, &mut report)
            .unwrap();

        // At minimum: it must not become "no hit".
        // Depending on QuantData internals, this may need adapting to your matrix API.
        assert!(
            !merged.is_empty(QuantData::EXONIC) || !merged.is_empty(QuantData::INTRONIC),
            "expected one gene or intron count from chr14/14 alias match"
        );
    }
}
