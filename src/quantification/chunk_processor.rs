//chunk_processor.rs
use anyhow::Result;
use gtf_splice_index::{MatchClass, MatchOptions, SpliceIndex};
use rayon::prelude::*;
use scdata::cell_data::GeneUmiHash;

use crate::quantification::cli::QuantMode;
use crate::quantification::job::Job;
use crate::quantification::snp::SnpSideChannel;
use crate::quantification::processor_options::ProcessorOptions;
use crate::results::QuantData;


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
        quant_mode: QuantMode,
        jobs: &[Job],
        merged: &mut QuantData,
    ) -> Result<()> {
        if jobs.is_empty() {
            return Ok(());
        }

        let threads = rayon::current_num_threads().max(1);
        let chunk_size = (jobs.len() / threads).max(10_000);

        let partials: Vec<QuantData> = jobs
            .par_chunks(chunk_size)
            .map(|chunk| self.process_partial(quant_mode, chunk))
            .collect();

        merged.report.stop_multi_processor_time();

        for partial in partials {
            merged.merge(&partial);
        }

        merged.report.stop_single_processor_time();

        Ok(())
    }

    fn process_partial(&self, quant_mode: QuantMode, jobs: &[Job]) -> QuantData {
        let mut out = QuantData::new();

        for job in jobs {
            match quant_mode {
                QuantMode::Gene => self.add_gene_hit(job, &mut out),
                QuantMode::Transcript => self.add_transcript_hit(job, &mut out),
            }

            self.add_snp_hits(job, &mut out);
        }

        out
    }

    fn add_gene_hit(&self, job: &Job, out: &mut QuantData) {
        let hits = self.idx.match_genes(&job.spliced, self.match_opts);

        if hits.is_empty() {
            out.report.report("no hit");
            return;
        }

        let hit = &hits[0];
        out.report.report(hit.best_hit.class.to_string());

        let feature_id = hit.gene_id as u64;
        let feature_umi = GeneUmiHash(feature_id, job.umi);

        if hit.best_hit.class == MatchClass::Intronic {
            out.intron
                .try_insert(&job.cell, feature_umi, 1.0, &mut out.report);
        } else {
            out.gene
                .try_insert(&job.cell, feature_umi, 1.0, &mut out.report);
        }
    }

    fn add_transcript_hit(&self, job: &Job, out: &mut QuantData) {
        let hits = self.idx.match_transcripts(&job.spliced, self.match_opts);

        if hits.is_empty() {
            out.report.report("no hit");
            return;
        }

        let hit = &hits[0];
        out.report.report(hit.hit.class.to_string());

        let feature_id = hit.transcript_id as u64;
        let feature_umi = GeneUmiHash(feature_id, job.umi);

        if hit.hit.class == MatchClass::Intronic {
            out.intron
                .try_insert(&job.cell, feature_umi, 1.0, &mut out.report);
        } else {
            out.gene
                .try_insert(&job.cell, feature_umi, 1.0, &mut out.report);
        }
    }

    fn add_snp_hits(&self, job: &Job, out: &mut QuantData) {
        let Some(snp) = self.snp else {
            return;
        };

        snp.add_hits(
            job.aligned.as_ref(),
            job.cell,
            job.umi,
            &mut out.snp_ref,
            &mut out.snp_alt,
            &mut out.report,
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
            vec![
                RefBlock::new(110, 150),
                RefBlock::new(200, 250),
            ],
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

        let mut merged = QuantData::new();

        processor
            .process_into(QuantMode::Gene, &[job], &mut merged)
            .unwrap();

        // At minimum: it must not become "no hit".
        // Depending on QuantData internals, this may need adapting to your matrix API.
        assert!(
            !merged.gene.is_empty() || !merged.intron.is_empty(),
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
            vec![
                RefBlock::new(110, 150),
                RefBlock::new(200, 250),
            ],
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

        let mut merged = QuantData::new();

        processor
            .process_into(QuantMode::Gene, &[job], &mut merged)
            .unwrap();

        // At minimum: it must not become "no hit".
        // Depending on QuantData internals, this may need adapting to your matrix API.
        assert!(
            !merged.gene.is_empty() || !merged.intron.is_empty(),
            "expected one gene or intron count from chr14/14 alias match"
        );
    }
}