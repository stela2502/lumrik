//chunk_processor.rs
use anyhow::Result;
use std::time::Instant;
use gtf_splice_index::{MatchClass, MatchOptions, QuantClass, SpliceIndex};
use rayon::prelude::*;
use run_stats::RunStats;
use scdata::cell_data::GeneUmiHash;

use crate::quantification::job::Job;
use crate::quantification::processor_options::ProcessorOptions;
use crate::quantification::snp::SnpSideChannel;
use scdata::{QuantData, SnpUmiEvidence};
use crate::quantification::new_transcript_data::NewTranscriptData;
use gtf_splice_index::{RefBlock, Transcript};

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
        novel: &mut NewTranscriptData,
        report: &mut RunStats,
    ) -> Result<()> {
        if jobs.is_empty() {
            return Ok(());
        }

        let threads = rayon::current_num_threads().max(1);
        let chunk_size = (jobs.len() / threads).max(10_000);

        report.start_counter();
        report.start_timer("bam_tide/multi_cpu/quantify_chunk");

        let partials: Vec<(QuantData, NewTranscriptData, RunStats)> = jobs
            .par_chunks(chunk_size)
            .map(|chunk| self.process_partial(chunk))
            .collect();

        report.stop_timer("bam_tide/multi_cpu/quantify_chunk");
        report.stop_multi_processor_time();

        report.start_timer("bam_tide/single_cpu/merge_quantification");
        for (partial_data, partial_novel, partial_report) in partials {
            novel.merge(partial_novel)?;
            let merge_report = merged.merge(&partial_data);
            report.merge(&partial_report);
            report.merge(&merge_report);
        }
        report.stop_timer("bam_tide/single_cpu/merge_quantification");
        report.stop_single_processor_time();

        Ok(())
    }

    fn process_partial(&self, jobs: &[Job]) -> (QuantData, NewTranscriptData, RunStats) {
        let mut out = QuantData::standard();
        let mut novel = NewTranscriptData::new(self.idx.bin_width);
        let mut report = RunStats::new(None, 0.0, jobs.len());

        for job in jobs {
            report.report("bam_tide/profile/jobs");
            let novel_id = self.add_hit(job, &mut out, &mut novel, &mut report);
            let started = Instant::now();
            let evidence = self.add_snp_hits(job, &mut out, &mut report);
            report.add_timer_duration("bam_tide/profile/snp", started.elapsed());
            if let Some(id) = novel_id {
                let started = Instant::now();
                novel.add_molecule(id, job.cell, job.umi, evidence);
                report.add_timer_duration("bam_tide/profile/novel_molecule", started.elapsed());
            }
        }

        (out, novel, report)
    }

    fn add_hit(&self, job: &Job, out: &mut QuantData, novel: &mut NewTranscriptData, report: &mut RunStats) -> Option<usize> {
        let started = Instant::now();
        let hits = self.idx.match_features(&job.spliced, self.match_opts);
        report.add_timer_duration("bam_tide/profile/match_features", started.elapsed());

        if hits.is_empty() {
            report.report("no hit");
            return None;
        }

        report.report("bam_tide/profile/matched");
        let hit = &hits[0];
        report.report(hit.hit.class.to_string());


        let feature_umi = GeneUmiHash(hit.feature_id, job.umi);
        // Only materialize a model when the alignment proves a new junction.
        // An incomplete read with a subset of reference junctions is NOT novel.
        let novel_id = if matches!(hit.hit.class, MatchClass::JunctionMismatch) {
            report.report("bam_tide/profile/novel_attempts");
            let started = Instant::now();
            let result = self.discover_novel(job, hit.transcript, novel, report);
            report.add_timer_duration("bam_tide/profile/novel_discovery", started.elapsed());
            if result.is_some() { report.report("bam_tide/profile/novel_insertions"); }
            result
        } else { None };
        if let Some(class) = hit.hit.class.quant_class() {
            let split_name = match class {
                QuantClass::Exonic => QuantData::EXONIC,
                QuantClass::Intronic => QuantData::INTRONIC,
            };
            out.try_insert(split_name, &job.cell, feature_umi, 1.0, report);
            // Canonical expression contains every gene-associated molecule,
            // independent of whether its evidence is exonic or intronic.
            let mut exprs_report = RunStats::new(None, 0.0, 0);
            out.try_insert(
                QuantData::EXPRS,
                &job.cell,
                feature_umi,
                1.0,
                &mut exprs_report,
            );
        }
        novel_id
    }

    /// Represent one new splice junction as one event, regardless of how many
    /// reference exons it bypasses. Only build a transcript on the first sighting.
    fn discover_novel(
        &self,
        job: &Job,
        reference: &Transcript,
        novel: &mut NewTranscriptData,
        report: &mut RunStats,
    ) -> Option<usize> {
        let exons = reference.exons();
        if exons.is_empty() { return None; }

        let known = reference.junctions();
        let mut event = None;
        // Avoid allocating a Vec of observed junctions on this hot path.
        for pair in job.spliced.blocks.windows(2) {
            let donor = pair[0].end;
            let acceptor = pair[1].start;
            if acceptor <= donor ||
                acceptor - donor <= self.match_opts.allowed_intronic_gap_size ||
                known.contains(&(donor, acceptor)) {
                continue;
            }
            if event.replace((donor, acceptor)).is_some() {
                return None; // Multiple novel junctions: not a single-event model.
            }
        }
        let (donor, acceptor) = event?;
        let chr_name = self.idx.chr_names.get(reference.chr_id)?;
        let gene_name = self.idx.genes.get(reference.gene_id)?.primary_name()?;
        // The junction is an event, not an isoform ID. Reuse either annotated
        // effect variant before doing any exon reconstruction.
        let base_name = format!("{gene_name}*{chr_name}:{donor}-{acceptor}*del");
        let shifted_name = format!("{base_name}_frameshift");
        if let Some(id) = novel.gtf_splice_index.novel_transcript_id(&base_name)
            .or_else(|| novel.gtf_splice_index.novel_transcript_id(&shifted_name)) {
            return Some(id);
        }
        let frameshift = reference.coding_bases_removed_by_junction(donor, acceptor)
            .is_some_and(|removed| removed % 3 != 0);
        let name = if frameshift { shifted_name } else { base_name };

        // A splice boundary may lie inside an annotated exon. Reconstruct only
        // the first model for this event; subsequent reads reuse its ID.
        let left = exons.iter().position(|e| e.start < donor && donor <= e.end)?;
        let right = exons.iter().position(|e| e.start <= acceptor && acceptor < e.end)?;
        if right < left { return None; }
        let mut revised = Vec::with_capacity(exons.len() + 1);
        revised.extend_from_slice(&exons[..left]);
        if exons[left].start < donor {
            revised.push(RefBlock::new(exons[left].start, donor));
        }
        if acceptor < exons[right].end {
            revised.push(RefBlock::new(acceptor, exons[right].end));
        }
        revised.extend_from_slice(&exons[right + 1..]);
        if revised.is_empty() || revised == exons { return None; }

        match novel.gtf_splice_index.insert_novel_transcript(
            gene_name, &name, chr_name, reference.strand, &revised,
        ) {
            Ok(id) => { report.report("novel splice event"); Some(id) }
            Err(err) => { eprintln!("[bam_tide] WARNING: novel splice event insertion: {err}"); None }
        }
    }

    fn add_snp_hits(&self, job: &Job, out: &mut QuantData, report: &mut RunStats) -> SnpUmiEvidence {
        let Some(snp) = self.snp else {
            return SnpUmiEvidence::default();
        };

        snp.add_hits(job.aligned.as_ref(), job.cell, job.umi, out, report)
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
    fn novel_exon_skip_is_inserted_and_umi_deduplicated() {
        let gtf = concat!(
            "chr14\tsrc\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\";\n",
            "chr14\tsrc\texon\t201\t250\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\";\n",
            "chr14\tsrc\texon\t301\t350\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\";\n",
        );
        let idx = SpliceIndex::new(100)
            .from_reader(Cursor::new(gtf.as_bytes()), IdNameKeys::default()).unwrap();
        let mut read = SplicedRead::new(idx.chr_id("chr14").unwrap(), Strand::Plus,
            vec![RefBlock::new(110, 150), RefBlock::new(300, 340)]);
        read.finalize();
        let job = Job { cell: 42, umi: 99, spliced: read, aligned: None };
        let processor = ChunkProcessor::new(&idx, None, MatchOptions::default(), ProcessorOptions::default());
        let mut novel = NewTranscriptData::new(idx.bin_width);
        let mut report = RunStats::new(None, 0.0, 2);
        let id = processor.discover_novel(&job, &idx.transcripts[0], &mut novel, &mut report).unwrap();
        assert_eq!(novel.gtf_splice_index.transcripts[id].primary_name(), Some("Gene1*chr14:150-300*del"));
        novel.add_molecule(id, 42, 99, SnpUmiEvidence::default());
        novel.add_molecule(id, 42, 99, SnpUmiEvidence::default());
        let mut observed = 0;
        novel.scdata.for_each_umi(|_, _, _| observed += 1);
        assert_eq!(observed, 1);
    }

    #[test]
    fn novel_name_marks_coding_frameshift_but_not_intron_or_noncoding_loss() {
        for (cds, expected) in [
            // Exon 2 is 50 coding bases, although the genomic gap is 150 bp.
            (true, "Gene1*chr14:150-300*del_frameshift"),
            (false, "Gene1*chr14:150-300*del"),
        ] {
            let mut gtf = concat!(
                "chr14\tsrc\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\";\n",
                "chr14\tsrc\texon\t201\t250\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\";\n",
                "chr14\tsrc\texon\t301\t350\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\";\n",
            ).to_owned();
            if cds {
                gtf.push_str("chr14\tsrc\tCDS\t101\t350\t.\t+\t0\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\";\n");
            }
            let idx = SpliceIndex::new(100).from_reader(Cursor::new(gtf.as_bytes()), IdNameKeys::default()).unwrap();
            let mut read = SplicedRead::new(idx.chr_id("chr14").unwrap(), Strand::Plus,
                vec![RefBlock::new(110, 150), RefBlock::new(300, 340)]);
            read.finalize();
            let job = Job { cell: 42, umi: 99, spliced: read, aligned: None };
            let processor = ChunkProcessor::new(&idx, None, MatchOptions::default(), ProcessorOptions::default());
            let mut novel = NewTranscriptData::new(idx.bin_width);
            let mut report = RunStats::new(None, 0.0, 2);
            let id = processor.discover_novel(&job, &idx.transcripts[0], &mut novel, &mut report).unwrap();
            assert_eq!(novel.gtf_splice_index.transcripts[id].primary_name(), Some(expected));
            assert_eq!(processor.discover_novel(&job, &idx.transcripts[0], &mut novel, &mut report), Some(id));
        }
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
        let mut novel = NewTranscriptData::new(idx.bin_width);
        let mut report = RunStats::new(None, 0.0, 1);

        processor
            .process_into(&[job], &mut merged, &mut novel, &mut report)
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
        let mut novel = NewTranscriptData::new(idx.bin_width);
        let mut report = RunStats::new(None, 0.0, 1);

        processor
            .process_into(&[job], &mut merged, &mut novel, &mut report)
            .unwrap();

        // At minimum: it must not become "no hit".
        // Depending on QuantData internals, this may need adapting to your matrix API.
        assert!(
            !merged.is_empty(QuantData::EXONIC) || !merged.is_empty(QuantData::INTRONIC),
            "expected one gene or intron count from chr14/14 alias match"
        );
    }
}
