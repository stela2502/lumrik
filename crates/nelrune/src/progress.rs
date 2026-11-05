use std::fmt::Display;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Mutex, RwLock};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use anyhow::{Context, Result};
use lumrik_status::{
    ServerContent, ServerSnapshot, StatusMetric, StatusSection, memory_status, snapshot_html,
};
use run_stats::RunStats;

#[derive(Debug)]
pub struct RunProgress {
    started: Instant,
    last_report: Instant,
    last_reads: usize,
    reads_seen: usize,
    processing_started: Option<Instant>,
    processing_start_reads: usize,
    report_every: usize,
    state: Arc<RwLock<RunStatus>>,
    log: Option<Mutex<File>>,
    metrics: Option<Mutex<File>>,
    run_stats: RunStats,
}

/// Cheap periodically-updated snapshot for the health server.
///
/// RunStats stays local to the hot processing path.  Only these few values
/// are copied when a normalizer finishes a chunk, so the server never needs to
/// lock RunStats itself.
#[derive(Debug, Clone)]
pub struct RunStatus {
    pub started_unix_ms: u128,
    pub finished_unix_ms: Option<u128>,
    pub stage: String,
    pub reads_processed: usize,
    pub reads_per_second: f64,
    pub average_reads_per_second: f64,
    pub steady_reads_per_second: f64,
    pub mapper_reads: usize,
    pub mapper_export_pct: f64,
    pub accepted_pairs: usize,
    pub failed_pairs: usize,
    pub candidate_pairs: usize,
    pub feature_tag_matches: usize,
    pub paired_r1_insert_found: usize,
    pub no_usable_paired_r1_insert: usize,
    pub zero_cassette: usize,
    pub one_cassette: usize,
    pub multi_cassette: usize,
    pub emitted_molecules: usize,
    pub single_cassette_emitted_molecules: usize,
    pub multi_cassette_emitted_molecules: usize,
    pub forward_molecules: usize,
    pub reverse_molecules: usize,
    pub barcode_crosscheck_primary_valid_other_failure: usize,
    pub barcode_crosscheck_v2_exact: usize,
    pub barcode_crosscheck_v3_exact: usize,
    pub barcode_crosscheck_v4_exact: usize,
    pub barcode_crosscheck_5p_exact: usize,
    pub barcode_crosscheck_arc_exact: usize,
    pub barcode_crosscheck_corrected: usize,
    pub barcode_crosscheck_ambiguous: usize,
    pub barcode_crosscheck_unmatched: usize,
    pub cassette_diag_anchor_seed_rejected: usize,
    pub cassette_diag_anchor_valid: usize,
    pub cassette_diag_barcode_truncated: usize,
    pub cassette_diag_barcode_rejected: usize,
    pub cassette_diag_barcode_valid: usize,
    pub cassette_diag_umi_truncated: usize,
    pub cassette_diag_polyt_missing: usize,
    pub cassette_diag_other_grammar: usize,
    pub cassette_diag_valid_structure_not_emitted: usize,

    pub low_complexity: usize,
    pub low_complexity_pct: f64,
    pub no_cell_umi: usize,
    pub no_cell_umi_pct: f64,
    pub duplicates: usize,
    pub unique_genomic: usize,
    pub unique_genomic_pct: f64,
    pub unique_feature: usize,
    pub unique_feature_pct: f64,
    pub duplicate_pct: f64,
    pub unique_yield_pct: f64,
    pub bam_records_seen: usize,
    pub quantified_bam_records: usize,
    pub quant_gex_provenance: usize,
    pub quant_pcr_duplicates: usize,
    pub compatible_bam_records: usize,
    pub unmapped_bam_records: usize,
    pub retained_cells: Option<usize>,
    pub observed_exonic_cells: usize,
    pub observed_intronic_cells: usize,
    pub observed_exonic_genes: usize,
    pub observed_intronic_genes: usize,
    pub exonic_umis: usize,
    pub intronic_umis: usize,
    pub match_exact_junction_chain: usize,
    pub match_compatible: usize,
    pub match_intronic: usize,
    pub match_incompatible: usize,
    pub match_junction_mismatch: usize,
    pub match_overhang_too_large: usize,
    pub process_rss_mib: f64,
    pub process_peak_rss_mib: f64,
    pub system_available_mib: f64,
    pub input_file: Option<String>,
    pub public_url: Option<String>,
}

impl Default for RunProgress {
    fn default() -> Self {
        Self::new()
    }
}

impl Default for RunStatus {
    fn default() -> Self {
        Self {
            started_unix_ms: 0,
            finished_unix_ms: None,
            stage: "startup".to_string(),
            reads_processed: 0,
            reads_per_second: 0.0,
            average_reads_per_second: 0.0,
            steady_reads_per_second: 0.0,
            mapper_reads: 0,
            mapper_export_pct: 0.0,
            accepted_pairs: 0,
            failed_pairs: 0,
            candidate_pairs: 0,
            feature_tag_matches: 0,
            paired_r1_insert_found: 0,
            no_usable_paired_r1_insert: 0,
            zero_cassette: 0,
            one_cassette: 0,
            multi_cassette: 0,
            emitted_molecules: 0,
            single_cassette_emitted_molecules: 0,
            multi_cassette_emitted_molecules: 0,
            forward_molecules: 0,
            reverse_molecules: 0,
            barcode_crosscheck_primary_valid_other_failure: 0,
            barcode_crosscheck_v2_exact: 0,
            barcode_crosscheck_v3_exact: 0,
            barcode_crosscheck_v4_exact: 0,
            barcode_crosscheck_5p_exact: 0,
            barcode_crosscheck_arc_exact: 0,
            barcode_crosscheck_corrected: 0,
            barcode_crosscheck_ambiguous: 0,
            barcode_crosscheck_unmatched: 0,
            cassette_diag_anchor_seed_rejected: 0,
            cassette_diag_anchor_valid: 0,
            cassette_diag_barcode_truncated: 0,
            cassette_diag_barcode_rejected: 0,
            cassette_diag_barcode_valid: 0,
            cassette_diag_umi_truncated: 0,
            cassette_diag_polyt_missing: 0,
            cassette_diag_other_grammar: 0,
            cassette_diag_valid_structure_not_emitted: 0,

            low_complexity: 0,
            low_complexity_pct: 0.0,
            no_cell_umi: 0,
            no_cell_umi_pct: 0.0,
            duplicates: 0,
            unique_genomic: 0,
            unique_genomic_pct: 0.0,
            unique_feature: 0,
            unique_feature_pct: 0.0,
            duplicate_pct: 0.0,
            unique_yield_pct: 0.0,
            bam_records_seen: 0,
            quantified_bam_records: 0,
            quant_gex_provenance: 0,
            quant_pcr_duplicates: 0,
            compatible_bam_records: 0,
            unmapped_bam_records: 0,
            retained_cells: None,
            observed_exonic_cells: 0,
            observed_intronic_cells: 0,
            observed_exonic_genes: 0,
            observed_intronic_genes: 0,
            exonic_umis: 0,
            intronic_umis: 0,
            match_exact_junction_chain: 0,
            match_compatible: 0,
            match_intronic: 0,
            match_incompatible: 0,
            match_junction_mismatch: 0,
            match_overhang_too_large: 0,
            process_rss_mib: 0.0,
            process_peak_rss_mib: 0.0,
            system_available_mib: 0.0,
            input_file: None,
            public_url: None,
        }
    }
}

impl RunProgress {
    pub fn new() -> Self {
        let now = Instant::now();
        let started_unix_ms = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_millis();
        let mut status = RunStatus::default();
        status.started_unix_ms = started_unix_ms;

        Self {
            started: now,
            last_report: now,
            last_reads: 0,
            reads_seen: 0,
            processing_started: None,
            processing_start_reads: 0,
            report_every: 100_000,
            state: Arc::new(RwLock::new(status)),
            log: None,
            metrics: None,
            run_stats: RunStats::new(None, 0.0, 0),
        }
    }

    pub fn start_timer(&mut self, name: impl Into<String>) {
        self.run_stats.start_timer(name);
    }

    pub fn stop_timer(&mut self, name: &str) {
        self.run_stats.stop_timer(name);
    }

    pub fn run_stats(&self) -> &RunStats {
        &self.run_stats
    }

    pub fn run_stats_mut(&mut self) -> &mut RunStats {
        &mut self.run_stats
    }

    /// Copy the few live counters needed by the health server from a normalizer
    /// RunStats snapshot.  This is called once per processed chunk, not once
    /// per read.
    pub fn update_from_run_stats(&mut self, info: &RunStats) {
        let reads = info.get_issue_count("reads_processed");
        let low_complexity = info.get_issue_count("low_complexity");
        let no_cell_umi = info.get_issue_count("no_cell_umi");
        let barcode_crosscheck_primary_valid_other_failure = info.get_issue_count("barcode_crosscheck_primary_valid_other_failure");
        let barcode_crosscheck_v2_exact = info.get_issue_count("barcode_crosscheck_v2_exact");
        let barcode_crosscheck_v3_exact = info.get_issue_count("barcode_crosscheck_v3_exact");
        let barcode_crosscheck_v4_exact = info.get_issue_count("barcode_crosscheck_v4_exact");
        let barcode_crosscheck_5p_exact = info.get_issue_count("barcode_crosscheck_5p_exact");
        let barcode_crosscheck_arc_exact = info.get_issue_count("barcode_crosscheck_arc_exact");
        let barcode_crosscheck_corrected = info.get_issue_count("barcode_crosscheck_corrected");
        let barcode_crosscheck_ambiguous = info.get_issue_count("barcode_crosscheck_ambiguous");
        let barcode_crosscheck_unmatched = info.get_issue_count("barcode_crosscheck_unmatched");
        let cassette_diag_anchor_seed_rejected = info.get_issue_count("cassette_diag_anchor_seed_rejected");
        let cassette_diag_anchor_valid = info.get_issue_count("cassette_diag_anchor_valid");
        let cassette_diag_barcode_truncated = info.get_issue_count("cassette_diag_barcode_truncated");
        let cassette_diag_barcode_rejected = info.get_issue_count("cassette_diag_barcode_rejected");
        let cassette_diag_barcode_valid = info.get_issue_count("cassette_diag_barcode_valid");
        let cassette_diag_umi_truncated = info.get_issue_count("cassette_diag_umi_truncated");
        let cassette_diag_polyt_missing = info.get_issue_count("cassette_diag_polyt_missing");
        let cassette_diag_other_grammar = info.get_issue_count("cassette_diag_other_grammar");
        let cassette_diag_valid_structure_not_emitted = info.get_issue_count("cassette_diag_valid_structure_not_emitted");


        let duplicates = info.get_issue_count("duplicate");
        let unique_genomic = info.get_issue_count("unique_genomic");
        let unique_feature = info.get_issue_count("unique_feature");
        let accepted_pairs = info.get_issue_count("accepted_pairs");
        let failed_pairs = info.get_issue_count("failed_pairs");
        let candidate_pairs = info.get_issue_count("candidate_pairs");
        let feature_tag_matches = info.get_issue_count("feature_tag_match");
        let paired_r1_insert_found = info.get_issue_count("paired_r1_insert_found");
        let no_usable_paired_r1_insert = info.get_issue_count("no_usable_paired_r1_insert");
        let zero_cassette = info.get_issue_count("zero_cassette");
        let one_cassette = info.get_issue_count("one_cassette");
        let multi_cassette = info.get_issue_count("multi_cassette");
        let emitted_molecules = info.get_issue_count("emitted_molecules");
        let single_cassette_emitted_molecules = info.get_issue_count("single_cassette_emitted_molecules");
        let multi_cassette_emitted_molecules = info.get_issue_count("multi_cassette_emitted_molecules");
        let forward_molecules = info.get_issue_count("forward_molecules");
        let reverse_molecules = info.get_issue_count("reverse_molecules");

        let pct = |n: usize| {
            if reads == 0 {
                0.0
            } else {
                100.0 * n as f64 / reads as f64
            }
        };
        let low_complexity_pct = pct(low_complexity);
        let no_cell_umi_pct = pct(no_cell_umi);
        let duplicate_pct = pct(duplicates);
        let unique_genomic_pct = pct(unique_genomic);
        let unique_feature_pct = pct(unique_feature);
        let mapper_export_pct = unique_genomic_pct;
        let unique_yield_pct = pct(unique_genomic.saturating_add(unique_feature));
        let memory = memory_status();
        let process_rss_mib = memory.process_rss_mib;
        let process_peak_rss_mib = memory.process_peak_rss_mib;
        let system_available_mib = memory.system_available_mib;

        self.reads_seen = reads;

        let now = Instant::now();
        if self.processing_started.is_none() && reads > 0 {
            self.processing_started = Some(now);
            self.processing_start_reads = reads;
        }
        let elapsed = now.duration_since(self.last_report).as_secs_f64();
        let delta = reads.saturating_sub(self.last_reads);
        let rate = if elapsed > 0.0 {
            delta as f64 / elapsed
        } else {
            0.0
        };

        if let Ok(mut state) = self.state.write() {
            state.reads_processed = reads;
            state.reads_per_second = rate;
            state.average_reads_per_second = self.average_reads_per_second();
            state.steady_reads_per_second = self.processing_reads_per_second();
            state.mapper_reads = unique_genomic;
            state.mapper_export_pct = mapper_export_pct;
            state.accepted_pairs = accepted_pairs;
            state.failed_pairs = failed_pairs;
            state.candidate_pairs = candidate_pairs;
            state.feature_tag_matches = feature_tag_matches;
            state.paired_r1_insert_found = paired_r1_insert_found;
            state.no_usable_paired_r1_insert = no_usable_paired_r1_insert;
            state.zero_cassette = zero_cassette;
            state.one_cassette = one_cassette;
            state.multi_cassette = multi_cassette;
            state.emitted_molecules = emitted_molecules;
            state.single_cassette_emitted_molecules = single_cassette_emitted_molecules;
            state.multi_cassette_emitted_molecules = multi_cassette_emitted_molecules;
            state.forward_molecules = forward_molecules;
            state.reverse_molecules = reverse_molecules;
            state.barcode_crosscheck_primary_valid_other_failure = barcode_crosscheck_primary_valid_other_failure;
            state.barcode_crosscheck_v2_exact = barcode_crosscheck_v2_exact;
            state.barcode_crosscheck_v3_exact = barcode_crosscheck_v3_exact;
            state.barcode_crosscheck_v4_exact = barcode_crosscheck_v4_exact;
            state.barcode_crosscheck_5p_exact = barcode_crosscheck_5p_exact;
            state.barcode_crosscheck_arc_exact = barcode_crosscheck_arc_exact;
            state.barcode_crosscheck_corrected = barcode_crosscheck_corrected;
            state.barcode_crosscheck_ambiguous = barcode_crosscheck_ambiguous;
            state.barcode_crosscheck_unmatched = barcode_crosscheck_unmatched;
            state.cassette_diag_anchor_seed_rejected = cassette_diag_anchor_seed_rejected;
            state.cassette_diag_anchor_valid = cassette_diag_anchor_valid;
            state.cassette_diag_barcode_truncated = cassette_diag_barcode_truncated;
            state.cassette_diag_barcode_rejected = cassette_diag_barcode_rejected;
            state.cassette_diag_barcode_valid = cassette_diag_barcode_valid;
            state.cassette_diag_umi_truncated = cassette_diag_umi_truncated;
            state.cassette_diag_polyt_missing = cassette_diag_polyt_missing;
            state.cassette_diag_other_grammar = cassette_diag_other_grammar;
            state.cassette_diag_valid_structure_not_emitted = cassette_diag_valid_structure_not_emitted;

            state.low_complexity = low_complexity;
            state.low_complexity_pct = low_complexity_pct;
            state.no_cell_umi = no_cell_umi;
            state.no_cell_umi_pct = no_cell_umi_pct;
            state.duplicates = duplicates;
            state.unique_genomic = unique_genomic;
            state.unique_genomic_pct = unique_genomic_pct;
            state.unique_feature = unique_feature;
            state.unique_feature_pct = unique_feature_pct;
            state.duplicate_pct = duplicate_pct;
            state.unique_yield_pct = unique_yield_pct;
            state.process_rss_mib = process_rss_mib;
            state.process_peak_rss_mib = process_peak_rss_mib;
            state.system_available_mib = system_available_mib;
        }

        self.metrics_line(&format!(
            "{:.3}\t{}\t{:.3}\t{:.3}\t{:.3}\t{}\t{:.3}\t{}\t{:.3}\t{}\t{:.3}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
            self.elapsed().as_secs_f64(), reads, rate, self.average_reads_per_second(),
            self.processing_reads_per_second(), low_complexity, low_complexity_pct, no_cell_umi, pct(no_cell_umi), duplicates, duplicate_pct,
            unique_genomic, unique_feature, unique_yield_pct,
            process_rss_mib, process_peak_rss_mib, system_available_mib,
        ));
        if delta >= self.report_every || self.last_reads == 0 || reads < self.last_reads {
            let message = format!(
                "{:>12} reads | {:>10.0} reads/s | low complexity {:>10} | no cell/UMI {:>10} | duplicate {:>10} | genomic {:>10} | feature {:>10}",
                reads, rate, low_complexity, no_cell_umi, duplicates, unique_genomic, unique_feature,
            );
            eprintln!("[nelrune] {message}");
            if barcode_crosscheck_primary_valid_other_failure + barcode_crosscheck_v2_exact + barcode_crosscheck_v3_exact + barcode_crosscheck_v4_exact + barcode_crosscheck_5p_exact + barcode_crosscheck_arc_exact + barcode_crosscheck_corrected + barcode_crosscheck_ambiguous + barcode_crosscheck_unmatched > 0 {
                eprintln!("[barcode] primary valid/other failure={} | alternative exact: v2={} v3={} v4={} 5p={} arc={} | corrected={} ambiguous={} unmatched={}", barcode_crosscheck_primary_valid_other_failure, barcode_crosscheck_v2_exact, barcode_crosscheck_v3_exact, barcode_crosscheck_v4_exact, barcode_crosscheck_5p_exact, barcode_crosscheck_arc_exact, barcode_crosscheck_corrected, barcode_crosscheck_ambiguous, barcode_crosscheck_unmatched);
            }
            if cassette_diag_anchor_valid + cassette_diag_anchor_seed_rejected > 0 {
                eprintln!("[cassette] anchor_seed_rejected={} | anchor_valid={} | barcode_truncated={} | barcode_rejected={} | barcode_valid={} | umi_truncated={} | polyt_missing={} | other_grammar={} | valid_structure_not_emitted={}", cassette_diag_anchor_seed_rejected, cassette_diag_anchor_valid, cassette_diag_barcode_truncated, cassette_diag_barcode_rejected, cassette_diag_barcode_valid, cassette_diag_umi_truncated, cassette_diag_polyt_missing, cassette_diag_other_grammar, cassette_diag_valid_structure_not_emitted);
            }
            self.log_line(&message);
        }

        self.last_report = now;
        self.last_reads = reads;
    }

    /// Print the accumulated broad Nelrune timings to stderr and nelrune.log.
    pub fn report_timings(&self) {
        self.report_block("Nelrune timings", &self.run_stats);
    }

    pub fn open_log(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        let file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(path)
            .with_context(|| format!("opening Nelrune log {}", path.display()))?;
        self.log = Some(Mutex::new(file));
        self.log_line("Nelrune started");

        let metrics_path = path.with_extension("metrics.tsv");
        let mut metrics = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(&metrics_path)
            .with_context(|| format!("opening Nelrune metrics {}", metrics_path.display()))?;
        writeln!(
            metrics,
            "elapsed_s\treads_processed\treads_per_second\taverage_reads_per_second\tsteady_reads_per_second\tlow_complexity\tlow_complexity_pct\tno_cell_umi\tno_cell_umi_pct\tduplicates\tduplicate_pct\tunique_genomic\tunique_feature\tunique_yield_pct\tprocess_rss_mib\tprocess_peak_rss_mib\tsystem_available_mib"
        )?;
        self.metrics = Some(Mutex::new(metrics));
        Ok(())
    }

    pub fn with_report_every(mut self, report_every: usize) -> Self {
        self.report_every = report_every.max(1);
        self
    }

    pub fn stage(&self, message: impl AsRef<str>) {
        let message = message.as_ref();
        eprintln!("[nelrune] {message}");
        self.log_line(&format!("stage: {message}"));

        if let Ok(mut state) = self.state.write() {
            state.stage = message.to_string();
        }
    }

    pub fn report_block(&self, title: &str, report: &impl Display) {
        let text = format!("{title}\n{}\n{report}", "-".repeat(title.len()));
        eprintln!("\n{text}");
        self.log_line(&text);
    }

    pub fn set_public_url(&self, url: impl Into<String>) {
        if let Ok(mut state) = self.state.write() {
            state.public_url = Some(url.into());
        }
    }

    pub fn state_handle(&self) -> Arc<RwLock<RunStatus>> {
        Arc::clone(&self.state)
    }

    pub fn reads_seen(&self) -> usize {
        self.reads_seen
    }

    pub fn elapsed(&self) -> Duration {
        self.started.elapsed()
    }

    pub fn average_reads_per_second(&self) -> f64 {
        let elapsed = self.elapsed().as_secs_f64();
        if elapsed <= 0.0 {
            0.0
        } else {
            self.reads_seen as f64 / elapsed
        }
    }

    /// Throughput after the first progress sample.  This deliberately excludes
    /// mapper/index warmup, which otherwise dominates short STAR runs.
    pub fn processing_reads_per_second(&self) -> f64 {
        let Some(started) = self.processing_started else {
            return 0.0;
        };
        let elapsed = started.elapsed().as_secs_f64();
        if elapsed <= 0.0 {
            return 0.0;
        }
        self.reads_seen.saturating_sub(self.processing_start_reads) as f64 / elapsed
    }

    /// Publish a live quantification snapshot after a BAM chunk has been processed.
    pub fn update_quantification_live(&self, data: &scdata::QuantData, report: &RunStats) {
        let memory = memory_status();
        if let Ok(mut state) = self.state.write() {
            state.bam_records_seen = report.get_issue_count("bam_records_seen");
            // Quantification has its own counters. Never overwrite the
            // normalization-stage candidate/duplicate/yield statistics.
            state.quant_gex_provenance = report.get_issue_count("primer provenance GEX");
            state.quant_pcr_duplicates = report.pcr_duplicates;
            state.quantified_bam_records = report.get_issue_count("cell_umi from BAM tags");
            state.compatible_bam_records = report.get_issue_count("Compatible");
            state.unmapped_bam_records = report.get_issue_count("unmapped");
            state.observed_exonic_cells = data.cell_ids(scdata::QuantData::EXONIC).len();
            state.observed_intronic_cells = data.cell_ids(scdata::QuantData::INTRONIC).len();
            state.observed_exonic_genes =
                data.observed_feature_ids(scdata::QuantData::EXONIC).len();
            state.observed_intronic_genes =
                data.observed_feature_ids(scdata::QuantData::INTRONIC).len();
            state.exonic_umis = data.total_umis(scdata::QuantData::EXONIC);
            state.intronic_umis = data.total_umis(scdata::QuantData::INTRONIC);
            state.match_exact_junction_chain = report.get_issue_count("ExactJunctionChain");
            state.match_compatible = report.get_issue_count("Compatible");
            state.match_intronic = report.get_issue_count("Intronic");
            state.match_incompatible = report.get_issue_count("Incompatible");
            state.match_junction_mismatch = report.get_issue_count("JunctionMismatch");
            state.match_overhang_too_large = report.get_issue_count("OverhangTooLarge");
            state.process_rss_mib = memory.process_rss_mib;
            state.process_peak_rss_mib = memory.process_peak_rss_mib;
            state.system_available_mib = memory.system_available_mib;
        }
    }

    /// Refresh memory while a non-streaming quant stage is running.
    pub fn update_memory(&self) {
        let memory = memory_status();
        if let Ok(mut state) = self.state.write() {
            state.process_rss_mib = memory.process_rss_mib;
            state.process_peak_rss_mib = memory.process_peak_rss_mib;
            state.system_available_mib = memory.system_available_mib;
        }
    }

    pub fn set_quantification_summary(
        &self,
        bam_records_seen: usize,
        quantified_bam_records: usize,
        compatible_bam_records: usize,
        unmapped_bam_records: usize,
        retained_cells: usize,
    ) {
        if let Ok(mut state) = self.state.write() {
            state.bam_records_seen = bam_records_seen;
            state.quantified_bam_records = quantified_bam_records;
            state.compatible_bam_records = compatible_bam_records;
            state.unmapped_bam_records = unmapped_bam_records;
            state.retained_cells = Some(retained_cells);
        }

        let message = format!(
            "quantification: {bam_records_seen} BAM records | {quantified_bam_records} with cell/UMI tags | {compatible_bam_records} compatible | {unmapped_bam_records} unmapped | {retained_cells} cells retained"
        );
        eprintln!("[nelrune] {message}");
        self.log_line(&message);
    }

    pub fn input_file(&self, file: impl AsRef<str>) {
        let file = file.as_ref();
        self.log_line(&format!("input: {file}"));
        let display = Path::new(file)
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or(file)
            .to_string();
        if let Ok(mut state) = self.state.write() {
            state.input_file = Some(display);
        }
    }

    pub fn clear_input_file(&self) {
        if let Ok(mut state) = self.state.write() {
            state.input_file = None;
        }
    }

    /// Persist the final health-server state so completed runs remain inspectable.
    pub fn write_final_status(&self, outdir: &Path, run_type: &str) -> Result<()> {
        let state = self
            .state
            .read()
            .map_err(|_| anyhow::anyhow!("Nelrune status lock poisoned"))?
            .clone();
        let jaml_name = format!("{run_type}log.jaml");
        let html_name = format!("{run_type}.log.html");
        let mut yaml = File::create(outdir.join(&jaml_name))
            .with_context(|| format!("creating {jaml_name}"))?;
        writeln!(yaml, "schema: lumrik-nelrune-run-summary-v1")?;
        writeln!(yaml, "stage: {:?}", state.stage)?;
        writeln!(yaml, "started_unix_ms: {}", state.started_unix_ms)?;
        match state.finished_unix_ms {
            Some(v) => writeln!(yaml, "finished_unix_ms: {v}")?,
            None => writeln!(yaml, "finished_unix_ms: null")?,
        }
        writeln!(yaml, "reads_processed: {}", state.reads_processed)?;
        writeln!(yaml, "zero_cassette: {}", state.zero_cassette)?;
        writeln!(yaml, "one_cassette: {}", state.one_cassette)?;
        writeln!(yaml, "multi_cassette: {}", state.multi_cassette)?;
        writeln!(yaml, "emitted_molecules: {}", state.emitted_molecules)?;
        writeln!(yaml, "single_cassette_emitted_molecules: {}", state.single_cassette_emitted_molecules)?;
        writeln!(yaml, "multi_cassette_emitted_molecules: {}", state.multi_cassette_emitted_molecules)?;
        writeln!(yaml, "mapper_reads: {}", state.mapper_reads)?;
        writeln!(yaml, "accepted_pairs: {}", state.accepted_pairs)?;
        writeln!(yaml, "failed_pairs: {}", state.failed_pairs)?;
        writeln!(yaml, "candidate_pairs: {}", state.candidate_pairs)?;
        writeln!(yaml, "feature_tag_matches: {}", state.feature_tag_matches)?;
        writeln!(yaml, "duplicates: {}", state.duplicates)?;
        writeln!(yaml, "unique_genomic: {}", state.unique_genomic)?;
        writeln!(yaml, "unique_feature: {}", state.unique_feature)?;
        writeln!(yaml, "bam_records_seen: {}", state.bam_records_seen)?;
        writeln!(
            yaml,
            "quantified_bam_records: {}",
            state.quantified_bam_records
        )?;
        match state.retained_cells {
            Some(v) => writeln!(yaml, "retained_cells: {v}")?,
            None => writeln!(yaml, "retained_cells: null")?,
        }
        writeln!(
            yaml,
            "observed_exonic_cells: {}",
            state.observed_exonic_cells
        )?;
        writeln!(
            yaml,
            "observed_intronic_cells: {}",
            state.observed_intronic_cells
        )?;
        writeln!(
            yaml,
            "observed_exonic_genes: {}",
            state.observed_exonic_genes
        )?;
        writeln!(
            yaml,
            "observed_intronic_genes: {}",
            state.observed_intronic_genes
        )?;
        writeln!(yaml, "exonic_umis: {}", state.exonic_umis)?;
        writeln!(yaml, "intronic_umis: {}", state.intronic_umis)?;
        writeln!(
            yaml,
            "match_exact_junction_chain: {}",
            state.match_exact_junction_chain
        )?;
        writeln!(yaml, "match_compatible: {}", state.match_compatible)?;
        writeln!(yaml, "match_intronic: {}", state.match_intronic)?;
        writeln!(yaml, "match_incompatible: {}", state.match_incompatible)?;
        writeln!(
            yaml,
            "match_junction_mismatch: {}",
            state.match_junction_mismatch
        )?;
        writeln!(
            yaml,
            "match_overhang_too_large: {}",
            state.match_overhang_too_large
        )?;
        writeln!(yaml, "process_rss_mib: {:.3}", state.process_rss_mib)?;
        writeln!(
            yaml,
            "process_peak_rss_mib: {:.3}",
            state.process_peak_rss_mib
        )?;
        writeln!(
            yaml,
            "system_available_mib: {:.3}",
            state.system_available_mib
        )?;
        yaml.flush()?;

        let html = snapshot_html(&state.server_snapshot());
        let mut report = File::create(outdir.join(&html_name))
            .with_context(|| format!("creating {html_name}"))?;
        report.write_all(html.as_bytes())?;
        report.flush()?;
        Ok(())
    }

    pub fn finish(&mut self) {
        let rate = self.average_reads_per_second();
        let finished_unix_ms = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_millis();
        if let Ok(mut state) = self.state.write() {
            state.reads_processed = self.reads_seen;
            state.reads_per_second = rate;
            state.average_reads_per_second = rate;
            state.steady_reads_per_second = self.processing_reads_per_second();
            state.finished_unix_ms = Some(finished_unix_ms);
            state.stage = "finished".to_string();
        }
        self.log_line(&format!(
            "finished: {} reads in {:.1}s ({:.0} reads/s average)",
            self.reads_seen,
            self.elapsed().as_secs_f64(),
            rate,
        ));
    }

    fn metrics_line(&self, text: &str) {
        let Some(metrics) = &self.metrics else {
            return;
        };
        if let Ok(mut file) = metrics.lock() {
            let _ = writeln!(file, "{text}");
            let _ = file.flush();
        }
    }

    fn log_line(&self, text: &str) {
        let Some(log) = &self.log else {
            return;
        };

        if let Ok(mut file) = log.lock() {
            let _ = writeln!(file, "[+{:>10.3}s] {text}", self.elapsed().as_secs_f64());
            let _ = file.flush();
        }
    }
}

impl ServerContent for RunStatus {
    fn server_snapshot(&self) -> ServerSnapshot {
        let pct = |n: usize| {
            if self.reads_processed == 0 {
                0.0
            } else {
                100.0 * n as f64 / self.reads_processed as f64
            }
        };
        let bam_pct = |n: usize| {
            if self.bam_records_seen == 0 {
                0.0
            } else {
                100.0 * n as f64 / self.bam_records_seen as f64
            }
        };

        let mut snapshot = ServerSnapshot {
            title: "Nelrune".to_string(),
            subtitle: "Live single-cell normalization, mapping and quantification".to_string(),
            started_unix_ms: self.started_unix_ms,
            finished_unix_ms: self.finished_unix_ms,
            stage: self.stage.clone(),
            public_url: self.public_url.clone(),
            sections: vec![
                StatusSection::new(
                    "Throughput",
                    vec![
                        StatusMetric::new("Reads processed", self.reads_processed.to_string()),
                        StatusMetric::new(
                            "Current reads / second",
                            format!("{:.0}", self.reads_per_second),
                        ),
                        StatusMetric::new(
                            "Steady reads / second",
                            format!("{:.0}", self.steady_reads_per_second),
                        ),
                        StatusMetric::new(
                            "Run-average reads / second",
                            format!("{:.0}", self.average_reads_per_second),
                        ),
                        StatusMetric::new(
                            "Current FASTQ",
                            self.input_file.clone().unwrap_or_else(|| "-".to_string()),
                        ),
                    ],
                ),
                StatusSection::new(
                    "ONT cassette detection (per input read)",
                    vec![
                        StatusMetric::new("Zero cassettes", format!("{} ({:.2}%)", self.zero_cassette, pct(self.zero_cassette))),
                        StatusMetric::new("One cassette", format!("{} ({:.2}%)", self.one_cassette, pct(self.one_cassette))),
                        StatusMetric::new("Multiple cassettes", format!("{} ({:.2}%)", self.multi_cassette, pct(self.multi_cassette))),
                        StatusMetric::new("Emitted genomic molecules", self.emitted_molecules.to_string()),
                        StatusMetric::new("  From single-cassette reads", self.single_cassette_emitted_molecules.to_string()),
                        StatusMetric::new("  From multi-cassette reads", self.multi_cassette_emitted_molecules.to_string()),
                    ],
                ),
                StatusSection::new(
                    "Read routing",
                    vec![
                        StatusMetric::new(
                            "Accepted genomic molecules",
                            format!("{} ({:.2}%)", self.accepted_pairs, pct(self.accepted_pairs)),
                        ),
                        StatusMetric::new(
                            "Exported to mapper",
                            format!("{} ({:.2}%)", self.mapper_reads, self.mapper_export_pct),
                        ),
                        StatusMetric::new(
                            "Feature-tag molecules",
                            format!("{} ({:.2}%)", self.unique_feature, self.unique_feature_pct),
                        ),
                        StatusMetric::new(
                            "Feature-tag matches",
                            format!(
                                "{} ({:.2}%)",
                                self.feature_tag_matches,
                                pct(self.feature_tag_matches)
                            ),
                        ),
                        StatusMetric::new(
                            "Rejected pairs",
                            format!("{} ({:.2}%)", self.failed_pairs, pct(self.failed_pairs)),
                        ),
                        StatusMetric::new(
                            "Low-complexity reads (filtered)",
                            format!("{} ({:.2}%)", self.low_complexity, self.low_complexity_pct),
                        ),
                        StatusMetric::new(
                            "Cell / UMI not detected",
                            format!("{} ({:.2}%)", self.no_cell_umi, self.no_cell_umi_pct),
                        ),
                    ],
                ),
                StatusSection::new(
                    "Molecules",
                    vec![
                        StatusMetric::new(
                            "Multi-cassette reads",
                            format!("{} ({:.2}%)", self.multi_cassette, pct(self.multi_cassette)),
                        ),
                        StatusMetric::new(
                            "Candidate molecules",
                            format!("{} ({:.2}%)", self.candidate_pairs, pct(self.candidate_pairs)),
                        ),
                        StatusMetric::new(
                            "Duplicates",
                            format!("{} ({:.2}%)", self.duplicates, self.duplicate_pct),
                        ),
                        StatusMetric::new(
                            "Unique molecule yield",
                            format!("{:.2}%", self.unique_yield_pct),
                        ),
                        StatusMetric::new(
                            "Paired R1 insert found",
                            format!(
                                "{} ({:.2}%)",
                                self.paired_r1_insert_found,
                                pct(self.paired_r1_insert_found)
                            ),
                        ),
                        StatusMetric::new(
                            "No usable paired R1 insert",
                            format!(
                                "{} ({:.2}%)",
                                self.no_usable_paired_r1_insert,
                                pct(self.no_usable_paired_r1_insert)
                            ),
                        ),
                        StatusMetric::new(
                            "Forward / reverse",
                            format!("{} / {}", self.forward_molecules, self.reverse_molecules),
                        ),
                    ],
                ),
                StatusSection::new(
                    "Mapper / quantification",
                    vec![
                        StatusMetric::new("BAM records seen", self.bam_records_seen.to_string()),
                        StatusMetric::new("GEX provenance alignments", self.quant_gex_provenance.to_string()),
                        StatusMetric::new("Quantification PCR duplicates", self.quant_pcr_duplicates.to_string()),
                        StatusMetric::new(
                            "BAM records with cell / UMI",
                            format!(
                                "{} ({:.2}%)",
                                self.quantified_bam_records,
                                bam_pct(self.quantified_bam_records)
                            ),
                        ),
                        StatusMetric::new(
                            "Compatible",
                            format!(
                                "{} ({:.2}%)",
                                self.compatible_bam_records,
                                bam_pct(self.compatible_bam_records)
                            ),
                        ),
                        StatusMetric::new(
                            "Unmapped",
                            format!(
                                "{} ({:.2}%)",
                                self.unmapped_bam_records,
                                bam_pct(self.unmapped_bam_records)
                            ),
                        ),
                        StatusMetric::new(
                            "Cells retained",
                            self.retained_cells
                                .map(|n| n.to_string())
                                .unwrap_or_else(|| "-".to_string()),
                        ),
                    ],
                ),
                StatusSection::new(
                    "Observed expression",
                    vec![
                        StatusMetric::new("Exonic cells", self.observed_exonic_cells.to_string()),
                        StatusMetric::new(
                            "Intronic cells",
                            self.observed_intronic_cells.to_string(),
                        ),
                        StatusMetric::new("Exonic genes", self.observed_exonic_genes.to_string()),
                        StatusMetric::new(
                            "Intronic genes",
                            self.observed_intronic_genes.to_string(),
                        ),
                        StatusMetric::new("Unique exonic UMIs", self.exonic_umis.to_string()),
                        StatusMetric::new("Unique intronic UMIs", self.intronic_umis.to_string()),
                    ],
                ),
                StatusSection::new(
                    "Transcript matching",
                    vec![
                        StatusMetric::new(
                            "Exact junction chain",
                            self.match_exact_junction_chain.to_string(),
                        ),
                        StatusMetric::new("Compatible", self.match_compatible.to_string()),
                        StatusMetric::new("Positive intronic", self.match_intronic.to_string()),
                        StatusMetric::new("Incompatible", self.match_incompatible.to_string()),
                        StatusMetric::new(
                            "Junction mismatch",
                            self.match_junction_mismatch.to_string(),
                        ),
                        StatusMetric::new(
                            "Overhang too large",
                            self.match_overhang_too_large.to_string(),
                        ),
                    ],
                ),
                StatusSection::new(
                    "Cassette failure diagnostics (candidate counts)",
                    vec![
                        StatusMetric::new("anchor seed rejected", self.cassette_diag_anchor_seed_rejected.to_string()),
                        StatusMetric::new("anchor valid", self.cassette_diag_anchor_valid.to_string()),
                        StatusMetric::new("barcode truncated", self.cassette_diag_barcode_truncated.to_string()),
                        StatusMetric::new("barcode rejected", self.cassette_diag_barcode_rejected.to_string()),
                        StatusMetric::new("barcode valid", self.cassette_diag_barcode_valid.to_string()),
                        StatusMetric::new("umi truncated", self.cassette_diag_umi_truncated.to_string()),
                        StatusMetric::new("polyt missing", self.cassette_diag_polyt_missing.to_string()),
                        StatusMetric::new("other grammar", self.cassette_diag_other_grammar.to_string()),
                        StatusMetric::new("valid structure not emitted", self.cassette_diag_valid_structure_not_emitted.to_string()),
                    ],
                ),
                StatusSection::new(
                    "Barcode crosscheck (diagnostic; no rescue)",
                    vec![
                        StatusMetric::new("Primary barcode valid; other failure", self.barcode_crosscheck_primary_valid_other_failure.to_string()),
                        StatusMetric::new("v2 exact", self.barcode_crosscheck_v2_exact.to_string()),
                        StatusMetric::new("v3 exact", self.barcode_crosscheck_v3_exact.to_string()),
                        StatusMetric::new("v4 exact", self.barcode_crosscheck_v4_exact.to_string()),
                        StatusMetric::new("5p exact", self.barcode_crosscheck_5p_exact.to_string()),
                        StatusMetric::new("arc exact", self.barcode_crosscheck_arc_exact.to_string()),
                        StatusMetric::new("corrected", self.barcode_crosscheck_corrected.to_string()),
                        StatusMetric::new("ambiguous", self.barcode_crosscheck_ambiguous.to_string()),
                        StatusMetric::new("unmatched", self.barcode_crosscheck_unmatched.to_string()),
                    ],
                ),
                StatusSection::new(
                    "Memory",
                    vec![
                        StatusMetric::new(
                            "Process RSS / peak",
                            format!(
                                "{:.0} / {:.0} MiB",
                                self.process_rss_mib, self.process_peak_rss_mib
                            ),
                        ),
                        StatusMetric::new(
                            "System memory available",
                            format!("{:.0} MiB", self.system_available_mib),
                        ),
                    ],
                ),
            ],
        };
        // A standalone quantification run does not execute the FASTQ/ONT
        // normalization stage. Hide its zero-valued counters instead of
        // presenting them as meaningful molecule statistics.
        if self.bam_records_seen > 0 && self.reads_processed == 0 {
            snapshot.sections.retain(|section| {
                !matches!(section.title.as_str(), "Read routing" | "Molecules")
            });
        }
        snapshot
    }
}
