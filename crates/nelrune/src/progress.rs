use std::fmt::Display;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Mutex, RwLock};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use anyhow::{Context, Result};
use lumrik_status::{memory_status, ServerContent, ServerSnapshot, StatusMetric, StatusSection};
use mapping_info::MappingInfo;

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
    mapping_info: MappingInfo,
}

/// Cheap periodically-updated snapshot for the health server.
///
/// MappingInfo stays local to the hot processing path.  Only these few values
/// are copied when a normalizer finishes a chunk, so the server never needs to
/// lock MappingInfo itself.
#[derive(Debug, Clone)]
pub struct RunStatus {
    pub started_unix_ms: u128,
    pub finished_unix_ms: Option<u128>,
    pub stage: String,
    pub reads_processed: usize,
    pub reads_per_second: f64,
    pub no_cell_umi: usize,
    pub duplicates: usize,
    pub unique_genomic: usize,
    pub unique_feature: usize,
    pub duplicate_pct: f64,
    pub unique_yield_pct: f64,
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
            no_cell_umi: 0,
            duplicates: 0,
            unique_genomic: 0,
            unique_feature: 0,
            duplicate_pct: 0.0,
            unique_yield_pct: 0.0,
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
            mapping_info: MappingInfo::new(None, 0.0, 0),
        }
    }

    pub fn start_timer(&mut self, name: impl Into<String>) {
        self.mapping_info.start_timer(name);
    }

    pub fn stop_timer(&mut self, name: &str) {
        self.mapping_info.stop_timer(name);
    }

    pub fn mapping_info(&self) -> &MappingInfo {
        &self.mapping_info
    }

    pub fn mapping_info_mut(&mut self) -> &mut MappingInfo {
        &mut self.mapping_info
    }

    /// Copy the few live counters needed by the health server from a normalizer
    /// MappingInfo snapshot.  This is called once per processed chunk, not once
    /// per read.
    pub fn update_from_mapping_info(&mut self, info: &MappingInfo) {
        let reads = info.get_issue_count("reads_processed");
        let no_cell_umi = info.get_issue_count("no_cell_umi");
        let duplicates = info.get_issue_count("duplicate");
        let unique_genomic = info.get_issue_count("unique_genomic");
        let unique_feature = info.get_issue_count("unique_feature");

        let pct = |n: usize| {
            if reads == 0 {
                0.0
            } else {
                100.0 * n as f64 / reads as f64
            }
        };
        let duplicate_pct = pct(duplicates);
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
            state.no_cell_umi = no_cell_umi;
            state.duplicates = duplicates;
            state.unique_genomic = unique_genomic;
            state.unique_feature = unique_feature;
            state.duplicate_pct = duplicate_pct;
            state.unique_yield_pct = unique_yield_pct;
            state.process_rss_mib = process_rss_mib;
            state.process_peak_rss_mib = process_peak_rss_mib;
            state.system_available_mib = system_available_mib;
        }

        self.metrics_line(&format!(
            "{:.3}\t{}\t{:.3}\t{:.3}\t{:.3}\t{}\t{:.3}\t{}\t{:.3}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
            self.elapsed().as_secs_f64(), reads, rate, self.average_reads_per_second(),
            self.processing_reads_per_second(), no_cell_umi, pct(no_cell_umi), duplicates, duplicate_pct,
            unique_genomic, unique_feature, unique_yield_pct,
            process_rss_mib, process_peak_rss_mib, system_available_mib,
        ));
        if delta >= self.report_every || self.last_reads == 0 || reads < self.last_reads {
            let message = format!(
                "{:>12} reads | {:>10.0} reads/s | no cell/UMI {:>10} | duplicate {:>10} | genomic {:>10} | feature {:>10}",
                reads, rate, no_cell_umi, duplicates, unique_genomic, unique_feature,
            );
            eprintln!("[nelrune] {message}");
            self.log_line(&message);
        }

        self.last_report = now;
        self.last_reads = reads;
    }

    /// Print the accumulated broad Nelrune timings to stderr and nelrune.log.
    pub fn report_timings(&self) {
        self.report_block("Nelrune timings", &self.mapping_info);
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
            "elapsed_s\treads_processed\treads_per_second\taverage_reads_per_second\tsteady_reads_per_second\tno_cell_umi\tno_cell_umi_pct\tduplicates\tduplicate_pct\tunique_genomic\tunique_feature\tunique_yield_pct\tprocess_rss_mib\tprocess_peak_rss_mib\tsystem_available_mib"
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

    pub fn input_file(&self, file: impl AsRef<str>) {
        let file = file.as_ref();
        self.log_line(&format!("input: {file}"));
        if let Ok(mut state) = self.state.write() {
            state.input_file = Some(file.to_string());
        }
    }

    pub fn clear_input_file(&self) {
        if let Ok(mut state) = self.state.write() {
            state.input_file = None;
        }
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
        ServerSnapshot {
            title: "Nelrune".to_string(),
            subtitle: "Live single-cell processing status".to_string(),
            started_unix_ms: self.started_unix_ms,
            finished_unix_ms: self.finished_unix_ms,
            stage: self.stage.clone(),
            public_url: self.public_url.clone(),
            sections: vec![
                StatusSection::new(
                    "Processing",
                    vec![
                        StatusMetric::new("Reads processed", self.reads_processed.to_string()),
                        StatusMetric::new(
                            "Reads / second",
                            format!("{:.0}", self.reads_per_second),
                        ),
                        StatusMetric::new(
                            "Current FASTQ",
                            self.input_file.clone().unwrap_or_else(|| "-".to_string()),
                        ),
                    ],
                ),
                StatusSection::new(
                    "Molecules",
                    vec![
                        StatusMetric::new("Cell / UMI not detected", self.no_cell_umi.to_string()),
                        StatusMetric::new("Duplicates", self.duplicates.to_string()),
                        StatusMetric::new("Unique genomic", self.unique_genomic.to_string()),
                        StatusMetric::new("Unique feature", self.unique_feature.to_string()),
                        StatusMetric::new(
                            "Duplicate fraction",
                            format!("{:.2}%", self.duplicate_pct),
                        ),
                        StatusMetric::new(
                            "Unique molecule yield",
                            format!("{:.2}%", self.unique_yield_pct),
                        ),
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
        }
    }
}
