use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use int_to_str::IntToStr;
use lumrik_status::{
    memory_status, public_hostname, snapshot_html, spawn_status_server, ServerContent, ServerSnapshot,
    StatusMetric, StatusSection,
};
use mapping_info::MappingInfo;
use sc_vdj::output::{write_mapping_info_report, ReportWriter};
use sc_primer::BdCellVersion;
use sc_vdj::{
    BamIngestProgress, CellEvidenceVdj, Chain, NelruneIdentityResolver, SegmentKind, VdjIndex,
    VdjIndexBuilder, VdjRunner, VdjRunnerConfig,
};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read as IoRead, Write as IoWrite};
use std::net::{IpAddr, Ipv4Addr, SocketAddr};
use std::path::{Path, PathBuf};
use std::sync::{Arc, RwLock};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    name = "nelrune-vdj",
    about = "Reconstruct single-cell V(D)J receptors from a retained Nelrune BAM",
    after_help = "LIVE STATUS\n  By default nelrune-vdj serves the Lumrik live dashboard on --health-port 8787.\n  Open the printed URL in a browser. On a cluster, use --health-hostname to control\n  the hostname shown in that URL, or forward the port over SSH.\n  Use --no-health-server to disable it."
)]
struct Cli {
    #[arg(long)]
    bam: PathBuf,
    #[arg(long)]
    index: Option<PathBuf>,
    #[arg(long)]
    gtf: Option<PathBuf>,
    #[arg(long)]
    genome: Option<PathBuf>,
    #[arg(long)]
    out: PathBuf,
    #[arg(long, default_value_t = 12)]
    min_sequence_overlap: usize,
    /// CPU threads for BAM decoding and independent V(D)J chain reconstruction.
    #[arg(long, default_value_t = 8)]
    threads: usize,
    /// Nelrune output directory or its final exonic MEX directory.
    /// When supplied, only barcodes present in the exonic output are admitted
    /// into the expensive V(D)J evidence pipeline.
    #[arg(long)]
    exonic: Option<PathBuf>,
    #[arg(long)]
    write_sequences: bool,
    /// BD/Rhapsody whitelist version used to emit the official positional
    /// Rustody cell id in AIRR output (v1, v2.96, or v2.384).
    #[arg(long)]
    bd_cell_version: Option<String>,
    /// Port for the live Lumrik V(D)J dashboard.
    #[arg(long, default_value_t = 8787)]
    health_port: u16,
    /// Hostname shown in the live-server URL. Useful on clusters.
    #[arg(long)]
    health_hostname: Option<String>,
    /// Disable the live Lumrik status server.
    #[arg(long, default_value_t = false)]
    no_health_server: bool,
}

#[derive(Debug, Clone)]
struct VdjRunStatus {
    started_unix_ms: u128,
    finished_unix_ms: Option<u128>,
    stage: String,
    public_url: Option<String>,
    threads: usize,
    phase: usize,
    phase_started: Instant,
    phase_times: [Option<Duration>; 5],
    reference_segments: usize,

    bam_records: usize,
    allowed_cell_records: usize,
    receptor_overlap_records: usize,
    unmapped_candidates: usize,
    unmapped_igh_admitted: usize,
    unmapped_igh_rescued_cells: usize,
    unmapped_igh_v_mappings: usize,
    unmapped_igh_d_mappings: usize,
    unmapped_igh_j_mappings: usize,
    unmapped_igh_c_mappings: usize,
    preliminary_cells: Option<usize>,
    evidence_cells: usize,
    compact_summaries: usize,
    physical_fragments: usize,
    cells_with_v: usize,
    cells_with_d: usize,
    cells_with_j: usize,
    cells_with_c: usize,
    j_constant_linked_cells: usize,
    j_constant_linked_fragments: usize,
    igh_vj_cells: usize,
    igh_dj_no_v_cells: usize,
    igh_j_no_vd_cells: usize,
    igh_intronic_c_cells: usize,
    igh_intronic_c_records: usize,

    calls_cells: usize,
    recombinations: usize,
    productive: usize,
    mean_calls_per_cell: f64,
    calls_by_chain: [usize; 7],
    knee_thresholds: [usize; 7],
    knee_evidence_cells: [usize; 7],
    knee_selected_cells: [usize; 7],

    rescan_records: usize,
    rescan_wanted_records: usize,
    rescan_batches: usize,
    receptor_rediscovery_hits: usize,
    constant_hits: usize,
    linked_fragments: usize,
    junction_support_reads: usize,
    junction_spanning_reads: usize,
    junction_conflicting_reads: usize,
    junction_refined_calls: usize,
    junction_refined_bases: usize,
    rescued_constants: usize,
}

impl VdjRunStatus {
    fn new(threads: usize) -> Self {
        Self {
            started_unix_ms: unix_ms(),
            finished_unix_ms: None,
            stage: "0/3 preparing V(D)J reference".to_string(),
            public_url: None,
            threads,
            phase: 0,
            phase_started: Instant::now(),
            phase_times: [None; 5],
            reference_segments: 0,
            bam_records: 0,
            allowed_cell_records: 0,
            receptor_overlap_records: 0,
            unmapped_candidates: 0,
            unmapped_igh_admitted: 0,
            unmapped_igh_rescued_cells: 0,
            unmapped_igh_v_mappings: 0,
            unmapped_igh_d_mappings: 0,
            unmapped_igh_j_mappings: 0,
            unmapped_igh_c_mappings: 0,
            preliminary_cells: None,
            evidence_cells: 0,
            compact_summaries: 0,
            physical_fragments: 0,
            cells_with_v: 0,
            cells_with_d: 0,
            cells_with_j: 0,
            cells_with_c: 0,
            j_constant_linked_cells: 0,
            j_constant_linked_fragments: 0,
            igh_vj_cells: 0,
            igh_dj_no_v_cells: 0,
            igh_j_no_vd_cells: 0,
            igh_intronic_c_cells: 0,
            igh_intronic_c_records: 0,
            calls_cells: 0,
            recombinations: 0,
            productive: 0,
            mean_calls_per_cell: 0.0,
            calls_by_chain: [0; 7],
            knee_thresholds: [0; 7],
            knee_evidence_cells: [0; 7],
            knee_selected_cells: [0; 7],
            rescan_records: 0,
            rescan_wanted_records: 0,
            rescan_batches: 0,
            receptor_rediscovery_hits: 0,
            constant_hits: 0,
            linked_fragments: 0,
            junction_support_reads: 0,
            junction_spanning_reads: 0,
            junction_conflicting_reads: 0,
            junction_refined_calls: 0,
            junction_refined_bases: 0,
            rescued_constants: 0,
        }
    }
}

impl VdjRunStatus {
    fn phase_time(&self, phase: usize) -> String {
        if let Some(elapsed) = self.phase_times.get(phase).and_then(|x| *x) {
            return format_duration(elapsed);
        }
        if self.phase == phase {
            return format_duration(self.phase_started.elapsed());
        }
        "pending".to_string()
    }

    fn run_elapsed(&self) -> Duration {
        let end_ms = self.finished_unix_ms.unwrap_or_else(unix_ms);
        Duration::from_millis(end_ms.saturating_sub(self.started_unix_ms) as u64)
    }
}

impl ServerContent for VdjRunStatus {
    fn server_snapshot(&self) -> ServerSnapshot {
        let memory = memory_status();
        let productive_pct = if self.recombinations == 0 {
            0.0
        } else {
            100.0 * self.productive as f64 / self.recombinations as f64
        };
        let rescan_complete = self.bam_records > 0 && self.rescan_records == self.bam_records;

        ServerSnapshot {
            title: "Lumrik V(D)J".to_string(),
            subtitle: "Evidence, receptor reconstruction, and splice-supported constant confirmation"
                .to_string(),
            started_unix_ms: self.started_unix_ms,
            finished_unix_ms: self.finished_unix_ms,
            stage: self.stage.clone(),
            public_url: self.public_url.clone(),
            sections: vec![
                StatusSection::new(
                    "Input",
                    vec![
                        StatusMetric::new("Stage time", self.phase_time(1)),
                        StatusMetric::new("BAM records", self.bam_records.to_string()),
                        StatusMetric::new(
                            "Allowed-cell records",
                            format_count_pct(self.allowed_cell_records, self.bam_records),
                        ),
                        StatusMetric::new(
                            "Receptor-overlap records",
                            self.receptor_overlap_records.to_string(),
                        ),
                        StatusMetric::new(
                            "Preliminary cells",
                            self.preliminary_cells
                                .map(|n| n.to_string())
                                .unwrap_or_else(|| "not supplied".to_string()),
                        ),
                    ],
                ),
                StatusSection::new(
                    "Receptor evidence",
                    vec![
                        StatusMetric::new("Cells with evidence", self.evidence_cells.to_string()),
                        StatusMetric::new(
                            "V / J evidence cells",
                            format!("{} / {}", self.cells_with_v, self.cells_with_j),
                        ),
                        StatusMetric::new(
                            "D evidence cells",
                            format_count_pct(self.cells_with_d, self.evidence_cells),
                        ),
                        StatusMetric::new(
                            "Retained exonic C cells",
                            format_count_pct(self.cells_with_c, self.evidence_cells),
                        ),
                        StatusMetric::new(
                            "J→C linked cells",
                            format_count_pct(self.j_constant_linked_cells, self.evidence_cells),
                        ),
                        StatusMetric::new(
                            "J→C splice-linked fragments",
                            self.j_constant_linked_fragments.to_string(),
                        ),
                        StatusMetric::new(
                            "Compact summaries / fragments",
                            format!("{} / {}", self.compact_summaries, self.physical_fragments),
                        ),
                    ],
                ),
                StatusSection::new(
                    "Unmapped IGH rescue",
                    vec![
                        StatusMetric::new(
                            "Candidates",
                            self.unmapped_candidates.to_string(),
                        ),
                        StatusMetric::new(
                            "JH-admitted records",
                            format_count_pct(
                                self.unmapped_igh_admitted,
                                self.unmapped_candidates,
                            ),
                        ),
                        StatusMetric::new(
                            "Rescued cells",
                            self.unmapped_igh_rescued_cells.to_string(),
                        ),
                        StatusMetric::new(
                            "V / D / J / C mappings",
                            format!(
                                "{} / {} / {} / {}",
                                self.unmapped_igh_v_mappings,
                                self.unmapped_igh_d_mappings,
                                self.unmapped_igh_j_mappings,
                                self.unmapped_igh_c_mappings,
                            ),
                        ),
                    ],
                ),
                StatusSection::new(
                    "IGH transcription / rearrangement state",
                    vec![
                        StatusMetric::new("V + J evidence cells", self.igh_vj_cells.to_string()),
                        StatusMetric::new("D + J / no V cells", self.igh_dj_no_v_cells.to_string()),
                        StatusMetric::new("J / no V,D cells", self.igh_j_no_vd_cells.to_string()),
                        StatusMetric::new(
                            "Intronic C cells / records",
                            format!("{} / {}", self.igh_intronic_c_cells, self.igh_intronic_c_records),
                        ),
                    ],
                ),
                StatusSection::new(
                    "Reconstruction",
                    vec![
                        StatusMetric::new("Stage time", self.phase_time(2)),
                        StatusMetric::new("IGH knee", format_knee(self.knee_thresholds[0], self.knee_selected_cells[0], self.knee_evidence_cells[0])),
                        StatusMetric::new("IGK knee", format_knee(self.knee_thresholds[1], self.knee_selected_cells[1], self.knee_evidence_cells[1])),
                        StatusMetric::new("IGL knee", format_knee(self.knee_thresholds[2], self.knee_selected_cells[2], self.knee_evidence_cells[2])),
                        StatusMetric::new("TRA knee", format_knee(self.knee_thresholds[3], self.knee_selected_cells[3], self.knee_evidence_cells[3])),
                        StatusMetric::new("TRB knee", format_knee(self.knee_thresholds[4], self.knee_selected_cells[4], self.knee_evidence_cells[4])),
                        StatusMetric::new("TRG / TRD knee", format!("{} / {}", format_knee(self.knee_thresholds[5], self.knee_selected_cells[5], self.knee_evidence_cells[5]), format_knee(self.knee_thresholds[6], self.knee_selected_cells[6], self.knee_evidence_cells[6]))),
                        StatusMetric::new("Called cells / recombinations", format!("{} / {}", self.calls_cells, self.recombinations)),
                        StatusMetric::new("Productive", format!("{} ({productive_pct:.1}%)", self.productive)),
                        StatusMetric::new("IGH / IGK / IGL", format!("{} / {} / {}", self.calls_by_chain[0], self.calls_by_chain[1], self.calls_by_chain[2])),
                        StatusMetric::new("TRA / TRB / TRG / TRD", format!("{} / {} / {} / {}", self.calls_by_chain[3], self.calls_by_chain[4], self.calls_by_chain[5], self.calls_by_chain[6])),
                    ],
                ),
                StatusSection::new(
                    "Confirmation rescan",
                    vec![
                        StatusMetric::new("Stage time", self.phase_time(3)),
                        StatusMetric::new(
                            "BAM scan progress",
                            format_count_pct(self.rescan_records, self.bam_records),
                        ),
                        StatusMetric::new(
                            "Full BAM pass",
                            if rescan_complete { "complete" } else { "in progress" },
                        ),
                        StatusMetric::new("Wanted-cell records", self.rescan_wanted_records.to_string()),
                        StatusMetric::new("Evidence batches", self.rescan_batches.to_string()),
                        StatusMetric::new("CDR3/receptor hits", self.receptor_rediscovery_hits.to_string()),
                        StatusMetric::new("Junction-support reads", self.junction_support_reads.to_string()),
                        StatusMetric::new("Full-junction spanning reads", self.junction_spanning_reads.to_string()),
                        StatusMetric::new("Junction-conflicting reads", self.junction_conflicting_reads.to_string()),
                        StatusMetric::new("Refined calls / bases", format!("{} / {}", self.junction_refined_calls, self.junction_refined_bases)),
                        StatusMetric::new("Constant-region hits", self.constant_hits.to_string()),
                        StatusMetric::new("Linked fragments", self.linked_fragments.to_string()),
                        StatusMetric::new("Constant calls rescued", self.rescued_constants.to_string()),
                    ],
                ),
                StatusSection::new(
                    "Performance",
                    vec![
                        StatusMetric::new("Total elapsed", format_duration(self.run_elapsed())),
                        StatusMetric::new("Reference preparation", self.phase_time(0)),
                        StatusMetric::new("Output writing", self.phase_time(4)),
                        StatusMetric::new("Worker threads", self.threads.to_string()),
                        StatusMetric::new("Reference segments", self.reference_segments.to_string()),
                        StatusMetric::new(
                            "Process RSS / peak",
                            format!("{:.0} / {:.0} MiB", memory.process_rss_mib, memory.process_peak_rss_mib),
                        ),
                        StatusMetric::new(
                            "System memory available",
                            format!("{:.0} MiB", memory.system_available_mib),
                        ),
                    ],
                ),
            ],
        }
    }
}

fn unix_ms() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_millis()
}

fn update_status<F>(status: &Arc<RwLock<VdjRunStatus>>, update: F)
where
    F: FnOnce(&mut VdjRunStatus),
{
    if let Ok(mut state) = status.write() {
        update(&mut state);
    }
}

fn advance_stage(status: &Arc<RwLock<VdjRunStatus>>, next_phase: usize, stage: impl Into<String>) {
    let stage = stage.into();
    update_status(status, |state| {
        if state.phase < state.phase_times.len() && state.phase_times[state.phase].is_none() {
            state.phase_times[state.phase] = Some(state.phase_started.elapsed());
        }
        state.phase = next_phase;
        state.phase_started = Instant::now();
        state.stage = stage;
    });
}

fn finish_status(status: &Arc<RwLock<VdjRunStatus>>, stage: impl Into<String>) {
    let stage = stage.into();
    update_status(status, |state| {
        if state.phase < state.phase_times.len() && state.phase_times[state.phase].is_none() {
            state.phase_times[state.phase] = Some(state.phase_started.elapsed());
        }
        state.stage = stage;
        state.finished_unix_ms = Some(unix_ms());
    });
}

fn format_duration(duration: Duration) -> String {
    let total = duration.as_secs();
    let days = total / 86_400;
    let hours = (total % 86_400) / 3_600;
    let minutes = (total % 3_600) / 60;
    let seconds = total % 60;
    if days > 0 {
        format!("{days}d {hours:02}h {minutes:02}m {seconds:02}s")
    } else if hours > 0 {
        format!("{hours}h {minutes:02}m {seconds:02}s")
    } else if minutes > 0 {
        format!("{minutes}m {seconds:02}s")
    } else {
        format!("{seconds}s")
    }
}

fn format_count_pct(count: usize, denominator: usize) -> String {
    let pct = if denominator == 0 {
        0.0
    } else {
        100.0 * count as f64 / denominator as f64
    };
    format!("{count} ({pct:.1}%)")
}

fn format_knee(threshold: usize, selected: usize, observed: usize) -> String {
    if observed == 0 {
        return "no evidence".to_string();
    }
    format!(">={threshold} reads · {selected}/{observed} cells")
}

fn find_file(dir: &Path, names: &[&str]) -> Result<PathBuf> {
    for name in names {
        let path = dir.join(name);
        if path.is_file() {
            return Ok(path);
        }
    }
    bail!("none of {} found in {}", names.join(", "), dir.display())
}

fn read_lines(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let input: Box<dyn IoRead> = if path.extension().and_then(|x| x.to_str()) == Some("gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(input)))
}

fn final_exonic_dir(path: &Path) -> Result<PathBuf> {
    if ["barcodes.tsv.gz", "barcodes.tsv"]
        .iter()
        .any(|name| path.join(name).is_file())
    {
        return Ok(path.to_path_buf());
    }

    let nested = path.join("exonic");
    if ["barcodes.tsv.gz", "barcodes.tsv"]
        .iter()
        .any(|name| nested.join(name).is_file())
    {
        return Ok(nested);
    }

    bail!(
        "{} is neither a Nelrune exonic MEX directory nor a Nelrune analysis directory containing exonic/",
        path.display()
    )
}

fn preliminary_cell_ids(exonic: Option<&Path>) -> Result<Option<HashSet<u64>>> {
    let Some(path) = exonic else {
        return Ok(None);
    };
    let dir = final_exonic_dir(path)?;
    let barcode_path = find_file(&dir, &["barcodes.tsv.gz", "barcodes.tsv"])?;
    let mut cells = HashSet::new();
    for line in read_lines(&barcode_path)?.lines() {
        let line = line?;
        let Some(barcode) = line.split('\t').next() else {
            continue;
        };
        if barcode.is_empty() {
            continue;
        }
        cells.insert(IntToStr::new(barcode.as_bytes()).into_u64());
    }
    eprintln!(
        "[nelrune-vdj] preliminary allowed-cell set: {} barcode(s) from {}",
        cells.len(),
        barcode_path.display()
    );
    Ok(Some(cells))
}

fn sync_evidence_status(
    status: &Arc<RwLock<VdjRunStatus>>,
    progress: BamIngestProgress,
    evidence: &CellEvidenceVdj,
    index: &VdjIndex,
) {
    let mut cells = 0usize;
    let mut summaries = 0usize;
    let mut fragments = 0usize;
    let mut cells_by_kind = [0usize; 4];
    let mut j_constant_linked_cells = 0usize;
    let mut j_constant_linked_fragments = 0usize;
    let mut igh_vj_cells = 0usize;
    let mut igh_dj_no_v_cells = 0usize;
    let mut igh_j_no_vd_cells = 0usize;
    let mut igh_intronic_c_cells = 0usize;
    let mut igh_intronic_c_records = 0usize;

    for (_cell_id, cell) in evidence.cells() {
        cells = cells.saturating_add(1);
        summaries = summaries.saturating_add(cell.summary_count());
        fragments = fragments.saturating_add(cell.physical_fragments());
        let linked = cell.j_constant_linked_fragments(index);
        j_constant_linked_fragments = j_constant_linked_fragments.saturating_add(linked);
        if linked > 0 {
            j_constant_linked_cells = j_constant_linked_cells.saturating_add(1);
        }

        let igh_v = cell.segment_mappings(Chain::Igh, SegmentKind::V) > 0;
        let igh_d = cell.segment_mappings(Chain::Igh, SegmentKind::D) > 0;
        let igh_j = cell.segment_mappings(Chain::Igh, SegmentKind::J) > 0;
        if igh_v && igh_j {
            igh_vj_cells = igh_vj_cells.saturating_add(1);
        } else if !igh_v && igh_d && igh_j {
            igh_dj_no_v_cells = igh_dj_no_v_cells.saturating_add(1);
        } else if !igh_v && !igh_d && igh_j {
            igh_j_no_vd_cells = igh_j_no_vd_cells.saturating_add(1);
        }
        let intronic = cell.intronic_constant_records(Chain::Igh);
        if intronic > 0 {
            igh_intronic_c_cells = igh_intronic_c_cells.saturating_add(1);
            igh_intronic_c_records = igh_intronic_c_records.saturating_add(intronic);
        }

        for (slot, kind) in [
            SegmentKind::V,
            SegmentKind::D,
            SegmentKind::J,
            SegmentKind::C,
        ]
        .into_iter()
        .enumerate()
        {
            if Chain::ALL
                .into_iter()
                .any(|chain| cell.segment_mappings(chain, kind) > 0)
            {
                cells_by_kind[slot] = cells_by_kind[slot].saturating_add(1);
            }
        }
    }

    update_status(status, |state| {
        state.bam_records = progress.bam_records;
        state.allowed_cell_records = progress.allowed_cell_records;
        state.receptor_overlap_records = progress.receptor_overlap_records;
        state.unmapped_candidates = progress.unmapped_candidates;
        state.unmapped_igh_admitted = progress.unmapped_igh_admitted;
        state.unmapped_igh_rescued_cells = progress.unmapped_igh_rescued_cells;
        state.unmapped_igh_v_mappings = progress.unmapped_igh_v_mappings;
        state.unmapped_igh_d_mappings = progress.unmapped_igh_d_mappings;
        state.unmapped_igh_j_mappings = progress.unmapped_igh_j_mappings;
        state.unmapped_igh_c_mappings = progress.unmapped_igh_c_mappings;
        state.evidence_cells = cells;
        state.compact_summaries = summaries;
        state.physical_fragments = fragments;
        state.cells_with_v = cells_by_kind[0];
        state.cells_with_d = cells_by_kind[1];
        state.cells_with_j = cells_by_kind[2];
        state.cells_with_c = cells_by_kind[3];
        state.j_constant_linked_cells = j_constant_linked_cells;
        state.j_constant_linked_fragments = j_constant_linked_fragments;
        state.igh_vj_cells = igh_vj_cells;
        state.igh_dj_no_v_cells = igh_dj_no_v_cells;
        state.igh_j_no_vd_cells = igh_j_no_vd_cells;
        state.igh_intronic_c_cells = igh_intronic_c_cells;
        state.igh_intronic_c_records = igh_intronic_c_records;
    });
}

fn chain_slot(chain: Chain) -> usize {
    match chain {
        Chain::Igh => 0,
        Chain::Igk => 1,
        Chain::Igl => 2,
        Chain::Tra => 3,
        Chain::Trb => 4,
        Chain::Trg => 5,
        Chain::Trd => 6,
    }
}

fn yaml_quote(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            _ => out.push(ch),
        }
    }
    out.push('"');
    out
}

fn duration_ms(value: Option<Duration>) -> u128 {
    value.map(|d| d.as_millis()).unwrap_or(0)
}

fn write_run_summary_yaml(path: &Path, state: &VdjRunStatus) -> Result<()> {
    let memory = memory_status();
    let mut w = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    writeln!(w, "schema: lumrik-vdj-run-summary-v1")?;
    writeln!(w, "run:")?;
    writeln!(w, "  stage: {}", yaml_quote(&state.stage))?;
    writeln!(w, "  started_unix_ms: {}", state.started_unix_ms)?;
    match state.finished_unix_ms {
        Some(value) => writeln!(w, "  finished_unix_ms: {value}")?,
        None => writeln!(w, "  finished_unix_ms: null")?,
    }
    writeln!(w, "  elapsed_ms: {}", state.run_elapsed().as_millis())?;
    writeln!(w, "  threads: {}", state.threads)?;
    writeln!(w, "  reference_segments: {}", state.reference_segments)?;

    writeln!(w, "input:")?;
    writeln!(w, "  bam_records: {}", state.bam_records)?;
    writeln!(w, "  allowed_cell_records: {}", state.allowed_cell_records)?;
    writeln!(w, "  receptor_overlap_records: {}", state.receptor_overlap_records)?;
    match state.preliminary_cells {
        Some(value) => writeln!(w, "  preliminary_cells: {value}")?,
        None => writeln!(w, "  preliminary_cells: null")?,
    }

    writeln!(w, "unmapped_igh_rescue:")?;
    writeln!(w, "  candidates: {}", state.unmapped_candidates)?;
    writeln!(w, "  admitted_records: {}", state.unmapped_igh_admitted)?;
    writeln!(w, "  rescued_cells: {}", state.unmapped_igh_rescued_cells)?;
    writeln!(w, "  v_mappings: {}", state.unmapped_igh_v_mappings)?;
    writeln!(w, "  d_mappings: {}", state.unmapped_igh_d_mappings)?;
    writeln!(w, "  j_mappings: {}", state.unmapped_igh_j_mappings)?;
    writeln!(w, "  c_mappings: {}", state.unmapped_igh_c_mappings)?;

    writeln!(w, "evidence:")?;
    writeln!(w, "  cells: {}", state.evidence_cells)?;
    writeln!(w, "  compact_summaries: {}", state.compact_summaries)?;
    writeln!(w, "  physical_fragments: {}", state.physical_fragments)?;
    writeln!(w, "  cells_with_v: {}", state.cells_with_v)?;
    writeln!(w, "  cells_with_d: {}", state.cells_with_d)?;
    writeln!(w, "  cells_with_j: {}", state.cells_with_j)?;
    writeln!(w, "  cells_with_retained_exonic_c: {}", state.cells_with_c)?;
    writeln!(w, "  j_constant_linked_cells: {}", state.j_constant_linked_cells)?;
    writeln!(w, "  j_constant_linked_fragments: {}", state.j_constant_linked_fragments)?;

    writeln!(w, "igh_state:")?;
    writeln!(w, "  v_j_cells: {}", state.igh_vj_cells)?;
    writeln!(w, "  d_j_no_v_cells: {}", state.igh_dj_no_v_cells)?;
    writeln!(w, "  j_no_v_d_cells: {}", state.igh_j_no_vd_cells)?;
    writeln!(w, "  intronic_c_cells: {}", state.igh_intronic_c_cells)?;
    writeln!(w, "  intronic_c_records: {}", state.igh_intronic_c_records)?;

    writeln!(w, "knees:")?;
    for (slot, chain) in Chain::ALL.into_iter().enumerate() {
        writeln!(w, "  {}:", chain.to_string())?;
        writeln!(w, "    threshold_records: {}", state.knee_thresholds[slot])?;
        writeln!(w, "    evidence_cells: {}", state.knee_evidence_cells[slot])?;
        writeln!(w, "    selected_cells: {}", state.knee_selected_cells[slot])?;
    }

    writeln!(w, "calls:")?;
    writeln!(w, "  called_cells: {}", state.calls_cells)?;
    writeln!(w, "  recombinations: {}", state.recombinations)?;
    writeln!(w, "  productive: {}", state.productive)?;
    writeln!(w, "  mean_per_called_cell: {:.6}", state.mean_calls_per_cell)?;
    for (slot, chain) in Chain::ALL.into_iter().enumerate() {
        writeln!(w, "  {}: {}", chain.to_string(), state.calls_by_chain[slot])?;
    }

    writeln!(w, "confirmation_rescan:")?;
    writeln!(w, "  bam_records_scanned: {}", state.rescan_records)?;
    writeln!(w, "  wanted_cell_records: {}", state.rescan_wanted_records)?;
    writeln!(w, "  evidence_batches: {}", state.rescan_batches)?;
    writeln!(w, "  receptor_rediscovery_hits: {}", state.receptor_rediscovery_hits)?;
    writeln!(w, "  constant_region_hits: {}", state.constant_hits)?;
    writeln!(w, "  linked_fragments: {}", state.linked_fragments)?;
    writeln!(w, "  junction_support_reads: {}", state.junction_support_reads)?;
    writeln!(w, "  junction_spanning_reads: {}", state.junction_spanning_reads)?;
    writeln!(w, "  junction_conflicting_reads: {}", state.junction_conflicting_reads)?;
    writeln!(w, "  junction_refined_calls: {}", state.junction_refined_calls)?;
    writeln!(w, "  junction_refined_bases: {}", state.junction_refined_bases)?;
    writeln!(w, "  constant_calls_rescued: {}", state.rescued_constants)?;
    writeln!(
        w,
        "  complete_bam_pass: {}",
        state.bam_records > 0 && state.rescan_records == state.bam_records
    )?;

    writeln!(w, "timings_ms:")?;
    writeln!(w, "  reference_preparation: {}", duration_ms(state.phase_times[0]))?;
    writeln!(w, "  evidence_collection: {}", duration_ms(state.phase_times[1]))?;
    writeln!(w, "  reconstruction: {}", duration_ms(state.phase_times[2]))?;
    writeln!(w, "  confirmation_rescan: {}", duration_ms(state.phase_times[3]))?;
    writeln!(w, "  output_writing: {}", duration_ms(state.phase_times[4]))?;
    writeln!(w, "memory_mib:")?;
    writeln!(w, "  process_rss: {:.3}", memory.process_rss_mib)?;
    writeln!(w, "  process_peak_rss: {:.3}", memory.process_peak_rss_mib)?;
    writeln!(w, "  system_available: {:.3}", memory.system_available_mib)?;
    w.flush()?;
    Ok(())
}

fn write_static_report(path: &Path, state: &VdjRunStatus) -> Result<()> {
    let html = snapshot_html(&state.server_snapshot());
    let mut w = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    w.write_all(html.as_bytes())?;
    w.flush()?;
    Ok(())
}

fn main() -> Result<()> {
    let c = Cli::parse();
    let status = Arc::new(RwLock::new(VdjRunStatus::new(c.threads.max(1))));
    let _status_server = if c.no_health_server {
        None
    } else {
        let addr = SocketAddr::new(IpAddr::V4(Ipv4Addr::UNSPECIFIED), c.health_port);
        let server = spawn_status_server(Arc::clone(&status), addr)?;
        let hostname = public_hostname(c.health_hostname.as_deref());
        let url = format!("http://{}:{}", hostname, server.addr().port());
        update_status(&status, |state| state.public_url = Some(url.clone()));
        eprintln!("[nelrune-vdj] live dashboard: {url}");
        eprintln!("[nelrune-vdj] status JSON: {url}/status   health probe: {url}/health");
        Some(server)
    };

    let mut mapping_info = MappingInfo::new(None, 0.0, 0);
    let loading_prebuilt_index = c.index.is_some();
    mapping_info.start_counter();
    mapping_info.start_timer("vdj.index_load_or_build");
    let index = match (&c.index, &c.gtf, &c.genome) {
        (Some(p), None, None) => {
            VdjIndex::load(p).with_context(|| format!("loading {}", p.display()))?
        }
        (None, Some(g), Some(f)) => VdjIndexBuilder::default().build(g, f)?,
        _ => bail!("provide either --index or --gtf + --genome"),
    };
    mapping_info.stop_timer("vdj.index_load_or_build");
    if loading_prebuilt_index && !index.has_precise_exon_blocks() {
        bail!(
            "VDJ index uses legacy span-only coordinates; rebuild it with the current vdj-index so constant-gene exons and introns can be distinguished"
        );
    }
    if loading_prebuilt_index {
        mapping_info.stop_file_io_time();
    } else {
        mapping_info.stop_single_processor_time();
    }
    update_status(&status, |state| state.reference_segments = index.len());

    let preliminary_cells = preliminary_cell_ids(c.exonic.as_deref())?;
    if preliminary_cells.is_none() {
        eprintln!(
            "[nelrune-vdj] no --exonic cell gate supplied; collecting receptor evidence from all BAM cell IDs"
        );
    }
    update_status(&status, |state| {
        state.preliminary_cells = preliminary_cells.as_ref().map(HashSet::len);
    });

    let mut runner = VdjRunner::new(
        index,
        VdjRunnerConfig {
            min_sequence_overlap: c.min_sequence_overlap,
        },
    );
    runner.set_threads(c.threads);

    advance_stage(&status, 1, "1/3 initial V(D)J evidence collection");
    mapping_info.start_counter();
    mapping_info.start_timer("vdj.bam_read");
    let status_for_read = Arc::clone(&status);
    let n = runner.read_bam_with_progress_for_cells(
        &c.bam,
        &NelruneIdentityResolver,
        preliminary_cells.as_ref(),
        move |progress, evidence, index| {
            sync_evidence_status(&status_for_read, progress, evidence, index);
        },
    )?;
    // The final callback already contains the exact BAM denominator counters.
    // `n` remains the public return value: receptor-overlap records ingested.
    let _ = n;
    mapping_info.stop_timer("vdj.bam_read");
    mapping_info.stop_file_io_time();

    let knee_selection = runner.receptor_knee_selection();
    update_status(&status, |state| {
        for selection in &knee_selection {
            let slot = chain_slot(selection.chain);
            state.knee_thresholds[slot] = selection.threshold_records;
            state.knee_evidence_cells[slot] = selection.evidence_cells;
            state.knee_selected_cells[slot] = selection.selected_cells;
        }
    });

    advance_stage(&status, 2, "2/3 reconstructing cell-specific receptors");
    mapping_info.start_counter();
    mapping_info.start_timer("vdj.recombination_calling");
    let mut calls = runner.identify_with_receptor_knee_selection(&knee_selection);
    for selection in &knee_selection {
        let prefix = format!(
            "vdj.knee.{}",
            selection.chain.to_string().to_ascii_lowercase()
        );
        mapping_info.report_n(
            format!("{prefix}.threshold_records"),
            selection.threshold_records,
        );
        mapping_info.report_n(format!("{prefix}.evidence_cells"), selection.evidence_cells);
        mapping_info.report_n(format!("{prefix}.selected_cells"), selection.selected_cells);
        eprintln!(
            "[nelrune-vdj] {} V/J knee: >= {} reads; {}/{} evidence cell(s) selected",
            selection.chain,
            selection.threshold_records,
            selection.selected_cells,
            selection.evidence_cells
        );
    }
    mapping_info.stop_timer("vdj.recombination_calling");
    mapping_info.stop_multi_processor_time();

    let mut nr = 0usize;
    let mut productive = 0usize;
    let mut calls_by_chain = [0usize; 7];
    let mut calls_cells = 0usize;
    for (_, cell_calls) in &calls {
        if !cell_calls.is_empty() {
            calls_cells += 1;
        }
        for call in cell_calls {
            nr += 1;
            if call.productive {
                productive += 1;
            }
            calls_by_chain[chain_slot(call.chain)] += 1;
        }
    }
    update_status(&status, |state| {
        state.calls_cells = calls_cells;
        state.recombinations = nr;
        state.productive = productive;
        state.calls_by_chain = calls_by_chain;
        state.mean_calls_per_cell = nr as f64 / calls_cells.max(1) as f64;
    });

    advance_stage(
        &status,
        3,
        "3/3 remapping CDR3s and confirming constant regions",
    );
    mapping_info.start_counter();
    mapping_info.start_timer("vdj.receptor_rediscovery");
    let status_for_rescan = Arc::clone(&status);
    let rescan = runner.rediscover_receptor_linkage_from_bam_with_report_and_progress(
        &c.bam,
        &NelruneIdentityResolver,
        &mut calls,
        move |progress| {
            update_status(&status_for_rescan, |state| {
                state.rescan_records = progress.bam_records_scanned;
                state.rescan_wanted_records = progress.wanted_cell_records;
                state.rescan_batches = progress.batches_completed;
                state.receptor_rediscovery_hits = progress.receptor_hit_records;
                state.constant_hits = progress.constant_hit_records;
                state.linked_fragments = progress.linked_fragments;
                state.junction_support_reads = progress.junction_support_reads;
                state.junction_spanning_reads = progress.junction_spanning_reads;
                state.junction_conflicting_reads = progress.junction_conflicting_reads;
            });
        },
    )?;
    mapping_info.stop_timer("vdj.receptor_rediscovery");
    mapping_info.stop_multi_processor_time();
    update_status(&status, |state| {
        state.rescan_records = rescan.bam_records_scanned;
        state.rescan_wanted_records = rescan.wanted_cell_records;
        state.receptor_rediscovery_hits = rescan.receptor_hit_records;
        state.constant_hits = rescan.constant_hit_records;
        state.linked_fragments = rescan.linked_fragments;
        state.junction_support_reads = rescan.junction_support_reads;
        state.junction_spanning_reads = rescan.junction_spanning_reads;
        state.junction_conflicting_reads = rescan.junction_conflicting_reads;
        state.junction_refined_calls = rescan.junction_refined_calls;
        state.junction_refined_bases = rescan.junction_refined_bases;
        state.rescued_constants = rescan.rescued;
    });

    advance_stage(&status, 4, "writing V(D)J summaries");
    mapping_info.start_counter();
    mapping_info.start_timer("vdj.output_write");
    let mut writer = ReportWriter::create(&c.out, c.write_sequences)?;
    if let Some(raw) = c.bd_cell_version.as_deref() {
        let version = BdCellVersion::parse(raw)
            .map_err(|e| anyhow::anyhow!("invalid --bd-cell-version {raw}: {e}"))?;
        writer = writer.with_bd_cell_version(version);
    }
    let mut by_chain = HashMap::<String, usize>::new();
    for (cell_id, rs) in &calls {
        let name = runner
            .cell_names
            .get(cell_id)
            .map(String::as_str)
            .unwrap_or("unknown");
        writer.write_cell(name, rs, &runner.index)?;
        for r in rs {
            *by_chain.entry(r.chain.to_string()).or_default() += 1;
        }
    }
    writer.finish()?;
    mapping_info.stop_timer("vdj.output_write");
    mapping_info.stop_file_io_time();

    mapping_info.total = n;
    mapping_info.report_n("vdj.threads", c.threads.max(1));
    mapping_info.report_n("vdj.receptor_overlap_records", n);
    let ingest_status = status
        .read()
        .expect("reading VDJ status after BAM ingestion")
        .clone();
    mapping_info.report_n("vdj.unmapped_candidates", ingest_status.unmapped_candidates);
    mapping_info.report_n("vdj.unmapped_igh_admitted", ingest_status.unmapped_igh_admitted);
    mapping_info.report_n("vdj.unmapped_igh_rescued_cells", ingest_status.unmapped_igh_rescued_cells);
    mapping_info.report_n("vdj.unmapped_igh_v_mappings", ingest_status.unmapped_igh_v_mappings);
    mapping_info.report_n("vdj.unmapped_igh_d_mappings", ingest_status.unmapped_igh_d_mappings);
    mapping_info.report_n("vdj.unmapped_igh_j_mappings", ingest_status.unmapped_igh_j_mappings);
    mapping_info.report_n("vdj.unmapped_igh_c_mappings", ingest_status.unmapped_igh_c_mappings);
    if let Some(cells) = preliminary_cells.as_ref() {
        mapping_info.report_n("vdj.preliminary_cells", cells.len());
    }
    mapping_info.report_n("vdj.cells_with_evidence", calls.len());
    mapping_info.report_n("vdj.recombinations", nr);
    mapping_info.report_n("vdj.receptor_rediscovery_constant_calls", rescan.rescued);
    mapping_info.report_n("vdj.junction_support_reads", rescan.junction_support_reads);
    mapping_info.report_n("vdj.junction_spanning_reads", rescan.junction_spanning_reads);
    mapping_info.report_n("vdj.junction_conflicting_reads", rescan.junction_conflicting_reads);
    mapping_info.report_n("vdj.junction_refined_calls", rescan.junction_refined_calls);
    mapping_info.report_n("vdj.junction_refined_bases", rescan.junction_refined_bases);
    for (chain, count) in &by_chain {
        mapping_info.report_n(format!("vdj.calls.{}", chain.to_ascii_lowercase()), *count);
    }
    write_mapping_info_report(c.out.join("vdj-mapping-info.txt"), &mapping_info)?;

    finish_status(&status, "finished — receptors reconstructed");

    let final_status = status
        .read()
        .map_err(|_| anyhow::anyhow!("VDJ run-status lock poisoned"))?
        .clone();
    write_run_summary_yaml(&c.out.join("vdj-run-summary.yaml"), &final_status)?;
    write_static_report(&c.out.join("vdj-report.html"), &final_status)?;

    eprintln!("{mapping_info}");
    eprintln!(
        "nelrune-vdj: {n} receptor-overlapping BAM records; {} evidence cell(s); {nr} recombination(s); {} constant call(s) added by CDR3/constant remapping; outputs in {}",
        runner.evidence.cell_count(),
        rescan.rescued,
        c.out.display()
    );
    Ok(())
}
