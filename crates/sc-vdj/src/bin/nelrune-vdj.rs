use anyhow::{bail, Context, Result};
use clap::Parser;
use lumrik_status::{
    memory_status, public_hostname, spawn_status_server, ServerContent, ServerSnapshot,
    StatusMetric, StatusSection,
};
use mapping_info::MappingInfo;
use sc_vdj::output::{write_mapping_info_report, ReportWriter};
use sc_vdj::{
    CellEvidenceVdj, Chain, NelruneIdentityResolver, SegmentKind, VdjIndex, VdjIndexBuilder,
    VdjRunner, VdjRunnerConfig,
};
use std::collections::{HashMap, HashSet};
use std::net::{IpAddr, Ipv4Addr, SocketAddr};
use std::path::PathBuf;
use std::sync::{Arc, RwLock};
use std::time::{SystemTime, UNIX_EPOCH};

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
    #[arg(long)]
    exonic: Option<PathBuf>,
    #[arg(long)]
    write_sequences: bool,
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
    reference_segments: usize,

    receptor_overlap_records: usize,
    evidence_cells: usize,
    compact_summaries: usize,
    physical_fragments: usize,
    mean_v: f64,
    mean_d: f64,
    mean_j: f64,
    mean_c: f64,

    calls_cells: usize,
    recombinations: usize,
    productive: usize,
    mean_calls_per_cell: f64,
    calls_by_chain: [usize; 7],

    rescan_records: usize,
    rescan_wanted_records: usize,
    rescan_batches: usize,
    receptor_rediscovery_hits: usize,
    constant_hits: usize,
    linked_fragments: usize,
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
            reference_segments: 0,
            receptor_overlap_records: 0,
            evidence_cells: 0,
            compact_summaries: 0,
            physical_fragments: 0,
            mean_v: 0.0,
            mean_d: 0.0,
            mean_j: 0.0,
            mean_c: 0.0,
            calls_cells: 0,
            recombinations: 0,
            productive: 0,
            mean_calls_per_cell: 0.0,
            calls_by_chain: [0; 7],
            rescan_records: 0,
            rescan_wanted_records: 0,
            rescan_batches: 0,
            receptor_rediscovery_hits: 0,
            constant_hits: 0,
            linked_fragments: 0,
            rescued_constants: 0,
        }
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
        ServerSnapshot {
            title: "Lumrik V(D)J".to_string(),
            subtitle: "Live receptor archaeology — evidence → receptors → CDR3 confirmation"
                .to_string(),
            started_unix_ms: self.started_unix_ms,
            finished_unix_ms: self.finished_unix_ms,
            stage: self.stage.clone(),
            public_url: self.public_url.clone(),
            sections: vec![
                StatusSection::new(
                    "1 · Initial evidence collection",
                    vec![
                        StatusMetric::new(
                            "Receptor-overlap records",
                            self.receptor_overlap_records.to_string(),
                        ),
                        StatusMetric::new("Distinct BAM cell IDs with receptor evidence", self.evidence_cells.to_string()),
                        StatusMetric::new("Persistent compact receptor summaries", self.compact_summaries.to_string()),
                        StatusMetric::new("Receptor-evidence fragments", self.physical_fragments.to_string()),
                        StatusMetric::new("Mean distinct V segments / evidence cell ID", format!("{:.2}", self.mean_v)),
                        StatusMetric::new("Mean distinct D segments / evidence cell ID", format!("{:.2}", self.mean_d)),
                        StatusMetric::new("Mean distinct J segments / evidence cell ID", format!("{:.2}", self.mean_j)),
                        StatusMetric::new("Mean distinct C segments / evidence cell ID", format!("{:.2}", self.mean_c)),
                    ],
                ),
                StatusSection::new(
                    "2 · Receptor reconstruction",
                    vec![
                        StatusMetric::new("Cells with reconstructed receptor", self.calls_cells.to_string()),
                        StatusMetric::new("Recombinations", self.recombinations.to_string()),
                        StatusMetric::new("Mean recombinations / called cell", format!("{:.2}", self.mean_calls_per_cell)),
                        StatusMetric::new("Productive", format!("{} ({productive_pct:.1}%)", self.productive)),
                        StatusMetric::new("IGH", self.calls_by_chain[0].to_string()),
                        StatusMetric::new("IGK", self.calls_by_chain[1].to_string()),
                        StatusMetric::new("IGL", self.calls_by_chain[2].to_string()),
                        StatusMetric::new("TRA", self.calls_by_chain[3].to_string()),
                        StatusMetric::new("TRB", self.calls_by_chain[4].to_string()),
                        StatusMetric::new("TRG / TRD", format!("{} / {}", self.calls_by_chain[5], self.calls_by_chain[6])),
                    ],
                ),
                StatusSection::new(
                    "3 · CDR3 / constant confirmation",
                    vec![
                        StatusMetric::new("BAM records rescanned", self.rescan_records.to_string()),
                        StatusMetric::new("Wanted-cell records", self.rescan_wanted_records.to_string()),
                        StatusMetric::new("Batches completed", self.rescan_batches.to_string()),
                        StatusMetric::new("CDR3/receptor rediscovery hits", self.receptor_rediscovery_hits.to_string()),
                        StatusMetric::new("Constant-region hits", self.constant_hits.to_string()),
                        StatusMetric::new("Linked fragments", self.linked_fragments.to_string()),
                        StatusMetric::new("Constant calls rescued", self.rescued_constants.to_string()),
                    ],
                ),
                StatusSection::new(
                    "Execution",
                    vec![
                        StatusMetric::new("Worker threads", self.threads.to_string()),
                        StatusMetric::new("Reference segments", self.reference_segments.to_string()),
                        StatusMetric::new(
                            "Process RSS / peak",
                            format!(
                                "{:.0} / {:.0} MiB",
                                memory.process_rss_mib, memory.process_peak_rss_mib
                            ),
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

fn set_stage(status: &Arc<RwLock<VdjRunStatus>>, stage: impl Into<String>) {
    let stage = stage.into();
    update_status(status, |state| state.stage = stage);
}

fn sync_evidence_status(
    status: &Arc<RwLock<VdjRunStatus>>,
    records: usize,
    evidence: &CellEvidenceVdj,
    index: &VdjIndex,
) {
    let cells = evidence.cell_count();
    let mut summaries = 0usize;
    let mut fragments = 0usize;
    let mut unique_elements = [0usize; 4];

    for (_, cell) in evidence.cells() {
        summaries = summaries.saturating_add(cell.summary_count());
        fragments = fragments.saturating_add(cell.physical_fragments());

        let mut by_kind: [HashSet<u16>; 4] = std::array::from_fn(|_| HashSet::new());
        for chain in cell.chains(index) {
            for summary in cell.summaries_for_chain(chain) {
                for segment_id in summary.segment_ids() {
                    let Some(segment) = index.segment(segment_id) else {
                        continue;
                    };
                    let slot = match segment.kind {
                        SegmentKind::V => 0,
                        SegmentKind::D => 1,
                        SegmentKind::J => 2,
                        SegmentKind::C => 3,
                    };
                    by_kind[slot].insert(segment_id);
                }
            }
        }
        for (slot, ids) in by_kind.into_iter().enumerate() {
            unique_elements[slot] = unique_elements[slot].saturating_add(ids.len());
        }
    }

    let denom = cells.max(1) as f64;
    update_status(status, |state| {
        state.receptor_overlap_records = records;
        state.evidence_cells = cells;
        state.compact_summaries = summaries;
        state.physical_fragments = fragments;
        state.mean_v = unique_elements[0] as f64 / denom;
        state.mean_d = unique_elements[1] as f64 / denom;
        state.mean_j = unique_elements[2] as f64 / denom;
        state.mean_c = unique_elements[3] as f64 / denom;
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
    if loading_prebuilt_index {
        mapping_info.stop_file_io_time();
    } else {
        mapping_info.stop_single_processor_time();
    }
    update_status(&status, |state| state.reference_segments = index.len());

    let mut runner = VdjRunner::new(
        index,
        VdjRunnerConfig {
            min_sequence_overlap: c.min_sequence_overlap,
        },
    );
    runner.set_threads(c.threads);

    set_stage(&status, "1/3 initial V(D)J evidence collection");
    mapping_info.start_counter();
    mapping_info.start_timer("vdj.bam_read");
    let status_for_read = Arc::clone(&status);
    let n = runner.read_bam_with_progress(&c.bam, &NelruneIdentityResolver, move |records, evidence, index| {
        sync_evidence_status(&status_for_read, records, evidence, index);
    })?;
    sync_evidence_status(&status, n, &runner.evidence, &runner.index);
    mapping_info.stop_timer("vdj.bam_read");
    mapping_info.stop_file_io_time();

    set_stage(&status, "2/3 reconstructing cell-specific receptors");
    mapping_info.start_counter();
    mapping_info.start_timer("vdj.recombination_calling");
    let mut calls = runner.identify();
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

    set_stage(&status, "3/3 remapping CDR3s and confirming constant regions");
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
        state.rescued_constants = rescan.rescued;
    });

    set_stage(&status, "writing V(D)J summaries");
    mapping_info.start_counter();
    mapping_info.start_timer("vdj.output_write");
    let mut writer = ReportWriter::create(&c.out, c.write_sequences)?;
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
    mapping_info.report_n("vdj.cells_with_evidence", calls.len());
    mapping_info.report_n("vdj.recombinations", nr);
    mapping_info.report_n("vdj.receptor_rediscovery_constant_calls", rescan.rescued);
    for (chain, count) in &by_chain {
        mapping_info.report_n(format!("vdj.calls.{}", chain.to_ascii_lowercase()), *count);
    }
    write_mapping_info_report(c.out.join("vdj-mapping-info.txt"), &mapping_info)?;

    update_status(&status, |state| {
        state.stage = "finished — receptors reconstructed".to_string();
        state.finished_unix_ms = Some(unix_ms());
    });

    eprintln!("{mapping_info}");
    let _ = &c.exonic;
    eprintln!(
        "nelrune-vdj: {n} receptor-overlapping BAM records; {} evidence cell(s); {nr} recombination(s); {} constant call(s) added by CDR3/constant remapping; outputs in {}",
        runner.evidence.cell_count(),
        rescan.rescued,
        c.out.display()
    );
    Ok(())
}
