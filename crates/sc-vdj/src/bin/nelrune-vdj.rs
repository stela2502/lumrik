use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use lumrik_status::{
    memory_status, public_hostname, spawn_status_server, ServerContent, ServerSnapshot, StatusMetric,
    StatusSection,
};
use mapping_info::MappingInfo;
use sc_vdj::{
    read_bam_receptor_evidence_with_progress, Chain, ExpressionMatrix, GermlineSegmentSupport,
    NelruneIdentityResolver, PosteriorAnalyzer, PosteriorConfig,
    RecombinationMeasurements, VdjMapper, VdjMapperConfig, VdjReference, VdjReferenceBuilder,
};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read as IoRead, Write};
use std::net::{IpAddr, Ipv4Addr, SocketAddr};
use std::path::{Path, PathBuf};
use std::sync::{Arc, RwLock};
use std::time::{SystemTime, UNIX_EPOCH};

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    name = "nelrune-vdj",
    about = "Scalable single-cell V(D)J reconstruction and structural clonotyping",
    override_usage = "nelrune-vdj --exonic <MEX_DIR> --bam <BAM> (--index <VDJIDX> | --gtf <GTF> --genome <FASTA>) --out <DIR> [OPTIONS]",
    after_help = "OUTPUTS
  vdj_calls.tsv         one biological rearrangement per row, including sequence artifacts
  vdj_receptors.tsv     one human-readable cell-level HC/LC receptor summary per row
  vdj-mapping-info.txt  performance counters and timings
  vdj_observed.fasta    optional (--write-sequences)
  vdj_naive.fasta       optional (--write-sequences)

LEGACY BD
  Older BD MEX barcodes may contain u64 A-padding. Use
  --cell-barcode-len 27 to match the original 27-base BAM cell barcode.

LIVE STATUS
  By default nelrune-vdj serves the Lumrik live dashboard on --health-port 8787.
  Use --no-health-server for batch environments that do not permit listening sockets."
)]
struct Cli {
    /// Nelrune exonic MEX directory.
    #[arg(long, value_name = "MEX_DIR")]
    exonic: PathBuf,

    /// Nelrune BAM with cell/UMI identity.
    #[arg(long, value_name = "BAM")]
    bam: PathBuf,

    /// Compiled V(D)J index (preferred reference input).
    #[arg(long, value_name = "VDJIDX")]
    index: Option<PathBuf>,

    /// Matching genome annotation; requires --genome when --index is absent.
    #[arg(long, value_name = "GTF")]
    gtf: Option<PathBuf>,

    /// Matching genome FASTA; requires --gtf when --index is absent.
    #[arg(long, value_name = "FASTA")]
    genome: Option<PathBuf>,

    /// Output directory.
    #[arg(long, value_name = "DIR")]
    out: PathBuf,

    /// Receptor-locus spatial bins.
    #[arg(long, default_value_t = 64, value_name = "N")]
    sterile_bins: usize,

    /// Normalize MEX barcodes to the first N bases (legacy BD: 27).
    #[arg(long, value_name = "N")]
    cell_barcode_len: Option<usize>,

    /// Also write observed and inferred-naive rearrangement FASTA files.
    #[arg(long)]
    write_sequences: bool,

    /// Port for the live Lumrik status server.
    #[arg(long, default_value_t = 8787)]
    health_port: u16,

    /// Hostname shown in the live-server URL. Useful on clusters.
    #[arg(long)]
    health_hostname: Option<String>,

    /// Disable the live Lumrik status server.
    #[arg(long, default_value_t = false)]
    no_health_server: bool,
}

impl Cli {
    fn validate(mut self) -> Result<Self> {
        if self.index.is_some() && (self.gtf.is_some() || self.genome.is_some()) {
            bail!("use either --index <VDJIDX> or --gtf <GTF> --genome <FASTA>, not both");
        }
        if self.index.is_none() && (self.gtf.is_none() || self.genome.is_none()) {
            bail!("missing reference: provide --index <VDJIDX> or both --gtf <GTF> and --genome <FASTA>");
        }
        if let Some(len) = self.cell_barcode_len {
            if !(1..=32).contains(&len) {
                bail!("--cell-barcode-len must be between 1 and 32");
            }
        }
        self.sterile_bins = self.sterile_bins.max(8);
        Ok(self)
    }
}

#[derive(Debug, Clone)]
struct VdjRunStatus {
    started_unix_ms: u128,
    finished_unix_ms: Option<u128>,
    stage: String,
    public_url: Option<String>,
    cells: usize,
    reference_segments: usize,
    receptor_chromosomes: usize,
    indexed_segment_intervals: usize,
    total_bam_records: usize,
    receptor_chromosome_records: usize,
    segment_overlap_records: usize,
    called_cell_records: usize,
    routed_evidence_records: usize,
    non_receptor_chromosome_records: usize,
    non_segment_overlap_records: usize,
    locus_records: [usize; 7],
    bam_scan_started_unix_ms: Option<u128>,
    bam_scan_finished_unix_ms: Option<u128>,
    reconstruction_started_unix_ms: Option<u128>,
    summaries_written: usize,
    umis: usize,
    umi_contigs: usize,
    seeded_contigs: usize,
    initial_local_alignments: usize,
    rearrangements: usize,
    heavy_ids: usize,
    light_ids: usize,
}

impl VdjRunStatus {
    fn new() -> Self {
        Self {
            started_unix_ms: unix_ms(),
            finished_unix_ms: None,
            stage: "startup".to_string(),
            public_url: None,
            cells: 0,
            reference_segments: 0,
            receptor_chromosomes: 0,
            indexed_segment_intervals: 0,
            total_bam_records: 0,
            receptor_chromosome_records: 0,
            segment_overlap_records: 0,
            called_cell_records: 0,
            routed_evidence_records: 0,
            non_receptor_chromosome_records: 0,
            non_segment_overlap_records: 0,
            locus_records: [0; 7],
            bam_scan_started_unix_ms: None,
            bam_scan_finished_unix_ms: None,
            reconstruction_started_unix_ms: None,
            summaries_written: 0,
            umis: 0,
            umi_contigs: 0,
            seeded_contigs: 0,
            initial_local_alignments: 0,
            rearrangements: 0,
            heavy_ids: 0,
            light_ids: 0,
        }
    }
}

impl ServerContent for VdjRunStatus {
    fn server_snapshot(&self) -> ServerSnapshot {
        let memory = memory_status();
        let now = unix_ms();
        let bam_rate = rate_per_second(self.total_bam_records, self.bam_scan_started_unix_ms, self.bam_scan_finished_unix_ms.unwrap_or(now));
        let cell_rate = rate_per_second(self.summaries_written, self.reconstruction_started_unix_ms, self.finished_unix_ms.unwrap_or(now));
        let segment_pct = if self.receptor_chromosome_records == 0 { 0.0 } else { 100.0 * self.segment_overlap_records as f64 / self.receptor_chromosome_records as f64 };
        ServerSnapshot {
            title: "Lumrik V(D)J".to_string(),
            subtitle: "Live receptor reconstruction internals".to_string(),
            started_unix_ms: self.started_unix_ms,
            finished_unix_ms: self.finished_unix_ms,
            stage: self.stage.clone(),
            public_url: self.public_url.clone(),
            sections: vec![
                StatusSection::new("BAM routing", vec![
                    StatusMetric::new("Records scanned", self.total_bam_records.to_string()),
                    StatusMetric::new("Records / sec", format!("{bam_rate:.0}")),
                    StatusMetric::new("On receptor chromosomes", self.receptor_chromosome_records.to_string()),
                    StatusMetric::new("Overlapping indexed segments", self.segment_overlap_records.to_string()),
                    StatusMetric::new("Called-cell VDJ records", self.called_cell_records.to_string()),
                    StatusMetric::new("Routed locus observations", self.routed_evidence_records.to_string()),
                    StatusMetric::new("Segment overlap / receptor chr", format!("{segment_pct:.2}%")),
                ]),
                StatusSection::new("Reconstruction", vec![
                    StatusMetric::new("Called cells", self.cells.to_string()),
                    StatusMetric::new("Cells completed", format!("{} / {}", self.summaries_written, self.cells)),
                    StatusMetric::new("Cells / sec", format!("{cell_rate:.2}")),
                    StatusMetric::new("UMIs examined", self.umis.to_string()),
                    StatusMetric::new("UMI contigs", self.umi_contigs.to_string()),
                    StatusMetric::new("Seeded contigs", self.seeded_contigs.to_string()),
                    StatusMetric::new("Local alignments", self.initial_local_alignments.to_string()),
                    StatusMetric::new("Rearrangements", self.rearrangements.to_string()),
                    StatusMetric::new("Cells with HC ID", self.heavy_ids.to_string()),
                    StatusMetric::new("Cells with LC ID", self.light_ids.to_string()),
                ]),
                StatusSection::new("Evidence by locus", vec![
                    StatusMetric::new("IGH", self.locus_records[0].to_string()),
                    StatusMetric::new("IGK", self.locus_records[1].to_string()),
                    StatusMetric::new("IGL", self.locus_records[2].to_string()),
                    StatusMetric::new("TRA", self.locus_records[3].to_string()),
                    StatusMetric::new("TRB", self.locus_records[4].to_string()),
                    StatusMetric::new("TRG", self.locus_records[5].to_string()),
                    StatusMetric::new("TRD", self.locus_records[6].to_string()),
                ]),
                StatusSection::new("Execution", vec![
                    StatusMetric::new("Reference segments", self.reference_segments.to_string()),
                    StatusMetric::new("Receptor chromosomes", self.receptor_chromosomes.to_string()),
                    StatusMetric::new("Indexed segment intervals", self.indexed_segment_intervals.to_string()),
                    StatusMetric::new("Non-receptor chromosome", self.non_receptor_chromosome_records.to_string()),
                    StatusMetric::new("Receptor chr, no segment", self.non_segment_overlap_records.to_string()),
                    StatusMetric::new("Process RSS / peak", format!("{:.0} / {:.0} MiB", memory.process_rss_mib, memory.process_peak_rss_mib)),
                    StatusMetric::new("System memory available", format!("{:.0} MiB", memory.system_available_mib)),
                ]),
            ],
        }
    }
}

fn rate_per_second(count: usize, started: Option<u128>, ended: u128) -> f64 {
    let Some(started) = started else { return 0.0; };
    let elapsed_ms = ended.saturating_sub(started);
    if elapsed_ms == 0 { 0.0 } else { count as f64 * 1000.0 / elapsed_ms as f64 }
}

fn sync_reconstruction_status(status: &Arc<RwLock<VdjRunStatus>>, mapping_info: &MappingInfo) {
    let counter = |name: &str| mapping_info.reads_log.get(name).copied().unwrap_or(0);
    update_status(status, |state| {
        state.umis = counter("vdj.umis");
        state.umi_contigs = counter("vdj.umi_contigs");
        state.seeded_contigs = counter("vdj.seeded_contigs");
        state.initial_local_alignments = counter("vdj.initial_local_alignments");
    });
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
    update_status(status, |state| state.stage = stage.into());
}

#[derive(Debug, Default)]
struct MexExpression {
    cells: Vec<String>,
    values: HashMap<(String, String), f64>,
}

impl ExpressionMatrix for MexExpression {
    fn expression(&self, cell: &str, gene: &str) -> f64 {
        *self
            .values
            .get(&(cell.to_string(), gene.to_ascii_uppercase()))
            .unwrap_or(&0.0)
    }
}

struct BatchWriter {
    calls: BufWriter<File>,
    receptors: BufWriter<File>,
    observed: Option<BufWriter<File>>,
    naive: Option<BufWriter<File>>,
}

impl BatchWriter {
    fn create(dir: &Path, write_sequences: bool) -> Result<Self> {
        fs::create_dir_all(dir)?;
        let mut calls = BufWriter::new(File::create(dir.join("vdj_calls.tsv"))?);
        writeln!(
            calls,
            "cell\trecombination_id\tchain\tstage\tv\td\tj\tc\tconstant_type\tsupport_umis\tv_score\tv_reads\tv_umis\tv_locus_fraction\tv_distance_bp\td_score\td_reads\td_umis\td_inferred_from_vj_junction\td_hypothesis_margin\tj_score\tj_reads\tj_umis\tc_score\tc_reads\tc_umis\tv_del_3\tp_v3_len\tp_v3\tn1_len\tn1\tp_d5_len\tp_d5\td_del_5\td_retained_len\td_del_3\tp_d3_len\tp_d3\tn2_len\tn2\tp_j5_len\tp_j5\tj_del_5\tpn_alternative\tgermline_v\tgermline_d\tgermline_j\tgermline_c\tobserved_v\tobserved_d\tobserved_j\tnaive_v\tnaive_d\tnaive_j\tobserved_rearrangement\tnaive_recombination"
        )?;
        let mut receptors = BufWriter::new(File::create(dir.join("vdj_receptors.tsv"))?);
        writeln!(
            receptors,
            "cell\theavy_recombination_id\tlight_recombination_id\theavy_chain\theavy_v\theavy_d\theavy_d_inferred_from_vj_junction\theavy_d_hypothesis_margin\theavy_j\theavy_c\theavy_observed_v\theavy_observed_d\theavy_observed_j\theavy_naive_recombination\tlight_chain\tlight_v\tlight_j\tlight_c\tlight_observed_v\tlight_observed_j\tlight_naive_recombination\trag1_expression\trag2_expression\tdntt_expression\trag_activity\trag_pair_detected\tlight_chain_status"
        )?;
        let observed = if write_sequences {
            Some(BufWriter::new(File::create(dir.join("vdj_observed.fasta"))?))
        } else {
            None
        };
        let naive = if write_sequences {
            Some(BufWriter::new(File::create(dir.join("vdj_naive.fasta"))?))
        } else {
            None
        };
        Ok(Self {
            calls,
            receptors,
            observed,
            naive,
        })
    }

    fn write_cells(&mut self, cells: &[sc_vdj::CellVdjSummary], reference: &VdjReference) -> Result<()> {
        for cell in cells {
            let heavy_call = strongest_call(&cell.rearrangements, |chain| chain.has_d());
            let light_call = strongest_call(&cell.rearrangements, |chain| !chain.has_d());
            let heavy_id = cell.strongest_heavy_id().map(|x| x.to_string()).unwrap_or_default();
            let light_id = cell.strongest_light_id().map(|x| x.to_string()).unwrap_or_default();
            let heavy_junction = heavy_call.and_then(|call| call.junction.as_ref());
            let light_junction = light_call.and_then(|call| call.junction.as_ref());
            let sequence = |value: Option<&[u8]>| {
                value.map(|x| String::from_utf8_lossy(x).into_owned()).unwrap_or_default()
            };
            let receptor_fields = vec![
                cell.cell.clone(),
                heavy_id,
                light_id,
                call_chain(heavy_call).to_string(),
                call_segment(heavy_call, |x| &x.v).to_string(),
                call_segment(heavy_call, |x| &x.d).to_string(),
                heavy_call.map(|x| x.d_inferred_from_vj_junction.to_string()).unwrap_or_default(),
                heavy_call.and_then(|x| x.d_hypothesis_margin).map(|x| x.to_string()).unwrap_or_default(),
                call_segment(heavy_call, |x| &x.j).to_string(),
                call_segment(heavy_call, |x| &x.c).to_string(),
                sequence(heavy_junction.map(|x| x.observed_v.as_slice())),
                sequence(heavy_junction.map(|x| x.observed_d.as_slice())),
                sequence(heavy_junction.map(|x| x.observed_j.as_slice())),
                sequence(heavy_junction.map(|x| x.inferred_naive_sequence.as_slice())),
                call_chain(light_call).to_string(),
                call_segment(light_call, |x| &x.v).to_string(),
                call_segment(light_call, |x| &x.j).to_string(),
                call_segment(light_call, |x| &x.c).to_string(),
                sequence(light_junction.map(|x| x.observed_v.as_slice())),
                sequence(light_junction.map(|x| x.observed_j.as_slice())),
                sequence(light_junction.map(|x| x.inferred_naive_sequence.as_slice())),
                format!("{:.6}", cell.recombination_activity.rag1_expression),
                format!("{:.6}", cell.recombination_activity.rag2_expression),
                format!("{:.6}", cell.recombination_activity.dntt_expression),
                format!("{:.6}", cell.recombination_activity.rag_activity),
                cell.recombination_activity.rag_pair_detected.to_string(),
                cell.light_chain_status().to_string(),
            ];
            writeln!(self.receptors, "{}", receptor_fields.join("\t"))?;

            for call in &cell.rearrangements {
                let measurements = RecombinationMeasurements::from_call(call);
                let rearrangement_id = sc_vdj::PackedRecombinationId::from_call(call)
                    .map(|x| x.to_string())
                    .unwrap_or_default();
                let m = measurements.as_ref();
                let junction = call.junction.as_ref();
                let seq = |value: Option<&[u8]>| {
                    value
                        .map(|x| String::from_utf8_lossy(x).into_owned())
                        .unwrap_or_default()
                };
                let fields = vec![
                    cell.cell.clone(),
                    rearrangement_id.clone(),
                    call.chain.to_string(),
                    format!("{:?}", call.stage),
                    segment_id(&call.v).to_string(),
                    segment_id(&call.d).to_string(),
                    segment_id(&call.j).to_string(),
                    segment_id(&call.c).to_string(),
                    segment_id(&call.c).to_string(),
                    call.total_supporting_umis.to_string(),
                    support_i32(&call.v, |x| x.local_alignment_score),
                    support_usize(&call.v, |x| x.supporting_reads),
                    support_usize(&call.v, |x| x.supporting_umis),
                    support_f64(&call.v, |x| x.locus_fraction),
                    support_u64(&call.v, |x| x.distance_to_recombination_center),
                    support_i32(&call.d, |x| x.local_alignment_score),
                    support_usize(&call.d, |x| x.supporting_reads),
                    support_usize(&call.d, |x| x.supporting_umis),
                    call.d_inferred_from_vj_junction.to_string(),
                    call.d_hypothesis_margin.map(|x| x.to_string()).unwrap_or_default(),
                    support_i32(&call.j, |x| x.local_alignment_score),
                    support_usize(&call.j, |x| x.supporting_reads),
                    support_usize(&call.j, |x| x.supporting_umis),
                    support_i32(&call.c, |x| x.local_alignment_score),
                    support_usize(&call.c, |x| x.supporting_reads),
                    support_usize(&call.c, |x| x.supporting_umis),
                    measurement(m.and_then(|x| x.v_del_3)),
                    measurement(m.and_then(|x| x.p_v3_len)),
                    seq(junction.map(|x| x.p_v3.as_slice())),
                    measurement(m.and_then(|x| x.n1_len)),
                    seq(junction.map(|x| x.n1.as_slice())),
                    measurement(m.and_then(|x| x.p_d5_len)),
                    seq(junction.map(|x| x.p_d5.as_slice())),
                    measurement(m.and_then(|x| x.d_del_5)),
                    measurement(m.and_then(|x| x.d_retained_len)),
                    measurement(m.and_then(|x| x.d_del_3)),
                    measurement(m.and_then(|x| x.p_d3_len)),
                    seq(junction.map(|x| x.p_d3.as_slice())),
                    measurement(m.and_then(|x| x.n2_len)),
                    seq(junction.map(|x| x.n2.as_slice())),
                    measurement(m.and_then(|x| x.p_j5_len)),
                    seq(junction.map(|x| x.p_j5.as_slice())),
                    measurement(m.and_then(|x| x.j_del_5)),
                    m.map(|x| x.pn_alternative.to_string()).unwrap_or_default(),
                    germline_sequence(reference, &call.v),
                    germline_sequence(reference, &call.d),
                    germline_sequence(reference, &call.j),
                    germline_sequence(reference, &call.c),
                    seq(junction.map(|x| x.observed_v.as_slice())),
                    seq(junction.map(|x| x.observed_d.as_slice())),
                    seq(junction.map(|x| x.observed_j.as_slice())),
                    seq(junction.map(|x| x.naive_v.as_slice())),
                    seq(junction.map(|x| x.naive_d.as_slice())),
                    seq(junction.map(|x| x.naive_j.as_slice())),
                    seq(junction.map(|x| x.observed_sequence.as_slice())),
                    seq(junction.map(|x| x.inferred_naive_sequence.as_slice())),
                ];
                writeln!(self.calls, "{}", fields.join("\t"))?;

                if let Some(junction) = junction {
                    let header = format!(
                        "{}|{}|{}",
                        fasta_token(&cell.cell),
                        call.chain,
                        fasta_token(&rearrangement_id)
                    );
                    if let Some(writer) = &mut self.observed {
                        writeln!(writer, ">{header}")?;
                        writeln!(writer, "{}", String::from_utf8_lossy(&junction.observed_sequence))?;
                    }
                    if let Some(writer) = &mut self.naive {
                        writeln!(writer, ">{header}")?;
                        writeln!(
                            writer,
                            "{}",
                            String::from_utf8_lossy(&junction.inferred_naive_sequence)
                        )?;
                    }
                }
            }
        }
        Ok(())
    }

    fn finish(mut self) -> Result<()> {
        self.calls.flush()?;
        self.receptors.flush()?;
        if let Some(mut writer) = self.observed {
            writer.flush()?;
        }
        if let Some(mut writer) = self.naive {
            writer.flush()?;
        }
        Ok(())
    }
}

fn strongest_call<F>(
    calls: &[sc_vdj::RearrangementCall],
    keep: F,
) -> Option<&sc_vdj::RearrangementCall>
where
    F: Fn(sc_vdj::Chain) -> bool,
{
    calls
        .iter()
        .filter(|call| keep(call.chain))
        .max_by(|a, b| {
            a.total_supporting_umis
                .cmp(&b.total_supporting_umis)
                .then_with(|| call_score(a).cmp(&call_score(b)))
                .then_with(|| a.chain.cmp(&b.chain))
                .then_with(|| a.notation.cmp(&b.notation))
        })
}

fn call_score(call: &sc_vdj::RearrangementCall) -> i32 {
    call.v.as_ref().map_or(0, |x| x.local_alignment_score)
        + call.d.as_ref().map_or(0, |x| x.local_alignment_score)
        + call.j.as_ref().map_or(0, |x| x.local_alignment_score)
}

fn call_chain(call: Option<&sc_vdj::RearrangementCall>) -> String {
    call.map(|x| x.chain.to_string()).unwrap_or_default()
}

fn call_segment<'a, F>(call: Option<&'a sc_vdj::RearrangementCall>, pick: F) -> &'a str
where
    F: Fn(&'a sc_vdj::RearrangementCall) -> &'a Option<GermlineSegmentSupport>,
{
    call.and_then(|x| pick(x).as_ref()).map(|x| x.id.as_str()).unwrap_or("")
}

fn germline_sequence(reference: &VdjReference, support: &Option<GermlineSegmentSupport>) -> String {
    support
        .as_ref()
        .and_then(|x| reference.segments.get(x.segment_index))
        .map(|x| String::from_utf8_lossy(&x.sequence).into_owned())
        .unwrap_or_default()
}

fn support_i32<F>(support: &Option<GermlineSegmentSupport>, get: F) -> String
where
    F: Fn(&GermlineSegmentSupport) -> i32,
{
    support.as_ref().map(|x| get(x).to_string()).unwrap_or_default()
}

fn support_usize<F>(support: &Option<GermlineSegmentSupport>, get: F) -> String
where
    F: Fn(&GermlineSegmentSupport) -> usize,
{
    support.as_ref().map(|x| get(x).to_string()).unwrap_or_default()
}

fn support_u64<F>(support: &Option<GermlineSegmentSupport>, get: F) -> String
where
    F: Fn(&GermlineSegmentSupport) -> u64,
{
    support.as_ref().map(|x| get(x).to_string()).unwrap_or_default()
}

fn support_f64<F>(support: &Option<GermlineSegmentSupport>, get: F) -> String
where
    F: Fn(&GermlineSegmentSupport) -> f64,
{
    support.as_ref().map(|x| format!("{:.6}", get(x))).unwrap_or_default()
}

fn main() -> Result<()> {
    let cli = Cli::parse().validate()?;
    let status = Arc::new(RwLock::new(VdjRunStatus::new()));
    let _status_server = if cli.no_health_server {
        None
    } else {
        let addr = SocketAddr::new(IpAddr::V4(Ipv4Addr::UNSPECIFIED), cli.health_port);
        let server = spawn_status_server(Arc::clone(&status), addr)?;
        let hostname = public_hostname(cli.health_hostname.as_deref());
        let url = format!("http://{}:{}", hostname, server.addr().port());
        update_status(&status, |state| state.public_url = Some(url.clone()));
        eprintln!("nelrune-vdj: live status: {url}");
        Some(server)
    };

    let mut mapping_info = MappingInfo::new(None, 0.0, 0);

    set_stage(&status, "loading V(D)J reference");
    mapping_info.start_timer("wall.reference_index");
    let mapper = if let Some(index) = &cli.index {
        VdjMapper::load_index(index)
            .with_context(|| format!("loading V(D)J index {}", index.display()))?
    } else {
        let reference = VdjReferenceBuilder::default()
            .build(
                cli.gtf.as_ref().expect("validated --gtf"),
                cli.genome.as_ref().expect("validated --genome"),
            )
            .context("building V(D)J reference")?;
        VdjMapper::new(reference, VdjMapperConfig::default())
    };
    mapping_info.stop_timer("wall.reference_index");
    let reference = mapper.reference();
    update_status(&status, |state| state.reference_segments = reference.len());

    set_stage(&status, "loading expression matrix");
    mapping_info.start_timer("wall.expression_matrix");
    let gex = read_mex(&cli.exonic, cli.cell_barcode_len)?;
    mapping_info.stop_timer("wall.expression_matrix");
    if gex.cells.is_empty() {
        bail!("MEX contains no cells");
    }

    let called_cells: HashSet<String> = gex.cells.iter().cloned().collect();
    update_status(&status, |state| state.cells = gex.cells.len());
    eprintln!(
        "nelrune-vdj: {} cells, one-pass in-memory BAM routing, {} reference segments",
        gex.cells.len(),
        reference.len()
    );

    let posterior_config = PosteriorConfig {
        sterile_bins: cli.sterile_bins,
        ..PosteriorConfig::default()
    };
    let analyzer = PosteriorAnalyzer::new(reference, &mapper, posterior_config);
    let mut writer = BatchWriter::create(&cli.out, cli.write_sequences)?;

    set_stage(&status, "routing BAM evidence by receptor locus");
    update_status(&status, |state| state.bam_scan_started_unix_ms = Some(unix_ms()));
    mapping_info.start_timer("wall.bam_route_in_memory");
    let progress_status = Arc::clone(&status);
    let (routed, stats) = read_bam_receptor_evidence_with_progress(
        &cli.bam,
        &NelruneIdentityResolver,
        |cell| called_cells.contains(cell),
        reference,
        move |stats| {
            update_status(&progress_status, |state| {
                state.receptor_chromosomes = stats.receptor_chromosomes;
                state.indexed_segment_intervals = stats.indexed_segment_intervals;
                state.total_bam_records = stats.total_records;
                state.receptor_chromosome_records = stats.receptor_chromosome_records;
                state.segment_overlap_records = stats.segment_overlap_records;
                state.called_cell_records = stats.called_cell_records;
                state.routed_evidence_records = stats.routed_evidence_records;
                state.non_receptor_chromosome_records = stats.non_receptor_chromosome_records;
                state.non_segment_overlap_records = stats.non_segment_overlap_records;
                state.locus_records = stats.locus_records;
            });
        },
    )?;
    mapping_info.stop_timer("wall.bam_route_in_memory");
    update_status(&status, |state| state.bam_scan_finished_unix_ms = Some(unix_ms()));

    mapping_info.total = stats.total_records;
    mapping_info.report_n("vdj.receptor_chromosome_records", stats.receptor_chromosome_records);
    mapping_info.report_n("vdj.segment_overlap_records", stats.segment_overlap_records);
    mapping_info.report_n("vdj.called_cell_vdj_records", stats.called_cell_records);
    mapping_info.report_n("vdj.routed_locus_observations", stats.routed_evidence_records);
    for (chain, count) in Chain::ALL.into_iter().zip(stats.locus_records) {
        mapping_info.report_n(&format!("vdj.routed.{}", chain.as_str().to_ascii_lowercase()), count);
    }
    eprintln!(
        "nelrune-vdj: BAM routed: total={} receptor_chr={} segment_overlap={} called_cell={} routed={}",
        stats.total_records,
        stats.receptor_chromosome_records,
        stats.segment_overlap_records,
        stats.called_cell_records,
        stats.routed_evidence_records
    );

    set_stage(&status, "reconstructing V(D)J receptors from locus-routed evidence");
    update_status(&status, |state| state.reconstruction_started_unix_ms = Some(unix_ms()));
    mapping_info.start_timer("wall.posterior");
    let reconstruction_status = Arc::clone(&status);
    let summaries = analyzer.analyze_routed_for_cells_with_progress(
        routed,
        &gex,
        gex.cells.iter().cloned(),
        &mut mapping_info,
        move |summary, local| {
            let counter = |name: &str| local.reads_log.get(name).copied().unwrap_or(0);
            update_status(&reconstruction_status, |state| {
                state.summaries_written += 1;
                state.umis += counter("vdj.umis");
                state.umi_contigs += counter("vdj.umi_contigs");
                state.seeded_contigs += counter("vdj.seeded_contigs");
                state.initial_local_alignments += counter("vdj.initial_local_alignments");
                state.rearrangements += summary.rearrangements.len();
                state.heavy_ids += if summary.strongest_heavy_id().is_some() { 1 } else { 0 };
                state.light_ids += if summary.strongest_light_id().is_some() { 1 } else { 0 };
            });
        },
    );
    mapping_info.stop_timer("wall.posterior");
    sync_reconstruction_status(&status, &mapping_info);
    let summary_count = summaries.len();
    let rearrangements: usize = summaries.iter().map(|cell| cell.rearrangements.len()).sum();
    let heavy_ids = summaries.iter().filter(|cell| cell.strongest_heavy_id().is_some()).count();
    let light_ids = summaries.iter().filter(|cell| cell.strongest_light_id().is_some()).count();
    update_status(&status, |state| {
        state.summaries_written = summary_count;
        state.rearrangements = rearrangements;
        state.heavy_ids = heavy_ids;
        state.light_ids = light_ids;
    });
    writer.write_cells(&summaries, reference)?;

    if summary_count == 0 {
        bail!("posterior analyzer produced no cell summaries");
    }
    set_stage(&status, "writing V(D)J outputs");
    writer.finish()?;
    let performance_path = cli.out.join("vdj-mapping-info.txt");
    fs::write(&performance_path, mapping_info.to_string())?;
    update_status(&status, |state| {
        state.finished_unix_ms = Some(unix_ms());
        state.stage = "finished".to_string();
    });
    eprintln!(
        "nelrune-vdj: wrote {} cell summaries to {}",
        summary_count,
        cli.out.display()
    );
    Ok(())
}

fn segment_id(segment: &Option<sc_vdj::GermlineSegmentSupport>) -> &str {
    segment.as_ref().map(|x| x.id.as_str()).unwrap_or("")
}

fn measurement(value: Option<u16>) -> String {
    value.map(|x| x.to_string()).unwrap_or_default()
}

fn fasta_token(value: &str) -> String {
    value
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() || matches!(c, '-' | '_' | '.' | ':') { c } else { '_' })
        .collect()
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

fn reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let input: Box<dyn IoRead> = if path.extension().and_then(|x| x.to_str()) == Some("gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(input)))
}

fn read_mex(dir: &Path, cell_barcode_len: Option<usize>) -> Result<MexExpression> {
    let barcode_path = find_file(dir, &["barcodes.tsv.gz", "barcodes.tsv"])?;
    let feature_path = find_file(
        dir,
        &["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"],
    )?;
    let matrix_path = find_file(dir, &["matrix.mtx.gz", "matrix.mtx"])?;

    let raw_cells: Vec<String> = reader(&barcode_path)?
        .lines()
        .collect::<std::io::Result<Vec<_>>>()?
        .into_iter()
        .filter_map(|line| line.split('\t').next().map(str::to_string))
        .filter(|x| !x.is_empty())
        .collect();

    let mut cells = Vec::with_capacity(raw_cells.len());
    let mut seen = HashSet::with_capacity(raw_cells.len());
    for raw in raw_cells {
        let cell = match cell_barcode_len {
            Some(len) => raw
                .get(..len)
                .with_context(|| format!("matrix barcode {raw} is shorter than --cell-barcode-len {len}"))?
                .to_string(),
            None => raw,
        };
        if !seen.insert(cell.clone()) {
            bail!("cell barcode normalization produced duplicate barcode {cell}");
        }
        cells.push(cell);
    }

    let genes: Vec<String> = reader(&feature_path)?
        .lines()
        .collect::<std::io::Result<Vec<_>>>()?
        .into_iter()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            fields
                .get(1)
                .copied()
                .unwrap_or(fields[0])
                .to_ascii_uppercase()
        })
        .collect();

    let mut input = reader(&matrix_path)?;
    let mut line = String::new();
    loop {
        line.clear();
        if input.read_line(&mut line)? == 0 {
            bail!("matrix ended before dimensions");
        }
        if !line.starts_with('%') {
            break;
        }
    }
    let dims: Vec<&str> = line.split_whitespace().collect();
    if dims.len() != 3 {
        bail!("invalid MatrixMarket dimensions line");
    }
    let rows: usize = dims[0].parse()?;
    let cols: usize = dims[1].parse()?;
    if rows != genes.len() || cols != cells.len() {
        bail!(
            "MEX dimensions {}x{} do not match {} features and {} cells",
            rows,
            cols,
            genes.len(),
            cells.len()
        );
    }

    let mut values = HashMap::new();
    for line in input.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split_whitespace();
        let row: usize = fields.next().context("matrix row")?.parse()?;
        let col: usize = fields.next().context("matrix col")?.parse()?;
        let value: f64 = fields.next().context("matrix value")?.parse()?;
        if row == 0 || col == 0 || row > genes.len() || col > cells.len() {
            bail!("matrix entry outside declared dimensions: row={row}, col={col}");
        }
        values.insert((cells[col - 1].clone(), genes[row - 1].clone()), value);
    }
    Ok(MexExpression { cells, values })
}
