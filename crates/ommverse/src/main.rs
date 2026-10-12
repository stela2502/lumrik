use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use lumrik_status::{
    ServerContent, ServerSnapshot, StatusMetric, StatusSection, memory_status, public_hostname,
    spawn_status_server,
};
use ommverse::{
    AaModelVault, ExactAaFeatureModel, InterProImportReport, Ommverse, ProteinFeatureKind,
    ProteinFeatureModel, TopologyEvaluation, TopologyState, ingest_interpro_many,
};
use std::net::{IpAddr, Ipv4Addr, SocketAddr};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::{Arc, RwLock};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

#[derive(Parser)]
#[command(
    name = "ommverse",
    about = "Build and explore Lumrik's integrated biological reference model"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}
#[derive(Subcommand)]
enum Command {
    /// Render Ommverse reference-model visualizations. Additional plot families can be added here.
    Plot {
        #[command(subcommand)]
        plot: PlotCommand,
    },
    Build {
        /// Existing downloaded reference root. Mutually exclusive with --assembly.
        #[arg(long, conflicts_with = "assembly")]
        reference: Option<PathBuf>,
        /// Fetch supported sources for this assembly into --cache and build them.
        #[arg(long, conflicts_with = "reference")]
        assembly: Option<String>,
        #[arg(long, default_value = ".ommverse-cache")]
        cache: PathBuf,
        /// Use this GTF exactly. Disables automatic GENCODE annotation selection.
        #[arg(long, conflicts_with_all = ["reference", "gencode_release"])]
        gtf: Option<PathBuf>,
        /// Pin a GENCODE release (for example 50 or M38). Defaults to latest when supported.
        #[arg(long, conflicts_with_all = ["reference", "gtf"])]
        gencode_release: Option<String>,
        /// Output Ommverse index. Defaults to <cache>/<assembly>/<assembly>.ommverse in --assembly mode.
        #[arg(long)]
        out: Option<PathBuf>,
        #[arg(long)]
        debug_failed_mappings: bool,
        #[arg(long, default_value_t = 8787)]
        health_port: u16,
        #[arg(long)]
        health_hostname: Option<String>,
        #[arg(long)]
        no_health_server: bool,
    },
    /// Stream InterPro protein2ipr annotations into an existing Ommverse index.
    IngestInterpro {
        /// One Ommverse index. Mutually exclusive with --index-dir.
        #[arg(long, conflicts_with = "index_dir")]
        index: Option<PathBuf>,
        /// Recursively enrich every base *.ommverse index below this directory in one InterPro pass.
        #[arg(long, conflicts_with = "index")]
        index_dir: Option<PathBuf>,
        #[arg(long)]
        protein2ipr: PathBuf,
        #[arg(long)]
        entry_list: PathBuf,
        /// Optional InterPro ParentChildTreeFile.txt; stored once in each enriched index.
        #[arg(long)]
        parent_child_tree: Option<PathBuf>,
        /// Output path for single-index mode. Directory mode writes <stem>.interpro.ommverse beside each input.
        #[arg(long, requires = "index")]
        out: Option<PathBuf>,
        #[arg(long, default_value_t = 8787)]
        health_port: u16,
        #[arg(long)]
        health_hostname: Option<String>,
        #[arg(long)]
        no_health_server: bool,
    },
    /// Import one or more MSigDB GMT files into a compact, queryable sidecar.
    /// Download licensed MSigDB collections separately; no redistribution is performed.
    IngestMsigdb {
        #[arg(long)]
        index: PathBuf,
        /// Optional local GMT files. If omitted, download Hallmark, C2 and C3.
        #[arg(long, num_args = 1..)]
        gmt: Vec<PathBuf>,
        /// MSigDB human release, pinned for reproducibility.
        #[arg(long, default_value = "2026.1.Hs")]
        release: String,
        /// Download cache directory (default: <index parent>/MSigDB).
        #[arg(long)]
        cache_dir: Option<PathBuf>,
    },
    /// Look up MSigDB sets containing a gene symbol, without loading the genome index.
    MsigdbGene {
        #[arg(long)]
        index: PathBuf,
        gene: String,
    },
    /// Retrieve members of a named MSigDB set, without loading the genome index.
    MsigdbSet {
        #[arg(long)]
        index: PathBuf,
        set: String,
    },
    /// Query the biological reference model with the Ommverse query language.
    Query {
        #[arg(long)]
        index: PathBuf,
        query: String,
    },
    /// Print the small capability schema without opening the Ommverse index.
    Schema {
        #[arg(long)]
        index: PathBuf,
        /// Create or replace the tiny schema sidecar for an existing index.
        /// This opens the Ommverse once; later schema reads and query validation do not.
        #[arg(long)]
        write: bool,
    },
    /// Inspect serialized core composition and external reference storage without changing either.
    Inspect {
        #[arg(long)]
        index: PathBuf,
    },
    /// Probe whether STRING contains interaction profiles that differ across InterPro feature architectures.
    StringProbe {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        info: PathBuf,
        #[arg(long)]
        links: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
    Protein {
        #[arg(long)]
        index: PathBuf,
        accession: String,
        #[arg(long)]
        sequence: bool,
    },
    Gene {
        #[arg(long)]
        index: PathBuf,
        symbol: String,
    },
    /// Build and validate a supervised protein-feature dataset from an Ommverse index.
    ///
    /// Complete protein sequences are reconstructed from the reference genome and
    /// transcript annotation. Curated feature intervals are projected onto those
    /// sequences and proteins with inconsistent annotations are rejected. The
    /// surviving proteins are split deterministically 50:50 at protein level.
    /// Residues outside the selected feature remain in the examples as background.
    /// This command validates/describes the corpus; it does not train a model.
    TrainingCorpus {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        feature: String,
        /// Context width used only to report collisions between adjacent feature flanks.
        #[arg(long, default_value_t = 10)]
        flank_width: usize,
    },

    /// Train and evaluate a feature HMM from exact amino-acid identity.
    ///
    /// The selected curated feature supplies supervised residue labels. Training
    /// uses only the deterministic training half; the other half is held out for
    /// evaluation. Observations are exact IntToProt residue identities: no chemistry,
    /// substitution matrix, or protein-language representation is supplied.
    ///
    /// The HMM has four biological states: BACKGROUND -> PRE_FEATURE -> FEATURE ->
    /// POST_FEATURE -> BACKGROUND. PRE/POST describe the immediately adjacent
    /// residues outside each curated feature. --flank-width controls that context
    /// window (default 10 aa). With --scan-unannotated the frozen model is then
    /// applied to proteins lacking the selected annotation; they never train it.
    TrainExactAa {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        feature: String,
        #[arg(long)]
        out: PathBuf,
        /// Number of residues immediately outside each curated feature assigned to PRE/POST states.
        #[arg(long, default_value_t = 10)]
        flank_width: usize,
        #[arg(long)]
        scan_unannotated: bool,
    },

    /// Train, inspect, and apply portable amino-acid feature model vaults.
    ModelVault {
        #[command(subcommand)]
        command: ModelVaultCommand,
    },

    /// Apply every frozen model in a vault to the matching curated feature in another Ommverse index.
    ///
    /// Nothing is fitted or updated: the vault is loaded exactly as trained and each
    /// model is evaluated against the target index's curated annotations. This is the
    /// cross-reference / cross-species transfer test.
    EvaluateModelVault {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        vault: PathBuf,
    },

    /// Apply a previously trained exact-AA feature model to unannotated proteins.
    ScanExactAa {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        model: PathBuf,
    },
 }

#[derive(Subcommand)]
enum ModelVaultCommand {
    /// Train all supported amino-acid feature HMMs and write <vault>.info beside the vault.
    Train {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, default_value_t = 10)]
        flank_width: usize,
        /// Maximum canonical InterPro feature length included in vault training.
        #[arg(long, default_value_t = 50)]
        interpro_max_len: usize,
    },
    /// Show provenance and held-out metrics stored with a model vault.
    Info {
        #[arg(long)]
        vault: PathBuf,
    },
    /// Apply every frozen model to arbitrary amino-acid FASTA sequences.
    ///
    /// Writes <out>.responses.tsv (every sequence x every model, no significance filtering)
    /// and <out>.features.tsv (0-based half-open predicted intervals). HMM raw scores are
    /// never labelled as p-values; calibrated_p_value remains empty until null calibration exists.
    Predict {
        #[arg(long)]
        vault: PathBuf,
        #[arg(long)]
        fasta: PathBuf,
        /// Output prefix for .responses.tsv and .features.tsv.
        #[arg(long)]
        out: PathBuf,
    },
}

fn optional_float(value: Option<f64>) -> String {
    value.map(|x| format!("{x:.6}")).unwrap_or_default()
}

fn read_aa_fasta(path: &Path) -> Result<Vec<(String, String)>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading amino-acid FASTA {}", path.display()))?;
    let mut records = Vec::new();
    let mut id: Option<String> = None;
    let mut sequence = String::new();
    for (line_no, raw) in text.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() { continue; }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = id.take() {
                if sequence.is_empty() { anyhow::bail!("FASTA record {previous} has no sequence"); }
                records.push((previous, std::mem::take(&mut sequence)));
            }
            let next = header.split_whitespace().next().unwrap_or("");
            if next.is_empty() { anyhow::bail!("empty FASTA identifier at line {}", line_no + 1); }
            id = Some(next.to_owned());
        } else {
            if id.is_none() { anyhow::bail!("FASTA sequence before first header at line {}", line_no + 1); }
            sequence.push_str(line);
        }
    }
    if let Some(previous) = id {
        if sequence.is_empty() { anyhow::bail!("FASTA record {previous} has no sequence"); }
        records.push((previous, sequence));
    }
    if records.is_empty() { anyhow::bail!("FASTA contains no sequences: {}", path.display()); }
    Ok(records)
}

fn format_model_vault_info(vault: &AaModelVault) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    writeln!(&mut out, "Ommverse AA model vault").unwrap();
    writeln!(&mut out, "  assembly:                  {}", vault.metadata.assembly).unwrap();
    writeln!(&mut out, "  source:                    {}", vault.metadata.source).unwrap();
    writeln!(&mut out, "  Ommverse format version:   {}", vault.metadata.ommverse_format_version).unwrap();
    writeln!(&mut out, "  vault format version:      {}", vault.metadata.vault_format_version).unwrap();
    writeln!(&mut out, "  flank width:               {} aa", vault.metadata.flank_width).unwrap();
    writeln!(&mut out, "  models:                    {}", vault.models.len() + vault.interpro_models.len()).unwrap();
    for entry in &vault.models {
        match entry {
            ProteinFeatureModel::ExactAa { model, evaluation } => writeln!(&mut out,
                "    {:?} / exact-AA: precision {:.4}, recall {:.4}, F1 {:.4}", model.feature,
                evaluation.precision(), evaluation.recall(), evaluation.f1()).unwrap(),
            ProteinFeatureModel::Chemistry { model, evaluation } => writeln!(&mut out,
                "    {:?} / chemistry: precision {:.4}, recall {:.4}, F1 {:.4}", model.feature,
                evaluation.precision(), evaluation.recall(), evaluation.f1()).unwrap(),
            ProteinFeatureModel::TopologyExactAa { evaluation, .. } => {
                writeln!(&mut out, "    topology / exact-AA: {} biological held-out state metrics", evaluation.states.len()).unwrap();
            }
            ProteinFeatureModel::TopologyChemistry { evaluation, .. } => {
                writeln!(&mut out, "    topology / chemistry: {} biological held-out state metrics", evaluation.states.len()).unwrap();
            }
        }
    }
    for entry in &vault.interpro_models {
        writeln!(&mut out, "    {} / exact-AA: precision {:.4}, recall {:.4}, F1 {:.4}",
            entry.model.feature, entry.evaluation.precision(), entry.evaluation.recall(), entry.evaluation.f1()).unwrap();
    }
    out
}

fn human_bytes(bytes: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];
    let mut value = bytes as f64;
    let mut unit = 0usize;
    while value >= 1024.0 && unit + 1 < UNITS.len() {
        value /= 1024.0;
        unit += 1;
    }
    if unit == 0 {
        format!("{bytes} B")
    } else {
        format!("{value:.2} {} ({bytes} B)", UNITS[unit])
    }
}

fn bincode_size<T: serde::Serialize>(value: &T) -> Result<u64> {
    bincode::serialized_size(value).context("measuring standalone bincode size")
}

fn directory_bytes(path: &Path, exclude: Option<&Path>) -> Result<(u64, usize)> {
    let mut bytes = 0u64;
    let mut files = 0usize;
    if !path.exists() {
        return Ok((0, 0));
    }
    let mut stack = vec![path.to_path_buf()];
    while let Some(dir) = stack.pop() {
        for entry in std::fs::read_dir(&dir).with_context(|| format!("reading {}", dir.display()))? {
            let entry = entry?;
            let entry_path = entry.path();
            if exclude.is_some_and(|excluded| entry_path == excluded) {
                continue;
            }
            let ty = entry.file_type()?;
            if ty.is_dir() {
                stack.push(entry_path);
            } else if ty.is_file() {
                bytes += entry.metadata()?.len();
                files += 1;
            }
        }
    }
    Ok((bytes, files))
}

fn print_inspection(index: &Path, omm: &Ommverse, load_seconds: f64) -> Result<()> {
    let index_bytes = std::fs::metadata(index)
        .with_context(|| format!("reading metadata for {}", index.display()))?
        .len();
    println!("Ommverse inspect — {}", omm.assembly);
    println!("  index:                     {}", index.display());
    println!("  format version:            {}", ommverse::OMMVERSE_FORMAT_VERSION);
    println!("  exact index bytes:         {}", human_bytes(index_bytes));
    println!("  deserialize + rebuild:     {load_seconds:.3} s");

    println!("\nIntegrated object counts (exact)");
    println!("  splice genes:              {}", omm.splice.genes.len());
    println!("  splice transcripts:        {}", omm.splice.transcripts.len());
    println!("  proteins:                  {}", omm.proteins.len());
    println!("  protein features:          {}", omm.proteins.iter().map(|p| p.features.len()).sum::<usize>());
    println!("  InterPro entries:          {}", omm.interpro_entries.len());
    println!("  chromatin elements:        {}", omm.chromatin.len());
    println!("  binding-union regions:     {}", omm.protein_binding.regions.len());
    println!("  source binding peaks:      {} (count only; raw peaks are not stored in v{})", omm.protein_binding.source_peak_count, ommverse::OMMVERSE_FORMAT_VERSION);
    println!("  CTCF anchors:              {}", omm.ctcf.anchors.len());
    println!("  CTCF domains:              {}", omm.ctcf.domains.len());
    println!("  observed loop anchors:     {}", omm.experimental_loops.anchors.len());
    println!("  observed loops:            {}", omm.experimental_loops.loops.len());

    let sections = [
        ("assembly", bincode_size(&omm.assembly)?),
        ("source_root", bincode_size(&omm.source_root)?),
        ("genome_twobit", bincode_size(&omm.genome_twobit)?),
        ("splice", bincode_size(&omm.splice)?),
        ("proteins", bincode_size(&omm.proteins)?),
        ("report", bincode_size(&omm.report)?),
        ("interpro_entries", bincode_size(&omm.interpro_entries)?),
        ("chromatin", bincode_size(&omm.chromatin)?),
        ("protein_binding", bincode_size(&omm.protein_binding)?),
        ("ctcf", bincode_size(&omm.ctcf)?),
        ("experimental_loops", bincode_size(&omm.experimental_loops)?),
        ("sources", bincode_size(&omm.sources)?),
    ];
    let diagnostic_sum: u64 = sections.iter().map(|(_, bytes)| *bytes).sum();
    println!("\nSerialized core diagnostics");
    println!("  These are exact bincode sizes when each field is serialized standalone.");
    println!("  They identify large payloads, but are not claimed as byte ranges within the file.");
    for (name, bytes) in sections {
        let pct = if index_bytes == 0 { 0.0 } else { 100.0 * bytes as f64 / index_bytes as f64 };
        println!("  {name:<25} {:>24}  {:>6.2}% of index", human_bytes(bytes), pct);
    }
    println!("  diagnostic field sum:      {}", human_bytes(diagnostic_sum));
    println!("  index framing/difference:  {} B", index_bytes as i128 - 8 - diagnostic_sum as i128);
    println!("  note: index has an 8-byte Ommverse magic/version header before the bincode payload");

    println!("\nCalculated owned payloads (not allocator/RAM measurements)");
    let protein_string_bytes: usize = omm.proteins.iter().map(|p| {
        p.accession.len() + p.entry_name.len() + p.name.len() + p.gene_symbol.len()
            + p.aliases.iter().map(String::len).sum::<usize>()
            + p.ensembl_gene.as_ref().map_or(0, String::len)
            + p.ensembl_protein.as_ref().map_or(0, String::len)
    }).sum();
    let feature_string_bytes: usize = omm.proteins.iter().flat_map(|p| &p.features).map(|f| {
        f.label.len() + f.description.len() + f.chromosome.len() + f.source_db.len()
    }).sum();
    println!("  protein UTF-8 contents:    {}", human_bytes(protein_string_bytes as u64));
    println!("  feature UTF-8 contents:    {}", human_bytes(feature_string_bytes as u64));
    println!("  protein transcript IDs:    {}", human_bytes((omm.proteins.iter().map(|p| p.transcript_ids.len()).sum::<usize>() * std::mem::size_of::<usize>()) as u64));
    println!("  binding regions structs:   {}", human_bytes((omm.protein_binding.regions.len() * std::mem::size_of::<ommverse::ProteinBindingRegion>()) as u64));
    println!("  CTCF anchor structs:       {}", human_bytes((omm.ctcf.anchors.len() * std::mem::size_of::<ommverse::CtcfAnchor>()) as u64));
    println!("  CTCF domain structs:       {}", human_bytes((omm.ctcf.domains.len() * std::mem::size_of::<ommverse::CtcfDomain>()) as u64));
    println!("  loop anchor structs:       {}", human_bytes((omm.experimental_loops.anchors.len() * std::mem::size_of::<ommverse::ExperimentalLoopAnchor>()) as u64));
    println!("  loop structs:              {}", human_bytes((omm.experimental_loops.loops.len() * std::mem::size_of::<ommverse::ExperimentalLoop>()) as u64));

    println!("\nExternal/source estate (exact filesystem bytes)");
    println!("  source root:               {}", omm.source_root.display());
    let index_canonical = std::fs::canonicalize(index).ok();
    let (root_bytes, root_files) = directory_bytes(&omm.source_root, index_canonical.as_deref())?;
    println!("  source tree excluding index: {} in {} files", human_bytes(root_bytes), root_files);
    if let Some(manifest) = &omm.sources {
        let mut manifest_total = 0u64;
        for resource in &manifest.resources {
            let Some(rel) = &resource.file else { continue; };
            let path = omm.source_root.join(rel);
            match std::fs::metadata(&path) {
                Ok(meta) if meta.is_file() => {
                    manifest_total += meta.len();
                    println!("  {:<27} {:<12} {:>24}  {}", resource.kind, resource.source, human_bytes(meta.len()), rel.display());
                }
                _ => println!("  {:<27} {:<12} {:>24}  {}", resource.kind, resource.source, "missing", rel.display()),
            }
        }
        println!("  manifest file total:       {}", human_bytes(manifest_total));
    } else {
        println!("  manifest:                  not embedded in this index");
    }
    Ok(())
}

fn percent(n: usize, d: usize) -> f64 {
    if d == 0 {
        0.0
    } else {
        100.0 * n as f64 / d as f64
    }
}

fn print_distribution(label: &str, s: &ommverse::DistributionSummary) {
    println!("  {label}:");
    println!("    n:                       {}", s.count);
    println!("    mean:                    {:.2} aa", s.mean);
    println!("    median:                  {:.2} aa", s.median);
    println!("    SD:                      {:.2} aa", s.sd);
    println!("    Q1 / Q3:                 {:.2} / {:.2} aa", s.q1, s.q3);
    println!("    min / max:               {} / {} aa", s.min, s.max);
}

#[derive(Debug, Clone)]
struct BuildStatus {
    started_unix_ms: u128,
    finished_unix_ms: Option<u128>,
    assembly: String,
    stage: String,
    stage_records: usize,
    public_url: Option<String>,
    transcripts: usize,
    proteins: usize,
    chromatin_elements: usize,
    binding_regions: usize,
    source_binding_peaks: usize,
}
impl BuildStatus {
    fn new(assembly: String) -> Self {
        Self {
            started_unix_ms: SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap_or_default()
                .as_millis(),
            finished_unix_ms: None,
            assembly,
            stage: "starting".to_owned(),
            stage_records: 0,
            public_url: None,
            transcripts: 0,
            proteins: 0,
            chromatin_elements: 0,
            binding_regions: 0,
            source_binding_peaks: 0,
        }
    }
    fn stage(&mut self, stage: &str, records: usize) {
        self.stage = stage.to_owned();
        self.stage_records = records;
    }
    fn model(&mut self, model: &Ommverse) {
        self.transcripts = model.report.transcripts;
        self.proteins = model.report.mapped_proteins;
        self.chromatin_elements = model.chromatin.len();
        self.binding_regions = model.protein_binding.regions.len();
        self.source_binding_peaks = model.protein_binding.source_peak_count;
    }
    fn finish(&mut self) {
        self.stage = "complete".to_owned();
        self.finished_unix_ms = Some(
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap_or_default()
                .as_millis(),
        );
    }
}
impl ServerContent for BuildStatus {
    fn server_snapshot(&self) -> ServerSnapshot {
        let mem = memory_status();
        ServerSnapshot {
            title: format!("Ommverse build — {}", self.assembly),
            subtitle: "Integrated biological reference build".to_owned(),
            started_unix_ms: self.started_unix_ms,
            finished_unix_ms: self.finished_unix_ms,
            stage: self.stage.clone(),
            public_url: self.public_url.clone(),
            sections: vec![
                StatusSection::new(
                    "Current stage",
                    vec![StatusMetric::new("Records", self.stage_records.to_string())],
                ),
                StatusSection::new(
                    "Annotations",
                    vec![
                        StatusMetric::new("Transcripts", self.transcripts.to_string()),
                        StatusMetric::new("Proteins", self.proteins.to_string()),
                        StatusMetric::new(
                            "Chromatin elements",
                            self.chromatin_elements.to_string(),
                        ),
                        StatusMetric::new(
                            "Binding-union regions",
                            self.binding_regions.to_string(),
                        ),
                        StatusMetric::new(
                            "Source binding peaks",
                            self.source_binding_peaks.to_string(),
                        ),
                    ],
                ),
                StatusSection::new(
                    "Memory",
                    vec![
                        StatusMetric::new("Process RSS", format!("{:.1} MiB", mem.process_rss_mib)),
                        StatusMetric::new(
                            "Peak RSS",
                            format!("{:.1} MiB", mem.process_peak_rss_mib),
                        ),
                        StatusMetric::new(
                            "System available",
                            format!("{:.1} MiB", mem.system_available_mib),
                        ),
                    ],
                ),
            ],
        }
    }
}

fn write_build_summary(
    path: &std::path::Path,
    model: &Ommverse,
    started_unix_ms: u128,
    finished_unix_ms: u128,
) -> Result<()> {
    let summary = format!(
        concat!(
            "assembly: {}\n",
            "started_unix_ms: {}\n",
            "finished_unix_ms: {}\n",
            "source_root: {}\n",
            "output: {}\n",
            "transcripts: {}\n",
            "proteins: {}\n",
            "linked_transcripts: {}\n",
            "protein_features: {}\n",
            "chromatin_elements: {}\n",
            "protein_binding_union_regions: {}\n",
            "source_binding_peaks: {}\n",
            "unlinked_mappings: {}\n",
            "orphan_features: {}\n"
        ),
        model.assembly,
        started_unix_ms,
        finished_unix_ms,
        model.source_root.display(),
        path.display(),
        model.report.transcripts,
        model.report.mapped_proteins,
        model.report.linked_transcripts,
        model.report.feature_records,
        model.chromatin.len(),
        model.protein_binding.regions.len(),
        model.protein_binding.source_peak_count,
        model.report.unlinked_mapping_records,
        model.report.feature_records_without_protein
    );
    let summary_path = path.with_extension("ommverse-build-summary.yaml");
    std::fs::write(&summary_path, summary)
        .with_context(|| format!("writing {}", summary_path.display()))?;
    println!("  build summary:           {}", summary_path.display());
    Ok(())
}

#[derive(Debug, Clone)]
struct InterProStatus {
    started_unix_ms: u128,
    finished_unix_ms: Option<u128>,
    stage: String,
    public_url: Option<String>,
    indices: Vec<String>,
    streamed: usize,
    records_per_second: f64,
    reports: Vec<InterProImportReport>,
}
impl InterProStatus {
    fn new(indices: Vec<String>) -> Self {
        Self {
            started_unix_ms: SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap_or_default()
                .as_millis(),
            finished_unix_ms: None,
            stage: "starting".to_owned(),
            public_url: None,
            indices,
            streamed: 0,
            records_per_second: 0.0,
            reports: Vec::new(),
        }
    }
    fn stage(&mut self, stage: &str) {
        self.stage = stage.to_owned();
    }
    fn update(&mut self, streamed: usize, reports: &[InterProImportReport], elapsed: f64) {
        self.streamed = streamed;
        self.records_per_second = if elapsed > 0.0 {
            streamed as f64 / elapsed
        } else {
            0.0
        };
        self.reports = reports.to_vec();
    }
    fn finish(&mut self) {
        self.stage = "complete".to_owned();
        self.finished_unix_ms = Some(
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap_or_default()
                .as_millis(),
        );
    }
}
impl ServerContent for InterProStatus {
    fn server_snapshot(&self) -> ServerSnapshot {
        let mut sections = vec![StatusSection::new(
            "InterPro stream",
            vec![
                StatusMetric::new("Records processed", self.streamed.to_string()),
                StatusMetric::new(
                    "Average records / second",
                    format!("{:.0}", self.records_per_second),
                ),
                StatusMetric::new("Indices", self.indices.len().to_string()),
            ],
        )];
        for (i, report) in self.reports.iter().enumerate() {
            let label = self.indices.get(i).map(String::as_str).unwrap_or("index");
            sections.push(StatusSection::new(
                label,
                vec![
                    StatusMetric::new("Matching records", report.matched_records.to_string()),
                    StatusMetric::new("Features added", report.features_added.to_string()),
                    StatusMetric::new("Duplicates", report.duplicate_features.to_string()),
                    StatusMetric::new("Annotations collapsed", report.annotations_collapsed.to_string()),
                    StatusMetric::new("Unknown entries", report.unknown_entries.to_string()),
                ],
            ));
        }
        ServerSnapshot {
            title: "Ommverse InterPro import".to_owned(),
            subtitle: "One-pass multi-index molecular annotation".to_owned(),
            started_unix_ms: self.started_unix_ms,
            finished_unix_ms: self.finished_unix_ms,
            stage: self.stage.clone(),
            public_url: self.public_url.clone(),
            sections,
        }
    }
}

fn discover_interpro_inputs(
    index: Option<&std::path::Path>,
    dir: Option<&std::path::Path>,
) -> Result<Vec<PathBuf>> {
    if let Some(index) = index {
        return Ok(vec![index.to_path_buf()]);
    }
    let dir = dir.context("provide either --index or --index-dir")?;
    let mut found = Vec::new();
    fn walk(dir: &std::path::Path, found: &mut Vec<PathBuf>) -> Result<()> {
        for entry in std::fs::read_dir(dir).with_context(|| format!("reading {}", dir.display()))? {
            let path = entry?.path();
            if path.is_dir() {
                walk(&path, found)?;
            } else if path.extension().is_some_and(|x| x == "ommverse")
                && !path
                    .file_name()
                    .and_then(|x| x.to_str())
                    .is_some_and(|x| x.ends_with(".interpro.ommverse"))
            {
                found.push(path);
            }
        }
        Ok(())
    }
    walk(dir, &mut found)?;
    found.sort();
    Ok(found)
}
fn interpro_output_path(input: &std::path::Path) -> Result<PathBuf> {
    let name = input
        .file_name()
        .and_then(|x| x.to_str())
        .context("Ommverse index has no filename")?;
    let stem = name
        .strip_suffix(".ommverse")
        .context("Ommverse index does not end in .ommverse")?;
    Ok(input.with_file_name(format!("{stem}.interpro.ommverse")))
}

fn validate_readable_file(path: &std::path::Path, label: &str) -> Result<()> {
    let metadata = std::fs::metadata(path).with_context(|| {
        format!(
            "{label} does not exist or is not accessible: {}",
            path.display()
        )
    })?;
    if !metadata.is_file() {
        anyhow::bail!("{label} is not a regular file: {}", path.display());
    }
    std::fs::File::open(path)
        .with_context(|| format!("{label} is not readable: {}", path.display()))?;
    Ok(())
}

#[derive(Subcommand)]
enum PlotCommand {
    /// Plot genomic transcript structures, CDSs and InterPro architecture for one or more genes.
    Genome {
        #[arg(long)]
        index: PathBuf,
        /// Output directory or filename prefix.
        #[arg(long)]
        out: PathBuf,
        /// Gene symbols to plot. The Ommverse reference is loaded only once.
        #[arg(required = true, num_args = 1..)]
        genes: Vec<String>,
    },
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Build {
            reference,
            assembly,
            cache,
            gtf,
            gencode_release,
            out,
            debug_failed_mappings,
            health_port,
            health_hostname,
            no_health_server,
        } => {
            if let Some(gtf) = &gtf {
                validate_readable_file(gtf, "GTF")?;
                if assembly.is_none() {
                    anyhow::bail!("--gtf requires --assembly so the remaining reference sources can be resolved");
                }
            }
            let assembly_name = assembly
                .clone()
                .or_else(|| {
                    reference
                        .as_ref()
                        .and_then(|p| p.file_name())
                        .and_then(|x| x.to_str())
                        .map(str::to_owned)
                })
                .unwrap_or_else(|| "unknown".to_owned());
            let state = Arc::new(RwLock::new(BuildStatus::new(assembly_name)));
            let _server = if no_health_server {
                None
            } else {
                let server = spawn_status_server(
                    state.clone(),
                    SocketAddr::new(IpAddr::V4(Ipv4Addr::UNSPECIFIED), health_port),
                )?;
                let url = format!(
                    "http://{}:{}",
                    public_hostname(health_hostname.as_deref()),
                    server.addr().port()
                );
                if let Ok(mut s) = state.write() {
                    s.public_url = Some(url.clone());
                }
                eprintln!("[ommverse] health server: {url}");
                Some(server)
            };
            let started_unix_ms = state.read().map(|s| s.started_unix_ms).unwrap_or_default();
            let update = |stage: &str, records: usize| {
                if let Ok(mut s) = state.write() {
                    if s.stage != stage {
                        eprintln!("[ommverse] stage: {stage}");
                    }
                    s.stage(stage, records);
                }
            };
            let out = match (&reference, &assembly, out) {
                (_, Some(assembly), None) => cache.join(assembly).join(format!("{assembly}.ommverse")),
                (_, _, Some(out)) => out,
                (Some(_), None, None) => anyhow::bail!("Build with --reference requires --out"),
                (None, None, None) => anyhow::bail!("Build requires either --reference or --assembly"),
            };
            let model = match (reference, assembly) {
                (Some(reference), None) => Ommverse::build_ucsc_with_debug_and_progress(
                    &reference,
                    debug_failed_mappings,
                    update,
                )?,
                (None, Some(assembly)) => {
                    let annotation = ommverse::sources::AnnotationOptions {
                        gencode_release,
                        explicit_gtf: gtf,
                    };
                    Ommverse::fetch_and_build_with_annotation_and_progress(
                        &assembly, &cache, annotation, update,
                    )?
                }
                (None, None) => anyhow::bail!("Build requires either --reference or --assembly"),
                (Some(_), Some(_)) => unreachable!("clap enforces conflicts"),
            };
            if let Ok(mut s) = state.write() {
                s.model(&model);
                s.stage("serializing Ommverse", model.protein_binding.regions.len());
            }
            eprintln!("[ommverse] stage: serializing Ommverse");
            model.save(&out)?;
            let finished_unix_ms = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap_or_default()
                .as_millis();
            write_build_summary(&out, &model, started_unix_ms, finished_unix_ms)?;
            if let Ok(mut s) = state.write() {
                s.finish();
            }
            println!("Ommverse {}", model.assembly);
            println!("  transcripts:             {}", model.report.transcripts);
            println!(
                "  proteins:                {}",
                model.report.mapped_proteins
            );
            println!(
                "  linked transcripts:      {}",
                model.report.linked_transcripts
            );
            println!(
                "  protein features:        {}",
                model.report.feature_records
            );
            println!("  chromatin elements:      {}", model.chromatin.len());
            println!(
                "  binding-union regions:   {}",
                model.protein_binding.regions.len()
            );
            println!(
                "  source binding peaks:    {}",
                model.protein_binding.source_peak_count
            );
            println!(
                "  loop experiments:        {}",
                model.experimental_loops.experiment_count
            );
            println!(
                "  recurrent loop anchors:  {}",
                model.experimental_loops.anchors.len()
            );
            println!(
                "  distinct observed loops: {}",
                model.experimental_loops.loops.len()
            );
            println!(
                "  unlinked mappings:       {}",
                model.report.unlinked_mapping_records
            );
            println!(
                "  orphan features:         {}",
                model.report.feature_records_without_protein
            );
            println!(
                "  transmembrane proteins:  {}",
                model
                    .proteins_with_feature(ProteinFeatureKind::Transmembrane)
                    .count()
            );
            for w in &model.report.warnings {
                eprintln!("warning: {w}");
            }
            println!("  saved: {}", out.display());
        }
        Command::IngestInterpro {
            index,
            index_dir,
            protein2ipr,
            entry_list,
            parent_child_tree,
            out,
            health_port,
            health_hostname,
            no_health_server,
        } => {
            let inputs = discover_interpro_inputs(index.as_deref(), index_dir.as_deref())?;
            if inputs.is_empty() {
                anyhow::bail!("no base .ommverse indices found");
            }
            if index.is_some() && out.is_none() {
                anyhow::bail!("--out is required with --index");
            }
            let outputs = if let Some(single_out) = out {
                vec![single_out]
            } else {
                inputs
                    .iter()
                    .map(|path| interpro_output_path(path.as_path()))
                    .collect::<Result<Vec<_>>>()?
            };

            // Start observability before any potentially expensive validation or model loading.
            let state = Arc::new(RwLock::new(InterProStatus::new(
                inputs.iter().map(|p| p.display().to_string()).collect(),
            )));
            let _server = if no_health_server {
                None
            } else {
                let server = spawn_status_server(
                    state.clone(),
                    SocketAddr::new(IpAddr::V4(Ipv4Addr::UNSPECIFIED), health_port),
                )?;
                let url = format!(
                    "http://{}:{}",
                    public_hostname(health_hostname.as_deref()),
                    server.addr().port()
                );
                if let Ok(mut s) = state.write() {
                    s.public_url = Some(url.clone());
                }
                eprintln!("[ommverse] health server: {url}");
                Some(server)
            };

            if let Ok(mut s) = state.write() {
                s.stage("validating inputs");
            }
            eprintln!("[ommverse] stage: validating inputs");
            validate_readable_file(&protein2ipr, "protein2ipr")?;
            validate_readable_file(&entry_list, "entry list")?;
            if let Some(tree) = &parent_child_tree {
                validate_readable_file(tree, "parent/child tree")?;
            }
            for path in &inputs {
                validate_readable_file(path, "Ommverse index")?;
            }

            println!("Ommverse InterPro import");
            println!("  indices:                  {}", inputs.len());
            println!("  protein2ipr:              {}", protein2ipr.display());
            println!("  entry list:               {}", entry_list.display());
            if let Some(tree) = &parent_child_tree {
                println!("  parent/child tree:        {}", tree.display());
            }

            if let Ok(mut s) = state.write() {
                s.stage("loading existing Ommverse models");
            }
            eprintln!("[ommverse] stage: loading existing Ommverse models");
            let mut ommverses = Vec::with_capacity(inputs.len());
            for path in &inputs {
                eprintln!("[ommverse] loading: {}", path.display());
                ommverses.push(Ommverse::load(path)?);
            }
            if let Ok(mut s) = state.write() {
                s.stage("streaming InterPro");
            }
            eprintln!("[ommverse] stage: streaming InterPro");
            let started = Instant::now();
            let reference_cache = if inputs.len() == 1 {
                let parent = outputs[0].parent().unwrap_or_else(|| std::path::Path::new("."));
                let cache = parent.join("interpro.reference.tsv.gz");
                let manifest = parent.join("interpro.sources.txt");
                if !manifest.is_file() {
                    let mut text = format!(
                        "protein2ipr={}\nentry_list={}\nreference_cache={}\n",
                        protein2ipr.display(),
                        entry_list.display(),
                        cache.display()
                    );
                    if let Some(tree) = &parent_child_tree {
                        text.push_str(&format!("parent_child_tree={}\n", tree.display()));
                    }
                    std::fs::write(&manifest, text)?;
                }
                if cache.is_file() {
                    eprintln!("[ommverse] using reference InterPro cache: {}", cache.display());
                } else {
                    eprintln!("[ommverse] reference InterPro cache will be written: {}", cache.display());
                }
                Some(cache)
            } else {
                None
            };
            let feature_source_register = if inputs.len() == 1 {
                let parent = outputs[0].parent().unwrap_or_else(|| std::path::Path::new("."));
                Some(parent.join("interpro.feature_sources.tsv"))
            } else {
                None
            };
            let reports = ingest_interpro_many(
                &mut ommverses,
                &protein2ipr,
                &entry_list,
                parent_child_tree.as_deref(),
                reference_cache.as_deref(),
                feature_source_register.as_deref(),
                |streamed, reports| {
                    if let Ok(mut s) = state.write() {
                        s.update(streamed, reports, started.elapsed().as_secs_f64());
                    }
                },
            )?;
            if let Ok(mut s) = state.write() {
                s.stage("saving enriched indices");
            }
            eprintln!("[ommverse] stage: saving enriched indices");
            for ((omm, output), report) in ommverses.iter().zip(&outputs).zip(&reports) {
                omm.save(output)?;
                println!("  {}", omm.assembly);
                println!("    records streamed:       {}", report.records_streamed);
                println!("    matching records:       {}", report.matched_records);
                println!("    features added:         {}", report.features_added);
                println!("    genomic projections:    {}", report.projected_features_added);
                println!("    duplicate features:     {}", report.duplicate_features);
                println!("    annotations collapsed:  {}", report.annotations_collapsed);
                println!("    annotations sanitized:  {}", report.annotations_sanitized);
                println!("    malformed records:      {}", report.malformed_records);
                println!("    unknown entries:        {}", report.unknown_entries);
                println!("    saved:                  {}", output.display());
            }
            if let Ok(mut s) = state.write() {
                s.finish();
            }
        }
        Command::Plot { plot } => match plot {
            PlotCommand::Genome { index, out, genes } => {
                eprintln!("[ommverse] loading reference for genome plots...");
                let omm = Ommverse::load(&index)?;
                eprintln!("[ommverse] reference ready; plotting {} gene(s)...", genes.len());
                for (i, gene) in genes.iter().enumerate() {
                    let plot_out = ommverse::plot::gene_plot_output_path(&out, gene);
                    if let Some(parent) = plot_out.parent() {
                        if !parent.as_os_str().is_empty() {
                            std::fs::create_dir_all(parent)?;
                        }
                    }
                    eprintln!(
                        "[ommverse] [{}/{}] {} -> {}",
                        i + 1,
                        genes.len(),
                        gene,
                        plot_out.display()
                    );
                    ommverse::plot::plot_gene_transcripts(&omm, gene, &plot_out)?;
                }
                eprintln!("[ommverse] wrote {} genome plot(s)", genes.len());
            }
        },
        Command::IngestMsigdb { index, gmt, release, cache_dir } => {
            let paths = if gmt.is_empty() {
                let cache = cache_dir.unwrap_or_else(|| index.parent().unwrap_or(std::path::Path::new(".")).join("MSigDB"));
                ommverse::msigdb::download_collections(&cache, &release)?
            } else {
                gmt
            };
            let count = ommverse::msigdb::ingest(&index, &paths, &release)?;
            println!("[ommverse] indexed {count} MSigDB gene sets in {}", ommverse::msigdb::sidecar(&index).display());
        }
        Command::MsigdbGene { index, gene } => {
            for row in ommverse::msigdb::find_gene(&index, &gene)? {
                println!("{}\t{}\t{}\t{}", row.name, row.description, row.source, row.release);
            }
        }
        Command::MsigdbSet { index, set } => {
            for row in ommverse::msigdb::find_set(&index, &set)? {
                println!("{}\t{}\t{}\t{}\t{}", row.name, row.description, row.source, row.release, row.genes.join(","));
            }
        }
        Command::Query { index, query } => {
            let query = ommverse::query::Query::parse(&query)?;
            let schema = ommverse::schema::OmmverseSchema::load_for_index(&index)?;
            schema.validate_query(&query)?;
            let omm = Ommverse::load(index)?;
            for row in ommverse::query::execute(&omm, &query)? {
                println!("{row}");
            }
        }
        Command::Schema { index, write } => {
            let schema = if write {
                let omm = Ommverse::load(&index)?;
                let schema = ommverse::schema::OmmverseSchema::for_ommverse(&omm);
                let path = schema.save_for_index(&index)?;
                eprintln!("[ommverse] wrote schema: {}", path.display());
                schema
            } else {
                ommverse::schema::OmmverseSchema::load_for_index(&index)?
            };
            print!("{}", serde_yaml::to_string(&schema)?);
        }
        Command::Inspect { index } => {
            let started = Instant::now();
            let omm = Ommverse::load(&index)?;
            let load_seconds = started.elapsed().as_secs_f64();
            print_inspection(&index, &omm, load_seconds)?;
        }
        Command::StringProbe { index, info, links, out } => {
            ommverse::string_probe::probe_string(&index, &info, &links, &out)?;
        }
        Command::Protein {
            index,
            accession,
            sequence,
        } => {
            let omm = Ommverse::load(index)?;
            let p = omm
                .protein(&accession)
                .with_context(|| format!("protein {accession} not found"))?;
            println!("{}  {}", p.accession, p.name);
            println!("  gene: {}", p.gene_symbol);
            println!("  entry: {}", p.entry_name);
            println!("  transcripts: {:?}", p.transcript_ids);
            println!("  features: {}", p.features.len());
            for f in &p.features {
                println!("    {:?}\t{:?}\t{}", f.kind, f.protein_range, f.description);
            }
            if sequence {
                println!("  sequence: {}", omm.protein_sequence(&accession)?);
            }
        }
        Command::Gene { index, symbol } => {
            let omm = Ommverse::load(index)?;
            let proteins: Vec<_> = omm.proteins_for_gene(&symbol).collect();
            if proteins.is_empty() {
                anyhow::bail!("gene {symbol} not found");
            }
            println!("{}", symbol);
            for p in proteins {
                println!(
                    "  {}\t{}\t{} features",
                    p.accession,
                    p.name,
                    p.features.len()
                );
            }
        }
        Command::TrainingCorpus {
            index,
            feature,
            flank_width,
        } => {
            let omm = Ommverse::load(index)?;
            let feature = ProteinFeatureKind::from_str(&feature)?;
            let corpus = omm.training_corpus_with_flank_diagnostic(feature, flank_width)?;
            let r = &corpus.report;
            println!("Ommverse supervised protein-feature corpus");
            println!("  feature:                   {:?}", feature);
            println!("  sequence:                  complete reconstructed proteins");
            println!("  labels:                    curated Ommverse feature intervals");
            println!("  background:                all residues outside the feature");
            println!("  split:                     deterministic protein-level 50:50");
            println!("  model training:            none (corpus validation only)");
            println!("  candidate proteins:        {}", r.candidate_proteins);
            println!("  reconstructed:             {}", r.reconstructed_proteins);
            println!("  rejected sequence:         {}", r.rejected_sequence);
            println!("  rejected missing range:    {}", r.rejected_missing_range);
            println!("  rejected out of bounds:    {}", r.rejected_out_of_bounds);
            println!("  training proteins:         {}", corpus.train.len());
            println!("  test proteins:             {}", corpus.test.len());
            println!("  total residues:            {}", r.total_residues);
            println!("  feature residues:          {}", r.feature_residues);
            println!("  feature segments:          {}", r.feature_segments);
            println!();
            println!("Feature geometry");
            print_distribution("feature length", &r.feature_lengths);
            println!();
            println!("Adjacent feature geometry");
            print_distribution("inter-feature gap", &r.inter_feature_gaps);
            println!(
                "  flank diagnostic width:    {} aa per side",
                r.diagnostic_flank_width
            );
            println!(
                "  overlapping flank pairs:   {} / {} ({:.2}%)",
                r.overlapping_flank_pairs,
                r.inter_feature_gaps.count,
                percent(r.overlapping_flank_pairs, r.inter_feature_gaps.count)
            );
        }
        Command::TrainExactAa {
            index,
            feature,
            out,
            flank_width,
            scan_unannotated,
        } => {
            let omm = Ommverse::load(index)?;
            let feature = ProteinFeatureKind::from_str(&feature)?;
            println!("Ommverse protein-feature training");
            println!("  feature:          {:?}", feature);
            println!("  observation:      exact amino-acid identity");
            println!("  states:           BACKGROUND / PRE_FEATURE / FEATURE / POST_FEATURE");
            println!("  flank width:      {} aa", flank_width);
            println!("  supervision:      curated Ommverse protein-feature coordinates");
            println!("  split:            deterministic protein-level 50:50");
            println!("  held-out labels:  never used for fitting");
            println!(
                "  unannotated scan: {}",
                if scan_unannotated {
                    "enabled; inference only"
                } else {
                    "disabled"
                }
            );
            println!();
            println!("Preparing supervised corpus...");
            let corpus = omm.training_corpus(feature)?;
            println!(
                "  accepted proteins: {} (train {}, test {})",
                corpus.train.len() + corpus.test.len(),
                corpus.train.len(),
                corpus.test.len()
            );
            println!("Training four-state exact-AA HMM...");
            let model = omm.train_exact_aa_feature_model(&corpus, flank_width)?;
            let evaluation = omm.evaluate_exact_aa_feature_model(&model, &corpus.test)?;
            model.save(&out)?;
            println!("Evaluating held-out proteins...");
            println!("Exact-AA {:?} model", feature);
            println!("  training proteins:       {}", corpus.train.len());
            println!("  held-out proteins:       {}", corpus.test.len());
            println!("  held-out residues:       {}", evaluation.residues);
            println!("  precision:               {:.4}", evaluation.precision());
            println!("  recall:                  {:.4}", evaluation.recall());
            println!("  specificity:             {:.4}", evaluation.specificity());
            println!("  F1:                      {:.4}", evaluation.f1());
            println!("  truth segments:          {}", evaluation.truth_segments);
            println!(
                "  predicted segments:      {}",
                evaluation.predicted_segments
            );
            println!(
                "  segments recovered >=50%:{} ({:.4})",
                evaluation.recovered_segments,
                evaluation.segment_recall()
            );
            println!(
                "  mean truth segment len:  {:.2} aa",
                ratio_usize(
                    evaluation.true_positive + evaluation.false_negative,
                    evaluation.truth_segments
                )
            );
            println!(
                "  mean predicted seg len:  {:.2} aa",
                ratio_usize(
                    evaluation.true_positive + evaluation.false_positive,
                    evaluation.predicted_segments
                )
            );
            println!("  initial:                 {:?}", model.initial);
            println!(
                "  transition rows:         BACKGROUND / PRE_FEATURE / FEATURE / POST_FEATURE"
            );
            for row in 0..4 {
                println!("    {:?}", &model.transition[row * 4..row * 4 + 4]);
            }
            println!("  learned AA emissions (BACKGROUND / PRE_FEATURE / FEATURE / POST_FEATURE):");
            for residue in b"ACDEFGHIKLMNPQRSTVWY" {
                // Canonical IntToProt codes are not alphabetical; ask IntToProt for the code.
                let aa = int_to_prot::IntToProt::new([*residue]).get(0).unwrap();
                let c = aa.code() as usize;
                println!(
                    "    {}  {:.6}  {:.6}  {:.6}  {:.6}",
                    *residue as char,
                    model.emission[0][c],
                    model.emission[1][c],
                    model.emission[2][c],
                    model.emission[3][c]
                );
            }
            println!("  saved: {}", out.display());
            if scan_unannotated {
                print_unannotated(&omm.scan_unannotated_exact_aa(&model)?);
            }
        }
        Command::ModelVault { command: ModelVaultCommand::Train {
            index,
            out,
            flank_width,
            interpro_max_len,
        } } => {
            let omm = Ommverse::load(&index)?;
            println!("Ommverse AA model-vault training");
            println!("  assembly:          {}", omm.assembly);
            println!("  observations:      exact-AA + intrinsic chemistry");
            println!(
                "  feature classes:   independent feature models + joint topology competition"
            );
            println!("  split:             deterministic protein-level 50:50 per feature");
            println!("  flank width:       {} aa", flank_width);
            println!("  held-out labels:   never used for fitting");
            println!("  InterPro max len:  {} aa", interpro_max_len);
            println!();
            let (vault, report) = omm.train_aa_model_vault(flank_width, interpro_max_len)?;
            vault.save(&out)?;
            println!("Model vault complete");
            println!(
                "  feature classes considered: {}",
                report.feature_classes_considered
            );
            println!(
                "  feature classes trained:    {}",
                report.feature_classes_trained
            );
            println!("  models trained:             {}", report.models_trained);
            println!("  models stored:              {}", vault.models.len() + vault.interpro_models.len());
            println!("  InterPro considered:        {}", report.interpro_features_considered);
            println!("  InterPro passed validation: {}", report.interpro_features_trained);
            println!(
                "  topology train proteins:    {}",
                report.topology_train_proteins
            );
            println!(
                "  topology test proteins:     {}",
                report.topology_test_proteins
            );
            println!(
                "  topology conflicts:         {} residues",
                report.topology_conflicting_residues
            );
            for (feature, reason) in &report.skipped {
                println!("  skipped {:?}: {}", feature, reason);
            }
            for (feature, reason) in &report.interpro_skipped {
                println!("  skipped InterPro {}: {}", feature, reason);
            }
            println!("  saved: {}", out.display());
            let info_path = PathBuf::from(format!("{}.info", out.display()));
            std::fs::write(&info_path, format_model_vault_info(&vault))?;
            println!("  info:  {}", info_path.display());
        }
        Command::ModelVault { command: ModelVaultCommand::Info { vault } } => {
            let vault = AaModelVault::load(vault)?;
            print!("{}", format_model_vault_info(&vault));
        }
        Command::ModelVault { command: ModelVaultCommand::Predict { vault, fasta, out } } => {
            let vault = AaModelVault::load(&vault)?;
            let sequences = read_aa_fasta(&fasta)?;
            let responses_path = PathBuf::from(format!("{}.responses.tsv", out.display()));
            let features_path = PathBuf::from(format!("{}.features.tsv", out.display()));
            let mut responses = String::from("sequence_id\tmodel_id\tmodel_name\tmodel_kind\traw_log_likelihood\theldout_precision\theldout_recall\theldout_f1\n");
            let mut features = String::from("sequence_id\taa_start_0based\taa_end_0based_exclusive\tmodel_id\tmodel_name\tmodel_kind\traw_log_likelihood\tcalibrated_p_value\theldout_precision\theldout_recall\theldout_f1\n");
            for (sequence_id, sequence) in sequences {
                let (raw, intervals) = vault.predict_sequence(sequence.as_bytes())?;
                for row in raw {
                    responses.push_str(&format!("{}\t{}\t{}\t{}\t{:.12}\t{}\t{}\t{}\n",
                        sequence_id, row.model_id, row.model_name, row.model_kind, row.raw_score,
                        optional_float(row.precision), optional_float(row.recall), optional_float(row.f1)));
                }
                for row in intervals {
                    features.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{:.12}\t{}\t{}\t{}\t{}\n",
                        sequence_id, row.start, row.end, row.model_id, row.model_name, row.model_kind,
                        row.raw_score, optional_float(row.calibrated_p_value), optional_float(row.precision),
                        optional_float(row.recall), optional_float(row.f1)));
                }
            }
            std::fs::write(&responses_path, responses)?;
            std::fs::write(&features_path, features)?;
            println!("Ommverse model-vault prediction");
            println!("  FASTA:     {}", fasta.display());
            println!("  vault:     {}", vault.metadata.source);
            println!("  raw:       {}", responses_path.display());
            println!("  features:  {}", features_path.display());
            println!("  coordinates: 0-based half-open [start, end)");
            println!("  significance: not calibrated; calibrated_p_value is intentionally empty");
        }
        Command::EvaluateModelVault { index, vault } => {
            let omm = Ommverse::load(&index)?;
            let vault = AaModelVault::load(&vault)?;
            println!("Ommverse model-vault transfer evaluation");
            println!("  trained assembly:          {}", vault.metadata.assembly);
            println!("  target assembly:           {}", omm.assembly);
            println!("  models:                    {}", vault.models.len() + vault.interpro_models.len());
            println!("  fitting on target:         none");
            println!("  target labels:             evaluation only");
            println!();
            println!(
                "{:<20} {:<10} {:>10} {:>10} {:>10} {:>10} {:>10}",
                "feature", "model", "source F1", "target F1", "precision", "recall", "seg.rec"
            );
            println!("{}", "-".repeat(86));

            let mut rows: Vec<(
                ProteinFeatureKind,
                &'static str,
                ommverse::FeatureEvaluation,
                ommverse::FeatureEvaluation,
            )> = Vec::new();
            let mut topology_diagnostics: Vec<(
                &'static str,
                TopologyEvaluation,
                TopologyEvaluation,
            )> = Vec::new();
            for entry in &vault.models {
                match entry {
                    ProteinFeatureModel::ExactAa {
                        model,
                        evaluation: source,
                    } => {
                        if matches!(
                            model.feature,
                            ProteinFeatureKind::Chain
                                | ProteinFeatureKind::Conflict
                                | ProteinFeatureKind::SignalPeptide
                        ) {
                            continue;
                        }
                        let corpus = omm.training_corpus(model.feature)?;
                        let target_examples: Vec<_> = corpus
                            .train
                            .iter()
                            .chain(corpus.test.iter())
                            .cloned()
                            .collect();
                        let target =
                            omm.evaluate_exact_aa_feature_model(model, &target_examples)?;
                        rows.push((model.feature, "exact-AA", source.clone(), target));
                    }
                    ProteinFeatureModel::Chemistry {
                        model,
                        evaluation: source,
                    } => {
                        if matches!(
                            model.feature,
                            ProteinFeatureKind::Chain
                                | ProteinFeatureKind::Conflict
                                | ProteinFeatureKind::SignalPeptide
                        ) {
                            continue;
                        }
                        let corpus = omm.training_corpus(model.feature)?;
                        let target_examples: Vec<_> = corpus
                            .train
                            .iter()
                            .chain(corpus.test.iter())
                            .cloned()
                            .collect();
                        let target =
                            omm.evaluate_chemistry_feature_model(model, &target_examples)?;
                        rows.push((model.feature, "chemistry", source.clone(), target));
                    }
                    ProteinFeatureModel::TopologyExactAa {
                        model,
                        evaluation: source,
                    } => {
                        let corpus = omm.topology_training_corpus()?;
                        let target_examples: Vec<_> = corpus
                            .train
                            .iter()
                            .chain(corpus.test.iter())
                            .cloned()
                            .collect();
                        let target =
                            omm.evaluate_exact_aa_topology_model(model, &target_examples)?;
                        for ((feature, source_state), target_state) in TopologyState::BIOLOGICAL
                            .iter()
                            .zip(&source.states)
                            .zip(&target.states)
                        {
                            rows.push((
                                *feature,
                                "joint-AA",
                                source_state.clone(),
                                target_state.clone(),
                            ));
                        }
                        topology_diagnostics.push(("joint-AA", source.clone(), target));
                    }
                    ProteinFeatureModel::TopologyChemistry {
                        model,
                        evaluation: source,
                    } => {
                        let corpus = omm.topology_training_corpus()?;
                        let target_examples: Vec<_> = corpus
                            .train
                            .iter()
                            .chain(corpus.test.iter())
                            .cloned()
                            .collect();
                        let target =
                            omm.evaluate_chemistry_topology_model(model, &target_examples)?;
                        for ((feature, source_state), target_state) in TopologyState::BIOLOGICAL
                            .iter()
                            .zip(&source.states)
                            .zip(&target.states)
                        {
                            rows.push((
                                *feature,
                                "joint-chem",
                                source_state.clone(),
                                target_state.clone(),
                            ));
                        }
                        topology_diagnostics.push(("joint-chem", source.clone(), target));
                    }
                }
            }
            let model_rank = |name: &str| match name {
                "exact-AA" => 0,
                "chemistry" => 1,
                "joint-AA" => 2,
                "joint-chem" => 3,
                _ => 4,
            };
            rows.sort_by_key(|(feature, model, _, _)| {
                (format!("{:?}", feature), model_rank(model))
            });
            for (feature, model, source, target) in rows {
                print_transfer_evaluation(feature, model, &source, &target);
            }
            if !topology_diagnostics.is_empty() {
                println!();
                println!("Joint topology latent-state diagnostics");
                println!("  counts are residues before CYTO/EXTRA substates are collapsed");
                for (model, source, target) in topology_diagnostics {
                    print_latent_topology_diagnostics(model, &source, &target);
                }
            }
        }
        Command::ScanExactAa { index, model } => {
            let omm = Ommverse::load(index)?;
            let model = ExactAaFeatureModel::load(model)?;
            print_unannotated(&omm.scan_unannotated_exact_aa(&model)?);
        }
    }
    Ok(())
}

fn print_unannotated(r: &ommverse::UnannotatedScanReport) {
    println!("Unannotated-protein scan");
    println!("  candidate proteins:      {}", r.candidate_proteins);
    println!("  reconstructed:           {}", r.reconstructed_proteins);
    println!("  rejected sequence:       {}", r.rejected_sequence);
    println!(
        "  predicted proteins:      {}",
        r.predicted_feature_proteins
    );
    println!("  predicted segments:      {}", r.predicted_segments);
    println!("  predicted residues:      {}", r.predicted_residues);
}

fn ratio_usize(num: usize, den: usize) -> f64 {
    if den == 0 {
        0.0
    } else {
        num as f64 / den as f64
    }
}

fn print_topology_info(label: &str, evaluation: &TopologyEvaluation) {
    println!("    {label}:");
    for (feature, e) in TopologyState::BIOLOGICAL.iter().zip(&evaluation.states) {
        println!(
            "      {:?}: precision {:.4}, recall {:.4}, F1 {:.4}, seg.rec {:.4}",
            feature,
            e.precision(),
            e.recall(),
            e.f1(),
            e.segment_recall()
        );
    }
}

fn print_latent_topology_diagnostics(
    model: &str,
    source: &TopologyEvaluation,
    target: &TopologyEvaluation,
) {
    println!();
    println!("  {model}");
    println!(
        "    {:<16} {:>12} {:>12} {:>12} {:>12}",
        "latent state", "src truth", "src pred", "tgt truth", "tgt pred"
    );
    for state in TopologyState::ALL {
        let i = state as usize;
        println!(
            "    {:<16} {:>12} {:>12} {:>12} {:>12}",
            state.label(),
            source.latent_truth_residues[i],
            source.latent_predicted_residues[i],
            target.latent_truth_residues[i],
            target.latent_predicted_residues[i]
        );
    }
    println!("    largest target latent confusions (truth -> predicted):");
    let mut confusions = Vec::new();
    for truth in 0..TopologyState::COUNT {
        for pred in 0..TopologyState::COUNT {
            if truth == pred {
                continue;
            }
            let n = target.latent_confusion[truth * TopologyState::COUNT + pred];
            if n > 0 {
                confusions.push((n, truth, pred));
            }
        }
    }
    confusions.sort_unstable_by(|a, b| b.0.cmp(&a.0));
    for (n, truth, pred) in confusions.into_iter().take(8) {
        println!(
            "      {:<16} -> {:<16} {:>12}",
            TopologyState::ALL[truth].label(),
            TopologyState::ALL[pred].label(),
            n
        );
    }
}

fn print_transfer_evaluation(
    feature: ProteinFeatureKind,
    observation: &str,
    source: &ommverse::FeatureEvaluation,
    target: &ommverse::FeatureEvaluation,
) {
    println!(
        "{:<20} {:<10} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>10.4}",
        format!("{:?}", feature),
        observation,
        source.f1(),
        target.f1(),
        target.precision(),
        target.recall(),
        target.segment_recall(),
    );
}
