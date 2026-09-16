//clonomap_batch.rs
use clap::Parser;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::Instant;

use clonomap::ClonoMap;

// ============================================================
// CLI
// ============================================================
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Batch runner for ClonoMap on Change-O clone files"
)]
struct Args {
    /// Input directory or single Change-O TSV file
    #[arg(short, long)]
    input: PathBuf,

    /// Filename pattern to match (substring match)
    #[arg(short, long, default_value = "ProductiveCloneDfined.tsv")]
    pattern: String,

    /// Minimum number of sequences per clone
    #[arg(short = 'n', long, default_value_t = 200)]
    min_size: usize,

    /// Number of PCA components
    #[arg(short, long, default_value_t = 30)]
    k: usize,

    /// Output directory
    #[arg(short, long)]
    outdir: Option<PathBuf>,

    /// Enable SVG plotting
    #[arg(long)]
    plots: bool,

    /// Number of threads for Rayon
    #[arg(long)]
    threads: Option<usize>,

    /// Make analysis amino acis based (give me the frame starting position (0 ,1 or 2)
    #[arg(long)]
    aa: Option<usize>,
}

// ============================================================
// Main
// ============================================================
fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to configure Rayon thread pool");
    }

    // --- determine output directory ---
    let outdir = match &args.outdir {
        Some(p) => p.clone(),
        None => {
            if args.input.is_dir() {
                args.input.clone()
            } else {
                args.input
                    .parent()
                    .unwrap_or_else(|| Path::new("."))
                    .to_path_buf()
            }
        }
    };

    fs::create_dir_all(&outdir)?;

    let files = collect_input_files(&args.input, &args.pattern)?;
    if files.is_empty() {
        anyhow::bail!("No input files found matching pattern");
    }

    let total_start = Instant::now();

    for file in files {
        process_file(&file, &args, &outdir)?;
    }

    println!(
        "\n============================================================\n\
         TOTAL RUNTIME: {:.2} seconds\n\
         ============================================================",
        total_start.elapsed().as_secs_f64()
    );

    Ok(())
}

// ============================================================
// File discovery
// ============================================================
fn collect_input_files(input: &Path, pattern: &str) -> anyhow::Result<Vec<PathBuf>> {
    if input.is_file() {
        return Ok(vec![input.to_path_buf()]);
    }

    let mut out = Vec::new();
    for entry in fs::read_dir(input)? {
        let p = entry?.path();
        if p.is_file() {
            if let Some(name) = p.file_name().and_then(|s| s.to_str()) {
                if name.contains(pattern) {
                    out.push(p);
                }
            }
        }
    }
    Ok(out)
}

// ============================================================
// Process one Change-O TSV file
// ============================================================
fn process_file(path: &Path, args: &Args, outdir: &PathBuf) -> anyhow::Result<()> {
    println!("\n------------------------------------------------------------");
    println!("Processing file: {}", path.display());

    let start = Instant::now();

    let clones = read_changeo_clones(path)?;
    let stem = path.file_stem().unwrap().to_string_lossy();

    let out_base = outdir.join(stem.as_ref());
    fs::create_dir_all(&out_base)?;

    for (clone_id, seqs) in clones {
        if seqs.len() < args.min_size {
            continue;
        }

        let clone_start = Instant::now();

        let clone_dir = out_base.join(format!("clone_{}", clone_id));
        fs::create_dir_all(&clone_dir)?;

        match run_clonomap_clone(&seqs, &clone_dir, args) {
            Ok(_) => {
                println!(
                    "   - clone {} (n={}) finished in {:.2} sec",
                    clone_id,
                    seqs.len(),
                    clone_start.elapsed().as_secs_f64()
                );
            }
            Err(e) => {
                eprintln!("❌ Clone {} failed: {}", clone_id, e);
            }
        }
    }

    println!("File completed in {:.2} sec", start.elapsed().as_secs_f64());

    Ok(())
}

// ============================================================
// Change-O TSV parsing
// ============================================================
fn read_changeo_clones(path: &Path) -> anyhow::Result<HashMap<String, Vec<String>>> {
    let f = File::open(path)?;
    let mut rdr = BufReader::new(f);

    let mut header = String::new();
    rdr.read_line(&mut header)?;
    let cols: Vec<&str> = header.trim_end().split('\t').collect();

    let clone_idx = cols
        .iter()
        .position(|c| *c == "clone_id")
        .ok_or_else(|| anyhow::anyhow!("Missing clone_id column"))?;

    let seq_idx = cols
        .iter()
        .position(|c| *c == "sequence_alignment")
        .ok_or_else(|| anyhow::anyhow!("Missing sequence_alignment column"))?;

    let mut clones: HashMap<String, Vec<String>> = HashMap::new();

    for line in rdr.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= seq_idx {
            continue;
        }

        clones
            .entry(fields[clone_idx].to_string())
            .or_default()
            .push(fields[seq_idx].to_string());
    }

    Ok(clones)
}

// ============================================================
// Per-clone execution (core logic reuse)
// ============================================================
fn run_clonomap_clone(seqs: &[String], outdir: &Path, args: &Args) -> anyhow::Result<(), String> {
    // === Core computation ===
    let model = match ClonoMap::new(seqs.to_vec(), args.k, false) {
        Ok(m) => m,
        Err(e) => {
            let msg = format!("ClonoMap failed: {}", e);
            return Err(msg.into());
        }
    };

    // === Outputs ===
    model
        .pca
        .to_tsv(&model.encoder.sequences, &outdir.join("coords.tsv"))
        .map_err(|e| format!("Failed to write PCA coords: {}", e))?;

    model
        .tree
        .to_tsv(&outdir.join("tree.tsv"))
        .map_err(|e| format!("Failed to write MST tree: {}", e))?;

    model
        .encoder
        .sequences
        .to_tsv(&outdir.join("rows.tsv"))
        .map_err(|e| format!("Failed to write AA rows: {}", e))?;

    // === Optional plots ===
        if args.plots {
        model
            .pca
            .plot_2d_clusters(&model.tree, outdir.join("pca.svg").to_str().unwrap())
            .map_err(|e| format!("Failed to write PCA plot: {}", e))?;

        model
            .tree
            .plot_2d(model.coords(), outdir.join("tree.svg").to_str().unwrap())
            .map_err(|e| format!("Failed to write tree plot: {}", e))?;
    }

    // === Rooting + Newick ===
    let root = model.pca.find_sparse_root(5);

    let rooted = model.tree.reroot(model.pca.coords.nrows(), root);

    rooted
        .to_tsv(&outdir.join("rooted_tree.tsv"))
        .map_err(|e| format!("Failed to write rooted tree: {}", e))?;

    let newick = rooted.to_newick(model.pca.coords.nrows(), root, &model.encoder.sequences);

    std::fs::write(outdir.join("tree.newick"), newick)
        .map_err(|e| format!("Failed to write Newick file: {}", e))?;

    Ok(())
}
