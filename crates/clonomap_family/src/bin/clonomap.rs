use clap::Parser;
use clonomap::ClonoMap;
use std::fs::read_to_string;
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(author, version, about)]
struct Args {
    /// Path to file with one DNA sequence per line
    input: String,

    /// Number of PCA components
    #[arg(short, long, default_value_t = 3)]
    k: usize,

    /// Output table with PCA coordinates
    #[arg(long)]
    coords: Option<String>,

    /// Output table with MST edges
    #[arg(long)]
    edges: Option<String>,

    /// amino acid representations + counts
    #[arg(long)]
    rows: Option<String>,

    /// Write PCA plot (SVG)
    #[arg(long)]
        plot_pca: Option<String>,

    /// Write PCA tree plot (SVG)
    #[arg(long)]
        plot_tree: Option<String>,

    /// Make analysis amino acis based (give me the frame starting position (0 ,1 or 2)
    #[arg(long)]
    aa: Option<usize>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse(); // ✅ parse ONCE

    let input = read_to_string(&args.input)?;

    let seqs: Vec<String> = match args.aa {
        Some(skip) if skip <= 2 => input
            .lines()
            .filter(|l| !l.trim().is_empty())
            .filter_map(|s| {
                if s.len() <= skip {
                    None
                } else {
                    Some(s[skip..].to_string())
                }
            })
            .collect(),

        Some(skip) => {
            return Err(format!("Invalid aa frame {} (must be 0, 1, or 2)", skip).into());
        }

        None => input
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|s| s.to_string())
            .collect(),
    };

    // ✅ NEW: constructor does all work
    let model = match args.aa {
        Some(_) => ClonoMap::new(seqs, args.k, true)?,
        None => ClonoMap::new(seqs, args.k, false)?,
    };

    // Determine PCA output path
    let coords_path: PathBuf = if let Some(user) = args.coords.as_ref() {
        PathBuf::from(user)
    } else {
        default_output_path(&args.input, "_pca.tsv")
    };

    // Determine MST edge output path
    let edges_path: PathBuf = if let Some(user) = args.edges.as_ref() {
        PathBuf::from(user)
    } else {
        default_output_path(&args.input, "_tree.tsv")
    };

    let rows_path: PathBuf = if let Some(user) = args.rows.as_ref() {
        PathBuf::from(user)
    } else {
        default_output_path(&args.input, "_rownames.tsv")
    };

    model.pca.to_tsv(&model.encoder.sequences, &coords_path)?;
    println!("Written PCA coords → {}", coords_path.display());

    model.tree.to_tsv(&edges_path)?;
    println!("Written MSt edges → {}", edges_path.display());

    model.encoder.sequences.to_tsv(&rows_path)?;
    println!("Written Amino Acid info → {}", rows_path.display());

        {
        if let Some(f) = args.plot_pca {
            model.pca.plot_2d_clusters(&model.tree, &f)?;
            eprintln!("✅ PCA plot written to {f}");
        }

        if let Some(f) = args.plot_tree {
            model.tree.plot_2d(&model.coords(), &f)?;
            eprintln!("✅ Tree plot written to {f}");
        }
    }

    // === Determine phylogenetic root from PCA trajectory ===
    let root = model.pca.find_sparse_root(5);
    println!("Ancestral root inferred at node {}", root);

    // === Produce rooted tree ===
    let rooted_tree = model.tree.reroot(model.pca.coords.nrows(), root);

    // === Write rooted tree to TSV ===
    let rooted_path = default_output_path(&args.input, "_rooted_tree.tsv");
    rooted_tree.to_tsv(&rooted_path)?;
    println!("Written rooted tree → {}", rooted_path.display());

    // === Write Newick ===
    let newick_path = default_output_path(&args.input, "_tree.newick");
    let newick = rooted_tree.to_newick(model.pca.coords.nrows(), root, &model.encoder.sequences);

    {
        use std::fs::write;
        write(&newick_path, newick)?;
        println!("Written Newick tree → {}", newick_path.display());
    }

    Ok(())
}

fn default_output_path(input: &str, suffix: &str) -> PathBuf {
    let path = Path::new(input);
    let stem = path.file_stem().unwrap_or_default();
    let parent = path.parent().unwrap_or_else(|| Path::new("."));

    let mut out = PathBuf::new();
    out.push(parent);
    out.push(format!("{}{}", stem.to_string_lossy(), suffix));

    out
}
