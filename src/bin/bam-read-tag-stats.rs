use anyhow::Result;
use bam_tide::read_tag_table::{
    PairStats, ReadTagTable, ReadTagTableCli, Summary,
};
use clap::Parser;

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Summarize external read-tag tables produced by bam-ont-normalizer or compatible tools."
)]

struct Cli {
    /// Shared read-tag table column options.
    ///
    /// For this stats command, --read-tag-table does not need to be supplied:
    /// the positional input path is used as the table path.
    #[command(flatten)]
    read_tags: ReadTagTableCli,

    /// Minimum number of observations for one cell + UMI pair.
    #[arg(long, default_value_t = 2)]
    min_pair_count: u64,

    /// Minimum number of surviving UMIs per cell after pair filtering.
    #[arg(long, default_value_t = 3)]
    min_cell_umis: u64,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let read_tags = cli.read_tags;

    let config = read_tags
        .to_config()
        .expect("positional input should provide read-tag-table config");

    let table = ReadTagTable::from_config(&config)?;

    let pre = table.summarize_pairs(1, 1);
    let post = table.summarize_pairs(cli.min_pair_count, cli.min_cell_umis);

    println!("Input: {:?}", read_tags.read_tag_table);
    println!();

    println!("Rows");
    println!("  accepted read-tag entries : {}", table.len());
    println!();

    println!("Filtering");
    println!("  min_pair_count           : {}", cli.min_pair_count);
    println!("  min_cell_umis            : {}", cli.min_cell_umis);
    println!();

    print_stats("Pre-filtered", &pre);
    println!();
    print_stats("Post-filtered", &post);

    Ok(())
}

fn print_stats(name: &str, stats: &PairStats) {
    println!("{name}");
    println!("  Cell entries              : {}", stats.cell_entries);
    println!(
        "  Unique cell + UMI         : {}",
        stats.unique_cell_umi_combos
    );
    println!(
        "  Total pair observations   : {}",
        stats.total_pair_observations
    );
    println!();

    print_summary("  UMIs per cell", &stats.umis_per_cell);
    println!();
    print_summary(
        "  Detections per cell + UMI",
        &stats.detections_per_cell_umi,
    );
}

fn print_summary(name: &str, s: &Summary) {
    println!("{name}");
    println!("    n      : {}", s.n);
    println!("    mean   : {:.3}", s.mean);
    println!("    median : {:.3}", s.median);
    println!("    min    : {}", s.min);
    println!("    max    : {}", s.max);
}