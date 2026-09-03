use std::collections::{HashMap, HashSet};
use std::collections::hash_map::DefaultHasher;
use std::fs::{self, File};
use std::hash::{Hash, Hasher};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Parser;
use int_to_str::int_to_str::IntToStr;
use mapping_info::MappingInfo;
use rust_htslib::bam::{self, Read, record::Aux};
use sc_te::TeIndex;
use scdata::FeatureIndex;

#[derive(Debug, Parser)]
#[command(
    name = "nelrune-te-tests",
    about = "Inspect a 10x BAM against a prebuilt TE splice index before TE rescue/EM"
)]
struct Cli {
    /// Coordinate-sorted Cell Ranger / 10x BAM.
    #[arg(long)]
    bam: PathBuf,

    /// Binary index produced by `gtf-splice-index build`.
    #[arg(long = "te-index")]
    te_index: PathBuf,

    /// Stop after this many BAM records. Omit to scan the complete BAM.
    #[arg(long)]
    max_reads: Option<usize>,

    #[arg(long, default_value = "CB")]
    cell_tag: String,

    #[arg(long, default_value = "UB")]
    umi_tag: String,

    /// Output directory for MappingInfo and diagnostic TSV files.
    #[arg(long)]
    out: PathBuf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct ReadKey {
    qname_hash: u64,
    cell_id: u64,
    umi_id: u64,
}

#[derive(Debug, Default)]
struct MultiRead {
    candidates: HashSet<u64>,
}

#[derive(Debug, Default, Clone, Copy)]
struct FeatureCounts {
    anchor_reads: usize,
    multimapper_reads: usize,
    primary_records: usize,
    secondary_records: usize,
}

#[derive(Debug, Default)]
struct AmbiguitySummary {
    same_gene_different_bins: usize,
    different_genes_same_bin: usize,
    different_genes_and_bins: usize,
    multiple_chromosomes: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    if cli.cell_tag.len() != 2 || cli.umi_tag.len() != 2 {
        bail!("--cell-tag and --umi-tag must be exactly two characters");
    }
    fs::create_dir_all(&cli.out)
        .with_context(|| format!("creating {}", cli.out.display()))?;

    let mut report = MappingInfo::new(None, 0.0, 0);

    report.start_timer("te_index.load");
    let mut index = TeIndex::load(&cli.te_index)?;
    report.stop_timer("te_index.load");
    if index.annotation_transcript_count() == 0 {
        bail!("TE splice index contains no transcripts/loci");
    }
    eprintln!(
        "[nelrune-te-tests] loaded TE splice index: {} genes, {} transcripts/loci, bin_width={} bp",
        index.annotation_gene_count(),
        index.annotation_transcript_count(),
        index.splice_index().bin_width,
    );

    report.start_timer("bam.open");
    let mut bam = bam::Reader::from_path(&cli.bam)
        .with_context(|| format!("opening {}", cli.bam.display()))?;
    let header = bam.header().to_owned();
    report.stop_timer("bam.open");

    let cell_tag: [u8; 2] = cli.cell_tag.as_bytes().try_into().unwrap();
    let umi_tag: [u8; 2] = cli.umi_tag.as_bytes().try_into().unwrap();

    let mut cells = HashSet::<u64>::new();
    let mut cell_umi = HashSet::<(u64, u64)>::new();
    let mut multimappers = HashMap::<ReadKey, MultiRead>::new();
    let mut features = Vec::<FeatureCounts>::new();

    let mut scanned = 0usize;
    report.start_timer("bam.scan");
    for record in bam.records() {
        if cli.max_reads.is_some_and(|max| scanned >= max) {
            break;
        }
        let rec = record?;
        scanned += 1;
        report.report("bam.records.total");

        if rec.is_secondary() {
            report.report("bam.records.secondary");
        } else if rec.is_supplementary() {
            report.report("bam.records.supplementary");
        } else {
            report.report("bam.records.primary");
        }
        if rec.is_duplicate() {
            report.report("bam.records.duplicate");
        }
        if rec.is_unmapped() {
            report.report("bam.records.unmapped");
        } else {
            report.report("bam.records.mapped");
        }

        let nh = aux_u32(&rec, b"NH");
        match nh {
            Some(1) => report.report("bam.records.nh1"),
            Some(n) if n > 1 => report.report("bam.records.nh_gt1"),
            Some(_) => report.report("bam.records.nh0"),
            None => report.report("bam.records.no_nh"),
        }

        let Some(cell_raw) = aux_string(&rec, &cell_tag) else {
            report.report("10x.records.no_cell_tag");
            continue;
        };
        let Some(umi) = aux_string(&rec, &umi_tag) else {
            report.report("10x.records.no_umi_tag");
            continue;
        };
        report.report("10x.records.with_cb_ub");

        let cell = cell_raw
            .split_once('-')
            .map_or(cell_raw.as_str(), |(barcode, _)| barcode);
        let cell_id = IntToStr::new(cell.as_bytes()).into_u64();
        let umi_id = IntToStr::new(umi.as_bytes()).into_u64();
        cells.insert(cell_id);
        cell_umi.insert((cell_id, umi_id));

        if rec.is_unmapped() {
            report.report("10x.records.unmapped_with_cb_ub");
            continue;
        }

        let overlaps = index.record_overlaps(&rec, &header)?;
        if features.len() < index.len() {
            features.resize(index.len(), FeatureCounts::default());
        }

        match overlaps.len() {
            0 => report.report("te.records.no_overlap"),
            1 => report.report("te.records.one_feature"),
            _ => report.report("te.records.multiple_features"),
        }
        if overlaps.is_empty() {
            continue;
        }

        for &feature_id in &overlaps {
            if rec.is_secondary() {
                features[feature_id as usize].secondary_records += 1;
            } else {
                features[feature_id as usize].primary_records += 1;
            }
        }

        let is_primary = !rec.is_secondary() && !rec.is_supplementary();
        if is_primary && nh == Some(1) {
            match overlaps.len() {
                1 => {
                    report.report("te.anchor_candidates.reads");
                    features[overlaps[0] as usize].anchor_reads += 1;
                }
                _ => report.report("te.nh1.multiple_features"),
            }
        }

        if nh.is_some_and(|n| n > 1) {
            let key = ReadKey {
                qname_hash: hash_qname(rec.qname()),
                cell_id,
                umi_id,
            };
            multimappers
                .entry(key)
                .or_default()
                .candidates
                .extend(overlaps);
        }
    }

    report.stop_timer("bam.scan");
    report.total = scanned;

    report.start_timer("te.postprocess");
    report.report_n("10x.cells.distinct", cells.len());
    report.report_n("10x.cell_umi_pairs.distinct", cell_umi.len());
    report.report_n("te.features.observed", index.len());
    report.report_n("te.multimapper.read_ids.with_te", multimappers.len());

    let mut candidate_histogram = HashMap::<usize, usize>::new();
    let mut ambiguity = AmbiguitySummary::default();
    for multi in multimappers.values() {
        let n = multi.candidates.len();
        *candidate_histogram.entry(n).or_insert(0) += 1;
        match n {
            0 => report.report("te.multimapper.no_candidate"),
            1 => report.report("te.multimapper.one_candidate"),
            _ => {
                report.report("te.multimapper.multiple_candidates");
                classify_ambiguity(&index, &multi.candidates, &mut ambiguity);
            }
        }
        for &feature_id in &multi.candidates {
            features[feature_id as usize].multimapper_reads += 1;
        }
    }

    let te_features_seen = features
        .iter()
        .filter(|x| x.anchor_reads > 0 || x.multimapper_reads > 0)
        .count();
    let te_features_with_anchor = features.iter().filter(|x| x.anchor_reads > 0).count();
    let te_features_only_multi = features
        .iter()
        .filter(|x| x.anchor_reads == 0 && x.multimapper_reads > 0)
        .count();

    report.report_n("te.features.seen_in_anchor_or_multimapper", te_features_seen);
    report.report_n("te.features.with_anchor", te_features_with_anchor);
    report.report_n("te.features.only_multimapper", te_features_only_multi);
    report.stop_timer("te.postprocess");

    report.start_timer("diagnostics.write");
    write_candidate_histogram(&cli.out, &candidate_histogram)?;
    write_feature_table(&cli.out, &index, &features)?;
    report.stop_timer("diagnostics.write");

    write_mapping_info(&cli.out, &report)?;
    write_diagnostic_report(
        &cli.out,
        &report,
        &index,
        scanned,
        &candidate_histogram,
        &ambiguity,
    )?;

    println!("{}", report);
    print_diagnostic_report(
        &mut std::io::stdout(),
        &report,
        &index,
        scanned,
        &candidate_histogram,
        &ambiguity,
    )?;
    eprintln!(
        "[nelrune-te-tests] scanned {} BAM records{}; observed {} spatial TE features; wrote {}",
        scanned,
        cli.max_reads
            .map(|n| format!(" (max-reads={n})"))
            .unwrap_or_default(),
        index.len(),
        cli.out.display()
    );
    Ok(())
}

fn classify_ambiguity(
    index: &TeIndex,
    candidates: &HashSet<u64>,
    summary: &mut AmbiguitySummary,
) {
    let mut chromosomes = HashSet::<usize>::new();
    let mut bins = HashSet::<(usize, usize)>::new();
    let mut genes = HashSet::<usize>::new();

    for &feature_id in candidates {
        let Some((chr_id, bin_id, gene_id)) = index.feature_key(feature_id) else {
            continue;
        };
        chromosomes.insert(chr_id);
        bins.insert((chr_id, bin_id));
        genes.insert(gene_id);
    }

    if chromosomes.len() > 1 {
        summary.multiple_chromosomes += 1;
    } else if genes.len() == 1 && bins.len() > 1 {
        summary.same_gene_different_bins += 1;
    } else if genes.len() > 1 && bins.len() == 1 {
        summary.different_genes_same_bin += 1;
    } else {
        summary.different_genes_and_bins += 1;
    }
}

fn pct(n: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        100.0 * n as f64 / denominator as f64
    }
}

fn write_count_line(
    writer: &mut impl Write,
    label: &str,
    count: usize,
    denominator: usize,
) -> Result<()> {
    writeln!(writer, "{label:<46} {count:>12}   ({:>6.2}%)", pct(count, denominator))?;
    Ok(())
}

fn write_diagnostic_report(
    out: &PathBuf,
    report: &MappingInfo,
    index: &TeIndex,
    scanned: usize,
    candidate_histogram: &HashMap<usize, usize>,
    ambiguity: &AmbiguitySummary,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(out.join("te_diagnostic_report.txt"))?);
    print_diagnostic_report(
        &mut writer,
        report,
        index,
        scanned,
        candidate_histogram,
        ambiguity,
    )
}

fn print_diagnostic_report(
    writer: &mut impl Write,
    report: &MappingInfo,
    index: &TeIndex,
    scanned: usize,
    candidate_histogram: &HashMap<usize, usize>,
    ambiguity: &AmbiguitySummary,
) -> Result<()> {
    let tagged = report.get_issue_count("10x.records.with_cb_ub");
    let cells = report.get_issue_count("10x.cells.distinct");
    let molecules = report.get_issue_count("10x.cell_umi_pairs.distinct");
    let nh1 = report.get_issue_count("bam.records.nh1");
    let nh_multi = report.get_issue_count("bam.records.nh_gt1");
    let no_overlap = report.get_issue_count("te.records.no_overlap");
    let one_feature = report.get_issue_count("te.records.one_feature");
    let multiple_features = report.get_issue_count("te.records.multiple_features");
    let anchors = report.get_issue_count("te.anchor_candidates.reads");
    let multi_te = report.get_issue_count("te.multimapper.read_ids.with_te");
    let multi_one = report.get_issue_count("te.multimapper.one_candidate");
    let multi_many = report.get_issue_count("te.multimapper.multiple_candidates");
    let with_anchor = report.get_issue_count("te.features.with_anchor");
    let only_multi = report.get_issue_count("te.features.only_multimapper");

    writeln!(writer, "\nNELRUNE-TE BAM DIAGNOSTIC")?;
    writeln!(writer, "=========================\n")?;

    writeln!(writer, "Input")?;
    writeln!(writer, "-----")?;
    writeln!(writer, "BAM records scanned:                         {scanned:>12}")?;
    write_count_line(writer, "Records with cell barcode + UMI:", tagged, scanned)?;
    writeln!(writer, "Cells observed:                              {cells:>12}")?;
    writeln!(writer, "Distinct cell/UMI pairs:                     {molecules:>12}\n")?;

    writeln!(writer, "Alignment structure")?;
    writeln!(writer, "-------------------")?;
    write_count_line(writer, "Unique alignments (NH=1):", nh1, scanned)?;
    write_count_line(writer, "Multimapping alignments (NH>1):", nh_multi, scanned)?;

    writeln!(writer, "\nTE overlap")?;
    writeln!(writer, "----------")?;
    writeln!(writer, "Percentages in this section use records carrying both CB and UB as denominator.")?;
    write_count_line(writer, "No spatial TE feature overlap:", no_overlap, tagged)?;
    write_count_line(writer, "Exactly one spatial TE feature:", one_feature, tagged)?;
    write_count_line(writer, "Multiple spatial TE features:", multiple_features, tagged)?;
    writeln!(writer)?;
    writeln!(writer, "A spatial TE feature is defined as:")?;
    writeln!(writer, "    chromosome + {} bp genomic bin + TE gene/subfamily", index.splice_index().bin_width)?;
    writeln!(writer, "Example:")?;
    writeln!(writer, "    chr1 + 100-101 Mb + L1M2")?;

    writeln!(writer, "\nTE anchor evidence")?;
    writeln!(writer, "------------------")?;
    write_count_line(writer, "Reads usable as unambiguous TE anchors:", anchors, tagged)?;
    writeln!(writer, "Spatial TE features observed:                 {:>12}", index.len())?;
    writeln!(writer, "Features supported by anchors:                {with_anchor:>12}")?;
    writeln!(writer, "Features seen only through multimappers:      {only_multi:>12}")?;

    writeln!(writer, "\nTE multimapping")?;
    writeln!(writer, "---------------")?;
    writeln!(writer, "Multimapping read IDs overlapping a TE:       {multi_te:>12}")?;
    write_count_line(writer, "Resolve to ONE spatial TE feature:", multi_one, multi_te)?;
    write_count_line(writer, "Remain ambiguous between features:", multi_many, multi_te)?;
    if multi_te > 0 {
        writeln!(writer)?;
        writeln!(writer, "Interpretation:")?;
        writeln!(
            writer,
            "  {:.2}% of TE-associated multimapping reads collapse to one spatial TE feature.",
            pct(multi_one, multi_te)
        )?;
        writeln!(
            writer,
            "  {:.2}% still have more than one candidate and require further resolution.",
            pct(multi_many, multi_te)
        )?;
    }

    writeln!(writer, "\nAmbiguity structure")?;
    writeln!(writer, "-------------------")?;
    writeln!(writer, "Among TE multimappers that still have >1 spatial candidate:")?;
    write_count_line(
        writer,
        "Same TE gene, different genomic bins:",
        ambiguity.same_gene_different_bins,
        multi_many,
    )?;
    write_count_line(
        writer,
        "Different TE genes, same genomic bin:",
        ambiguity.different_genes_same_bin,
        multi_many,
    )?;
    write_count_line(
        writer,
        "Different TE genes and genomic bins:",
        ambiguity.different_genes_and_bins,
        multi_many,
    )?;
    write_count_line(
        writer,
        "Candidates on multiple chromosomes:",
        ambiguity.multiple_chromosomes,
        multi_many,
    )?;

    writeln!(writer, "\nCandidate features per TE multimapping read")?;
    writeln!(writer, "-------------------------------------------")?;
    let mut rows: Vec<_> = candidate_histogram.iter().map(|(&n, &count)| (n, count)).collect();
    rows.sort_unstable_by_key(|x| x.0);
    for (n, count) in rows {
        writeln!(writer, "{n:>4} candidate{:1}: {count:>12}   ({:>6.2}%)", if n == 1 { " " } else { "s" }, pct(count, multi_te))?;
    }

    writeln!(writer, "\nTiming")?;
    writeln!(writer, "------")?;
    for (name, label) in [
        ("te_index.load", "TE index loading"),
        ("bam.open", "BAM opening"),
        ("bam.scan", "BAM scanning"),
        ("te.postprocess", "Post-processing"),
        ("diagnostics.write", "Diagnostic TSV writing"),
    ] {
        writeln!(
            writer,
            "{label:<28} {:>10.3} s",
            report.timer_duration(name).as_secs_f64()
        )?;
    }
    let scan_seconds = report.timer_duration("bam.scan").as_secs_f64();
    if scan_seconds > 0.0 {
        writeln!(writer, "BAM scan throughput          {:>10.0} records/s", scanned as f64 / scan_seconds)?;
    }

    Ok(())
}

fn aux_string(rec: &bam::Record, tag: &[u8; 2]) -> Option<String> {
    match rec.aux(tag).ok()? {
        Aux::String(s) => Some(s.to_string()),
        _ => None,
    }
}

fn aux_u32(rec: &bam::Record, tag: &[u8; 2]) -> Option<u32> {
    match rec.aux(tag).ok()? {
        Aux::U8(x) => Some(x as u32),
        Aux::I8(x) if x >= 0 => Some(x as u32),
        Aux::U16(x) => Some(x as u32),
        Aux::I16(x) if x >= 0 => Some(x as u32),
        Aux::U32(x) => Some(x),
        Aux::I32(x) if x >= 0 => Some(x as u32),
        _ => None,
    }
}

fn hash_qname(qname: &[u8]) -> u64 {
    let mut hasher = DefaultHasher::new();
    qname.hash(&mut hasher);
    hasher.finish()
}

fn write_mapping_info(out: &PathBuf, report: &MappingInfo) -> Result<()> {
    let mut writer = BufWriter::new(File::create(out.join("mapping_info.txt"))?);
    write!(writer, "{}", report)?;
    Ok(())
}

fn write_candidate_histogram(out: &PathBuf, histogram: &HashMap<usize, usize>) -> Result<()> {
    let mut rows: Vec<_> = histogram.iter().map(|(&n, &count)| (n, count)).collect();
    rows.sort_unstable_by_key(|x| x.0);
    let mut writer = BufWriter::new(File::create(out.join("candidate_set_histogram.tsv"))?);
    writeln!(writer, "te_candidates\tread_ids")?;
    for (n, count) in rows {
        writeln!(writer, "{n}\t{count}")?;
    }
    Ok(())
}

fn write_feature_table(out: &PathBuf, index: &TeIndex, counts: &[FeatureCounts]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(out.join("te_features.tsv"))?);
    writeln!(
        writer,
        "feature_id\tfeature\tgene\tchrom\tbin\tstart\tend\tanchor_reads\tmultimapper_reads\tprimary_records\tsecondary_records"
    )?;

    for (id, count) in counts.iter().enumerate() {
        if count.anchor_reads == 0
            && count.multimapper_reads == 0
            && count.primary_records == 0
            && count.secondary_records == 0
        {
            continue;
        }

        let Some((chr_id, bin_id, _gene_id)) = index.feature_key(id as u64) else {
            continue;
        };
        let Some((chrom, start, end, gene)) = index.feature_coordinates(id as u64) else {
            continue;
        };

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            id,
            index.feature_name(id as u64),
            gene,
            chrom,
            bin_id,
            start,
            end,
            count.anchor_reads,
            count.multimapper_reads,
            count.primary_records,
            count.secondary_records,
        )?;

        debug_assert!(chr_id < index.splice_index().chr_names.len());
    }
    Ok(())
}
