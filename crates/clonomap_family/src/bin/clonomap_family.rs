use anyhow::{Context, Result, bail};
use clap::Parser;
use clonomap::{
    AlignedCell, CellReceptor, ClonoMap, Family, FamilyConfig, MutationMeasurement, Receptor,
    ReferenceModels, align_fragment, rooted_categorical_hex, rooted_continuous_hex,
};
use ndarray::Array2;
use sc_primer::{Chemistry, PrimerDetector};
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::{File, create_dir_all};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Validate ClonoMap receptor families from a Lumrik VDJ output directory without plotting"
)]
struct Args {
    /// Lumrik/nelrune-vdj output directory containing airr_rearrangements.tsv and vdj_calls.tsv
    #[arg(long, num_args = 1..)]
    vdj_out: Vec<PathBuf>,
    /// Output directory for family and rejection reports
    #[arg(long)]
    out: PathBuf,
    /// Single-cell chemistry/chemistries used to identify and clip known technical sequence.
    #[arg(long, value_enum, num_args = 1.., default_values_t = [Chemistry::BdV2_384, Chemistry::BdV2_384Vdj])]
    chemistry: Vec<Chemistry>,
    /// Maximum nucleotide CDR3 edit distance during initial HC family collection
    #[arg(long, default_value_t = 3)]
    max_cdr3_distance: usize,
    /// Minimum final HC family size that is interesting enough to plot/report regardless of Pearson significance
    #[arg(long, default_value_t = 100)]
    min_family_size: usize,
    /// Maximum nominal Pearson p-value that is interesting enough to plot/report regardless of family size
    #[arg(long, default_value_t = 0.05)]
    max_pearson_p: f64,
    /// Minimum identity for the first HC mutation-alignment pass
    #[arg(long, default_value_t = 0.50)]
    min_alignment_identity: f64,
    /// Hard identity required when an ejected HC tries another family
    #[arg(long, default_value_t = 0.75)]
    hard_alignment_identity: f64,
    /// Restore SVG family plots for final families at or above --min-family-size.
    #[arg(long)]
    plots: bool,
    /// Draw rooted family plots radially instead of the default layered layout.
    #[arg(long)]
    radial_layout: bool,
    /// Persistent reference_curator store. Discovered sequences keep stable identities across runs.
    #[arg(long)]
    reference_curator: Option<PathBuf>,
}

#[derive(Clone)]
struct CallRow {
    cell: String,
    source: String,
    receptor: Receptor,
    productive: bool,
}

#[derive(Debug)]
struct UnassignedCell {
    cell: CellReceptor,
    from_family: String,
}

fn main() -> Result<()> {
    let args = Args::parse();
    if !(0.0..=1.0).contains(&args.max_pearson_p) {
        bail!("--max-pearson-p must be between 0 and 1");
    }
    create_dir_all(&args.out).with_context(|| format!("create {}", args.out.display()))?;
    let cfg = FamilyConfig {
        max_cdr3_distance: args.max_cdr3_distance,
        min_alignment_identity: args.min_alignment_identity,
        hard_alignment_identity: args.hard_alignment_identity,
        ..FamilyConfig::default()
    };

    let sources = input_sources(&args.vdj_out)?;
    let mut rows = Vec::new();
    for (source, dir) in &sources {
        let mut input = read_vdj_output(dir, source)?;
        println!(
            "Loaded {} receptor calls from {} ({})",
            input.len(),
            source,
            dir.display()
        );
        rows.append(&mut input);
    }
    // ClonoMap discovers reference-incompatible fragments; reference_curator
    // owns their persistent identity, resolution state and external evidence.
    let curator_store = args
        .reference_curator
        .clone()
        .unwrap_or_else(|| args.out.join("reference_curator.bin"));
    let mut novel_v = ReferenceModels::open(Some(&curator_store), &args.chemistry)
        .map_err(anyhow::Error::msg)
        .context("opening reference curator")?;

    // Learn the ordinary distance-to-supplied-germline background from the full
    // repertoire, before any family-size or CDR3 grouping can bias it. HC and LC
    // have separate backgrounds. Receptors beyond mean + 2 SD receive a
    // hypothetical V model through reference_curator, then the normal ClonoMap
    // family machinery starts from scratch with those augmented assignments.
    let refinement = refine_global_reference_outliers(&mut rows, &mut novel_v);
    println!("\nGlobal reference refinement (max distance from supplied germline)");
    print_refinement_line("HC", &refinement.hc);
    print_refinement_line("LC", &refinement.lc);

    let mut novel_receptors = refinement.hc.hypothesized + refinement.lc.hypothesized;
    for row in &mut rows {
        if row.productive
            && !row.receptor.v.starts_with("RC-")
            && novel_v.resolve_receptor_for_sample(
                &mut row.receptor,
                cfg.hard_alignment_identity,
                Some(&row.source),
            )
        {
            novel_receptors += 1;
        }
    }
    let mut by_cell: BTreeMap<String, Vec<CallRow>> = BTreeMap::new();
    for row in rows {
        by_cell.entry(row.cell.clone()).or_default().push(row);
    }

    let examined = by_cell.len();
    let mut no_hc = 0usize;
    let mut multiple_hc = 0usize;
    let mut qualified = Vec::new();
    for (cell_id, rows) in by_cell {
        let hcs: Vec<&CallRow> = rows
            .iter()
            .filter(|r| r.productive && r.receptor.chain == "IGH" && receptor_complete(&r.receptor))
            .collect();
        match hcs.as_slice() {
            [] => no_hc += 1,
            [hc] => {
                let lc = rows
                    .iter()
                    .filter(|r| {
                        r.productive
                            && matches!(r.receptor.chain.as_str(), "IGK" | "IGL")
                            && receptor_complete(&r.receptor)
                    })
                    .map(|r| r.receptor.clone())
                    .collect();
                qualified.push(CellReceptor {
                    cell_id,
                    hc: hc.receptor.clone(),
                    lc,
                });
            }
            _ => multiple_hc += 1,
        }
    }

    println!("HC qualification");
    println!("  cells examined:                 {examined}");
    println!("  accepted: exactly one valid HC: {}", qualified.len());
    println!("  skipped: no valid HC:           {no_hc}");
    println!("  skipped: multiple valid HCs:    {multiple_hc}");
    println!("  NOTE: apparent duplicate HCs are intentionally not collapsed; investigate later.");

    let mut families = build_initial_families(qualified, &cfg);
    let provisional_cells: usize = families.iter().map(Family::provisional_len).sum();
    println!(
        "\nInitial HC families: {} ({provisional_cells} cells)",
        families.len()
    );

    // Family owns only accepted members. align() returns bare CellReceptors;
    // reassignment policy and unresolved cells belong to this outer layer.
    let mut rejected = Vec::<(String, CellReceptor)>::new();
    for family in &mut families {
        let from = family.name.clone();
        rejected.extend(
            family
                .align(&cfg, &novel_v)
                .into_iter()
                .map(|cell| (from.clone(), cell)),
        );
    }
    let failed_original = rejected.len();
    let unassigned = reassign_cells(&mut families, rejected, &cfg, &novel_v);
    let reassigned = failed_original.saturating_sub(unassigned.len());
    let final_cells: usize = families.iter().map(|f| f.members.len()).sum();

    println!("\nHC post-cluster mutation validation");
    println!("  provisional cells:             {provisional_cells}");
    println!("  passed/reassigned final cells: {final_cells}");
    println!("  failed original family:        {failed_original}");
    println!("  reassigned to another family:  {reassigned}");
    println!("  remained unassigned:           {}", unassigned.len());

    // LC clone ownership starts only after HC membership is final.
    for family in &mut families {
        family.split_light_chains(&cfg, &novel_v);
    }

    write_family_report(&args.out.join("families.tsv"), &families)?;
    let cell_id_detectors: Vec<PrimerDetector> = args
        .chemistry
        .iter()
        .copied()
        .filter_map(|chemistry| PrimerDetector::from_chemistry(chemistry).ok())
        .collect();
    write_cell_report(&args.out.join("cells.tsv"), &families, &cell_id_detectors)?;
    write_overlap_reports(&args.out, &families)?;
    write_reference_candidate_report(&args.out.join("reference_candidates.tsv"), &novel_v)?;
    novel_v
        .save(&curator_store)
        .map_err(anyhow::Error::msg)
        .with_context(|| format!("saving reference curator {}", curator_store.display()))?;
    write_unassigned(&args.out.join("hc_unassigned.tsv"), &unassigned)?;
    if args.plots {
        write_family_plots(
            &args.out.join("plots"),
            &families,
            args.min_family_size,
            args.max_pearson_p,
            args.radial_layout,
        )?;
    }

    println!(
        "\nReference-incompatible models: {} ({} receptors assigned; evidence in reference_candidates.tsv)",
        novel_v.total(),
        novel_receptors
    );

    let mut selected: Vec<_> = families
        .iter()
        .filter(|f| family_is_selected(f, args.min_family_size, args.max_pearson_p))
        .collect();
    selected.sort_by_key(|f| std::cmp::Reverse(f.members.len()));
    println!(
        "\nFinal HC families selected by size >= {} OR Pearson p <= {}: {}",
        args.min_family_size,
        args.max_pearson_p,
        selected.len()
    );
    println!(
        "  {:<72} {:>6}  {:>6}  {:>9}  {:>6}  {:>8}  {:>8}  {:>16}",
        "HC family", "cells", "HC max", "LC clones", "LC n", "LC mean", "LC SD", "LC top 3"
    );
    println!("  {}", "-".repeat(145));
    for f in selected {
        let r = f.mutation_report();
        println!(
            "  {:<72} {:>6}  {:>6}  {:>9}  {:>6}  {:>8}  {:>8}  {:>16}",
            r.family,
            r.cells,
            fmt_usize_option(r.max_hc_mutations),
            r.lc_clones,
            r.lc_mutations.n,
            fmt_float_option(r.lc_mutations.mean),
            fmt_float_option(r.lc_mutations.sd),
            fmt_max3(&r.lc_mutations.max3),
        );
    }
    println!("\nReports written to {}", args.out.display());
    Ok(())
}

#[derive(Default)]
struct PlotState {
    cells: BTreeSet<String>,
    isotypes: BTreeMap<String, usize>,
    lc_name: String,
    hc_depth: usize,
    lc_depth: Option<usize>,
}

#[derive(Clone, Copy, Default)]
struct PearsonStat {
    n: usize,
    r: Option<f64>,
    p: Option<f64>,
}

#[derive(Default)]
struct FamilyPlotAnalysis {
    abundance_hc: PearsonStat,
    abundance_lc: PearsonStat,
    abundance_paired: PearsonStat,
    hc_lc: PearsonStat,
    isotype_hex: String,
    light_chain_hex: String,
    hc_depth_hex: String,
    lc_depth_hex: String,
    paired_depth_hex: String,
    hclc: BTreeMap<String, HclcAnalysis>,
}

#[derive(Clone, Copy, Default)]
struct HclcAnalysis {
    abundance_hc: PearsonStat,
    abundance_lc: PearsonStat,
    abundance_paired: PearsonStat,
    hc_lc: PearsonStat,
}

fn mutation_distance(m: &MutationMeasurement) -> usize {
    m.substitutions + m.indels.iter().map(|x| x.len).sum::<usize>()
}

#[derive(Clone, Copy, Default)]
struct RefinementStats {
    n: usize,
    mean: Option<f64>,
    sd: Option<f64>,
    threshold: Option<f64>,
    outliers: usize,
    hypothesized: usize,
}

#[derive(Clone, Copy, Default)]
struct GlobalRefinement {
    hc: RefinementStats,
    lc: RefinementStats,
}

fn distance_stats(values: &[usize]) -> RefinementStats {
    if values.is_empty() {
        return RefinementStats::default();
    }
    let n = values.len();
    let mean = values.iter().map(|&x| x as f64).sum::<f64>() / n as f64;
    let variance = values
        .iter()
        .map(|&x| {
            let d = x as f64 - mean;
            d * d
        })
        .sum::<f64>()
        / n as f64;
    let sd = variance.sqrt();
    RefinementStats {
        n,
        mean: Some(mean),
        sd: Some(sd),
        threshold: Some(mean + 2.0 * sd),
        ..Default::default()
    }
}

fn supplied_germline_distance(row: &CallRow) -> Option<usize> {
    if !row.productive || !receptor_complete(&row.receptor) {
        return None;
    }
    align_fragment(&row.receptor.naive, &row.receptor.observed).map(|m| mutation_distance(&m))
}

fn refine_global_reference_outliers(
    rows: &mut [CallRow],
    models: &mut ReferenceModels,
) -> GlobalRefinement {
    let hc_values: Vec<usize> = rows
        .iter()
        .filter(|r| r.receptor.chain == "IGH")
        .filter_map(supplied_germline_distance)
        .collect();
    let lc_values: Vec<usize> = rows
        .iter()
        .filter(|r| matches!(r.receptor.chain.as_str(), "IGK" | "IGL"))
        .filter_map(supplied_germline_distance)
        .collect();
    let mut result = GlobalRefinement {
        hc: distance_stats(&hc_values),
        lc: distance_stats(&lc_values),
    };

    for row in rows {
        let stats = match row.receptor.chain.as_str() {
            "IGH" => &mut result.hc,
            "IGK" | "IGL" => &mut result.lc,
            _ => continue,
        };
        let Some(threshold) = stats.threshold else {
            continue;
        };
        let Some(distance) = supplied_germline_distance(row) else {
            continue;
        };
        if distance as f64 <= threshold {
            continue;
        }
        stats.outliers += 1;
        if models.hypothesize_receptor_for_sample(&mut row.receptor, Some(&row.source)) {
            stats.hypothesized += 1;
        }
    }
    result
}

fn print_refinement_line(label: &str, stats: &RefinementStats) {
    match (stats.mean, stats.sd, stats.threshold) {
        (Some(mean), Some(sd), Some(threshold)) => println!(
            "  {label}: n={} mean={mean:.2} SD={sd:.2} threshold={threshold:.2} outliers={} hypothetical models assigned={}",
            stats.n, stats.outliers, stats.hypothesized
        ),
        _ => println!("  {label}: no evaluable receptors"),
    }
}

fn pearson(x: &[f64], y: &[f64]) -> PearsonStat {
    let n = x.len().min(y.len());
    if n < 3 {
        return PearsonStat {
            n,
            ..Default::default()
        };
    }
    let (x, y) = (&x[..n], &y[..n]);
    let mx = x.iter().sum::<f64>() / n as f64;
    let my = y.iter().sum::<f64>() / n as f64;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    let mut sxy = 0.0;
    for (&a, &b) in x.iter().zip(y) {
        let dx = a - mx;
        let dy = b - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    if sxx == 0.0 || syy == 0.0 {
        return PearsonStat {
            n,
            ..Default::default()
        };
    }
    let r = (sxy / (sxx * syy).sqrt()).clamp(-1.0, 1.0);
    let p = if r.abs() >= 1.0 {
        0.0
    } else {
        let t = r.abs() * (((n - 2) as f64) / (1.0 - r * r)).sqrt();
        let dist = StudentsT::new(0.0, 1.0, (n - 2) as f64).expect("valid Student t distribution");
        2.0 * (1.0 - dist.cdf(t))
    };
    PearsonStat {
        n,
        r: Some(r),
        p: Some(p),
    }
}

fn family_plot_analysis(family: &Family) -> FamilyPlotAnalysis {
    let hc_by_cell: HashMap<&str, (&str, usize)> = family
        .members
        .iter()
        .map(|m| {
            (
                m.cell.cell_id.as_str(),
                (m.cell.hc.c.as_str(), m.hc_mutations.mutation_events()),
            )
        })
        .collect();
    let mut paired_cells = BTreeSet::new();
    let mut states: BTreeMap<(usize, String, Option<usize>), PlotState> = BTreeMap::new();
    for lc in &family.light_clones {
        for member in &lc.members {
            let Some(&(isotype, hc_depth)) = hc_by_cell.get(member.cell.as_str()) else {
                continue;
            };
            paired_cells.insert(member.cell.clone());
            let lc_depth = member.mutations.mutation_events();
            let state = states
                .entry((hc_depth, lc.name.clone(), Some(lc_depth)))
                .or_default();
            state.cells.insert(member.cell.clone());
            *state
                .isotypes
                .entry(
                    if isotype.is_empty() {
                        "unknown"
                    } else {
                        isotype
                    }
                    .to_string(),
                )
                .or_default() += 1;
            state.lc_name = lc.name.clone();
            state.hc_depth = hc_depth;
            state.lc_depth = Some(lc_depth);
        }
    }
    for member in &family.members {
        if paired_cells.contains(&member.cell.cell_id) {
            continue;
        }
        let hc_depth = member.hc_mutations.mutation_events();
        let state = states
            .entry((hc_depth, "unpaired".to_string(), None))
            .or_default();
        state.cells.insert(member.cell.cell_id.clone());
        let iso = if member.cell.hc.c.is_empty() {
            "unknown"
        } else {
            member.cell.hc.c.as_str()
        };
        *state.isotypes.entry(iso.to_string()).or_default() += 1;
        state.lc_name = "unpaired".to_string();
        state.hc_depth = hc_depth;
    }
    let sv: Vec<PlotState> = states.into_values().collect();
    let paired: Vec<&PlotState> = sv.iter().filter(|s| s.lc_depth.is_some()).collect();
    let abundance: Vec<f64> = paired.iter().map(|s| s.cells.len() as f64).collect();
    let hc: Vec<f64> = paired.iter().map(|s| s.hc_depth as f64).collect();
    let lc: Vec<f64> = paired.iter().map(|s| s.lc_depth.unwrap() as f64).collect();
    let total: Vec<f64> = paired
        .iter()
        .map(|s| (s.hc_depth + s.lc_depth.unwrap()) as f64)
        .collect();

    let mut iso_cat = vec!["HC NAIVE".to_string()];
    let mut lc_cat = vec!["HC NAIVE".to_string()];
    let mut hd = vec![None];
    let mut ld = vec![None];
    let mut pd = vec![None];
    for s in &sv {
        iso_cat.push(
            s.isotypes
                .iter()
                .max_by(|(ka, va), (kb, vb)| va.cmp(vb).then_with(|| kb.cmp(ka)))
                .map(|(k, _)| k.clone())
                .unwrap_or_else(|| "unknown".into()),
        );
        lc_cat.push(s.lc_name.clone());
        hd.push(Some(s.hc_depth as f32));
        ld.push(s.lc_depth.map(|x| x as f32));
        pd.push(Some((s.hc_depth + s.lc_depth.unwrap_or(0)) as f32));
    }
    let mut lc_order: Vec<String> = family.light_clones.iter().map(|x| x.name.clone()).collect();
    lc_order.sort_by_key(|name| {
        std::cmp::Reverse(
            family
                .light_clones
                .iter()
                .find(|x| &x.name == name)
                .map(|x| x.members.len())
                .unwrap_or(0),
        )
    });
    lc_order.dedup();
    if lc_cat.iter().any(|x| x == "unpaired") {
        lc_order.push("unpaired".into());
    }
    let mut hclc = BTreeMap::new();
    for lc_name in family.light_clones.iter().map(|x| x.name.clone()) {
        let subset: Vec<&PlotState> = paired
            .iter()
            .copied()
            .filter(|s| s.lc_name == lc_name)
            .collect();
        let abundance: Vec<f64> = subset.iter().map(|s| s.cells.len() as f64).collect();
        let hc: Vec<f64> = subset.iter().map(|s| s.hc_depth as f64).collect();
        let lc: Vec<f64> = subset.iter().map(|s| s.lc_depth.unwrap() as f64).collect();
        let total: Vec<f64> = subset
            .iter()
            .map(|s| (s.hc_depth + s.lc_depth.unwrap()) as f64)
            .collect();
        hclc.insert(
            lc_name,
            HclcAnalysis {
                abundance_hc: pearson(&abundance, &hc),
                abundance_lc: pearson(&abundance, &lc),
                abundance_paired: pearson(&abundance, &total),
                hc_lc: pearson(&hc, &lc),
            },
        );
    }
    let iso_colors = rooted_categorical_hex(&iso_cat, None, 0);
    let lc_colors = rooted_categorical_hex(&lc_cat, Some(&lc_order), 0);
    let hd_colors = rooted_continuous_hex(&hd, 0);
    let ld_colors = rooted_continuous_hex(&ld, 0);
    let pd_colors = rooted_continuous_hex(&pd, 0);
    let encode = |colors: &[String]| -> String {
        sv.iter()
            .enumerate()
            .flat_map(|(i, s)| {
                s.cells
                    .iter()
                    .map(move |cell| format!("{}={}", cell, colors[i + 1]))
            })
            .collect::<Vec<_>>()
            .join(";")
    };
    FamilyPlotAnalysis {
        abundance_hc: pearson(&abundance, &hc),
        abundance_lc: pearson(&abundance, &lc),
        abundance_paired: pearson(&abundance, &total),
        hc_lc: pearson(&hc, &lc),
        isotype_hex: encode(&iso_colors),
        light_chain_hex: encode(&lc_colors),
        hc_depth_hex: encode(&hd_colors),
        lc_depth_hex: encode(&ld_colors),
        paired_depth_hex: encode(&pd_colors),
        hclc,
    }
}

fn stats_significant(stats: &[PearsonStat], max_p: f64) -> bool {
    stats.iter().any(|s| s.p.is_some_and(|p| p <= max_p))
}

fn family_is_significant(a: &FamilyPlotAnalysis, max_p: f64) -> bool {
    stats_significant(
        &[a.abundance_hc, a.abundance_lc, a.abundance_paired, a.hc_lc],
        max_p,
    ) || a.hclc.values().any(|x| {
        stats_significant(
            &[x.abundance_hc, x.abundance_lc, x.abundance_paired, x.hc_lc],
            max_p,
        )
    })
}

fn family_is_selected(family: &Family, min_size: usize, max_p: f64) -> bool {
    family.members.len() >= min_size || family_is_significant(&family_plot_analysis(family), max_p)
}

fn write_family_plots(
    out: &Path,
    families: &[Family],
    min_size: usize,
    max_p: f64,
    radial_layout: bool,
) -> Result<()> {
    create_dir_all(out).with_context(|| format!("create {}", out.display()))?;
    for (plot_index, family) in families
        .iter()
        .filter(|f| family_is_selected(f, min_size, max_p))
        .enumerate()
    {
        write_one_family_plot(out, family, plot_index + 1, min_size, max_p, radial_layout)?;
    }
    Ok(())
}

fn write_one_family_plot(
    out: &Path,
    family: &Family,
    plot_index: usize,
    _min_size: usize,
    _max_p: f64,
    radial_layout: bool,
) -> Result<()> {
    // Plotting consumes the frozen Family.  It deliberately uses the mutation
    // measurements already accepted by HC/LC family logic; it never realigns a
    // receptor or changes membership.  Until span-aware positional mutation
    // coordinates are available, topology uses the honest comparable quantities
    // we have for every partial receptor: HC and LC mutation-event depth.
    let hc_by_cell: HashMap<&str, (&str, usize)> = family
        .members
        .iter()
        .map(|m| {
            (
                m.cell.cell_id.as_str(),
                (m.cell.hc.c.as_str(), m.hc_mutations.mutation_events()),
            )
        })
        .collect();
    let mut paired_cells = BTreeSet::new();
    let mut states: BTreeMap<(usize, String, Option<usize>), PlotState> = BTreeMap::new();

    for lc in &family.light_clones {
        for member in &lc.members {
            let Some(&(isotype, hc_depth)) = hc_by_cell.get(member.cell.as_str()) else {
                continue;
            };
            paired_cells.insert(member.cell.clone());
            let lc_depth = member.mutations.mutation_events();
            let key = (hc_depth, lc.name.clone(), Some(lc_depth));
            let state = states.entry(key).or_default();
            state.cells.insert(member.cell.clone());
            *state
                .isotypes
                .entry(
                    if isotype.is_empty() {
                        "unknown"
                    } else {
                        isotype
                    }
                    .to_string(),
                )
                .or_default() += 1;
            state.lc_name = lc.name.clone();
            state.hc_depth = hc_depth;
            state.lc_depth = Some(lc_depth);
        }
    }
    for member in &family.members {
        if paired_cells.contains(&member.cell.cell_id) {
            continue;
        }
        let hc_depth = member.hc_mutations.mutation_events();
        let key = (hc_depth, "unpaired".to_string(), None);
        let state = states.entry(key).or_default();
        state.cells.insert(member.cell.cell_id.clone());
        let iso = if member.cell.hc.c.is_empty() {
            "unknown"
        } else {
            member.cell.hc.c.as_str()
        };
        *state.isotypes.entry(iso.to_string()).or_default() += 1;
        state.lc_name = "unpaired".to_string();
        state.hc_depth = hc_depth;
    }
    if states.len() < 2 {
        return Ok(());
    }

    const ROOT: &str = "HC NAIVE";
    let state_values: Vec<PlotState> = states.into_values().collect();
    let n = state_values.len() + 1;
    let mut labels = Vec::with_capacity(n);
    labels.push(ROOT.to_string());
    labels.extend(state_values.iter().enumerate().map(|(i, s)| {
        format!(
            "state{}|HC{}|{}|LC{}",
            i + 1,
            s.hc_depth,
            s.lc_name,
            s.lc_depth
                .map(|x| x.to_string())
                .unwrap_or_else(|| "NA".into())
        )
    }));
    let mut features = Array2::<f32>::zeros((n, 2));
    let mut groups = Vec::with_capacity(n);
    groups.push(ROOT.to_string());
    for (i, s) in state_values.iter().enumerate() {
        features[[i + 1, 0]] = s.hc_depth as f32;
        features[[i + 1, 1]] = s.lc_depth.unwrap_or(0) as f32;
        groups.push(s.lc_name.clone());
    }
    let model = ClonoMap::from_grouped_feature_matrix(labels, features, groups, 0, 1, 2)
        .map_err(|e| anyhow::anyhow!("ClonoMap plotting model failed for {}: {e}", family.name))?;

    let mut abundance = vec![1usize];
    let mut iso_cat = vec![ROOT.to_string()];
    let mut iso_mixed = vec![false];
    let mut lc_cat = vec![ROOT.to_string()];
    let mut lc_mixed = vec![false];
    let mut hc_depth = vec![None];
    let mut lc_depth = vec![None];
    let mut paired_depth = vec![None];
    for s in &state_values {
        abundance.push(s.cells.len());
        let dominant_iso = s
            .isotypes
            .iter()
            .max_by(|(ka, va), (kb, vb)| va.cmp(vb).then_with(|| kb.cmp(ka)))
            .map(|(k, _)| k.clone())
            .unwrap_or_else(|| "unknown".to_string());
        iso_cat.push(dominant_iso);
        iso_mixed.push(s.isotypes.len() > 1);
        lc_cat.push(s.lc_name.clone());
        lc_mixed.push(false);
        hc_depth.push(Some(s.hc_depth as f32));
        lc_depth.push(s.lc_depth.map(|x| x as f32));
        paired_depth.push(Some((s.hc_depth + s.lc_depth.unwrap_or(0)) as f32));
    }
    let mut lc_order: Vec<String> = family.light_clones.iter().map(|x| x.name.clone()).collect();
    lc_order.sort_by_key(|name| {
        std::cmp::Reverse(
            family
                .light_clones
                .iter()
                .find(|x| &x.name == name)
                .map(|x| x.members.len())
                .unwrap_or(0),
        )
    });
    lc_order.dedup();
    if lc_cat.iter().any(|x| x == "unpaired") {
        lc_order.push("unpaired".to_string());
    }

    // Keep filesystem names deliberately short. Family/LC identities and the
    // selection statistics belong in the plot title and tabular reports; putting
    // them into a directory component can exceed Linux NAME_MAX (typically 255).
    let dir = out.join(format!("family_{plot_index:04}"));
    create_dir_all(&dir).with_context(|| format!("create {}", dir.display()))?;
    let title = format!(
        "{} | {} cells | {} LC clones",
        family.name,
        family.members.len(),
        family.light_clones.len()
    );
    let path = |name: &str| dir.join(name).to_string_lossy().into_owned();
    model
        .tree
        .plot_rooted_annotated_cached(
            n,
            0,
            &iso_cat,
            None,
            &iso_mixed,
            &abundance,
            "IGH constant class",
            &title,
            &path("mst_rooted_isotype.svg"),
            radial_layout,
        )
        .map_err(|e| anyhow::anyhow!("isotype SVG failed for {}: {e}", family.name))?;
    model
        .tree
        .plot_rooted_annotated_cached(
            n,
            0,
            &lc_cat,
            Some(&lc_order),
            &lc_mixed,
            &abundance,
            "Light-chain clone",
            &title,
            &path("mst_rooted_light_chain.svg"),
            radial_layout,
        )
        .map_err(|e| anyhow::anyhow!("LC SVG failed for {}: {e}", family.name))?;
    model
        .tree
        .plot_rooted_continuous_cached(
            n,
            0,
            &hc_depth,
            &abundance,
            "HC mutational depth (events)",
            &title,
            &path("mst_rooted_hc_depth.svg"),
            radial_layout,
        )
        .map_err(|e| anyhow::anyhow!("HC-depth SVG failed for {}: {e}", family.name))?;
    model
        .tree
        .plot_rooted_continuous_cached(
            n,
            0,
            &lc_depth,
            &abundance,
            "LC mutational depth (events)",
            &title,
            &path("mst_rooted_lc_depth.svg"),
            radial_layout,
        )
        .map_err(|e| anyhow::anyhow!("LC-depth SVG failed for {}: {e}", family.name))?;
    model
        .tree
        .plot_rooted_continuous_cached(
            n,
            0,
            &paired_depth,
            &abundance,
            "Paired HC+LC mutational depth (events)",
            &title,
            &path("mst_rooted_paired_depth.svg"),
            radial_layout,
        )
        .map_err(|e| anyhow::anyhow!("paired-depth SVG failed for {}: {e}", family.name))?;
    Ok(())
}

fn receptor_complete(r: &Receptor) -> bool {
    !r.v.is_empty()
        && !r.j.is_empty()
        && !r.cdr3_nt.is_empty()
        && !r.naive.is_empty()
        && !r.observed.is_empty()
}

fn build_initial_families(mut cells: Vec<CellReceptor>, cfg: &FamilyConfig) -> Vec<Family> {
    cells.sort_by(|a, b| {
        (&a.hc.v, &a.hc.j, &a.hc.cdr3_nt, &a.cell_id).cmp(&(
            &b.hc.v,
            &b.hc.j,
            &b.hc.cdr3_nt,
            &b.cell_id,
        ))
    });
    let mut families: Vec<Family> = Vec::new();
    for cell in cells {
        let mut pending = Some(cell);
        for family in &mut families {
            let c = pending.take().unwrap();
            match family.add_candidate(c, cfg) {
                Ok(()) => {
                    pending = None;
                    break;
                }
                Err(c) => pending = Some(c),
            }
        }
        if let Some(seed) = pending {
            let name = format!("HC:{}:{}:CDR3:{}", seed.hc.v, seed.hc.j, seed.hc.cdr3_nt);
            families.push(Family::new(name, seed));
        }
    }
    families
}

fn reassign_cells(
    families: &mut [Family],
    rejected: Vec<(String, CellReceptor)>,
    cfg: &FamilyConfig,
    novel_v: &ReferenceModels,
) -> Vec<UnassignedCell> {
    let mut unassigned = Vec::new();
    for (from_family, rejected_cell) in rejected {
        let mut cell = Some(rejected_cell);
        for family in families.iter_mut() {
            if family.name == from_family {
                continue;
            }
            let candidate = cell.take().unwrap();
            match family.try_integrate(candidate, cfg, novel_v) {
                Ok(()) => {
                    cell = None;
                    break;
                }
                Err(candidate) => cell = Some(candidate),
            }
        }
        if let Some(cell) = cell {
            unassigned.push(UnassignedCell { cell, from_family });
        }
    }
    unassigned
}

const CELL_SCOPE_SEPARATOR: char = '\u{1f}';

fn scoped_cell_id(source: &str, cell: &str) -> String {
    format!("{source}{CELL_SCOPE_SEPARATOR}{cell}")
}

fn split_scoped_cell_id(cell: &str) -> (&str, &str) {
    cell.split_once(CELL_SCOPE_SEPARATOR)
        .unwrap_or(("input", cell))
}

fn source_label(dir: &Path) -> String {
    let name = dir.file_name().and_then(|x| x.to_str()).unwrap_or("input");
    if name == "vdj_out" {
        dir.parent()
            .and_then(Path::parent)
            .and_then(Path::file_name)
            .and_then(|x| x.to_str())
            .or_else(|| {
                dir.parent()
                    .and_then(Path::file_name)
                    .and_then(|x| x.to_str())
            })
            .unwrap_or("input")
            .to_string()
    } else {
        name.to_string()
    }
}

fn input_sources(paths: &[PathBuf]) -> Result<Vec<(String, PathBuf)>> {
    let mut seen = BTreeSet::new();
    let mut out = Vec::with_capacity(paths.len());
    for path in paths {
        let source = source_label(path);
        if !seen.insert(source.clone()) {
            bail!(
                "multiple --vdj-out inputs resolve to source label {source:?}; use distinct sample directories"
            );
        }
        out.push((source, path.clone()));
    }
    Ok(out)
}

fn read_vdj_output(dir: &Path, source: &str) -> Result<Vec<CallRow>> {
    let airr_path = dir.join("airr_rearrangements.tsv");
    let calls_path = dir.join("vdj_calls.tsv");
    if !airr_path.is_file() {
        bail!("missing {}", airr_path.display());
    }
    if !calls_path.is_file() {
        bail!("missing {}", calls_path.display());
    }

    // naive_recombination is a per-cell reconstruction span. The compact
    // recombination ID may legitimately be shared by hundreds of cells, so the
    // join key is (cell, recombination_id), never recombination_id alone.
    let reconstruction_by_cell_and_id = read_reconstructions(&calls_path)?;

    let f = File::open(&airr_path).with_context(|| format!("open {}", airr_path.display()))?;
    let mut lines = BufReader::new(f).lines();
    let header = lines.next().context("empty airr_rearrangements.tsv")??;
    let names: Vec<&str> = header.split('\t').collect();
    let index: HashMap<&str, usize> = names.iter().enumerate().map(|(i, n)| (*n, i)).collect();
    let req = |name: &str| -> Result<usize> {
        index.get(name).copied().with_context(|| {
            format!(
                "{}: missing required AIRR column {name}",
                airr_path.display()
            )
        })
    };
    let cell_i = req("cell_id")?;
    let id_i = req("lumrik_recombination_id")?;
    let chain_i = req("locus")?;
    let v_i = req("v_call")?;
    let d_i = req("d_call")?;
    let j_i = req("j_call")?;
    let c_i = req("c_call")?;
    let cdr3_i = index
        .get("cdr3")
        .copied()
        .or_else(|| index.get("junction").copied())
        .with_context(|| format!("{}: missing cdr3/junction", airr_path.display()))?;
    let observed_i = req("sequence")?;
    let productive_i = req("productive")?;
    let max_i = [
        cell_i,
        id_i,
        chain_i,
        v_i,
        d_i,
        j_i,
        c_i,
        cdr3_i,
        observed_i,
        productive_i,
    ]
    .into_iter()
    .max()
    .unwrap();

    let mut out = Vec::new();
    let mut missing_naive = 0usize;
    for (line_no, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= max_i {
            bail!(
                "{}:{} has {} fields, need at least {}",
                airr_path.display(),
                line_no + 2,
                fields.len(),
                max_i + 1
            );
        }
        let cell = fields[cell_i].trim();
        let id = fields[id_i].trim();
        let key = (cell.to_string(), id.to_string());
        let Some(reconstruction) = reconstruction_by_cell_and_id.get(&key) else {
            missing_naive += 1;
            continue;
        };
        out.push(CallRow {
            cell: scoped_cell_id(source, cell),
            source: source.to_string(),
            productive: productive_value(fields[productive_i]),
            receptor: Receptor {
                id: id.to_string(),
                chain: fields[chain_i].to_string(),
                v: fields[v_i].to_string(),
                d: fields[d_i].to_string(),
                j: fields[j_i].to_string(),
                c: fields[c_i].to_string(),
                cdr3_nt: fields[cdr3_i].to_string(),
                naive: reconstruction.naive.clone(),
                observed: reconstruction.observed_receptor.clone(),
                alternative_reconstructions: reconstruction.alternatives,
                evidence_reads: reconstruction.evidence_reads,
            },
        });
    }
    if missing_naive > 0 {
        eprintln!(
            "ClonoMap: skipped {missing_naive} AIRR rows without matching (cell, recombination_id) naive_recombination in vdj_calls.tsv"
        );
    }
    Ok(out)
}

#[derive(Debug, Clone)]
struct SelectedReconstruction {
    naive: String,
    observed_receptor: String,
    alternatives: usize,
    evidence_reads: usize,
}

#[derive(Debug, Clone)]
struct ReconstructionCandidate {
    naive: String,
    observed_rearrangement: String,
    observed_receptor: String,
    productive: bool,
    mutation_events: usize,
    support_features: usize,
    junction_support_reads: usize,
    junction_spanning_reads: usize,
    receptor_rediscovery_reads: usize,
}

/// Collapse raw sc-vdj rows to one receptor reconstruction per
/// (cell, recombination_id) before CellReceptor construction.
///
/// A compact recombination ID describes the recombination structure, but the
/// same cell can carry competing nucleotide reconstructions for that structure.
/// We deliberately resolve that ambiguity here, outside the Family model:
///   1. productive beats non-productive;
///   2. lower observed-vs-naive mutation burden wins;
///   3. stronger junction/feature/rediscovery evidence breaks ties;
///   4. longer informative alignment, then lexical sequence order make the
///      result deterministic.
fn read_reconstructions(path: &Path) -> Result<HashMap<(String, String), SelectedReconstruction>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut lines = BufReader::new(f).lines();
    let header = lines.next().context("empty vdj_calls.tsv")??;
    let names: Vec<&str> = header.split('\t').collect();
    let index: HashMap<&str, usize> = names.iter().enumerate().map(|(i, n)| (*n, i)).collect();
    let req = |name: &str| -> Result<usize> {
        index
            .get(name)
            .copied()
            .with_context(|| format!("{}: missing required column {name}", path.display()))
    };
    let cell_i = req("cell")?;
    let id_i = req("recombination_id")?;
    let naive_i = req("naive_recombination")?;
    let rearr_i = req("observed_rearrangement")?;
    let receptor_i = req("observed_receptor_sequence")?;
    let productivity_i = req("productivity_status")?;
    let support_i = req("support_features")?;
    let junction_support_i = req("junction_support_reads")?;
    let junction_spanning_i = req("junction_spanning_reads")?;
    let rediscovery_i = req("receptor_rediscovery_reads")?;
    let max_i = [
        cell_i,
        id_i,
        naive_i,
        rearr_i,
        receptor_i,
        productivity_i,
        support_i,
        junction_support_i,
        junction_spanning_i,
        rediscovery_i,
    ]
    .into_iter()
    .max()
    .unwrap();

    let mut candidates: HashMap<(String, String), Vec<ReconstructionCandidate>> = HashMap::new();
    for (line_no, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= max_i {
            bail!(
                "{}:{} has {} fields, need at least {}",
                path.display(),
                line_no + 2,
                fields.len(),
                max_i + 1
            );
        }
        let cell = fields[cell_i].trim();
        let id = fields[id_i].trim();
        let naive = fields[naive_i].trim();
        let observed_rearrangement = fields[rearr_i].trim();
        let observed_receptor = fields[receptor_i].trim();
        if cell.is_empty()
            || id.is_empty()
            || naive.is_empty()
            || observed_rearrangement.is_empty()
            || observed_receptor.is_empty()
        {
            continue;
        }
        let measurement = align_fragment(naive, observed_rearrangement);
        let mutation_events = measurement
            .as_ref()
            .map_or(usize::MAX, |m| m.mutation_events());
        candidates
            .entry((cell.to_string(), id.to_string()))
            .or_default()
            .push(ReconstructionCandidate {
                naive: naive.to_string(),
                observed_rearrangement: observed_rearrangement.to_string(),
                observed_receptor: observed_receptor.to_string(),
                productive: productive_value(fields[productivity_i]),
                mutation_events,
                support_features: parse_count(fields[support_i]),
                junction_support_reads: parse_count(fields[junction_support_i]),
                junction_spanning_reads: parse_count(fields[junction_spanning_i]),
                receptor_rediscovery_reads: parse_count(fields[rediscovery_i]),
            });
    }

    let mut out = HashMap::new();
    let mut ambiguous_keys = 0usize;
    let mut discarded_rows = 0usize;
    for (key, mut rows) in candidates {
        let alternatives = rows.len();
        if alternatives > 1 {
            ambiguous_keys += 1;
            discarded_rows += alternatives - 1;
        }
        rows.sort_by(|a, b| reconstruction_rank(a).cmp(&reconstruction_rank(b)));
        let best = rows.remove(0);
        let evidence_reads = best
            .receptor_rediscovery_reads
            .max(best.junction_spanning_reads)
            .max(best.junction_support_reads);
        out.insert(
            key,
            SelectedReconstruction {
                naive: best.naive,
                observed_receptor: best.observed_receptor,
                alternatives,
                evidence_reads,
            },
        );
    }
    eprintln!(
        "ClonoMap receptor normalization: {ambiguous_keys} cell/recombination calls had alternative reconstructions; selected one and discarded {discarded_rows} alternatives"
    );
    Ok(out)
}

fn reconstruction_rank(
    c: &ReconstructionCandidate,
) -> (
    bool,
    usize,
    std::cmp::Reverse<usize>,
    std::cmp::Reverse<usize>,
    std::cmp::Reverse<usize>,
    std::cmp::Reverse<usize>,
    std::cmp::Reverse<usize>,
    String,
    String,
) {
    let informative =
        align_fragment(&c.naive, &c.observed_rearrangement).map_or(0, |m| m.informative_pairs);
    (
        !c.productive,
        c.mutation_events,
        std::cmp::Reverse(c.junction_support_reads),
        std::cmp::Reverse(c.junction_spanning_reads),
        std::cmp::Reverse(c.support_features),
        std::cmp::Reverse(c.receptor_rediscovery_reads),
        std::cmp::Reverse(informative),
        c.naive.clone(),
        c.observed_receptor.clone(),
    )
}

fn parse_count(s: &str) -> usize {
    s.trim().parse().unwrap_or(0)
}

fn productive_value(s: &str) -> bool {
    matches!(
        s.trim().to_ascii_lowercase().as_str(),
        "productive" | "true" | "t" | "yes" | "1"
    )
}

fn fmt_usize_option(value: Option<usize>) -> String {
    value.map_or_else(|| "-".to_string(), |v| v.to_string())
}

fn fmt_float_option(value: Option<f64>) -> String {
    value.map_or_else(|| "-".into(), |x| format!("{x:.2}"))
}

fn fmt_max3(values: &[usize]) -> String {
    if values.is_empty() {
        return "[]".into();
    }
    format!(
        "[{}]",
        values
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn write_reference_candidate_report(path: &Path, models: &ReferenceModels) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "candidate_id\tchain\toriginal_v\tobservations_this_run\tevidence_reads_this_run\tresolved\tsequence"
    )?;
    for entry in models.entries() {
        let candidate = models
            .curator()
            .candidate(&entry.id)
            .expect("session candidate must exist in curator");
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            entry.id,
            entry.chain,
            entry.original_v,
            entry.observations,
            entry.evidence_reads,
            candidate.resolved.is_some(),
            String::from_utf8_lossy(&candidate.sequence)
        )?;
    }
    Ok(())
}

fn write_family_report(path: &Path, families: &[Family]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "family\tcells\tmax_hc_mutations\tlc_clones\tlc_n\tmean_lc_mutations\tsd_lc_mutations\tlc_max_1\tlc_max_2\tlc_max_3\tpearson_abundance_hc_n\tpearson_abundance_hc_r\tpearson_abundance_hc_p\tpearson_abundance_lc_n\tpearson_abundance_lc_r\tpearson_abundance_lc_p\tpearson_abundance_paired_n\tpearson_abundance_paired_r\tpearson_abundance_paired_p\tpearson_hc_lc_n\tpearson_hc_lc_r\tpearson_hc_lc_p\tmst_rooted_isotype_svg_hex\tmst_rooted_light_chain_svg_hex\tmst_rooted_hc_depth_svg_hex\tmst_rooted_lc_depth_svg_hex\tmst_rooted_paired_depth_svg_hex"
    )?;
    let mut rows: Vec<_> = families.iter().collect();
    rows.sort_by_key(|f| std::cmp::Reverse(f.members.len()));
    for f in rows {
        let r = f.mutation_report();
        let a = family_plot_analysis(f);
        let mut max3 = r.lc_mutations.max3.iter().rev();
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.family,
            r.cells,
            r.max_hc_mutations.map_or("NA".into(), |x| x.to_string()),
            r.lc_clones,
            r.lc_mutations.n,
            r.lc_mutations
                .mean
                .map_or("NA".into(), |x| format!("{x:.4}")),
            r.lc_mutations.sd.map_or("NA".into(), |x| format!("{x:.4}")),
            max3.next().map_or("NA".into(), |x| x.to_string()),
            max3.next().map_or("NA".into(), |x| x.to_string()),
            max3.next().map_or("NA".into(), |x| x.to_string()),
            a.abundance_hc.n,
            fmt_stat(a.abundance_hc.r),
            fmt_stat(a.abundance_hc.p),
            a.abundance_lc.n,
            fmt_stat(a.abundance_lc.r),
            fmt_stat(a.abundance_lc.p),
            a.abundance_paired.n,
            fmt_stat(a.abundance_paired.r),
            fmt_stat(a.abundance_paired.p),
            a.hc_lc.n,
            fmt_stat(a.hc_lc.r),
            fmt_stat(a.hc_lc.p),
            a.isotype_hex,
            a.light_chain_hex,
            a.hc_depth_hex,
            a.lc_depth_hex,
            a.paired_depth_hex,
        )?;
    }
    Ok(())
}

fn color_for(encoded: &str, cell: &str) -> String {
    encoded
        .split(';')
        .find_map(|x| {
            let (id, color) = x.split_once('=')?;
            (id == cell).then(|| color.to_string())
        })
        .unwrap_or_else(|| "NA".into())
}

fn write_cell_report(
    path: &Path,
    families: &[Family],
    cell_id_detectors: &[PrimerDetector],
) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "source\tcell\tbd_cell_id\tfamily\tlc_clone\thc_mutation_count\tlc_mutation_count\ttotal_mutation_count\thc_distance\tlc_distance\ttotal_distance\tfamily_pearson_abundance_hc_n\tfamily_pearson_abundance_hc_r\tfamily_pearson_abundance_hc_p\tfamily_pearson_abundance_lc_n\tfamily_pearson_abundance_lc_r\tfamily_pearson_abundance_lc_p\tfamily_pearson_abundance_paired_n\tfamily_pearson_abundance_paired_r\tfamily_pearson_abundance_paired_p\tfamily_pearson_hc_lc_n\tfamily_pearson_hc_lc_r\tfamily_pearson_hc_lc_p\thclc_pearson_abundance_hc_n\thclc_pearson_abundance_hc_r\thclc_pearson_abundance_hc_p\thclc_pearson_abundance_lc_n\thclc_pearson_abundance_lc_r\thclc_pearson_abundance_lc_p\thclc_pearson_abundance_paired_n\thclc_pearson_abundance_paired_r\thclc_pearson_abundance_paired_p\thclc_pearson_hc_lc_n\thclc_pearson_hc_lc_r\thclc_pearson_hc_lc_p\tmst_rooted_isotype_svg_hex\tmst_rooted_light_chain_svg_hex\tmst_rooted_hc_depth_svg_hex\tmst_rooted_lc_depth_svg_hex\tmst_rooted_paired_depth_svg_hex"
    )?;
    for family in families {
        let a = family_plot_analysis(family);
        let hc_by_cell: HashMap<&str, &AlignedCell> = family
            .members
            .iter()
            .map(|m| (m.cell.cell_id.as_str(), m))
            .collect();
        let mut emitted = BTreeSet::new();
        for lc in &family.light_clones {
            let hs = a.hclc.get(&lc.name).copied().unwrap_or_default();
            for lm in &lc.members {
                let Some(hm) = hc_by_cell.get(lm.cell.as_str()) else {
                    continue;
                };
                emitted.insert(lm.cell.clone());
                let hc_count = hm.hc_mutations.mutation_events();
                let lc_count = lm.mutations.mutation_events();
                let hc_dist = mutation_distance(&hm.hc_mutations);
                let lc_dist = mutation_distance(&lm.mutations);
                let (source, cell) = split_scoped_cell_id(&lm.cell);
                let bd_cell_id = bd_cell_id_for_seq(cell, cell_id_detectors);
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    source,
                    cell,
                    bd_cell_id,
                    family.name,
                    lc.name,
                    hc_count,
                    lc_count,
                    hc_count + lc_count,
                    hc_dist,
                    lc_dist,
                    hc_dist + lc_dist,
                    a.abundance_hc.n,
                    fmt_stat(a.abundance_hc.r),
                    fmt_stat(a.abundance_hc.p),
                    a.abundance_lc.n,
                    fmt_stat(a.abundance_lc.r),
                    fmt_stat(a.abundance_lc.p),
                    a.abundance_paired.n,
                    fmt_stat(a.abundance_paired.r),
                    fmt_stat(a.abundance_paired.p),
                    a.hc_lc.n,
                    fmt_stat(a.hc_lc.r),
                    fmt_stat(a.hc_lc.p),
                    hs.abundance_hc.n,
                    fmt_stat(hs.abundance_hc.r),
                    fmt_stat(hs.abundance_hc.p),
                    hs.abundance_lc.n,
                    fmt_stat(hs.abundance_lc.r),
                    fmt_stat(hs.abundance_lc.p),
                    hs.abundance_paired.n,
                    fmt_stat(hs.abundance_paired.r),
                    fmt_stat(hs.abundance_paired.p),
                    hs.hc_lc.n,
                    fmt_stat(hs.hc_lc.r),
                    fmt_stat(hs.hc_lc.p),
                    color_for(&a.isotype_hex, &lm.cell),
                    color_for(&a.light_chain_hex, &lm.cell),
                    color_for(&a.hc_depth_hex, &lm.cell),
                    color_for(&a.lc_depth_hex, &lm.cell),
                    color_for(&a.paired_depth_hex, &lm.cell)
                )?;
            }
        }
        for hm in &family.members {
            if emitted.contains(&hm.cell.cell_id) {
                continue;
            }
            let hc_count = hm.hc_mutations.mutation_events();
            let hc_dist = mutation_distance(&hm.hc_mutations);
            let (source, cell) = split_scoped_cell_id(&hm.cell.cell_id);
            let bd_cell_id = bd_cell_id_for_seq(cell, cell_id_detectors);
            writeln!(
                w,
                "{}\t{}\t{}\t{}\tunpaired\t{}\tNA\t{}\t{}\tNA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\tNA\tNA\t0\tNA\tNA\t0\tNA\tNA\t0\tNA\tNA\t{}\t{}\t{}\t{}\t{}",
                source,
                cell,
                bd_cell_id,
                family.name,
                hc_count,
                hc_count,
                hc_dist,
                hc_dist,
                a.abundance_hc.n,
                fmt_stat(a.abundance_hc.r),
                fmt_stat(a.abundance_hc.p),
                a.abundance_lc.n,
                fmt_stat(a.abundance_lc.r),
                fmt_stat(a.abundance_lc.p),
                a.abundance_paired.n,
                fmt_stat(a.abundance_paired.r),
                fmt_stat(a.abundance_paired.p),
                a.hc_lc.n,
                fmt_stat(a.hc_lc.r),
                fmt_stat(a.hc_lc.p),
                color_for(&a.isotype_hex, &hm.cell.cell_id),
                color_for(&a.light_chain_hex, &hm.cell.cell_id),
                color_for(&a.hc_depth_hex, &hm.cell.cell_id),
                color_for(&a.lc_depth_hex, &hm.cell.cell_id),
                color_for(&a.paired_depth_hex, &hm.cell.cell_id)
            )?;
        }
    }
    Ok(())
}

fn mutation_signature(m: &MutationMeasurement) -> String {
    let substitutions = m
        .substitution_positions
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let indels = m
        .indels
        .iter()
        .map(|x| {
            format!(
                "{}{}:{}",
                if x.inserted { "+" } else { "-" },
                x.naive_pos,
                x.len
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("S[{substitutions}]I[{indels}]")
}

fn source_counts<'a>(cells: impl Iterator<Item = &'a str>) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for cell in cells {
        let (source, _) = split_scoped_cell_id(cell);
        *counts.entry(source.to_string()).or_default() += 1;
    }
    counts
}

fn counts_text(counts: &BTreeMap<String, usize>) -> String {
    counts
        .iter()
        .map(|(source, n)| format!("{source}={n}"))
        .collect::<Vec<_>>()
        .join(";")
}

fn write_overlap_reports(out: &Path, families: &[Family]) -> Result<()> {
    let mut events: Vec<(String, String, String, String, BTreeMap<String, usize>)> = Vec::new();
    for family in families {
        let hc_counts = source_counts(family.members.iter().map(|m| m.cell.cell_id.as_str()));
        if hc_counts.len() > 1 {
            events.push((
                "HC".into(),
                family.name.clone(),
                "NA".into(),
                "NA".into(),
                hc_counts,
            ));
        }
        let hc_by_cell: HashMap<&str, &AlignedCell> = family
            .members
            .iter()
            .map(|m| (m.cell.cell_id.as_str(), m))
            .collect();
        for lc in &family.light_clones {
            let lc_counts = source_counts(lc.members.iter().map(|m| m.cell.as_str()));
            if lc_counts.len() > 1 {
                events.push((
                    "HC_LC".into(),
                    family.name.clone(),
                    lc.name.clone(),
                    "NA".into(),
                    lc_counts,
                ));
            }
            let mut signatures: BTreeMap<String, Vec<&str>> = BTreeMap::new();
            for lm in &lc.members {
                let Some(hm) = hc_by_cell.get(lm.cell.as_str()) else {
                    continue;
                };
                let signature = format!(
                    "HC:{}|LC:{}",
                    mutation_signature(&hm.hc_mutations),
                    mutation_signature(&lm.mutations)
                );
                signatures
                    .entry(signature)
                    .or_default()
                    .push(lm.cell.as_str());
            }
            for (signature, cells) in signatures {
                let counts = source_counts(cells.into_iter());
                if counts.len() > 1 {
                    events.push((
                        "HC_LC_MUTATION_SET".into(),
                        family.name.clone(),
                        lc.name.clone(),
                        signature,
                        counts,
                    ));
                }
            }
        }
    }
    if events.is_empty() {
        return Ok(());
    }

    let mut tsv = BufWriter::new(File::create(out.join("overlap_events.tsv"))?);
    writeln!(
        tsv,
        "event_type\tfamily\tlc_clone\tmutation_set\tsources\tsource_counts\ttotal_cells"
    )?;
    for (kind, family, lc, mutations, counts) in &events {
        writeln!(
            tsv,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            kind,
            family,
            lc,
            mutations,
            counts.keys().cloned().collect::<Vec<_>>().join(";"),
            counts_text(counts),
            counts.values().sum::<usize>()
        )?;
    }

    let mut md = BufWriter::new(File::create(out.join("overlap_events.md"))?);
    writeln!(md, "# Cross-source ClonoMap overlaps\n")?;
    writeln!(
        md,
        "Joint family construction found **{}** cross-source overlap events. These are derived summaries; `cells.tsv` remains the cell-level authoritative output.\n",
        events.len()
    )?;
    for (kind, family, lc, mutations, counts) in &events {
        writeln!(
            md,
            "- **{}** — `{}`{}: {} cells ({}){}",
            kind,
            family,
            if lc == "NA" {
                String::new()
            } else {
                format!(" / `{lc}`")
            },
            counts.values().sum::<usize>(),
            counts_text(counts),
            if mutations == "NA" {
                String::new()
            } else {
                format!("; mutation set `{mutations}`")
            }
        )?;
    }
    Ok(())
}

fn bd_cell_id_for_seq(cell: &str, detectors: &[PrimerDetector]) -> String {
    detectors
        .iter()
        .find_map(|detector| detector.cell_id_for_seq(cell.as_bytes()))
        .map_or_else(|| "NA".to_string(), |id| id.to_string())
}

fn fmt_stat(x: Option<f64>) -> String {
    x.map_or("NA".into(), |v| format!("{v:.6}"))
}

fn write_unassigned(path: &Path, rows: &[UnassignedCell]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "source\tcell\treceptor_id\tfrom_family")?;
    for r in rows {
        let (source, cell) = split_scoped_cell_id(&r.cell.cell_id);
        writeln!(
            w,
            "{}\t{}\t{}\t{}",
            source, cell, r.cell.hc.id, r.from_family
        )?;
    }
    Ok(())
}
