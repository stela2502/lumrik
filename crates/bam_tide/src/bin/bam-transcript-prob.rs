use anyhow::{bail, Context, Result};
use bam_tide::core::ref_block::record_to_blocks;
use clap::Parser;
use gtf_splice_index::{MatchClass, MatchOptions, SpliceIndex, SplicedRead, Strand};
use rayon::prelude::*;
use ommverse::Ommverse;
use rust_htslib::bam::{self, Read};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
#[command(author, version, about = "Estimate transcript usage from paired genomic BAM fragments")]
struct Cli {
    #[arg(long, required = true)]
    bam: Vec<PathBuf>,
    #[arg(long)]
    index: PathBuf,
    #[arg(long)]
    gene: String,
    #[arg(long, default_value = "transcript_probability.tsv")]
    out: PathBuf,
    /// Optional fragment-by-transcript evidence table. Not written unless requested.
    #[arg(long)]
    evidence_out: Option<PathBuf>,
    /// Optional exhaustive fragment x transcript geometry audit. This writes
    /// every transcript of the requested gene, including rejected ones.
    #[arg(long)]
    audit_out: Option<PathBuf>,
    #[arg(long, default_value_t = 0)]
    min_mapq: u8,
    #[arg(long, default_value_t = 200)]
    em_iterations: usize,
    #[arg(long, default_value_t = 1e-8)]
    em_epsilon: f64,
    /// Number of fragment bootstrap replicates used to estimate per-sample
    /// transcript-fraction uncertainty. Set to 0 to disable bootstrapping.
    #[arg(long, default_value_t = 100)]
    bootstrap_replicates: usize,
    /// Deterministic seed for fragment bootstrap resampling.
    #[arg(long, default_value_t = 0x4d45495332_u64)]
    bootstrap_seed: u64,
    /// Number of randomized EM starts used as a practical transcript-
    /// identifiability diagnostic. Set to 0 to disable. Unlike the fragment
    /// bootstrap, this keeps the observed data fixed and asks whether different
    /// valid starting mixtures converge to the same transcript solution.
    #[arg(long, default_value_t = 100)]
    identifiability_starts: usize,
    /// Deterministic seed for randomized EM starts.
    #[arg(long, default_value_t = 0x4944454e54_u64)]
    identifiability_seed: u64,
    /// Define exactly two comparison groups as LABEL=SUBSTRING. Samples are
    /// assigned by matching SUBSTRING against the BAM file name. When two
    /// groups are supplied, a summary table, two SVG figures and a Markdown
    /// report are written next to --out.
    #[arg(long = "group", value_name = "LABEL=SUBSTRING")]
    groups: Vec<String>,
    /// Minimum mean transcript fraction in either group for inclusion in plots.
    #[arg(long, default_value_t = 0.02)]
    plot_min_fraction: f64,
    /// Optional InterPro-enriched Ommverse index. When supplied in two-group
    /// report mode, a genomic transcript-structure SVG is produced with
    /// transcript-specific InterPro protein features projected through the CDS
    /// back onto the coding exons.
    #[arg(long)]
    ommverse: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
enum EvidenceKind {
    #[default]
    ExonCompatible,
    SpliceCompatible,
    ExactJunctionChain,
}

impl EvidenceKind {
    const fn as_str(self) -> &'static str {
        match self {
            Self::ExonCompatible => "exon_compatible",
            Self::SpliceCompatible => "splice_compatible",
            Self::ExactJunctionChain => "exact_junction_chain",
        }
    }

    const fn is_splice(self) -> bool {
        matches!(self, Self::SpliceCompatible | Self::ExactJunctionChain)
    }
}


#[derive(Clone, Debug)]
struct SampleEstimate {
    sample: String,
    tx_id: usize,
    fraction: f64,
    bootstrap: BootstrapSummary,
    identifiability: IdentifiabilitySummary,
    unique_fragments: usize,
    splice_informative: usize,
    exon_only: usize,
    total_fragments: usize,
}

#[derive(Clone, Debug)]
struct GroupSpec {
    label: String,
    pattern: String,
}

#[derive(Debug)]
struct ReadEvidence {
    tx: HashMap<usize, EvidenceKind>,
    blocks: Vec<gtf_splice_index::RefBlock>,
}

#[derive(Debug, Default)]
struct PendingFragment {
    mates: Vec<ReadEvidence>,
}

#[derive(Debug)]
struct FragmentEvidence {
    read_id: Vec<u8>,
    tx: Vec<(usize, EvidenceKind)>,
    mate_blocks: Vec<Vec<gtf_splice_index::RefBlock>>,
}

fn main() -> Result<()> {
    let args = Cli::parse();
    let idx = SpliceIndex::load(&args.index)
        .with_context(|| format!("reading splice index {}", args.index.display()))?;

    let gene = idx.genes.iter().find(|g| g.names.iter().any(|n| n == &args.gene))
        .with_context(|| format!("gene {} not found in splice index", args.gene))?;
    let gene_tx: HashSet<usize> = gene.transcript_ids().iter().copied().collect();
    if gene_tx.is_empty() { bail!("gene {} has no transcripts", args.gene); }

    let mut span: Option<(usize, u32, u32)> = None;
    for &tx_id in &gene_tx {
        let tx = &idx.transcripts[tx_id];
        let Some(first) = tx.exons().first() else { continue };
        let Some(last) = tx.exons().last() else { continue };
        span = Some(match span {
            None => (tx.chr_id, first.start, last.end),
            Some((chr, start, end)) => {
                if chr != tx.chr_id { bail!("gene {} spans multiple chromosomes", args.gene); }
                (chr, start.min(first.start), end.max(last.end))
            }
        });
    }
    let (chr_id, start, end) = span.context("gene has no exon span")?;
    let chr_name = idx.chr_names.get(chr_id).context("invalid chromosome id in splice index")?;

    let mut out = BufWriter::new(File::create(&args.out)?);
    let mut evidence_out = args.evidence_out.as_ref()
        .map(File::create)
        .transpose()?
        .map(BufWriter::new);
    let mut audit_out = args.audit_out.as_ref()
        .map(File::create)
        .transpose()?
        .map(BufWriter::new);
    writeln!(out, "sample\tgene\ttranscript\ttranscript_aliases\texpected_fragments\tfraction\tbootstrap_mean_fraction\tbootstrap_sd_fraction\tbootstrap_ci025\tbootstrap_ci975\tbootstrap_nonzero_fraction\tmultistart_mean_fraction\tmultistart_sd_fraction\tmultistart_p025\tmultistart_p975\tmultistart_range\tunique_fragments\tsplice_informative\texon_only\ttotal_fragments")?;
    if let Some(evidence) = evidence_out.as_mut() {
        writeln!(evidence, "sample\tgene\tread_id\ttranscript\ttranscript_aliases\tevidence\tcompatible_transcripts\tposterior")?;
    }
    if let Some(audit) = audit_out.as_mut() {
        writeln!(audit, "sample\tgene\tread_id\ttranscript\ttranscript_aliases\tfragment_compatible\tmate\tread_start\tread_end\tread_blocks\tread_junctions\ttranscript_start\ttranscript_end\ttranscript_read_start\ttranscript_read_end\tclass\texonic_bases\tintronic_bases\tmatched_junctions\tunmatched_junctions")?;
    }

    let mut sample_estimates: Vec<SampleEstimate> = Vec::new();

    for bam_path in &args.bam {
        let fragments = collect_fragments(&idx, &gene_tx, chr_name, start, end, bam_path, args.min_mapq)?;
        let estimates = em(&fragments, &gene_tx, args.em_iterations, args.em_epsilon);
        let bootstrap = bootstrap_fractions(
            &fragments,
            &gene_tx,
            args.em_iterations,
            args.em_epsilon,
            args.bootstrap_replicates,
            args.bootstrap_seed,
        );
        let identifiability = multistart_identifiability(
            &fragments,
            &gene_tx,
            args.em_iterations,
            args.em_epsilon,
            args.identifiability_starts,
            args.identifiability_seed,
        );
        let total = fragments.len();
        let sample = bam_path.file_name().and_then(|x| x.to_str()).unwrap_or("sample");

        for &tx_id in gene.transcript_ids() {
            let expected = estimates.get(&tx_id).copied().unwrap_or(0.0);
            let fraction = if total > 0 { expected / total as f64 } else { 0.0 };
            let unique = fragments.iter().filter(|f| f.tx.len() == 1 && f.tx[0].0 == tx_id).count();
            let splice = fragments.iter().filter(|f| f.tx.iter().any(|(t, k)| *t == tx_id && k.is_splice())).count();
            let exon = fragments.iter().filter(|f| f.tx.iter().any(|(t, k)| *t == tx_id && *k == EvidenceKind::ExonCompatible)).count();
            let tx_name = idx.transcript_name(tx_id).unwrap_or("NA");
            let aliases = idx.transcripts[tx_id].names.join(";");
            let b = bootstrap.get(&tx_id).copied().unwrap_or_default();
            let ident = identifiability.get(&tx_id).copied().unwrap_or_default();
            writeln!(out, "{sample}\t{}\t{tx_name}\t{aliases}\t{expected:.6}\t{fraction:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.6}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{unique}\t{splice}\t{exon}\t{total}",
                args.gene, b.mean, b.sd, b.ci025, b.ci975, b.nonzero_fraction,
                ident.mean, ident.sd, ident.p025, ident.p975, ident.range)?;
            sample_estimates.push(SampleEstimate {
                sample: sample.to_string(), tx_id, fraction, bootstrap: b, identifiability: ident,
                unique_fragments: unique, splice_informative: splice,
                exon_only: exon, total_fragments: total,
            });
        }

        if let Some(evidence) = evidence_out.as_mut() {
            for f in &fragments {
                let denom: f64 = f.tx.iter().map(|(t, _)| estimates.get(t).copied().unwrap_or(0.0)).sum();
                let read_id = String::from_utf8_lossy(&f.read_id);
                for &(tx_id, kind) in &f.tx {
                    let posterior = if denom > 0.0 {
                        estimates.get(&tx_id).copied().unwrap_or(0.0) / denom
                    } else {
                        0.0
                    };
                    let tx_name = idx.transcript_name(tx_id).unwrap_or("NA");
                    let aliases = idx.transcripts[tx_id].names.join(";");
                    writeln!(evidence, "{sample}\t{}\t{read_id}\t{tx_name}\t{aliases}\t{}\t{}\t{posterior:.8}",
                        args.gene, kind.as_str(), f.tx.len())?;
                }
            }
        }

        if let Some(audit) = audit_out.as_mut() {
            write_audit(audit, &idx, gene.transcript_ids(), chr_id, sample, &args.gene, &fragments)?;
        }
    }
    out.flush()?;
    if !args.groups.is_empty() {
        let groups = parse_groups(&args.groups)?;
        write_comparison_outputs(&args, &idx, &sample_estimates, &groups)?;
    }
    Ok(())
}


fn parse_groups(raw: &[String]) -> Result<[GroupSpec; 2]> {
    if raw.len() != 2 {
        bail!("--group must be supplied exactly twice, e.g. --group LACZ=LACZ --group MEIS2-g3=MEIS2-g3");
    }
    let parse = |value: &str| -> Result<GroupSpec> {
        let (label, pattern) = value.split_once('=')
            .with_context(|| format!("invalid --group {value:?}; expected LABEL=SUBSTRING"))?;
        if label.is_empty() || pattern.is_empty() { bail!("group label and substring must not be empty"); }
        Ok(GroupSpec { label: label.to_string(), pattern: pattern.to_string() })
    };
    Ok([parse(&raw[0])?, parse(&raw[1])?])
}

fn output_prefix(out: &Path) -> PathBuf {
    let parent = out.parent().unwrap_or_else(|| Path::new("."));
    let name = out.file_name().and_then(|x| x.to_str()).unwrap_or("transcript_probability.tsv");
    let stem = name.strip_suffix(".tsv").unwrap_or(name);
    parent.join(stem)
}

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
        .replace('"', "&quot;").replace('\'', "&apos;")
}

fn sample_name(raw: &str) -> String {
    raw.strip_suffix(".markdup.sorted.bam").unwrap_or(raw).to_string()
}

fn write_comparison_outputs(
    args: &Cli,
    idx: &SpliceIndex,
    rows: &[SampleEstimate],
    groups: &[GroupSpec; 2],
) -> Result<()> {
    let mut assigned: HashMap<String, usize> = HashMap::new();
    for row in rows {
        if assigned.contains_key(&row.sample) { continue; }
        let hits: Vec<usize> = groups.iter().enumerate()
            .filter_map(|(i, g)| row.sample.contains(&g.pattern).then_some(i)).collect();
        if hits.len() != 1 {
            bail!("sample {} matched {} comparison groups; each sample must match exactly one --group substring", row.sample, hits.len());
        }
        assigned.insert(row.sample.clone(), hits[0]);
    }
    let mut group_samples: [Vec<String>; 2] = [Vec::new(), Vec::new()];
    for (sample, &g) in &assigned { group_samples[g].push(sample.clone()); }
    for samples in &mut group_samples { samples.sort(); }
    if group_samples.iter().any(|x| x.is_empty()) { bail!("both comparison groups must contain at least one sample"); }

    let prefix = output_prefix(&args.out);
    let table_path = PathBuf::from(format!("{}.comparison.tsv", prefix.display()));
    let bootstrap_svg = PathBuf::from(format!("{}.usage_bootstrap.svg", prefix.display()));
    let shift_svg = PathBuf::from(format!("{}.usage_shift.svg", prefix.display()));
    let structure_svg = PathBuf::from(format!("{}.transcript_structure.svg", prefix.display()));
    let md_path = PathBuf::from(format!("{}.report.md", prefix.display()));

    let mut tx_ids: Vec<usize> = rows.iter().map(|r| r.tx_id).collect();
    tx_ids.sort_unstable(); tx_ids.dedup();
    tx_ids.sort_by(|a,b| {
        let max_mean = |tx: usize| -> f64 {
            (0..2).map(|g| {
                let vals: Vec<f64> = rows.iter().filter(|r| r.tx_id == tx && assigned.get(&r.sample) == Some(&g)).map(|r| r.fraction).collect();
                if vals.is_empty() { 0.0 } else { vals.iter().sum::<f64>() / vals.len() as f64 }
            }).fold(0.0, f64::max)
        };
        max_mean(*b).total_cmp(&max_mean(*a))
    });
    let plotted: Vec<usize> = tx_ids.iter().copied().filter(|&tx| {
        (0..2).any(|g| {
            let vals: Vec<f64> = rows.iter().filter(|r| r.tx_id == tx && assigned.get(&r.sample) == Some(&g)).map(|r| r.fraction).collect();
            !vals.is_empty() && vals.iter().sum::<f64>() / vals.len() as f64 >= args.plot_min_fraction
        })
    }).collect();

    let mut table = BufWriter::new(File::create(&table_path)?);
    writeln!(table, "gene\ttranscript\ttranscript_aliases\tgroup\tn\tmean_fraction\tsd_between_samples\tmean_bootstrap_sd\tmean_multistart_sd\tmean_multistart_range\tmean_unique_fragments\tmean_splice_informative\tmean_exon_only")?;
    for &tx in &tx_ids {
        for g in 0..2 {
            let rs: Vec<&SampleEstimate> = rows.iter().filter(|r| r.tx_id == tx && assigned.get(&r.sample) == Some(&g)).collect();
            let vals: Vec<f64> = rs.iter().map(|r| r.fraction).collect();
            let mean = vals.iter().sum::<f64>() / vals.len() as f64;
            let sd = if vals.len()>1 { (vals.iter().map(|x|(x-mean).powi(2)).sum::<f64>()/(vals.len()-1) as f64).sqrt() } else { 0.0 };
            let avg = |f: fn(&SampleEstimate)->f64| rs.iter().map(|r|f(r)).sum::<f64>()/rs.len() as f64;
            let name = idx.transcript_name(tx).unwrap_or("NA");
            let aliases = idx.transcripts[tx].names.join(";");
            writeln!(table, "{}\t{name}\t{aliases}\t{}\t{}\t{mean:.8}\t{sd:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.3}\t{:.3}\t{:.3}", args.gene, groups[g].label, rs.len(),
                avg(|r|r.bootstrap.sd), avg(|r|r.identifiability.sd), avg(|r|r.identifiability.range), avg(|r|r.unique_fragments as f64), avg(|r|r.splice_informative as f64), avg(|r|r.exon_only as f64))?;
        }
    }

    write_bootstrap_svg(&bootstrap_svg, args, idx, rows, &assigned, groups, &group_samples, &plotted)?;
    write_shift_svg(&shift_svg, args, idx, rows, &assigned, groups, &group_samples, &plotted)?;
    let structure_written = if let Some(ommverse_path) = args.ommverse.as_deref() {
        let omm = Ommverse::load(ommverse_path)
            .with_context(|| format!("loading Ommverse {}", ommverse_path.display()))?;
        write_structure_svg(&structure_svg, args, idx, rows, &assigned, groups, &plotted, &omm)?;
        true
    } else { false };
    write_report(&md_path, args, idx, rows, &assigned, groups, &group_samples, &plotted, &table_path, &bootstrap_svg, &shift_svg, structure_written.then_some(structure_svg.as_path()))?;
    eprintln!("comparison outputs: {} {} {} {}{}", table_path.display(), bootstrap_svg.display(), shift_svg.display(), md_path.display(), if structure_written { format!(" {}", structure_svg.display()) } else { String::new() });
    Ok(())
}

const BLUE: &str = "#2f6fbb";
const PURPLES: [&str; 3] = ["#5b3c88", "#7b5aa6", "#a184c2"];
const ORANGES: [&str; 3] = ["#d97706", "#f59e0b", "#fbbf24"];

fn sample_colour(group: usize, sample_index: usize) -> &'static str {
    let palette = if group == 0 { &PURPLES } else { &ORANGES };
    palette[sample_index % palette.len()]
}

fn svg_header(w: f64, h: f64) -> String {
    format!(r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial,sans-serif;fill:#222}} .small{{font-size:11px}} .label{{font-size:12px}} .title{{font-size:18px;font-weight:bold}} .grid{{stroke:#ddd;stroke-width:1}} </style>"#)
}

fn write_bootstrap_svg(path: &Path, args: &Cli, idx: &SpliceIndex, rows: &[SampleEstimate], _assigned: &HashMap<String,usize>, groups: &[GroupSpec;2], group_samples: &[Vec<String>;2], plotted: &[usize]) -> Result<()> {
    let w=1180.0; let left=250.0; let right=270.0; let top=80.0; let row_h=78.0; let bottom=75.0;
    let h=top+row_h*plotted.len() as f64+bottom;
    let maxv=plotted.iter().flat_map(|&tx| rows.iter().filter(move |r| r.tx_id==tx).map(|r| r.bootstrap.ci975.max(r.identifiability.p975).max(r.fraction))).fold(0.0,f64::max).max(0.01)*1.08;
    let plot_w=w-left-right; let x=|v:f64| left+plot_w*(v/maxv);
    let mut s=svg_header(w,h);
    s.push_str(&format!(r#"<text x="{left}" y="32" class="title">{} transcript usage by sample</text><text x="{left}" y="52" class="label">Dots are biological samples; coloured bars are 95% fragment-bootstrap intervals; grey bars show randomized-EM identifiability spread</text>"#, xml_escape(&args.gene)));
    for (ri,&tx) in plotted.iter().enumerate() {
        let y=top+ri as f64*row_h+row_h/2.0;
        let tx_name=idx.transcript_name(tx).unwrap_or("NA");
        let alias=idx.transcripts[tx].names.iter().find(|x|x.starts_with("ENST")).map(String::as_str).unwrap_or("");
        s.push_str(&format!(r#"<line x1="{left}" y1="{y:.1}" x2="{:.1}" y2="{y:.1}" class="grid"/><text x="{:.1}" y="{:.1}" text-anchor="end" class="label">{}</text><text x="{:.1}" y="{:.1}" text-anchor="end" class="small">{}</text>"#, w-right, left-12.0,y-4.0,xml_escape(tx_name),left-12.0,y+13.0,xml_escape(alias)));
        for g in 0..2 {
            let gy=y + if g==0 {-12.0} else {12.0};
            for (si,sample) in group_samples[g].iter().enumerate() {
                if let Some(r)=rows.iter().find(|r|r.tx_id==tx && &r.sample==sample) {
                    let dy=(si as f64-(group_samples[g].len() as f64-1.0)/2.0)*5.0;
                    let cy=gy+dy; let lo=x(r.bootstrap.ci025); let hi=x(r.bootstrap.ci975); let cx=x(r.fraction);
                    let ilo=x(r.identifiability.p025); let ihi=x(r.identifiability.p975);
                    let c=sample_colour(g,si);
                    s.push_str(&format!(r##"<line x1="{ilo:.1}" y1="{cy:.1}" x2="{ihi:.1}" y2="{cy:.1}" stroke="#999" stroke-width="5" opacity="0.35"/><line x1="{lo:.1}" y1="{cy:.1}" x2="{hi:.1}" y2="{cy:.1}" stroke="{c}" stroke-width="1.5"/><line x1="{lo:.1}" y1="{:.1}" x2="{lo:.1}" y2="{:.1}" stroke="{c}"/><line x1="{hi:.1}" y1="{:.1}" x2="{hi:.1}" y2="{:.1}" stroke="{c}"/><circle cx="{cx:.1}" cy="{cy:.1}" r="4.5" fill="{c}"/>"##,cy-4.0,cy+4.0,cy-4.0,cy+4.0));
                }
            }
        }
    }
    for i in 0..=4 { let v=maxv*i as f64/4.0; let xx=x(v); s.push_str(&format!(r#"<text x="{xx}" y="{:.1}" text-anchor="middle" class="small">{:.1}%</text>"#,h-35.0,v*100.0)); }
    let lx=w-right+30.0; let mut ly=top;
    s.push_str(&format!(r#"<text x="{lx}" y="{ly}" class="label" font-weight="bold">Biological sample</text>"#)); ly+=22.0;
    for g in 0..2 { for (si,sample) in group_samples[g].iter().enumerate() { let c=sample_colour(g,si); s.push_str(&format!(r#"<circle cx="{lx}" cy="{ly}" r="5" fill="{c}"/><text x="{:.1}" y="{:.1}" class="small">{}</text>"#,lx+12.0,ly+4.0,xml_escape(&sample_name(sample)))); ly+=20.0; } ly+=8.0; }
    ly+=4.0;
    s.push_str(&format!(r#"<text x="{lx}" y="{ly}" class="label" font-weight="bold">What the marks mean</text>"#)); ly+=24.0;
    s.push_str(&format!(r##"<line x1="{lx}" y1="{ly}" x2="{:.1}" y2="{ly}" stroke="#999" stroke-width="5" opacity="0.35"/><text x="{:.1}" y="{:.1}" class="small" font-weight="bold">Model identifiability</text><text x="{lx}" y="{:.1}" class="small">Spread across randomized EM starts.</text><text x="{lx}" y="{:.1}" class="small">Wide = transcript fraction is not uniquely</text><text x="{lx}" y="{:.1}" class="small">constrained by the observed read classes.</text>"##,lx+32.0,lx+42.0,ly+4.0,ly+19.0,ly+34.0,ly+49.0)); ly+=72.0;
    let example_colour=sample_colour(0,0);
    s.push_str(&format!(r#"<line x1="{lx}" y1="{ly}" x2="{:.1}" y2="{ly}" stroke="{example_colour}" stroke-width="1.5"/><line x1="{lx}" y1="{:.1}" x2="{lx}" y2="{:.1}" stroke="{example_colour}"/><line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{example_colour}"/><text x="{:.1}" y="{:.1}" class="small" font-weight="bold">Sampling uncertainty</text><text x="{lx}" y="{:.1}" class="small">95% fragment-bootstrap interval.</text><text x="{lx}" y="{:.1}" class="small">Wide = estimate depends strongly on which</text><text x="{lx}" y="{:.1}" class="small">fragments were sampled.</text>"#,lx+32.0,ly-4.0,ly+4.0,lx+32.0,ly-4.0,lx+32.0,ly+4.0,lx+42.0,ly+4.0,ly+19.0,ly+34.0,ly+49.0)); ly+=72.0;
    s.push_str(&format!(r#"<circle cx="{lx}" cy="{ly}" r="5" fill="{example_colour}"/><text x="{:.1}" y="{:.1}" class="small" font-weight="bold">Estimated transcript usage</text><text x="{lx}" y="{:.1}" class="small">Dot position is the fitted EM fraction.</text>"#,lx+12.0,ly+4.0,ly+19.0)); ly+=46.0;
    s.push_str(&format!(r#"<text x="{lx}" y="{ly}" class="small" font-weight="bold">Best supported:</text><text x="{lx}" y="{:.1}" class="small">narrow grey + narrow coloured interval.</text>"#,ly+16.0));
    s.push_str("</svg>"); std::fs::write(path,s)?; Ok(())
}

fn write_shift_svg(path: &Path, args: &Cli, idx: &SpliceIndex, rows: &[SampleEstimate], assigned: &HashMap<String,usize>, groups: &[GroupSpec;2], group_samples: &[Vec<String>;2], plotted: &[usize]) -> Result<()> {
    let w=1180.0; let left=210.0; let right=280.0; let top=80.0; let row_h=42.0; let bottom=75.0; let h=top+row_h*plotted.len() as f64+bottom;
    let maxv=plotted.iter().flat_map(|&tx|rows.iter().filter(move|r|r.tx_id==tx).map(|r|r.fraction)).fold(0.0,f64::max).max(0.01)*1.12;
    let plot_w=w-left-right; let x=|v:f64| left+plot_w*(v/maxv);
    let mut s=svg_header(w,h); s.push_str(&format!(r#"<text x="{left}" y="32" class="title">{} transcript usage shift</text><text x="{left}" y="52" class="label">Blue points and lines = group means; coloured points = biological samples</text>"#,xml_escape(&args.gene)));
    for (ri,&tx) in plotted.iter().enumerate() { let y=top+ri as f64*row_h+row_h/2.0; let name=idx.transcript_name(tx).unwrap_or("NA"); s.push_str(&format!(r#"<text x="{:.1}" y="{:.1}" text-anchor="end" class="label">{}</text><line x1="{left}" y1="{:.1}" x2="{:.1}" y2="{:.1}" class="grid"/>"#,left-10.0,y+4.0,xml_escape(name),y,w-right,y));
        let mut means=[0.0;2]; for g in 0..2 { let rs:Vec<&SampleEstimate>=rows.iter().filter(|r|r.tx_id==tx&&assigned.get(&r.sample)==Some(&g)).collect(); means[g]=rs.iter().map(|r|r.fraction).sum::<f64>()/rs.len() as f64; }
        s.push_str(&format!(r#"<line x1="{:.1}" y1="{y}" x2="{:.1}" y2="{y}" stroke="{BLUE}" stroke-width="2"/><circle cx="{:.1}" cy="{y}" r="5" fill="{BLUE}"/><circle cx="{:.1}" cy="{y}" r="5" fill="{BLUE}"/>"#,x(means[0]),x(means[1]),x(means[0]),x(means[1])));
        for g in 0..2 { let dy=if g==0{-8.0}else{8.0}; for (si,sample) in group_samples[g].iter().enumerate() { if let Some(r)=rows.iter().find(|r|r.tx_id==tx&&&r.sample==sample) { let c=sample_colour(g,si); s.push_str(&format!(r#"<circle cx="{:.1}" cy="{:.1}" r="3.5" fill="{c}"/>"#,x(r.fraction),y+dy)); } } }
    }
    for i in 0..=4 { let v=maxv*i as f64/4.0; let xx=x(v); s.push_str(&format!(r#"<text x="{xx}" y="{:.1}" text-anchor="middle" class="small">{:.1}%</text>"#,h-35.0,v*100.0)); }
    let lx=w-right+30.0; let mut ly=top;
    s.push_str(&format!(r#"<text x="{lx}" y="{ly}" class="label" font-weight="bold">What the marks mean</text>"#)); ly+=24.0;
    s.push_str(&format!(r#"<line x1="{lx}" y1="{ly}" x2="{:.1}" y2="{ly}" stroke="{BLUE}" stroke-width="2"/><circle cx="{:.1}" cy="{ly}" r="5" fill="{BLUE}"/><text x="{:.1}" y="{:.1}" class="small" font-weight="bold">Group mean</text><text x="{lx}" y="{:.1}" class="small">Blue line connects the two group means.</text>"#,lx+24.0,lx+12.0,lx+34.0,ly+4.0,ly+19.0)); ly+=46.0;
    s.push_str(&format!(r#"<text x="{lx}" y="{ly}" class="label" font-weight="bold">Biological sample</text>"#)); ly+=22.0;
    for g in 0..2 { for (si,sample) in group_samples[g].iter().enumerate() { let c=sample_colour(g,si); s.push_str(&format!(r#"<circle cx="{lx}" cy="{ly}" r="5" fill="{c}"/><text x="{:.1}" y="{:.1}" class="small">{}</text>"#,lx+12.0,ly+4.0,xml_escape(&sample_name(sample)))); ly+=20.0; } ly+=8.0; }
    ly+=4.0;
    s.push_str(&format!(r#"<text x="{lx}" y="{ly}" class="small" font-weight="bold">Uncertainty is intentionally omitted here.</text><text x="{lx}" y="{:.1}" class="small">See the transcript-usage plot for sampling</text><text x="{lx}" y="{:.1}" class="small">uncertainty and model identifiability.</text>"#,ly+16.0,ly+31.0));
    s.push_str("</svg>"); std::fs::write(path,s)?; Ok(())
}

fn write_structure_svg(
    path: &Path,
    args: &Cli,
    idx: &SpliceIndex,
    rows: &[SampleEstimate],
    assigned: &HashMap<String, usize>,
    groups: &[GroupSpec; 2],
    plotted: &[usize],
    omm: &Ommverse,
) -> Result<()> {
    if plotted.is_empty() { return Ok(()); }
    let mut gstart = u32::MAX;
    let mut gend = 0u32;
    for &tx_id in plotted {
        for b in idx.transcripts[tx_id].exons() {
            gstart = gstart.min(b.start);
            gend = gend.max(b.end);
        }
    }
    if gstart >= gend { return Ok(()); }

    let w = 1500.0;
    let left = 250.0;
    let right = 360.0;
    let top = 95.0;
    let row_h = 112.0;
    let bottom = 70.0;
    let h = top + row_h * plotted.len() as f64 + bottom;
    let plot_w = w - left - right;
    let gx = |p: u32| left + plot_w * ((p.saturating_sub(gstart)) as f64 / (gend - gstart) as f64);
    let mut s = svg_header(w, h);
    s.push_str(&format!(r##"<text x="{left}" y="30" class="title">{} transcript structure and InterPro architecture</text><text x="{left}" y="51" class="label">InterPro amino-acid features are projected through each transcript CDS onto the same genomic coordinate axis; introns are never painted as protein.</text><text x="{left}" y="70" class="small">Genomic span: {}–{} (0-based half-open); exon geometry is drawn to genomic scale.</text>"##, xml_escape(&args.gene), gstart, gend));

    for (ri, &tx_id) in plotted.iter().enumerate() {
        let tx = &idx.transcripts[tx_id];
        let y = top + ri as f64 * row_h + 24.0;
        let name = idx.transcript_name(tx_id).unwrap_or("NA");
        let enst = tx.names.iter().find(|x| x.starts_with("ENST")).map(String::as_str).unwrap_or("");
        s.push_str(&format!(r##"<text x="{:.1}" y="{:.1}" text-anchor="end" class="label" font-weight="bold">{}</text><text x="{:.1}" y="{:.1}" text-anchor="end" class="small">{}</text>"##, left-12.0, y+4.0, xml_escape(name), left-12.0, y+19.0, xml_escape(enst)));
        s.push_str(&format!(r##"<line x1="{:.1}" y1="{y:.1}" x2="{:.1}" y2="{y:.1}" stroke="#777" stroke-width="1"/>"##, gx(tx.exons().first().map(|b|b.start).unwrap_or(gstart)), gx(tx.exons().last().map(|b|b.end).unwrap_or(gend))));
        for b in tx.exons() {
            let x1=gx(b.start); let x2=gx(b.end);
            s.push_str(&format!(r##"<rect x="{x1:.1}" y="{:.1}" width="{:.1}" height="12" fill="#c7c7c7" stroke="#666" stroke-width="0.7"/>"##, y-6.0, (x2-x1).max(1.0)));
        }
        if let Some((cs,ce))=tx.cds_span() {
            for exon in tx.exons() {
                let a=exon.start.max(cs); let b=exon.end.min(ce);
                if a<b { let x1=gx(a); let x2=gx(b); s.push_str(&format!(r##"<rect x="{x1:.1}" y="{:.1}" width="{:.1}" height="8" fill="#444"/>"##, y-4.0,(x2-x1).max(1.0))); }
            }
        }

        let mut proteins = Vec::new();
        let mut seen_proteins = HashSet::new();
        for alias in &tx.names {
            for protein in omm.proteins_for_transcript(alias) {
                if seen_proteins.insert(protein.accession.clone()) { proteins.push(protein); }
            }
        }
        let mut feature_y = y + 25.0;
        let mut seen_features = HashSet::new();
        for protein in proteins {
            for feature in protein.features.iter().filter(|f| f.source_db == "InterPro") {
                let Some((aa0,aa1)) = feature.protein_range else { continue; };
                let blocks=tx.protein_range_to_genomic_blocks(aa0,aa1);
                if blocks.is_empty() { continue; }
                let key=(protein.accession.clone(),feature.label.clone(),aa0,aa1);
                if !seen_features.insert(key) { continue; }
                let first_x=gx(blocks.first().unwrap().start); let last_x=gx(blocks.last().unwrap().end);
                if blocks.len()>1 { s.push_str(&format!(r##"<line x1="{first_x:.1}" y1="{feature_y:.1}" x2="{last_x:.1}" y2="{feature_y:.1}" stroke="#777" stroke-width="1" stroke-dasharray="3,3"/>"##)); }
                for b in &blocks { let x1=gx(b.start); let x2=gx(b.end); s.push_str(&format!(r##"<rect x="{x1:.1}" y="{:.1}" width="{:.1}" height="8" rx="2" fill="#2f6fbb"/>"##,feature_y-4.0,(x2-x1).max(1.5))); }
                let label = if feature.description.is_empty() { feature.label.clone() } else { format!("{} · {}", feature.label, feature.description) };
                s.push_str(&format!(r##"<text x="{:.1}" y="{:.1}" class="small">{}</text>"##, last_x+6.0, feature_y+4.0, xml_escape(&label)));
                feature_y += 14.0;
                if feature_y > y + 52.0 { break; }
            }
            if feature_y > y + 52.0 { break; }
        }
        if seen_features.is_empty() {
            s.push_str(&format!(r##"<text x="{left}" y="{:.1}" class="small">No transcript-linked InterPro feature with protein coordinates</text>"##,y+31.0));
        }

        let mut means=[0.0;2];
        for g in 0..2 {
            let vals:Vec<f64>=rows.iter().filter(|r|r.tx_id==tx_id&&assigned.get(&r.sample)==Some(&g)).map(|r|r.fraction).collect();
            if !vals.is_empty() { means[g]=vals.iter().sum::<f64>()/vals.len() as f64; }
        }
        let rx=w-right+25.0;
        s.push_str(&format!(r##"<text x="{rx}" y="{:.1}" class="small"><tspan font-weight="bold">{}</tspan> {:.1}%</text><text x="{rx}" y="{:.1}" class="small"><tspan font-weight="bold">{}</tspan> {:.1}%</text>"##,y-3.0,xml_escape(&groups[0].label),means[0]*100.0,y+14.0,xml_escape(&groups[1].label),means[1]*100.0));
    }
    s.push_str(&format!(r##"<text x="{left}" y="{:.1}" class="small">Grey exon = transcribed exon · dark inset = CDS · blue blocks = InterPro protein feature projected onto the coding genomic bases · dashed connector = one protein feature spanning splice junction(s)</text>"##,h-28.0));
    s.push_str("</svg>");
    std::fs::write(path,s)?;
    Ok(())
}

fn write_report(path:&Path,args:&Cli,idx:&SpliceIndex,rows:&[SampleEstimate],assigned:&HashMap<String,usize>,groups:&[GroupSpec;2],group_samples:&[Vec<String>;2],plotted:&[usize],table:&Path,boot:&Path,shift:&Path,structure:Option<&Path>)->Result<()> {
    let mut w=BufWriter::new(File::create(path)?);
    writeln!(w,"# {} transcript-usage comparison\n",args.gene)?;
    writeln!(w,"## Question\n\nDoes relative transcript usage within **{}** differ between **{}** and **{}**? This analysis compares exactly two sample groups by design.\n",args.gene,groups[0].label,groups[1].label)?;
    writeln!(w,"## What was done\n\nPaired BAM fragments overlapping the gene were matched against immutable transcript models from `gtf_splice_index`. A fragment can remain compatible with more than one transcript; `ExactJunctionChain` and splice-compatible matches are retained at equal selection rank so a shorter annotated transcript does not exclude a longer transcript that supports the same observed splice. Transcript fractions were estimated with an EM model over transcript-compatibility equivalence classes. Per-sample sampling uncertainty was estimated by {} fragment-bootstrap replicates with replacement, rerunning the EM for each replicate. Practical transcript identifiability was tested independently with {} randomized EM starting mixtures while keeping the observed fragments fixed. If many distinct starts converge to the same solution, the transcript mixture is well constrained by the compatibility-class geometry; if they retain substantially different solutions, individual transcript fractions are not identifiable even when the fragment-bootstrap interval is narrow.\n",args.bootstrap_replicates,args.identifiability_starts)?;
    writeln!(w,"The figures show only transcripts whose mean inferred fraction is at least {:.1}% in either group. The complete numeric results remain in `{}`.\n",args.plot_min_fraction*100.0,args.out.display())?;
    writeln!(w,"## Outputs\n\n- `{}` — group-level table supporting the figures.\n- `{}` — biological-sample dot plot with within-sample 95% bootstrap intervals.\n- `{}` — compact transcript-usage shift plot; blue points/lines are group means and coloured points are biological samples.",table.display(),boot.display(),shift.display())?;
    if let Some(structure) = structure {
        writeln!(w,"- `{}` — genomic transcript structures with transcript-specific InterPro protein features projected through each CDS back onto the coding exons. Protein feature blocks therefore mirror exon geometry and never fill introns.",structure.display())?;
    }
    writeln!(w)?;
    writeln!(w,"## Samples\n")?; for g in 0..2 { writeln!(w,"- **{}:** {}",groups[g].label,group_samples[g].iter().map(|x|sample_name(x)).collect::<Vec<_>>().join(", "))?; }
    writeln!(w,"\n## Interpretation guide\n\nThe coloured bootstrap interval answers how stable a transcript estimate is **within one BAM** under resampling of observed fragments. The grey randomized-EM interval asks a different question: whether the same fixed fragment evidence admits materially different transcript mixtures from different valid starting points. Separation between biological samples is a different source of variation. A narrow bootstrap interval therefore does not prove that a transcript is identifiable. A wide grey interval is a warning that the individual transcript abundance is underdetermined by the observed compatibility classes. Conversely, a narrow grey interval indicates practical convergence of the EM solution, but still does not prove that the annotation/model itself is complete or correct. Transcripts with zero unique fragments can still be constrained if overlapping compatibility classes jointly identify them. The comparison table therefore retains unique-fragment and splice-informative counts alongside the inferred fractions.\n")?;
    writeln!(w,"## Important limitation\n\nThis is a transcript-usage model, not a direct transcript-specific molecule count. Ambiguous fragments contribute probabilistically through the EM model. The SVGs are intended to make the inferred shift and its uncertainty explicit; raw splice inspection (for example a sashimi plot) remains useful as a validation/debugging view but is not the primary quantitative result.\n")?;
    writeln!(w,"## Plotted transcripts\n")?; for &tx in plotted { writeln!(w,"- {} ({})",idx.transcript_name(tx).unwrap_or("NA"),idx.transcripts[tx].names.join(";"))?; }
    Ok(())
}

fn collect_fragments(
    idx: &SpliceIndex,
    gene_tx: &HashSet<usize>,
    chr_name: &str,
    start: u32,
    end: u32,
    bam_path: &Path,
    min_mapq: u8,
) -> Result<Vec<FragmentEvidence>> {
    let mut reader = bam::IndexedReader::from_path(bam_path)
        .with_context(|| format!("opening indexed BAM {}", bam_path.display()))?;
    let header = reader.header().clone();
    let tid = find_tid(&header, chr_name)
        .with_context(|| format!("chromosome {chr_name} absent from {}", bam_path.display()))?;
    reader.fetch((tid, start as i64, end as i64))
        .with_context(|| format!("fetching {chr_name}:{start}-{end} from {}", bam_path.display()))?;

    let opts = MatchOptions::default();
    let mut pending: HashMap<Vec<u8>, PendingFragment> = HashMap::new();

    for rec in reader.records() {
        let rec = rec?;
        if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() || rec.mapq() < min_mapq { continue; }
        let blocks = record_to_blocks(&rec);
        if blocks.is_empty() { continue; }
        let mut spliced = SplicedRead::new(idx.chr_id(chr_name).context("splice chromosome disappeared")?, Strand::Unknown, blocks);
        spliced.finalize();
        let has_junction = !spliced.junctions().is_empty();
        let tx: HashMap<usize, EvidenceKind> = idx.match_transcripts(&spliced, opts)
            .into_iter()
            .filter(|m| gene_tx.contains(&m.transcript_id))
            .filter_map(|m| {
                let kind = match m.hit.class {
                    MatchClass::ExactJunctionChain => EvidenceKind::ExactJunctionChain,
                    MatchClass::Compatible if has_junction => EvidenceKind::SpliceCompatible,
                    MatchClass::Compatible => EvidenceKind::ExonCompatible,
                    _ => return None,
                };
                Some((m.transcript_id, kind))
            })
            .collect();
        if tx.is_empty() { continue; }
        pending.entry(rec.qname().to_vec()).or_default().mates.push(ReadEvidence { tx, blocks: spliced.blocks.clone() });
    }

    let mut fragments = Vec::with_capacity(pending.len());
    for (read_id, p) in pending {
        let mut it = p.mates.into_iter();
        let Some(first) = it.next() else { continue };
        let mut tx = first.tx;
        let mut mate_blocks = vec![first.blocks];
        for mate in it {
            tx.retain(|id, kind| {
                if let Some(mate_kind) = mate.tx.get(id) {
                    *kind = (*kind).max(*mate_kind);
                    true
                } else {
                    false
                }
            });
            mate_blocks.push(mate.blocks);
        }
        if tx.is_empty() { continue; }
        let mut tx: Vec<(usize, EvidenceKind)> = tx.into_iter().collect();
        tx.sort_unstable_by_key(|(id, _)| *id);
        fragments.push(FragmentEvidence { read_id, tx, mate_blocks });
    }
    Ok(fragments)
}

fn format_blocks(blocks: &[gtf_splice_index::RefBlock]) -> String {
    blocks.iter().map(|b| format!("{}-{}", b.start, b.end)).collect::<Vec<_>>().join(",")
}

fn format_junctions(junctions: &[(u32, u32)]) -> String {
    junctions.iter().map(|(a, b)| format!("{}-{}", a, b)).collect::<Vec<_>>().join(",")
}

fn opt_usize(v: Option<usize>) -> String {
    v.map(|x| x.to_string()).unwrap_or_else(|| "NA".to_string())
}

fn write_audit<W: Write>(
    out: &mut W,
    idx: &SpliceIndex,
    transcript_ids: &[usize],
    chr_id: usize,
    sample: &str,
    gene: &str,
    fragments: &[FragmentEvidence],
) -> Result<()> {
    let opts = MatchOptions::default();
    for fragment in fragments {
        let read_id = String::from_utf8_lossy(&fragment.read_id);
        let compatible: HashSet<usize> = fragment.tx.iter().map(|(t, _)| *t).collect();
        for &tx_id in transcript_ids {
            let tx = &idx.transcripts[tx_id];
            let tx_name = idx.transcript_name(tx_id).unwrap_or("NA");
            let aliases = tx.names.join(";");
            let fragment_compatible = compatible.contains(&tx_id);
            for (mate_i, blocks) in fragment.mate_blocks.iter().enumerate() {
                let mut read = SplicedRead::new(chr_id, Strand::Unknown, blocks.clone());
                read.finalize();
                let a = tx.placement_audit(&read, opts);
                writeln!(out,
                    "{sample}\t{gene}\t{read_id}\t{tx_name}\t{aliases}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    fragment_compatible,
                    mate_i + 1,
                    a.read_start,
                    a.read_end,
                    format_blocks(&read.blocks),
                    format_junctions(&a.read_junctions),
                    a.transcript_start,
                    a.transcript_end,
                    opt_usize(a.transcript_read_start),
                    opt_usize(a.transcript_read_end),
                    a.hit.class,
                    a.exonic_bases,
                    a.intronic_bases,
                    a.matched_junctions,
                    a.unmatched_junctions,
                )?;
            }
        }
    }
    Ok(())
}

fn find_tid(header: &bam::HeaderView, chr: &str) -> Option<u32> {
    let aliases = if let Some(rest) = chr.strip_prefix("chr") {
        vec![chr.to_string(), rest.to_string()]
    } else {
        vec![chr.to_string(), format!("chr{chr}")]
    };
    aliases.into_iter().find_map(|name| header.tid(name.as_bytes()))
}

#[derive(Clone, Debug)]
struct CompatibilityClass {
    tx: Vec<usize>,
    count: usize,
}

fn compatibility_classes(fragments: &[FragmentEvidence]) -> Vec<CompatibilityClass> {
    let mut grouped: HashMap<Vec<usize>, usize> = HashMap::new();
    for fragment in fragments {
        let mut tx: Vec<usize> = fragment.tx.iter().map(|(tx_id, _)| *tx_id).collect();
        tx.sort_unstable();
        tx.dedup();
        if !tx.is_empty() {
            *grouped.entry(tx).or_insert(0) += 1;
        }
    }
    grouped.into_iter()
        .map(|(tx, count)| CompatibilityClass { tx, count })
        .collect()
}

fn em(fragments: &[FragmentEvidence], gene_tx: &HashSet<usize>, max_iter: usize, epsilon: f64) -> HashMap<usize, f64> {
    let classes = compatibility_classes(fragments);
    em_weighted(&classes, gene_tx, max_iter, epsilon)
}

fn em_weighted(classes: &[CompatibilityClass], gene_tx: &HashSet<usize>, max_iter: usize, epsilon: f64) -> HashMap<usize, f64> {
    if classes.is_empty() { return HashMap::new(); }
    let init = 1.0 / gene_tx.len() as f64;
    let theta: HashMap<usize, f64> = gene_tx.iter().map(|&t| (t, init)).collect();
    em_weighted_from_theta(classes, gene_tx, max_iter, epsilon, theta)
}

fn em_weighted_from_theta(classes: &[CompatibilityClass], gene_tx: &HashSet<usize>, max_iter: usize, epsilon: f64, mut theta: HashMap<usize, f64>) -> HashMap<usize, f64> {
    if classes.is_empty() { return HashMap::new(); }
    let mut counts = HashMap::new();
    for _ in 0..max_iter {
        counts = gene_tx.iter().map(|&t| (t, 0.0)).collect();
        for class in classes {
            let denom: f64 = class.tx.iter().map(|t| theta.get(t).copied().unwrap_or(0.0)).sum();
            if denom <= 0.0 { continue; }
            let weight = class.count as f64;
            for t in &class.tx {
                *counts.entry(*t).or_insert(0.0) +=
                    weight * theta.get(t).copied().unwrap_or(0.0) / denom;
            }
        }
        let total: f64 = counts.values().sum();
        if total <= 0.0 { break; }
        let mut delta: f64 = 0.0;
        for &t in gene_tx {
            let next = counts.get(&t).copied().unwrap_or(0.0) / total;
            delta = delta.max((next - theta.get(&t).copied().unwrap_or(0.0)).abs());
            theta.insert(t, next);
        }
        if delta < epsilon { break; }
    }
    counts
}


#[derive(Clone, Copy, Debug, Default)]
struct IdentifiabilitySummary {
    mean: f64,
    sd: f64,
    p025: f64,
    p975: f64,
    range: f64,
}

fn multistart_identifiability(
    fragments: &[FragmentEvidence],
    gene_tx: &HashSet<usize>,
    max_iter: usize,
    epsilon: f64,
    starts: usize,
    seed: u64,
) -> HashMap<usize, IdentifiabilitySummary> {
    if fragments.is_empty() || starts == 0 { return HashMap::new(); }
    let classes = compatibility_classes(fragments);
    let n = fragments.len() as f64;
    let tx_ids: Vec<usize> = gene_tx.iter().copied().collect();
    let runs: Vec<HashMap<usize, f64>> = (0..starts).into_par_iter().map(|start| {
        let mut rng = BootstrapRng::new(seed ^ (start as u64).wrapping_mul(0xD1B54A32D192ED03));
        let mut init = HashMap::with_capacity(tx_ids.len());
        let mut sum = 0.0;
        for &tx in &tx_ids {
            // Positive, strongly variable starts.  -ln(U) gives an exponential
            // draw and therefore a Dirichlet(1) mixture after normalization.
            let u = ((rng.next_u64() >> 11) as f64 + 1.0) / ((1u64 << 53) as f64 + 1.0);
            let v = -u.ln();
            init.insert(tx, v); sum += v;
        }
        for v in init.values_mut() { *v /= sum; }
        let counts = em_weighted_from_theta(&classes, gene_tx, max_iter, epsilon, init);
        tx_ids.iter().map(|&tx| (tx, counts.get(&tx).copied().unwrap_or(0.0) / n)).collect()
    }).collect();
    tx_ids.into_iter().map(|tx| {
        let mut values: Vec<f64> = runs.iter().map(|r| r.get(&tx).copied().unwrap_or(0.0)).collect();
        values.sort_by(|a,b| a.total_cmp(b));
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let sd = if values.len()>1 {
            (values.iter().map(|x|(x-mean).powi(2)).sum::<f64>()/(values.len()-1) as f64).sqrt()
        } else { 0.0 };
        let min = values.first().copied().unwrap_or(0.0);
        let max = values.last().copied().unwrap_or(0.0);
        (tx, IdentifiabilitySummary { mean, sd, p025: percentile(&values,0.025), p975: percentile(&values,0.975), range: max-min })
    }).collect()
}

#[derive(Clone, Copy, Debug, Default)]
struct BootstrapSummary {
    mean: f64,
    sd: f64,
    ci025: f64,
    ci975: f64,
    nonzero_fraction: f64,
}

fn bootstrap_fractions(
    fragments: &[FragmentEvidence],
    gene_tx: &HashSet<usize>,
    max_iter: usize,
    epsilon: f64,
    replicates: usize,
    seed: u64,
) -> HashMap<usize, BootstrapSummary> {
    if fragments.is_empty() || replicates == 0 {
        return HashMap::new();
    }

    let n = fragments.len();
    // Collapse identical transcript-compatibility sets once.  The EM only
    // depends on these equivalence classes and their multiplicities, not on
    // fragment identity.  Bootstrap resampling therefore only changes class
    // counts; each EM iteration scales with the number of distinct classes.
    let base_classes = compatibility_classes(fragments);
    let class_lookup: Vec<usize> = base_classes.iter().enumerate()
        .flat_map(|(class_id, class)| std::iter::repeat(class_id).take(class.count))
        .collect();

    let runs: Vec<HashMap<usize, f64>> = (0..replicates)
        .into_par_iter()
        .map(|replicate| {
            let mut rng = BootstrapRng::new(seed ^ (replicate as u64).wrapping_mul(0x9E3779B97F4A7C15));
            let mut sampled_counts = vec![0usize; base_classes.len()];
            for _ in 0..n {
                sampled_counts[class_lookup[rng.index(n)]] += 1;
            }
            let sampled_classes: Vec<CompatibilityClass> = base_classes.iter()
                .zip(sampled_counts)
                .filter_map(|(class, count)| (count > 0).then(|| CompatibilityClass {
                    tx: class.tx.clone(),
                    count,
                }))
                .collect();
            let counts = em_weighted(&sampled_classes, gene_tx, max_iter, epsilon);
            gene_tx.iter().map(|&tx_id| {
                let fraction = counts.get(&tx_id).copied().unwrap_or(0.0) / n as f64;
                (tx_id, fraction)
            }).collect()
        })
        .collect();

    gene_tx.iter().map(|&tx_id| {
        let mut values: Vec<f64> = runs.iter()
            .map(|run| run.get(&tx_id).copied().unwrap_or(0.0))
            .collect();
        values.sort_by(|a, b| a.total_cmp(b));
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let sd = if values.len() > 1 {
            let ss = values.iter().map(|x| (x - mean).powi(2)).sum::<f64>();
            (ss / (values.len() - 1) as f64).sqrt()
        } else { 0.0 };
        let nonzero = values.iter().filter(|&&x| x > 0.0).count() as f64 / values.len() as f64;
        let summary = BootstrapSummary {
            mean,
            sd,
            ci025: percentile(&values, 0.025),
            ci975: percentile(&values, 0.975),
            nonzero_fraction: nonzero,
        };
        (tx_id, summary)
    }).collect()
}

fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() { return 0.0; }
    if sorted.len() == 1 { return sorted[0]; }
    let pos = p.clamp(0.0, 1.0) * (sorted.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    if lo == hi { sorted[lo] } else {
        let w = pos - lo as f64;
        sorted[lo] * (1.0 - w) + sorted[hi] * w
    }
}

struct BootstrapRng(u64);

impl BootstrapRng {
    fn new(seed: u64) -> Self {
        Self(if seed == 0 { 0xA0761D6478BD642F } else { seed })
    }

    fn next_u64(&mut self) -> u64 {
        // SplitMix64: tiny, deterministic and sufficient for bootstrap resampling.
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }

    fn index(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}
