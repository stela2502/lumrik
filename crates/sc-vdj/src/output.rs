use crate::posterior::CellVdjSummary;
use crate::types::Chain;
use anyhow::Result;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn write_reports<P: AsRef<Path>>(dir: P, cells: &[CellVdjSummary]) -> Result<()> {
    let dir = dir.as_ref();
    fs::create_dir_all(dir)?;
    let mut summary = BufWriter::new(File::create(dir.join("vdj_cell_summary.tsv"))?);
    writeln!(
        summary,
        "cell\tprogram\tdevelopment_score\tevidence_code\trag_activity\trearrangements"
    )?;
    for c in cells {
        let calls = c
            .rearrangements
            .iter()
            .map(|x| x.notation.as_str())
            .collect::<Vec<_>>()
            .join(";");
        writeln!(
            summary,
            "{}\t{:?}\t{:.6}\t{}\t{:.6}\t{}",
            c.cell,
            c.development.program,
            c.development.probability_like_score,
            c.development.evidence_code,
            c.development.rag_activity,
            calls
        )?;
    }

    let mut calls = BufWriter::new(File::create(dir.join("vdj_rearrangements.tsv"))?);
    writeln!(
        calls,
        "cell\tchain\tstage\tnotation\tv\td\tj\tc\tumis\tv_locus_fraction\tv_distance_to_center"
    )?;
    for cell in cells {
        for r in &cell.rearrangements {
            writeln!(
                calls,
                "{}\t{}\t{:?}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                cell.cell,
                r.chain,
                r.stage,
                r.notation,
                id(&r.v),
                id(&r.d),
                id(&r.j),
                id(&r.c),
                r.total_supporting_umis,
                num(&r.v, |x| x.locus_fraction),
                num(&r.v, |x| x.distance_to_recombination_center as f64)
            )?;
        }
    }

    let mut sterile = BufWriter::new(File::create(dir.join("vdj_sterile_spatial.tsv"))?);
    writeln!(sterile,"cell\tchain\tbin\tstart_fraction\tend_fraction\tunique_umis\treads\tbreadth\tcentroid\tproximal_fraction\tdistal_fraction")?;
    for c in cells {
        for p in &c.sterile {
            for b in &p.bins {
                writeln!(
                    sterile,
                    "{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                    c.cell,
                    p.chain,
                    b.bin,
                    b.start_fraction,
                    b.end_fraction,
                    b.unique_umis,
                    b.reads,
                    p.breadth,
                    p.centroid,
                    p.proximal_fraction,
                    p.distal_fraction
                )?;
            }
        }
    }

    let mut intervals = BufWriter::new(File::create(dir.join("vdj_sterile_intervals.tsv"))?);
    writeln!(intervals, "cell\tchain\tstart\tend\tunique_umis\treads")?;
    for c in cells {
        for p in &c.sterile {
            for i in &p.supported_intervals {
                writeln!(
                    intervals,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    c.cell, p.chain, i.start, i.end, i.unique_umis, i.reads
                )?;
            }
        }
    }

    let mut rationale = BufWriter::new(File::create(dir.join("vdj_development_rationale.tsv"))?);
    writeln!(
        rationale,
        "cell\tprogram\tscore\tevidence_code\tgene\texpression\tweight\tcontribution"
    )?;
    for c in cells {
        for x in &c.development.contributions {
            writeln!(
                rationale,
                "{}\t{:?}\t{:.6}\t{}\t{}\t{:.6}\t{:.3}\t{:.6}",
                c.cell,
                c.development.program,
                c.development.probability_like_score,
                c.development.evidence_code,
                x.gene,
                x.expression,
                x.weight,
                x.contribution
            )?;
        }
    }

    // Always emit all seven loci, including the common all-zero case.
    let mut sample = BufWriter::new(File::create(dir.join("vdj_sample_summary.tsv"))?);
    writeln!(sample,"chain\tcells\tcells_with_sterile\tfraction_with_sterile\tmean_sterile_umis\tmean_breadth\tmean_centroid\tmean_distal_fraction\tmax_distal_fraction")?;
    for chain in Chain::ALL {
        let mut present = 0usize;
        let mut umis = 0usize;
        let mut breadth = 0.0;
        let mut centroid = 0.0;
        let mut distal = 0.0;
        let mut max_distal = 0.0;
        for c in cells {
            if let Some(p) = c.sterile.iter().find(|p| p.chain == chain) {
                if p.total_unique_umis > 0 {
                    present += 1
                }
                umis += p.total_unique_umis;
                breadth += p.breadth;
                centroid += p.centroid;
                distal += p.distal_fraction;
                if p.distal_fraction > max_distal {
                    max_distal = p.distal_fraction;
                }
            }
        }
        let n = cells.len().max(1) as f64;
        writeln!(
            sample,
            "{}\t{}\t{}\t{:.8}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            chain,
            cells.len(),
            present,
            present as f64 / n,
            umis as f64 / n,
            breadth / n,
            centroid / n,
            distal / n,
            max_distal
        )?;
    }
    Ok(())
}

fn id(x: &Option<crate::posterior::GermlineSegmentSupport>) -> &str {
    x.as_ref().map(|x| x.id.as_str()).unwrap_or("")
}
fn num<F: FnOnce(&crate::posterior::GermlineSegmentSupport) -> f64>(
    x: &Option<crate::posterior::GermlineSegmentSupport>,
    f: F,
) -> String {
    x.as_ref()
        .map(|x| format!("{:.6}", f(x)))
        .unwrap_or_default()
}
