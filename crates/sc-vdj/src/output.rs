use crate::posterior::CellVdjSummary;
use crate::types::Chain;
use anyhow::Result;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Default, Clone)]
struct SampleAccumulator {
    present: usize,
    umis: usize,
    breadth: f64,
    centroid: f64,
    distal: f64,
    max_distal: f64,
}

pub struct ReportWriter {
    summary: BufWriter<File>,
    calls: BufWriter<File>,
    sterile: BufWriter<File>,
    intervals: BufWriter<File>,
    rationale: BufWriter<File>,
    sample_path: std::path::PathBuf,
    sample: HashMap<Chain, SampleAccumulator>,
    total_cells: usize,
}

impl ReportWriter {
    pub fn create<P: AsRef<Path>>(dir: P) -> Result<Self> {
        let dir = dir.as_ref();
        fs::create_dir_all(dir)?;
        let mut summary = BufWriter::new(File::create(dir.join("vdj_cell_summary.tsv"))?);
        writeln!(
            summary,
            "cell\tprogram\tdevelopment_score\tevidence_code\trag_activity\trearrangements"
        )?;
        let mut calls = BufWriter::new(File::create(dir.join("vdj_rearrangements.tsv"))?);
        writeln!(
            calls,
            "cell\tchain\tstage\tnotation\tv\td\tj\tc\tumis\tv_locus_fraction\tv_distance_to_center"
        )?;
        let mut sterile = BufWriter::new(File::create(dir.join("vdj_sterile_spatial.tsv"))?);
        writeln!(sterile,"cell\tchain\tbin\tstart_fraction\tend_fraction\tunique_umis\treads\tbreadth\tcentroid\tproximal_fraction\tdistal_fraction")?;
        let mut intervals = BufWriter::new(File::create(dir.join("vdj_sterile_intervals.tsv"))?);
        writeln!(intervals, "cell\tchain\tstart\tend\tunique_umis\treads")?;
        let mut rationale = BufWriter::new(File::create(dir.join("vdj_development_rationale.tsv"))?);
        writeln!(
            rationale,
            "cell\tprogram\tscore\tevidence_code\tgene\texpression\tweight\tcontribution"
        )?;
        Ok(Self {
            summary,
            calls,
            sterile,
            intervals,
            rationale,
            sample_path: dir.join("vdj_sample_summary.tsv"),
            sample: HashMap::new(),
            total_cells: 0,
        })
    }

    pub fn write_cells(&mut self, cells: &[CellVdjSummary]) -> Result<()> {
        for c in cells {
            self.total_cells += 1;
            let calls = c
                .rearrangements
                .iter()
                .map(|x| x.notation.as_str())
                .collect::<Vec<_>>()
                .join(";");
            writeln!(
                self.summary,
                "{}\t{:?}\t{:.6}\t{}\t{:.6}\t{}",
                c.cell,
                c.development.program,
                c.development.probability_like_score,
                c.development.evidence_code,
                c.development.rag_activity,
                calls
            )?;

            for r in &c.rearrangements {
                writeln!(
                    self.calls,
                    "{}\t{}\t{:?}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    c.cell,
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

            for p in &c.sterile {
                for b in &p.bins {
                    writeln!(
                        self.sterile,
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
                for i in &p.supported_intervals {
                    writeln!(
                        self.intervals,
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        c.cell, p.chain, i.start, i.end, i.unique_umis, i.reads
                    )?;
                }

                let acc = self.sample.entry(p.chain).or_default();
                if p.total_unique_umis > 0 {
                    acc.present += 1;
                }
                acc.umis += p.total_unique_umis;
                acc.breadth += p.breadth;
                acc.centroid += p.centroid;
                acc.distal += p.distal_fraction;
                acc.max_distal = acc.max_distal.max(p.distal_fraction);
            }

            for x in &c.development.contributions {
                writeln!(
                    self.rationale,
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
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        self.summary.flush()?;
        self.calls.flush()?;
        self.sterile.flush()?;
        self.intervals.flush()?;
        self.rationale.flush()?;

        let mut sample = BufWriter::new(File::create(&self.sample_path)?);
        writeln!(sample,"chain\tcells\tcells_with_sterile\tfraction_with_sterile\tmean_sterile_umis\tmean_breadth\tmean_centroid\tmean_distal_fraction\tmax_distal_fraction")?;
        for chain in Chain::ALL {
            let acc = self.sample.get(&chain).cloned().unwrap_or_default();
            let n = self.total_cells.max(1) as f64;
            writeln!(
                sample,
                "{}\t{}\t{}\t{:.8}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                chain,
                self.total_cells,
                acc.present,
                acc.present as f64 / n,
                acc.umis as f64 / n,
                acc.breadth / n,
                acc.centroid / n,
                acc.distal / n,
                acc.max_distal
            )?;
        }
        sample.flush()?;
        Ok(())
    }
}

pub fn write_reports<P: AsRef<Path>>(dir: P, cells: &[CellVdjSummary]) -> Result<()> {
    let mut writer = ReportWriter::create(dir)?;
    writer.write_cells(cells)?;
    writer.finish()
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
