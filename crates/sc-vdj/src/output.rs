use crate::identity::{PackedRecombinationId, RecombinationMeasurements};
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
            "cell\trag1_expression\trag2_expression\tdntt_expression\trag_activity\trag_pair_detected\tlight_chain_status\theavy_recombination_id\tlight_recombination_id\trearrangements"
        )?;
        let mut calls = BufWriter::new(File::create(dir.join("vdj_rearrangements.tsv"))?);
        writeln!(
            calls,
            "cell\trecombination_id\tchain\tstage\tnotation\tv\td\tj\tc\tumis\tv_locus_fraction\tv_distance_to_center\tv_del_3\tp_v3_len\tp_v3\tn1_len\tn1\tp_d5_len\tp_d5\td_del_5\td_retained_len\td_del_3\td_inferred_from_vj_junction\td_hypothesis_margin\tp_d3_len\tp_d3\tn2_len\tn2\tp_j5_len\tp_j5\tj_del_5\tpn_alternative\tobserved_v\tobserved_d\tobserved_j\tnaive_v\tnaive_d\tnaive_j\tobserved_rearrangement\tnaive_recombination"
        )?;
        let mut sterile = BufWriter::new(File::create(dir.join("vdj_sterile_spatial.tsv"))?);
        writeln!(sterile,"cell\tchain\tbin\tstart_fraction\tend_fraction\tunique_umis\treads\tbreadth\tcentroid\tproximal_fraction\tdistal_fraction")?;
        let mut intervals = BufWriter::new(File::create(dir.join("vdj_sterile_intervals.tsv"))?);
        writeln!(intervals, "cell\tchain\tstart\tend\tunique_umis\treads")?;
        let mut rationale = BufWriter::new(File::create(dir.join("vdj_recombination_activity.tsv"))?);
        writeln!(rationale, "cell\tgene\texpression\tactivity")?;
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
                "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}",
                c.cell,
                c.recombination_activity.rag1_expression,
                c.recombination_activity.rag2_expression,
                c.recombination_activity.dntt_expression,
                c.recombination_activity.rag_activity,
                c.recombination_activity.rag_pair_detected,
                c.light_chain_status(),
                c.strongest_heavy_id().map(|x| x.to_string()).unwrap_or_default(),
                c.strongest_light_id().map(|x| x.to_string()).unwrap_or_default(),
                calls,
            )?;

            for r in &c.rearrangements {
                let measurements = RecombinationMeasurements::from_call(r);
                let m = measurements.as_ref();
                let junction = r.junction.as_ref();
                let seq = |value: Option<&[u8]>| {
                    value
                        .map(|x| String::from_utf8_lossy(x).into_owned())
                        .unwrap_or_default()
                };
                let fields = vec![
                    c.cell.clone(),
                    PackedRecombinationId::from_call(r).map(|x| x.to_string()).unwrap_or_default(),
                    r.chain.to_string(),
                    format!("{:?}", r.stage),
                    r.notation.clone(),
                    id(&r.v).to_string(),
                    id(&r.d).to_string(),
                    id(&r.j).to_string(),
                    id(&r.c).to_string(),
                    r.total_supporting_umis.to_string(),
                    num(&r.v, |x| x.locus_fraction),
                    num(&r.v, |x| x.distance_to_recombination_center as f64),
                    opt_u16(m.and_then(|x| x.v_del_3)),
                    opt_u16(m.and_then(|x| x.p_v3_len)),
                    seq(junction.map(|x| x.p_v3.as_slice())),
                    opt_u16(m.and_then(|x| x.n1_len)),
                    seq(junction.map(|x| x.n1.as_slice())),
                    opt_u16(m.and_then(|x| x.p_d5_len)),
                    seq(junction.map(|x| x.p_d5.as_slice())),
                    opt_u16(m.and_then(|x| x.d_del_5)),
                    opt_u16(m.and_then(|x| x.d_retained_len)),
                    opt_u16(m.and_then(|x| x.d_del_3)),
                    r.d_inferred_from_vj_junction.to_string(),
                    r.d_hypothesis_margin.map(|x| x.to_string()).unwrap_or_default(),
                    opt_u16(m.and_then(|x| x.p_d3_len)),
                    seq(junction.map(|x| x.p_d3.as_slice())),
                    opt_u16(m.and_then(|x| x.n2_len)),
                    seq(junction.map(|x| x.n2.as_slice())),
                    opt_u16(m.and_then(|x| x.p_j5_len)),
                    seq(junction.map(|x| x.p_j5.as_slice())),
                    opt_u16(m.and_then(|x| x.j_del_5)),
                    m.map(|x| x.pn_alternative.to_string()).unwrap_or_default(),
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

            for x in &c.recombination_activity.contributions {
                writeln!(
                    self.rationale,
                    "{}\t{}\t{:.6}\t{:.6}",
                    c.cell, x.gene, x.expression, x.activity
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

fn opt_u16(value: Option<u16>) -> String {
    value.map(|x| x.to_string()).unwrap_or_default()
}
