use crate::{Recombination, VdjIndex};
use anyhow::{Context, Result};
use mapping_info::MappingInfo;
use sc_primer::{BdCellVersion, RhapsodyWhitelist};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub struct ReportWriter {
    calls: BufWriter<File>,
    airr: BufWriter<File>,
    receptors: BufWriter<File>,
    observed: Option<BufWriter<File>>,
    naive: Option<BufWriter<File>>,
    rhapsody: Option<RhapsodyWhitelist>,
}

impl ReportWriter {
    pub fn create<P: AsRef<Path>>(dir: P, write_sequences: bool) -> Result<Self> {
        let dir = dir.as_ref();
        fs::create_dir_all(dir)?;

        let mut calls = writer(dir.join("vdj_calls.tsv"))?;
        writeln!(calls, "cell\trustody_cell_id\trecombination_id\tchain\tstage\tv\td\tj\tc\tproductivity_status\tsupport_features\treceptor_rediscovery_reads\tjunction_support_reads\tjunction_spanning_reads\tjunction_conflicting_reads\tjunction_refined_bases\tconstant_link_fragments\tconstant_spanning_reads\tconstant_link_call\tv_del_3\tp_v3_len\tp_v3\tn1_len\tn1\tp_d5_len\tp_d5\td_del_5\td_retained_len\td_retained\td_del_3\tp_d3_len\tp_d3\tn2_len\tn2\tp_j5_len\tp_j5\tj_del_5\tpn_alternative\tobserved_rearrangement\tnaive_recombination\tobserved_receptor_sequence")?;

        let mut airr = writer(dir.join("airr_rearrangements.tsv"))?;
        writeln!(airr, "sequence_id\tsequence\tproductive\tvj_in_frame\tstop_codon\tcomplete_vdj\tlocus\tv_call\td_call\tj_call\tc_call\tjunction\tjunction_aa\tcdr3\tcdr3_aa\tnp1\tnp2\tnp1_length\tnp2_length\tcell_id\tlumrik_rustody_cell_id\tlumrik_productivity_status\tlumrik_recombination_id\tlumrik_supporting_features\tlumrik_receptor_rediscovery_reads\tlumrik_junction_support_reads\tlumrik_junction_spanning_reads\tlumrik_junction_conflicting_reads\tlumrik_junction_refined_bases\tlumrik_constant_link_fragments\tlumrik_constant_spanning_reads\tlumrik_constant_link_call")?;

        let mut receptors = writer(dir.join("vdj_receptors.tsv"))?;
        writeln!(receptors, "cell\theavy_recombination_id\tlight_recombination_id\theavy_chain\theavy_v\theavy_d\theavy_j\theavy_c\theavy_support_features\theavy_receptor_rediscovery_reads\theavy_constant_link_fragments\theavy_constant_spanning_reads\theavy_naive_recombination\tlight_chain\tlight_v\tlight_j\tlight_c\tlight_support_features\tlight_receptor_rediscovery_reads\tlight_constant_link_fragments\tlight_constant_spanning_reads\tlight_naive_recombination")?;

        let observed = write_sequences
            .then(|| writer(dir.join("vdj_observed.fasta")))
            .transpose()?;
        let naive = write_sequences
            .then(|| writer(dir.join("vdj_naive.fasta")))
            .transpose()?;
        Ok(Self {
            calls,
            airr,
            receptors,
            observed,
            naive,
            rhapsody: None,
        })
    }

    pub fn with_bd_cell_version(mut self, version: BdCellVersion) -> Self {
        self.rhapsody = Some(RhapsodyWhitelist::builtin(version));
        self
    }

    pub fn write_cell(
        &mut self,
        cell: &str,
        recombinations: &[Recombination],
        index: &VdjIndex,
    ) -> Result<()> {
        let heavy = strongest(recombinations, true);
        let light = strongest(recombinations, false);
        writeln!(
            self.receptors,
            "{}",
            receptor_row(cell, heavy, light, index).join("\t")
        )?;

        for r in recombinations {
            let v = seg(index, Some(r.v));
            let d = seg(index, r.d);
            let j = seg(index, Some(r.j));
            let c = r
                .constant
                .as_ref()
                .map(|x| seg(index, Some(x.segment)))
                .unwrap_or_default();
            let stage = if r.chain.has_d() { "Vdj" } else { "Vj" };
            let x = &r.junction;
            let fields = vec![
                cell.to_string(),
                self.rhapsody
                    .as_ref()
                    .and_then(|wl| wl.cell_id_for_seq(cell.as_bytes()))
                    .map(|x| x.to_string())
                    .unwrap_or_default(),
                r.stable_id.to_string(),
                r.chain.to_string(),
                stage.into(),
                v,
                d,
                j,
                c,
                r.productivity_status.as_str().to_string(),
                r.supporting_features.to_string(),
                r.receptor_linkage.rediscovery_reads.to_string(),
                r.receptor_linkage.junction_support_reads.to_string(),
                r.receptor_linkage.junction_spanning_reads.to_string(),
                r.receptor_linkage.junction_conflicting_reads.to_string(),
                r.receptor_linkage.junction_refined_bases.to_string(),
                r.receptor_linkage.constant_link_fragments.to_string(),
                r.receptor_linkage.constant_spanning_reads.to_string(),
                r.receptor_linkage
                    .constant_segment
                    .map(|segment| seg(index, Some(segment)))
                    .unwrap_or_default(),
                x.v_del_3.to_string(),
                x.p_v3_len().to_string(),
                dna(&x.p_v3),
                x.n1_len().to_string(),
                dna(&x.n1),
                x.p_d5_len().to_string(),
                dna(&x.p_d5),
                opt(x.d_del_5),
                x.d_retained_len().to_string(),
                dna(&x.d_retained),
                opt(x.d_del_3),
                x.p_d3_len().to_string(),
                dna(&x.p_d3),
                x.n2_len().to_string(),
                dna(&x.n2),
                x.p_j5_len().to_string(),
                dna(&x.p_j5),
                x.j_del_5.to_string(),
                x.pn_alternative.to_string(),
                dna(&r.observed_rearrangement),
                dna(&r.naive_recombination),
                dna(&r.observed_receptor_sequence),
            ];
            writeln!(self.calls, "{}", fields.join("\t"))?;

            let np1 = if r.chain.has_d() {
                [x.p_v3.as_slice(), x.n1.as_slice(), x.p_d5.as_slice()].concat()
            } else {
                [x.p_v3.as_slice(), x.n1.as_slice(), x.p_j5.as_slice()].concat()
            };
            let np2 = if r.chain.has_d() {
                [x.p_d3.as_slice(), x.n2.as_slice(), x.p_j5.as_slice()].concat()
            } else {
                Vec::new()
            };
            let airr = vec![
                format!("{}|{}", cell, r.stable_id),
                dna(&r.observed_receptor_sequence),
                airr_tf(r.productive, r.productivity_status.is_unknown()),
                airr_tf(r.in_frame, r.productivity_status.is_unknown()),
                airr_tf(r.stop_codon, r.productivity_status.is_unknown()),
                "T".into(),
                r.chain.to_string(),
                seg(index, Some(r.v)),
                seg(index, r.d),
                seg(index, Some(r.j)),
                r.constant
                    .as_ref()
                    .map(|q| seg(index, Some(q.segment)))
                    .unwrap_or_default(),
                dna(&r.airr_junction),
                String::from_utf8_lossy(&r.airr_junction_aa).into_owned(),
                dna(&r.cdr3),
                String::from_utf8_lossy(&r.cdr3_aa).into_owned(),
                dna(&np1),
                dna(&np2),
                np1.len().to_string(),
                np2.len().to_string(),
                cell.to_string(),
                self.rhapsody
                    .as_ref()
                    .and_then(|wl| wl.cell_id_for_seq(cell.as_bytes()))
                    .map(|x| x.to_string())
                    .unwrap_or_default(),
                r.productivity_status.as_str().to_string(),
                r.stable_id.to_string(),
                r.supporting_features.to_string(),
                r.receptor_linkage.rediscovery_reads.to_string(),
                r.receptor_linkage.junction_support_reads.to_string(),
                r.receptor_linkage.junction_spanning_reads.to_string(),
                r.receptor_linkage.junction_conflicting_reads.to_string(),
                r.receptor_linkage.junction_refined_bases.to_string(),
                r.receptor_linkage.constant_link_fragments.to_string(),
                r.receptor_linkage.constant_spanning_reads.to_string(),
                r.receptor_linkage
                    .constant_segment
                    .map(|segment| seg(index, Some(segment)))
                    .unwrap_or_default(),
            ];
            writeln!(self.airr, "{}", airr.join("\t"))?;

            let header = format!(
                "{}|{}|{}",
                fasta_token(cell),
                r.chain,
                fasta_token(&r.stable_id.to_string())
            );
            if let Some(w) = &mut self.observed {
                writeln!(w, ">{header}")?;
                writeln!(w, "{}", dna(&r.observed_receptor_sequence))?;
            }
            if let Some(w) = &mut self.naive {
                writeln!(w, ">{header}")?;
                writeln!(w, "{}", dna(&r.naive_recombination))?;
            }
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        self.calls.flush()?;
        self.airr.flush()?;
        self.receptors.flush()?;
        if let Some(w) = &mut self.observed {
            w.flush()?;
        }
        if let Some(w) = &mut self.naive {
            w.flush()?;
        }
        Ok(())
    }
}

pub fn write_mapping_info<P: AsRef<Path>>(
    path: P,
    receptor_records: usize,
    cells: usize,
    recombinations: usize,
    by_chain: &HashMap<String, usize>,
) -> Result<()> {
    let mut info = MappingInfo::new(None, 0.0, 0);
    info.total = receptor_records;
    info.report_n("vdj.receptor_overlap_records", receptor_records);
    info.report_n("vdj.cells_with_evidence", cells);
    info.report_n("vdj.recombinations", recombinations);
    for (chain, count) in by_chain {
        info.report_n(format!("vdj.calls.{}", chain.to_ascii_lowercase()), *count);
    }
    write_mapping_info_report(path, &info)
}

pub fn write_mapping_info_report<P: AsRef<Path>>(path: P, info: &MappingInfo) -> Result<()> {
    let mut w = writer(path.as_ref().to_path_buf())?;
    write!(w, "{info}")?;
    w.flush()?;
    Ok(())
}

fn strongest<'a>(rs: &'a [Recombination], heavy: bool) -> Option<&'a Recombination> {
    rs.iter()
        .filter(|r| r.chain.has_d() == heavy)
        .max_by_key(|r| {
            (
                r.supporting_features,
                std::cmp::Reverse(r.v),
                std::cmp::Reverse(r.j),
            )
        })
}
fn receptor_row(
    cell: &str,
    heavy: Option<&Recombination>,
    light: Option<&Recombination>,
    index: &VdjIndex,
) -> Vec<String> {
    let heavy_id = heavy.map(|r| r.stable_id.to_string()).unwrap_or_default();
    let light_id = light.map(|r| r.stable_id.to_string()).unwrap_or_default();
    let mut v = vec![cell.to_string(), heavy_id, light_id];
    v.extend(role_details(heavy, index, true));
    v.extend(role_details(light, index, false));
    v
}
fn role_details(r: Option<&Recombination>, index: &VdjIndex, heavy: bool) -> Vec<String> {
    match r {
        None => vec![String::new(); if heavy { 10 } else { 9 }],
        Some(r) => {
            let c = r
                .constant
                .as_ref()
                .map(|x| seg(index, Some(x.segment)))
                .unwrap_or_default();
            if heavy {
                vec![
                    r.chain.to_string(),
                    seg(index, Some(r.v)),
                    seg(index, r.d),
                    seg(index, Some(r.j)),
                    c,
                    r.supporting_features.to_string(),
                    r.receptor_linkage.rediscovery_reads.to_string(),
                    r.receptor_linkage.constant_link_fragments.to_string(),
                    r.receptor_linkage.constant_spanning_reads.to_string(),
                    dna(&r.naive_recombination),
                ]
            } else {
                vec![
                    r.chain.to_string(),
                    seg(index, Some(r.v)),
                    seg(index, Some(r.j)),
                    c,
                    r.supporting_features.to_string(),
                    r.receptor_linkage.rediscovery_reads.to_string(),
                    r.receptor_linkage.constant_link_fragments.to_string(),
                    r.receptor_linkage.constant_spanning_reads.to_string(),
                    dna(&r.naive_recombination),
                ]
            }
        }
    }
}
fn writer(path: PathBuf) -> Result<BufWriter<File>> {
    if let Some(p) = path.parent() {
        fs::create_dir_all(p)?
    }
    Ok(BufWriter::new(
        File::create(&path).with_context(|| format!("creating {}", path.display()))?,
    ))
}
fn seg(index: &VdjIndex, id: Option<crate::SegmentId>) -> String {
    id.and_then(|x| index.segment(x))
        .map(|x| x.name.clone())
        .unwrap_or_default()
}
fn dna(x: &[u8]) -> String {
    String::from_utf8_lossy(x).into_owned()
}
fn opt(x: Option<u16>) -> String {
    x.map(|x| x.to_string()).unwrap_or_default()
}
fn airr_tf(x: bool, unknown: bool) -> String {
    if unknown {
        String::new()
    } else if x {
        "T".into()
    } else {
        "F".into()
    }
}
fn fasta_token(x: &str) -> String {
    x.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '-' | '_' | '.' | ':') {
                c
            } else {
                '_'
            }
        })
        .collect()
}
