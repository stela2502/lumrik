use crate::Ommverse;
use anyhow::{Context, Result, bail};
use gtf_splice_index::types::RefBlock;
use std::collections::{HashMap, HashSet};
use std::path::Path;

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

fn svg_header(w: f64, h: f64) -> String {
    format!(r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial,sans-serif;fill:#222}} .small{{font-size:11px}} .label{{font-size:12px}} .title{{font-size:18px;font-weight:bold}}</style>"#)
}

#[derive(Clone, Debug, Default)]
pub struct GeneTranscriptPlotOptions {
    /// Restrict the plot to these transcript names/accessions. Empty means all.
    pub transcripts: Vec<String>,
    /// Optional free-form notes rendered at the right of a transcript row.
    pub row_notes: HashMap<String, Vec<String>>,
}

/// Plot every annotated transcript of one gene together with the genomic
/// InterPro feature models retained by each transcript CDS.
pub fn plot_gene_transcripts(omm: &Ommverse, gene_name: &str, path: &Path) -> Result<()> {
    plot_gene_transcripts_with_options(omm, gene_name, path, &GeneTranscriptPlotOptions::default())
}

pub fn plot_gene_transcripts_with_options(
    omm: &Ommverse,
    gene_name: &str,
    path: &Path,
    options: &GeneTranscriptPlotOptions,
) -> Result<()> {
    let gene = omm
        .gene(gene_name)
        .with_context(|| format!("gene {gene_name} not found in Ommverse"))?;
    let requested: HashSet<&str> = options.transcripts.iter().map(String::as_str).collect();
    let plotted: Vec<usize> = gene
        .transcript_ids()
        .iter()
        .copied()
        .filter(|&tx_id| {
            requested.is_empty()
                || omm.splice.transcripts[tx_id]
                    .names
                    .iter()
                    .any(|name| requested.contains(name.as_str()))
        })
        .collect();
    if plotted.is_empty() {
        bail!("gene {gene_name} has no transcripts");
    }

    let mut gstart = u32::MAX;
    let mut gend = 0u32;
    for &tx_id in &plotted {
        let tx = &omm.splice.transcripts[tx_id];
        for b in tx.exons() {
            gstart = gstart.min(b.start);
            gend = gend.max(b.end);
        }
    }
    if gstart >= gend {
        bail!("gene {gene_name} has no plottable genomic span");
    }

    let chr_id = omm.splice.transcripts[plotted[0]].chr_id;
    let chr_name = omm
        .splice
        .chr_names
        .get(chr_id)
        .context("gene chromosome is missing from Ommverse splice index")?;
    let feature_ids = omm
        .protein_features
        .candidates_for_region(chr_id, gstart, gend);

    let w = 1500.0;
    let left = 250.0;
    let right = 360.0;
    let top = 95.0;
    let row_h = (70.0 + 14.0 * feature_ids.len() as f64).max(112.0);
    let bottom = 70.0;
    let h = top + row_h * plotted.len() as f64 + bottom;
    let plot_w = w - left - right;
    let gx = |p: u32| {
        left + plot_w * (p.saturating_sub(gstart) as f64 / (gend - gstart) as f64)
    };

    let mut s = svg_header(w, h);
    s.push_str(&format!(
        r##"<text x="{left}" y="30" class="title">{} transcript structure and InterPro architecture</text><text x="{left}" y="51" class="label">InterPro features are stored as spliced genomic models and intersected with each transcript CDS.</text><text x="{left}" y="70" class="small">{}:{}–{} (0-based half-open); exon geometry is drawn to genomic scale.</text>"##,
        xml_escape(gene_name), xml_escape(chr_name), gstart, gend
    ));

    for (ri, &tx_id) in plotted.iter().enumerate() {
        let tx = &omm.splice.transcripts[tx_id];
        let y = top + ri as f64 * row_h + 24.0;
        let name = omm.splice.transcript_name(tx_id).unwrap_or("NA");
        let enst = tx
            .names
            .iter()
            .find(|x| x.starts_with("ENST"))
            .map(String::as_str)
            .unwrap_or("");
        s.push_str(&format!(
            r##"<text x="{:.1}" y="{:.1}" text-anchor="end" class="label" font-weight="bold">{}</text><text x="{:.1}" y="{:.1}" text-anchor="end" class="small">{}</text>"##,
            left - 12.0, y + 4.0, xml_escape(name), left - 12.0, y + 19.0, xml_escape(enst)
        ));
        s.push_str(&format!(
            r##"<line x1="{:.1}" y1="{y:.1}" x2="{:.1}" y2="{y:.1}" stroke="#777" stroke-width="1"/>"##,
            gx(tx.exons().first().map(|b| b.start).unwrap_or(gstart)),
            gx(tx.exons().last().map(|b| b.end).unwrap_or(gend))
        ));
        for b in tx.exons() {
            let x1 = gx(b.start);
            let x2 = gx(b.end);
            s.push_str(&format!(
                r##"<rect x="{x1:.1}" y="{:.1}" width="{:.1}" height="12" fill="#c7c7c7" stroke="#666" stroke-width="0.7"/>"##,
                y - 6.0, (x2 - x1).max(1.0)
            ));
        }
        if let Some((cs, ce)) = tx.cds_span() {
            for exon in tx.exons() {
                let a = exon.start.max(cs);
                let b = exon.end.min(ce);
                if a < b {
                    let x1 = gx(a);
                    let x2 = gx(b);
                    s.push_str(&format!(
                        r##"<rect x="{x1:.1}" y="{:.1}" width="{:.1}" height="8" fill="#444"/>"##,
                        y - 4.0, (x2 - x1).max(1.0)
                    ));
                }
            }
        }

        let cds_blocks: Vec<RefBlock> = if let Some((cs, ce)) = tx.cds_span() {
            tx.exons()
                .iter()
                .filter_map(|exon| {
                    let start = exon.start.max(cs);
                    let end = exon.end.min(ce);
                    (start < end).then(|| RefBlock::new(start, end))
                })
                .collect()
        } else {
            Vec::new()
        };

        let mut feature_y = y + 25.0;
        let mut seen_features = HashSet::new();
        for &feature_id in &feature_ids {
            let Some(feature) = omm.protein_features.features.get(feature_id) else {
                continue;
            };
            if feature.chr_id != chr_id {
                continue;
            }

            let mut blocks = Vec::new();
            let mut feature_bases = 0u32;
            let mut retained_bases = 0u32;
            for fb in &feature.blocks {
                feature_bases += fb.end.saturating_sub(fb.start);
                for cb in &cds_blocks {
                    let start = fb.start.max(cb.start);
                    let end = fb.end.min(cb.end);
                    if start < end {
                        retained_bases += end - start;
                        blocks.push(RefBlock::new(start, end));
                    }
                }
            }
            if blocks.is_empty() {
                continue;
            }

            let key = (feature.label.clone(), feature.signature.clone(), blocks.clone());
            if !seen_features.insert(key) {
                continue;
            }
            let complete = retained_bases == feature_bases;
            let first_x = gx(blocks.first().unwrap().start);
            let last_x = gx(blocks.last().unwrap().end);
            if blocks.len() > 1 {
                s.push_str(&format!(
                    r##"<line x1="{first_x:.1}" y1="{feature_y:.1}" x2="{last_x:.1}" y2="{feature_y:.1}" stroke="#777" stroke-width="1" stroke-dasharray="3,3"/>"##
                ));
            }
            for b in &blocks {
                let x1 = gx(b.start);
                let x2 = gx(b.end);
                s.push_str(&format!(
                    r##"<rect x="{x1:.1}" y="{:.1}" width="{:.1}" height="8" rx="2" fill="#2f6fbb"{} />"##,
                    feature_y - 4.0,
                    (x2 - x1).max(1.5),
                    if complete { "" } else { r#" opacity="0.45""# }
                ));
            }
            let base_label = if feature.description.is_empty() {
                feature.label.clone()
            } else {
                format!("{} · {}", feature.label, feature.description)
            };
            let label = if complete {
                base_label
            } else {
                format!("{base_label} · partial {retained_bases}/{feature_bases} bp")
            };
            s.push_str(&format!(
                r##"<text x="{:.1}" y="{:.1}" class="small">{}</text>"##,
                last_x + 6.0,
                feature_y + 4.0,
                xml_escape(&label)
            ));
            feature_y += 14.0;
        }
        if seen_features.is_empty() {
            s.push_str(&format!(
                r##"<text x="{left}" y="{:.1}" class="small">No genomic InterPro feature overlaps this transcript CDS</text>"##,
                y + 31.0
            ));
        }

        let row_key = tx.names.iter().find(|x| x.starts_with("ENST")).map(String::as_str).unwrap_or(name);
        if let Some(notes) = options.row_notes.get(row_key) {
            let rx = w - right + 25.0;
            for (ni, note) in notes.iter().enumerate() {
                s.push_str(&format!(
                    r##"<text x="{rx}" y="{:.1}" class="small">{}</text>"##,
                    y - 3.0 + ni as f64 * 17.0,
                    xml_escape(note)
                ));
            }
        }
    }

    s.push_str(&format!(
        r##"<text x="{left}" y="{:.1}" class="small">Grey exon = transcribed exon · dark inset = CDS · blue blocks = genomic InterPro feature retained by this transcript · faded blue = partial feature · dashed connector = one genomic feature spanning splice junction(s)</text>"##,
        h - 28.0
    ));
    s.push_str("</svg>");
    std::fs::write(path, s)?;
    Ok(())
}
