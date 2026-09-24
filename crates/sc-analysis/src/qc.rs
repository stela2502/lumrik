use crate::data::SingleCellData;

#[derive(Debug, Clone)]
pub(crate) struct QcMetrics {
    pub total: Vec<f64>,
    pub detected: Vec<usize>,
    pub mito: Vec<f64>,
    pub ribo: Vec<f64>,
    pub vdj: Vec<f64>,
    pub surviving: Vec<f64>,
    pub keep_gene: Vec<bool>,
}

pub(crate) fn is_vdj_segment_gene(name: &str) -> bool {
    let u = name.to_ascii_uppercase();
    const PREFIXES: &[&str] = &[
        "IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ", "TRAV", "TRAJ", "TRBV", "TRBD",
        "TRBJ", "TRGV", "TRGJ", "TRDV", "TRDD", "TRDJ",
    ];
    PREFIXES.iter().any(|prefix| u.starts_with(prefix))
}

pub(crate) fn classify_gene(name: &str) -> (bool, bool, bool) {
    let u = name.to_ascii_uppercase();
    let mito = u.starts_with("MT-");
    let ribo = u.starts_with("RPS") || u.starts_with("RPL");
    let vdj = is_vdj_segment_gene(name);
    (mito, ribo, vdj)
}

pub(crate) fn qc(data: &SingleCellData, exclude_vdj: bool) -> QcMetrics {
    let n = data.n_cells();
    let mut total = vec![0.0; n];
    let mut detected = vec![0usize; n];
    let mut mito = vec![0.0; n];
    let mut ribo = vec![0.0; n];
    let mut vdj = vec![0.0; n];
    let mut keep_gene = vec![true; data.n_features()];
    for (g, row) in data.matrix.outer_iterator().enumerate() {
        let (is_mito, is_ribo, is_vdj) = classify_gene(&data.features[g]);
        keep_gene[g] = !(is_mito || is_ribo || (exclude_vdj && is_vdj));
        for (c, v) in row.iter() {
            if *v > 0.0 {
                let x = *v as f64;
                total[c] += x;
                detected[c] += 1;
                if is_mito {
                    mito[c] += x;
                }
                if is_ribo {
                    ribo[c] += x;
                }
                if is_vdj {
                    vdj[c] += x;
                }
            }
        }
    }
    let surviving = (0..n)
        .map(|c| {
            let excluded_vdj = if exclude_vdj { vdj[c] } else { 0.0 };
            (total[c] - mito[c] - ribo[c] - excluded_vdj).max(0.0)
        })
        .collect();
    QcMetrics {
        total,
        detected,
        mito,
        ribo,
        vdj,
        surviving,
        keep_gene,
    }
}

#[cfg(test)]
mod tests {
    use super::is_vdj_segment_gene;

    #[test]
    fn receptor_vdj_segments_are_classified_without_constants() {
        for gene in [
            "Ighv1-39", "Ighd1-1", "Ighj3", "Igkv8-24", "Igkj2", "Iglv1", "Iglj1", "Trav1",
            "Traj1", "Trbv1", "Trbd1", "Trbj1", "Trgv1", "Trgj1", "Trdv1", "Trdd1", "Trdj1",
            "IGHV3-23", "TRBJ2-7",
        ] {
            assert!(
                is_vdj_segment_gene(gene),
                "{gene} should be a V/D/J segment"
            );
        }
        for gene in [
            "Ighm", "Ighg1", "Igha", "Igkc", "Iglc1", "Trac", "Trbc1", "IGHM", "IGKC", "TRAC",
            "Cd79a", "Ms4a1", "Jchain",
        ] {
            assert!(
                !is_vdj_segment_gene(gene),
                "{gene} should remain ordinary expression"
            );
        }
    }
}
