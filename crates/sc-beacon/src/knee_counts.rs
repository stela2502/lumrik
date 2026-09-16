use anyhow::{bail, Result};

#[derive(Debug, Clone)]
pub struct KneeCountFit {
    pub knee_rank: usize,
    pub umi_cutoff: u32,
    pub candidate_barcodes: usize,
    pub informative_barcodes: usize,
    pub score: f64,
}

/// Detect the knee of the descending barcode-rank curve in log10(rank) x log10(UMI) space.
///
/// The terminal one-UMI plateau is excluded from knee fitting because it contains no rank-curve
/// shape information. The knee is the point with the greatest perpendicular distance from the
/// chord joining the first and last informative points. All barcodes with UMI count >= the UMI
/// count at that knee are called.
pub fn fit_knee_counts(counts: &[u32]) -> Result<KneeCountFit> {
    let mut sorted: Vec<u32> = counts.iter().copied().filter(|&n| n > 0).collect();
    if sorted.len() < 3 {
        bail!("need at least three positive barcode counts for knee calling");
    }
    sorted.sort_unstable_by(|a, b| b.cmp(a));

    let informative = sorted.iter().take_while(|&&n| n > 1).count();
    if informative < 3 {
        bail!("need at least three barcodes above the one-UMI plateau for knee calling");
    }

    let x0 = 1.0_f64.log10();
    let y0 = (sorted[0] as f64).log10();
    let x1 = (informative as f64).log10();
    let y1 = (sorted[informative - 1] as f64).log10();
    let dx = x1 - x0;
    let dy = y1 - y0;
    let denom = (dx * dx + dy * dy).sqrt();
    if denom <= f64::EPSILON {
        bail!("barcode-rank curve has no usable extent for knee calling");
    }

    let mut best_idx = 1usize;
    let mut best_distance = f64::NEG_INFINITY;
    for (idx, &count) in sorted[..informative - 1].iter().enumerate().skip(1) {
        let x = ((idx + 1) as f64).log10();
        let y = (count as f64).log10();
        // Signed distance from the endpoint chord. For a descending barcode-rank curve the
        // biologically useful knee lies above/right of the chord; absolute distance makes the
        // method insensitive to orientation while retaining the same extremum.
        let distance = (dy * (x - x0) - dx * (y - y0)).abs() / denom;
        if distance > best_distance {
            best_distance = distance;
            best_idx = idx;
        }
    }

    let umi_cutoff = sorted[best_idx];
    let knee_rank = sorted.partition_point(|&n| n >= umi_cutoff);

    Ok(KneeCountFit {
        knee_rank,
        umi_cutoff,
        candidate_barcodes: sorted.len(),
        informative_barcodes: informative,
        score: best_distance,
    })
}
