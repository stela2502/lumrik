use anyhow::{Result, bail};

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

    fn chord_knee(sorted: &[u32], start: usize, end: usize) -> Result<(usize, f64)> {
        if end.saturating_sub(start) < 2 {
            bail!("barcode-rank segment is too short for knee calling");
        }
        let x0 = ((start + 1) as f64).log10();
        let y0 = (sorted[start] as f64).log10();
        let x1 = (end as f64).log10();
        let y1 = (sorted[end - 1] as f64).log10();
        let dx = x1 - x0;
        let dy = y1 - y0;
        let denom = (dx * dx + dy * dy).sqrt();
        if denom <= f64::EPSILON {
            bail!("barcode-rank curve has no usable extent for knee calling");
        }

        let mut best_idx = start + 1;
        let mut best_distance = f64::NEG_INFINITY;
        for idx in (start + 1)..(end - 1) {
            let x = ((idx + 1) as f64).log10();
            let y = (sorted[idx] as f64).log10();
            let distance = (dy * (x - x0) - dx * (y - y0)).abs() / denom;
            if distance > best_distance {
                best_distance = distance;
                best_idx = idx;
            }
        }
        Ok((best_idx, best_distance))
    }

    // A heterogeneous sample can contain more than one genuine bend. The global chord
    // preferentially finds the high-count bend and can therefore discard an entire lower-RNA
    // cell population. After finding that bend, recursively inspect the lower-count segment.
    // If it contains a well-resolved second bend, use the terminal bend as the cell/background
    // boundary. Requiring substantial rank span on both sides prevents tiny tail wiggles from
    // replacing the primary knee.
    let (primary_idx, primary_score) = chord_knee(&sorted, 0, informative)?;
    let mut best_idx = primary_idx;
    let mut best_distance = primary_score;
    let remaining = informative.saturating_sub(primary_idx);
    if primary_idx >= 10 && remaining >= 20 {
        let (secondary_idx, secondary_score) = chord_knee(&sorted, primary_idx, informative)?;
        let left = secondary_idx.saturating_sub(primary_idx);
        let right = informative.saturating_sub(secondary_idx);
        if left >= 10 && right >= 10 && secondary_score >= 0.05 {
            best_idx = secondary_idx;
            best_distance = secondary_score;
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
