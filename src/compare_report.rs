use std::fmt::{Display, Formatter, Result};

#[derive(Debug, Clone, Default)]
pub struct MaxDiffBin {
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub x: f64,
    pub y: f64,
    pub abs_diff: f64,
}

impl Display for MaxDiffBin {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(
            f,
            "chr:start-end\tpython\trust\tabs_diff\n{}:{}-{}\t{:.6}\t{:.6}\t{:.6}",
            self.chr, self.start, self.end, self.x, self.y, self.abs_diff
        )
    }
}

#[derive(Debug, Default, Clone)]
pub struct CompareReport {
    /// Total number of bins compared
    pub n: usize,

    /// Number of bins where |x - y| > eps
    pub n_over_eps: usize,

    pub max_abs: f64,

    /// Sum of absolute differences
    pub sum_abs: f64,

    /// Sum of squared differences ( (x - y)^2 )
    pub sum_sq: f64,

    // --- for Pearson correlation ---
    sum_x: f64,
    sum_y: f64,
    sum_x2: f64,
    sum_y2: f64,
    sum_xy: f64,

    // NEW:
    pub max_bin: Option<MaxDiffBin>,
}

impl CompareReport {
    /// Update this report from two equal-length bin vectors.
    pub fn update_from_bins(
        &mut self,
        chr: &str,
        chr_len: u64,
        bin_width: u64,
        x: &[f64],
        y: &[f64],
        eps: f64,
    ) {
        debug_assert_eq!(x.len(), y.len(), "bin vectors must have equal length");

        for (i, (&xi, &yi)) in x.iter().zip(y.iter()).enumerate() {
            let start = (i as u64) * bin_width;
            if start >= chr_len {
                break;
            }
            let end = (start + bin_width).min(chr_len);

            self.update_one(chr, start, end, xi, yi, eps);
        }
    }

    /// Update this report with a single pair of values.
    #[inline]
    pub fn update_one(&mut self, chr: &str, start: u64, end: u64, x: f64, y: f64, eps: f64) {
        let d = x - y;
        let ad = d.abs();

        self.n += 1;

        if ad > eps {
            self.n_over_eps += 1;
        }

        if ad > self.max_abs {
            self.max_abs = ad;
            self.max_bin = Some(MaxDiffBin {
                chr: chr.to_string(),
                start,
                end,
                x,
                y,
                abs_diff: ad,
            });
        }

        self.sum_abs += ad;
        self.sum_sq += d * d;

        self.sum_x += x;
        self.sum_y += y;
        self.sum_x2 += x * x;
        self.sum_y2 += y * y;
        self.sum_xy += x * y;
    }

    /// Merge another report into this one (for TOTAL aggregation).
    pub fn merge(&mut self, other: &CompareReport) {
        self.n += other.n;
        self.n_over_eps += other.n_over_eps;
        self.max_abs = self.max_abs.max(other.max_abs);
        self.sum_abs += other.sum_abs;
        self.sum_sq += other.sum_sq;

        self.sum_x += other.sum_x;
        self.sum_y += other.sum_y;
        self.sum_x2 += other.sum_x2;
        self.sum_y2 += other.sum_y2;
        self.sum_xy += other.sum_xy;

        self.max_bin = match (self.max_bin.take(), other.max_bin.clone()) {
            (None, None) => None,
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (Some(a), Some(b)) => {
                if a.abs_diff >= b.abs_diff {
                    Some(a)
                } else {
                    Some(b)
                }
            }
        };
    }

    /// Finalize metrics:
    /// (
    ///   n_over_eps,
    ///   frac_n_over_eps,
    ///   mean_abs,
    ///   var_abs,
    ///   rmse,
    ///   max_abs
    /// )
    pub fn finish(&self) -> (usize, f64, f64, f64, f64, f64) {
        if self.n == 0 {
            (0, 0.0, 0.0, 0.0, 0.0, 0.0)
        } else {
            let n = self.n as f64;

            let n_over_eps = self.n_over_eps;
            let frac_n_over_eps = (self.n_over_eps as f64) / n;

            let mean_abs = self.sum_abs / n;

            // var(|d|) = E[d^2] - (E[|d|])^2, and d^2 == |d|^2
            let var_abs = (self.sum_sq / n) - (mean_abs * mean_abs);

            let rmse = (self.sum_sq / n).sqrt();

            (
                n_over_eps,
                frac_n_over_eps,
                mean_abs,
                var_abs.max(0.0),
                rmse,
                self.max_abs,
            )
        }
    }

    pub fn finish_names() -> (
        &'static str,
        &'static str,
        &'static str,
        &'static str,
        &'static str,
        &'static str,
    ) {
        (
            "n_over_eps",
            "frac_n_over_eps",
            "mean_abs",
            "var_abs",
            "rmse",
            "max_abs",
        )
    }

    /// Pearson correlation between X and Y.
    pub fn pearson(&self) -> f64 {
        if self.n < 2 {
            return 0.0;
        }
        let n = self.n as f64;

        let num = n * self.sum_xy - self.sum_x * self.sum_y;
        let den_x = n * self.sum_x2 - self.sum_x * self.sum_x;
        let den_y = n * self.sum_y2 - self.sum_y * self.sum_y;

        if den_x <= 0.0 || den_y <= 0.0 {
            return 0.0; // degenerate case (constant vector)
        }

        let rho = num / (den_x.sqrt() * den_y.sqrt());
        rho.clamp(-1.0, 1.0)
    }
}
