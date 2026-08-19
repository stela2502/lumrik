use std::collections::BTreeMap;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};

use crate::utils::{quantile, percent};

use crate::cell_guide_assignments::CellGuideAssignments;

#[derive(Debug, Clone)]
pub struct MultiGuideGapStats {
    n_called_guides: Option<usize>,
    n_cells: usize,

    mean: f64,

    min: f64,
    q1: f64,
    median: f64,
    q3: f64,
    p90: f64,
    p95: f64,
    max: f64,
}

impl MultiGuideGapStats {
    pub fn new(
        n_called_guides: Option<usize>,
        values: &[f64],
    ) -> Self {
        assert!(
            !values.is_empty(),
            "cannot calculate statistics from an empty value set"
        );

        let mut sorted = values.to_vec();
        sorted.sort_unstable_by(f64::total_cmp);

        let n_cells = sorted.len();
        let mean = sorted.iter().sum::<f64>() / n_cells as f64;

        Self {
            n_called_guides,
            n_cells,
            mean,
            min: sorted[0],
            q1: quantile(&sorted, 0.25),
            median: quantile(&sorted, 0.50),
            q3: quantile(&sorted, 0.75),
            p90: quantile(&sorted, 0.90),
            p95: quantile(&sorted, 0.95),
            max: *sorted.last().expect("non-empty values"),
        }
    }

    pub fn collect(
        assignments: &CellGuideAssignments,
    ) -> Vec<Self> {
        let mut by_multiplicity:
            BTreeMap<usize, Vec<f64>> =
            BTreeMap::new();

        let mut all_multi = Vec::new();

        for (n_called, log_odds_gap)
            in assignments.multi_guide_gaps()
        {
            by_multiplicity
                .entry(n_called)
                .or_default()
                .push(log_odds_gap);

            all_multi.push(log_odds_gap);
        }

        let mut result =
            Vec::with_capacity(by_multiplicity.len() + 1);

        for (n_called_guides, values)
            in by_multiplicity
        {
            result.push(
                Self::new(
                    Some(n_called_guides),
                    &values,
                )
            );
        }

        if !all_multi.is_empty() {
            result.push(
                Self::new(
                    None,
                    &all_multi,
                )
            );
        }

        result
    }

    pub fn n_called_guides(&self) -> Option<usize> {
        self.n_called_guides
    }

    pub fn n_cells(&self) -> usize {
        self.n_cells
    }

    pub fn mean(&self) -> f64 {
        self.mean
    }

    pub fn min(&self) -> f64 {
        self.min
    }

    pub fn q1(&self) -> f64 {
        self.q1
    }

    pub fn median(&self) -> f64 {
        self.median
    }

    pub fn q3(&self) -> f64 {
        self.q3
    }

    pub fn p90(&self) -> f64 {
        self.p90
    }

    pub fn p95(&self) -> f64 {
        self.p95
    }

    pub fn max(&self) -> f64 {
        self.max
    }

    pub fn group_name(&self) -> String {
        self.n_called_guides
            .map(|n| n.to_string())
            .unwrap_or_else(|| "ALL".to_string())
    }

    pub fn header() -> &'static str {
        concat!(
            "called_guides",
            "\tcells",
            "\tmean_log_odds_gap",
            "\tmin",
            "\tq1",
            "\tmedian",
            "\tq3",
            "\tp90",
            "\tp95",
            "\tmax"
        )
    }
}

impl Display for MultiGuideGapStats {
    fn fmt(
        &self,
        f: &mut Formatter<'_>,
    ) -> FmtResult {
        write!(
            f,
            "{}\t{}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}",
            self.group_name(),
            self.n_cells,
            self.mean,
            self.min,
            self.q1,
            self.median,
            self.q3,
            self.p90,
            self.p95,
            self.max,
        )
    }
}


pub trait MultiGuideGapStatsTable {
    fn print_table(&self);

    fn print_assignment_summary(
        &self,
        assignments: &CellGuideAssignments,
    ) -> String;

    fn primary_guide_counts(
        &self,
        assignments: &CellGuideAssignments,
        minimum_odds_ratio: f64,
    ) -> String;

    fn write_table(
        &self,
        out: &PathBuf,
    ) -> Result<()>;
}

impl MultiGuideGapStatsTable for [MultiGuideGapStats] {
    fn print_table(&self) {
        println!("Multi-guide log-odds-gap statistics");
        println!("{}", MultiGuideGapStats::header());

        for stat in self {
            println!("{stat}");
        }
    }

    fn write_table(
        &self,
        out: &PathBuf,
    ) -> Result<()> {
        let path =
            out.join("multi_guide_log_odds_gap_stats.tsv");

        let file = File::create(&path)
            .with_context(|| format!("creating {}", path.display()))?;
        let mut writer = BufWriter::new(file);

        writeln!(
            writer,
            "{}",
            MultiGuideGapStats::header()
        )?;

        for stat in self {
            writeln!(writer, "{stat}")?;
        }

        writer
            .flush()
            .with_context(|| format!("flushing {}", path.display()))?;

        Ok(())
    }

    fn primary_guide_counts(
	    &self,
	    assignments: &CellGuideAssignments,
	    minimum_odds_ratio: f64,
	) -> String {
	    use std::collections::BTreeMap;
	    use std::fmt::Write;

	    let minimum_log_odds_gap =
	        minimum_odds_ratio.ln();

	    let mut best_guide_counts:
	        BTreeMap<&str, usize> =
	        BTreeMap::new();

	    for assignment in &assignments.rows {
	        if assignment.n_called_guides == 0 {
	            continue;
	        }

	        if assignment.log_odds_gap
	            < minimum_log_odds_gap
	        {
	            continue;
	        }

	        *best_guide_counts
	            .entry(
	                assignment.best_guide.as_str()
	            )
	            .or_default() += 1;
	    }

	    let total_confident: usize =
	        best_guide_counts
	            .values()
	            .sum();

	    let mut out =
	        String::new();

	    writeln!(
	        out,
	        "\nHigh-confidence primary-guide counts"
	    ).unwrap();

	    writeln!(
	        out,
	        "------------------------------------"
	    ).unwrap();

	    writeln!(
	        out,
	        "{:<32} {:>8}",
	        "guide",
	        "cells",
	    ).unwrap();

	    for (guide, count) in &best_guide_counts {
	        writeln!(
	            out,
	            "{:<32} {:>8}",
	            guide,
	            count,
	        ).unwrap();
	    }

	    writeln!(
	        out,
	        "{:<32} {:>8}",
	        "TOTAL",
	        total_confident,
	    ).unwrap();

	    out
	}

	fn print_assignment_summary(
	    &self,
	    assignments: &CellGuideAssignments,
	) -> String {
	    use std::fmt::Write;

	    const LOG_ODDS_10: f64 = std::f64::consts::LN_10;
	    const LOG_ODDS_100: f64 =
	        2.0 * std::f64::consts::LN_10;

	    let total = assignments.rows.len();

	    let mut no_guide = 0usize;
	    let mut single = 0usize;
	    let mut multi = 0usize;

	    let mut ambiguous = 0usize;
	    let mut moderate = 0usize;
	    let mut clear = 0usize;

	    let mut min_gap = f64::INFINITY;

	    for row in &assignments.rows {
	        match row.n_called_guides {
	            0 => {
	                no_guide += 1;
	            }

	            1 => {
	                single += 1;
	            }

	            _ => {
	                multi += 1;

	                let gap = row.log_odds_gap;

	                if !gap.is_finite() {
	                    continue;
	                }

	                min_gap = min_gap.min(gap);

	                if gap < LOG_ODDS_10 {
	                    ambiguous += 1;
	                } else if gap < LOG_ODDS_100 {
	                    moderate += 1;
	                } else {
	                    clear += 1;
	                }
	            }
	        }
	    }

	    let mut out = String::new();

	    writeln!(out).unwrap();
	    writeln!(out, "Guide assignment summary").unwrap();
	    writeln!(out, "------------------------").unwrap();

	    writeln!(
	        out,
	        "{:<32} {:>8}  ({:>5.1}%)",
	        "Total cells:",
	        total,
	        100.0,
	    )
	    .unwrap();

	    writeln!(
	        out,
	        "{:<32} {:>8}  ({:>5.1}%)",
	        "No guide assigned:",
	        no_guide,
	        percent(no_guide, total),
	    )
	    .unwrap();

	    writeln!(
	        out,
	        "{:<32} {:>8}  ({:>5.1}%)",
	        "Single guide called:",
	        single,
	        percent(single, total),
	    )
	    .unwrap();

	    writeln!(
	        out,
	        "{:<32} {:>8}  ({:>5.1}%)",
	        "Multiple guides called:",
	        multi,
	        percent(multi, total),
	    )
	    .unwrap();

	    if multi == 0 {
	        return out;
	    }

	    writeln!(out).unwrap();
	    writeln!(out, "Multi-guide resolution").unwrap();
	    writeln!(out, "----------------------").unwrap();

	    writeln!(
	        out,
	        "{:<32} {:>8}  ({:>5.1}%)",
	        "Clear best guide (>100:1):",
	        clear,
	        percent(clear, multi),
	    )
	    .unwrap();

	    writeln!(
	        out,
	        "{:<32} {:>8}  ({:>5.1}%)",
	        "Moderate separation (10-100:1):",
	        moderate,
	        percent(moderate, multi),
	    )
	    .unwrap();

	    writeln!(
	        out,
	        "{:<32} {:>8}  ({:>5.1}%)",
	        "Ambiguous best guide (<10:1):",
	        ambiguous,
	        percent(ambiguous, multi),
	    )
	    .unwrap();

	    /*
	     * `self` already contains the distribution statistics.
	     *
	     * The ALL row is the combined multi-guide population.
	     */
	    if let Some(all) = self
	        .iter()
	        .find(|stats| stats.n_called_guides().is_none())
	    {
	        writeln!(out).unwrap();
	        writeln!(out, "Best-vs-second evidence").unwrap();
	        writeln!(out, "-----------------------").unwrap();

	        writeln!(
	            out,
	            "{:<32} {:>12.2}",
	            "Median log-odds difference:",
	            all.median(),
	        )
	        .unwrap();

	        writeln!(
	            out,
	            "{:<32} {:>12.2}",
	            "Q1 log-odds difference:",
	            all.q1(),
	        )
	        .unwrap();

	        writeln!(
	            out,
	            "{:<32} {:>12.2}",
	            "Smallest log-odds difference:",
	            all.min(),
	        )
	        .unwrap();
	    }

	    if min_gap.is_finite() && min_gap < 700.0 {
	        writeln!(
	            out,
	            "{:<32} {:>12.1}:1",
	            "Weakest odds advantage:",
	            min_gap.exp(),
	        )
	        .unwrap();
	    }

	    writeln!(out).unwrap();

	    out
	}
}


