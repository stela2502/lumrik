use std::fmt;
use std::fs;
use std::path::Path;
use std::time::Duration;

use anyhow::{Context, Result};

use bam_tide::results::QuantData;
use mapping_info::MappingInfo;

use crate::fast_features::FastFeatures;
use crate::progress::RunProgress;

#[derive(Debug, Default)]
pub struct RunSummary {
    pub elapsed: Duration,

    pub input: InputSummary,
    pub mapping: MappingSummary,
    pub quantification: QuantificationSummary,

    pub fast_features: Option<FastFeatureSummary>,
    pub beacon: Option<BeaconSummary>,
}

#[derive(Debug, Default)]
pub struct InputSummary {
    pub reads_processed: usize,
    pub average_reads_per_second: f64,
}

#[derive(Debug, Default)]
pub struct MappingSummary {
    pub mapped: usize,
    pub unmapped: usize,
    pub duplicates: usize,
    pub filtered: usize,
}

#[derive(Debug, Default)]
pub struct QuantificationSummary {
    pub called_cells: usize,

    pub exonic_features: usize,
    pub exonic_cells: usize,
    pub exonic_entries: usize,

    pub intronic_features: usize,
    pub intronic_cells: usize,
    pub intronic_entries: usize,
}


#[derive(Debug, Default)]
pub struct FastFeatureSummary {
    pub reference_features: usize,

    pub droplets_with_signal: usize,
    pub called_cells_with_signal: usize,

    pub observed_feature_ids: usize,
}

#[derive(Debug, Default)]
pub struct BeaconSummary {
    pub background_droplets: usize,
    pub ambient_umis: usize,

    pub iterations: usize,
    pub mathematically_converged: bool,

    pub observed_pairs: usize,
    pub called_pairs: usize,

    pub cells_no_feature: usize,
    pub cells_single_feature: usize,
    pub cells_multi_feature: usize,
}




impl RunSummary {
    pub fn from_run(
        progress: &RunProgress,
        data: &QuantData,
        called_cells: usize,
        fast_features: Option<&FastFeatures>,
        beacon: Option<&sc_beacon::BeaconResult>,
    ) -> Self {
        let input =
            InputSummary {
                reads_processed:
                    progress.reads_seen(),

                average_reads_per_second:
                    progress
                        .average_reads_per_second(),
            };

        let mapping =
            MappingSummary::from_report(
                &data.report,
            );

        let (
            exonic_features,
            exonic_cells,
            exonic_entries,
        ) = data.gene.dimensions();

        let (
            intronic_features,
            intronic_cells,
            intronic_entries,
        ) = data.intron.dimensions();

        let quantification = QuantificationSummary {
            called_cells,

            exonic_features,
            exonic_cells,
            exonic_entries,

            intronic_features,
            intronic_cells,
            intronic_entries,
        };

        let fast_features =
            fast_features.map(
                FastFeatureSummary::from_features
            );

        let beacon =
            beacon.map(
                BeaconSummary::from_beacon
            );

        Self {
            elapsed:
                progress.elapsed(),

            input,
            mapping,
            quantification,
            fast_features,
            beacon,
        }
    }
}

impl MappingSummary {
    pub fn from_report(
        report: &MappingInfo,
    ) -> Self {
        Self {
            mapped:
                report.ok_reads,

            duplicates:
                report.pcr_duplicates,

            // Replace with actual MappingInfo fields
            // if they exist.
            unmapped: 0,
            filtered: 0,
        }
    }
}

impl FastFeatureSummary {
    pub fn from_features(
        features: &FastFeatures,
    ) -> Self {
        Self {
            reference_features:
                features.mapper
                    .feature_count(),

            droplets_with_signal:
                features.data.len(),

            called_cells_with_signal:
                features.data.len(),

            observed_feature_ids:
                features.data
                    .observed_feature_ids()
                    .len(),
        }
    }
}


impl BeaconSummary {
    pub fn from_beacon(
        beacon: &sc_beacon::BeaconResult,
    ) -> Self {
        let observed_pairs =
            beacon.calls.flat.len();

        let called_pairs =
            beacon.calls
                .flat
                .iter()
                .filter(|call| call.called)
                .count();

        let mut cells_no_feature = 0usize;
        let mut cells_single_feature = 0usize;
        let mut cells_multi_feature = 0usize;

        for row in &beacon.assignments.rows {
            match row.n_called_guides {
                0 => {
                    cells_no_feature += 1;
                }

                1 => {
                    cells_single_feature += 1;
                }

                _ => {
                    cells_multi_feature += 1;
                }
            }
        }

        Self {
            background_droplets:
                beacon
                    .ambient
                    .background_droplets,

            ambient_umis:
                beacon
                    .ambient
                    .total_umis
                    as usize,

            iterations:
                beacon
                    .fitted
                    .iterations,

            mathematically_converged:
                beacon
                    .fitted
                    .mathematical_converged,

            observed_pairs,
            called_pairs,

            cells_no_feature,
            cells_single_feature,
            cells_multi_feature,
        }
    }
}

fn format_duration(
    duration: Duration,
) -> String {
    let secs =
        duration.as_secs();

    let hours =
        secs / 3600;

    let minutes =
        (secs % 3600) / 60;

    let seconds =
        secs % 60;

    if hours > 0 {
        format!(
            "{hours}h {minutes:02}m {seconds:02}s"
        )
    } else if minutes > 0 {
        format!(
            "{minutes}m {seconds:02}s"
        )
    } else {
        format!("{seconds}s")
    }
}

impl fmt::Display for RunSummary {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        writeln!(
            f,
            "Nelrune run complete"
        )?;

        writeln!(
            f,
            "===================="
        )?;

        writeln!(f)?;

        writeln!(
            f,
            "Input"
        )?;

        writeln!(
            f,
            "-----"
        )?;

        writeln!(
            f,
            "Reads processed:          {:>12}",
            self.input.reads_processed,
        )?;

        writeln!(
            f,
            "Elapsed time:             {:>12}",
            format_duration(self.elapsed),
        )?;

        writeln!(
            f,
            "Average throughput:       {:>12.0} reads/s",
            self.input.average_reads_per_second,
        )?;

        writeln!(f)?;

        writeln!(
            f,
            "Mapping"
        )?;

        writeln!(
            f,
            "-------"
        )?;

        writeln!(
            f,
            "Mapped/accepted reads:    {:>12}",
            self.mapping.mapped,
        )?;

        writeln!(
            f,
            "PCR duplicates:           {:>12}",
            self.mapping.duplicates,
        )?;

        writeln!(f)?;

        writeln!(
            f,
            "Quantification"
        )?;

        writeln!(
            f,
            "--------------"
        )?;

        writeln!(
            f,
            "Called cells:             {:>12}",
            self.quantification.called_cells,
        )?;

        writeln!(
            f,
            "Observed exonic features: {:>12}",
            self.quantification.exonic_features,
        )?;

        writeln!(
            f,
            "Observed intron features: {:>12}",
            self.quantification.intronic_features,
        )?;

        writeln!(
            f,
            "Exonic sparse entries:    {:>12}",
            self.quantification.exonic_entries,
        )?;

        writeln!(
            f,
            "Intronic sparse entries:  {:>12}",
            self.quantification.intronic_entries,
        )?;

        if let Some(features) =
            &self.fast_features
        {
            writeln!(f)?;

            writeln!(
                f,
                "Fast features"
            )?;

            writeln!(
                f,
                "-------------"
            )?;

            writeln!(
                f,
                "Reference features:      {:>12}",
                features.reference_features,
            )?;

            writeln!(
                f,
                "Cells with signal:       {:>12}",
                features.called_cells_with_signal,
            )?;

            writeln!(
                f,
                "Observed feature IDs:    {:>12}",
                features.observed_feature_ids,
            )?;
        }

        if let Some(beacon) =
            &self.beacon
        {
            writeln!(f)?;

            writeln!(
                f,
                "Beacon"
            )?;

            writeln!(
                f,
                "------"
            )?;

            writeln!(
                f,
                "Background droplets:     {:>12}",
                beacon.background_droplets,
            )?;

            writeln!(
                f,
                "Ambient feature UMIs:    {:>12}",
                beacon.ambient_umis,
            )?;

            writeln!(
                f,
                "Iterations:              {:>12}",
                beacon.iterations,
            )?;

            writeln!(
                f,
                "Mathematical convergence:{:>12}",
                if beacon.mathematically_converged {
                    "yes"
                } else {
                    "no"
                },
            )?;

            writeln!(
                f,
                "Observed cell/features:  {:>12}",
                beacon.observed_pairs,
            )?;

            writeln!(
                f,
                "Called genuine pairs:    {:>12}",
                beacon.called_pairs,
            )?;

            writeln!(
                f,
                "Cells with no feature:   {:>12}",
                beacon.cells_no_feature,
            )?;

            writeln!(
                f,
                "Cells with one feature:  {:>12}",
                beacon.cells_single_feature,
            )?;

            writeln!(
                f,
                "Cells with >1 feature:   {:>12}",
                beacon.cells_multi_feature,
            )?;
        }

        Ok(())
    }
}

impl RunSummary {
    pub fn write<P: AsRef<Path>>(
        &self,
        path: P,
    ) -> Result<()> {
        let path =
            path.as_ref();

        fs::write(
            path,
            self.to_string(),
        )
        .with_context(|| {
            format!(
                "writing Nelrune summary {}",
                path.display()
            )
        })
    }
}