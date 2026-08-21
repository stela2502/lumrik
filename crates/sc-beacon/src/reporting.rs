use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

use scdata::FeatureIndex;
use int_to_str::IntToStr;

use crate::background::AmbientModel;
use crate::caller::GuideCalls;
use crate::model::FittedModel;


impl AmbientModel {
    fn write_tsv<W, I>(
        &self,
        writer: &mut W,
        index: &I,
    ) -> Result<()>
    where
        W: Write,
        I: FeatureIndex,
    {
        writeln!(
            writer,
            "guide_id\tguide_name\tambient_umis\tp_g"
        )?;

        for guide_id in 0..self.guide_umis.len() {
            let feature_id = guide_id as u64;

            writeln!(
                writer,
                "{}\t{}\t{}\t{:.12}",
                feature_id,
                index.feature_name(feature_id),
                self.guide_umis[guide_id],
                self.guide_probability[guide_id],
            )?;
        }

        Ok(())
    }

    pub fn print_table<I>(
        &self,
        index: &I,
    ) -> Result<()>
    where
        I: FeatureIndex,
    {
        let stdout = std::io::stdout();
        let mut writer = stdout.lock();

        self.write_tsv(
            &mut writer,
            index,
        )
    }

    pub fn write_table<P, I>(
        &self,
        out: P,
        index: &I,
    ) -> Result<()>
    where
        P: AsRef<Path>,
        I: FeatureIndex,
    {
        let out = out.as_ref();

        let path =
            out.join("ambient_guides.tsv");

        let file =
            File::create(&path)
                .with_context(|| {
                    format!(
                        "creating {}",
                        path.display()
                    )
                })?;

        let mut writer =
            BufWriter::new(file);

        self.write_tsv(
            &mut writer,
            index,
        )?;

        writer
            .flush()
            .with_context(|| {
                format!(
                    "flushing {}",
                    path.display()
                )
            })?;

        Ok(())
    }
}


impl FittedModel {
    fn write_tsv<W, I>(
        &self,
        writer: &mut W,
        index: &I,
    ) -> Result<()>
    where
        W: Write,
        I: FeatureIndex,
    {
        writeln!(
            writer,
            "guide_id\tguide_name\tprior_real\ttrue_mean\ttheta"
        )?;

        for (guide_id, model) in
            self.guides.iter().enumerate()
        {
            let feature_id =
                guide_id as u64;

            writeln!(
                writer,
                "{}\t{}\t{:.8}\t{:.8}\t{:.8}",
                feature_id,
                index.feature_name(feature_id),
                model.prior_real,
                model.mean,
                model.theta,
            )?;
        }

        Ok(())
    }

    pub fn print_table<I>(
        &self,
        index: &I,
    ) -> Result<()>
    where
        I: FeatureIndex,
    {
        let stdout =
            std::io::stdout();

        let mut writer =
            stdout.lock();

        self.write_tsv(
            &mut writer,
            index,
        )
    }

    pub fn write_table<P, I>(
        &self,
        out: P,
        index: &I,
    ) -> Result<()>
    where
        P: AsRef<Path>,
        I: FeatureIndex,
    {
        let out =
            out.as_ref();

        let path =
            out.join("guide_models.tsv");

        let file =
            File::create(&path)
                .with_context(|| {
                    format!(
                        "creating {}",
                        path.display()
                    )
                })?;

        let mut writer =
            BufWriter::new(file);

        self.write_tsv(
            &mut writer,
            index,
        )?;

        writer
            .flush()
            .with_context(|| {
                format!(
                    "flushing {}",
                    path.display()
                )
            })?;

        Ok(())
    }
}


impl GuideCalls {
    fn write_tsv<W, I>(
        &self,
        writer: &mut W,
        index: &I,
        cell_len: usize,
    ) -> Result<()>
    where
        W: Write,
        I: FeatureIndex,
    {
        writeln!(
            writer,
            concat!(
                "barcode",
                "\tguide_id",
                "\tguide_name",
                "\tumi_count",
                "\tlambda_c",
                "\tp_g",
                "\texpected_ambient",
                "\tposterior",
                "\tlog_odds",
                "\tambient_p",
                "\tq_value",
                "\tcalled"
            )
        )?;

        for call in &self.flat {
            let feature_id =
                call.guide_id as u64;

            let barcode =
            IntToStr::from_u64(call.cell_id)
                .to_string(cell_len);

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{:.8}\t{:.12}\t{:.8}\t{:.8}\t{:.8}\t{:.4e}\t{:.4e}\t{}",
                barcode,
                feature_id,
                index.feature_name(feature_id),
                call.count,
                call.lambda_cell,
                call.ambient_probability,
                call.expected_ambient,
                call.posterior.probability,
                call.posterior.log_odds,
                call.ambient_p_value,
                call.q_value,
                call.called,
            )?;
        }

        Ok(())
    }

    pub fn print_table<I>(
        &self,
        index: &I,
        data: usize,
    ) -> Result<()>
    where
        I: FeatureIndex,
    {
        let stdout =
            std::io::stdout();

        let mut writer =
            stdout.lock();

        self.write_tsv(
            &mut writer,
            index,
            data,
        )
    }

    pub fn write_table<P, I>(
        &self,
        out: P,
        index: &I,
        data: usize,
    ) -> Result<()>
    where
        P: AsRef<Path>,
        I: FeatureIndex,
    {
        let out =
            out.as_ref();

        let path =
            out.join("guide_calls.tsv");

        let file =
            File::create(&path)
                .with_context(|| {
                    format!(
                        "creating {}",
                        path.display()
                    )
                })?;

        let mut writer =
            BufWriter::new(file);

        self.write_tsv(
            &mut writer,
            index,
            data,
        )?;

        writer
            .flush()
            .with_context(|| {
                format!(
                    "flushing {}",
                    path.display()
                )
            })?;

        Ok(())
    }
}