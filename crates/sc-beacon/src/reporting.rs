use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};

use crate::background::AmbientModel;
use crate::caller::GuideCalls;
use crate::dataset::GuideDataset;
use crate::model::FittedModel;
use crate::tenx::GuideFeatureIndex;

impl AmbientModel {
    fn write_tsv<W: Write>(
        &self,
        writer: &mut W,
        index: &GuideFeatureIndex,
    ) -> Result<()> {
        writeln!(writer, "guide_id\tguide_name\tambient_umis\tp_g")?;

        for (guide_id, feature) in index.guides().iter().enumerate() {
            writeln!(
                writer,
                "{}\t{}\t{}\t{:.12}",
                feature.id,
                feature.name,
                self.guide_umis[guide_id],
                self.guide_probability[guide_id],
            )?;
        }

        Ok(())
    }

    pub fn print_table(
        &self,
        index: &GuideFeatureIndex,
    ) -> Result<()> {
        let stdout = std::io::stdout();
        let mut writer = stdout.lock();
        self.write_tsv(&mut writer, index)
    }

    pub fn write_table(
        &self,
        out: &PathBuf,
        index: &GuideFeatureIndex,
    ) -> Result<()> {
        let path = out.join("ambient_guides.tsv");
        let file = File::create(&path)
            .with_context(|| format!("creating {}", path.display()))?;
        let mut writer = BufWriter::new(file);

        self.write_tsv(&mut writer, index)?;
        writer
            .flush()
            .with_context(|| format!("flushing {}", path.display()))?;

        Ok(())
    }
}


impl FittedModel {
    fn write_tsv<W: Write>(
        &self,
        writer: &mut W,
        index: &GuideFeatureIndex,
    ) -> Result<()> {
        writeln!(
            writer,
            "guide_id\tguide_name\tprior_real\ttrue_mean\ttheta"
        )?;

        for (guide_id, model) in self.guides.iter().enumerate() {
            let feature = &index.guides()[guide_id];

            writeln!(
                writer,
                "{}\t{}\t{:.8}\t{:.8}\t{:.8}",
                feature.id,
                feature.name,
                model.prior_real,
                model.mean,
                model.theta,
            )?;
        }

        Ok(())
    }

    pub fn print_table(
        &self,
        index: &GuideFeatureIndex,
    ) -> Result<()> {
        let stdout = std::io::stdout();
        let mut writer = stdout.lock();
        self.write_tsv(&mut writer, index)
    }

    pub fn write_table(
        &self,
        out: &PathBuf,
        index: &GuideFeatureIndex,
    ) -> Result<()> {
        let path = out.join("guide_models.tsv");
        let file = File::create(&path)
            .with_context(|| format!("creating {}", path.display()))?;
        let mut writer = BufWriter::new(file);

        self.write_tsv(&mut writer, index)?;
        writer
            .flush()
            .with_context(|| format!("flushing {}", path.display()))?;

        Ok(())
    }
}


impl GuideCalls {
    fn write_tsv<W: Write>(
        &self,
        writer: &mut W,
        index: &GuideFeatureIndex,
        data: &GuideDataset,
    ) -> Result<()> {
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
            let feature = &index.guides()[call.guide_id as usize];

            let barcode = data
                .barcode_by_id
                .get(&call.cell_id)
                .map(String::as_str)
                .unwrap_or("UNKNOWN");

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{:.8}\t{:.12}\t{:.8}\t{:.8}\t{:.8}\t{:.4e}\t{:.4e}\t{}",
                barcode,
                feature.id,
                feature.name,
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

    pub fn print_table(
        &self,
        index: &GuideFeatureIndex,
        data: &GuideDataset,
    ) -> Result<()> {
        let stdout = std::io::stdout();
        let mut writer = stdout.lock();
        self.write_tsv(&mut writer, index, data)
    }

    pub fn write_table(
        &self,
        out: &PathBuf,
        index: &GuideFeatureIndex,
        data: &GuideDataset,
    ) -> Result<()> {
        let path = out.join("guide_calls.tsv");
        let file = File::create(&path)
            .with_context(|| format!("creating {}", path.display()))?;
        let mut writer = BufWriter::new(file);

        self.write_tsv(&mut writer, index, data)?;
        writer
            .flush()
            .with_context(|| format!("flushing {}", path.display()))?;

        Ok(())
    }
}
