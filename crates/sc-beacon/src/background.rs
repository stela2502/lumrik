use anyhow::{bail, Result};

use crate::dataset::GuideDataset;

#[derive(Debug, Clone)]
pub struct BackgroundConfig {
    /// Symmetric Dirichlet pseudocount for guide frequencies.
    pub alpha: f64,
}

impl Default for BackgroundConfig {
    fn default() -> Self {
        Self { alpha: 0.5 }
    }
}

#[derive(Debug, Clone)]
pub struct AmbientModel {
    /// p_g: fraction of ambient guide molecules belonging to guide g.
    pub guide_probability: Vec<f64>,
    pub guide_umis: Vec<u64>,
    pub total_umis: u64,
    pub background_droplets: usize,
}

impl AmbientModel {
    pub fn fit(data: &GuideDataset, cfg: &BackgroundConfig) -> Result<Self> {
        if data.n_guides() == 0 {
            bail!("Cannot fit ambient model without guides");
        }
        if data.n_cells() == 0 {
            bail!("Cannot fit ambient model without background droplets");
        }
        if cfg.alpha < 0.0 {
            bail!("Dirichlet alpha must be >= 0");
        }

        let guide_umis: Vec<u64> = data
            .by_guide
            .iter()
            .map(|obs| obs.iter().map(|x| x.count as u64).sum())
            .collect();

        let total_umis: u64 = guide_umis.iter().sum();
        if total_umis == 0 {
            bail!("No guide UMIs were observed in background droplets");
        }

        let denom = total_umis as f64 + cfg.alpha * guide_umis.len() as f64;
        let guide_probability = guide_umis
            .iter()
            .map(|&n| (n as f64 + cfg.alpha) / denom)
            .collect();

        Ok(Self {
            guide_probability,
            guide_umis,
            total_umis,
            background_droplets: data.n_cells(),
        })
    }

    pub fn p(&self, guide_id: u32) -> f64 {
        self.guide_probability[guide_id as usize]
    }
}
