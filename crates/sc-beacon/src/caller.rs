use std::collections::HashMap;

use crate::model::FittedModel;
use crate::stats::{benjamini_hochberg, poisson_upper_tail};
use crate::stats::PosteriorEvidence;

#[derive(Debug, Clone)]
pub struct CallConfig {
    pub minimum_posterior: f64,
    pub maximum_fdr: f64,
}

impl Default for CallConfig {
    fn default() -> Self {
        Self {
            minimum_posterior: 0.95,
            maximum_fdr: 0.05,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GuideCall {
    pub cell_id: u64,
    pub guide_id: u32,
    pub count: u32,
    pub lambda_cell: f64,
    pub ambient_probability: f64,
    pub expected_ambient: f64,
    pub posterior: PosteriorEvidence,
    pub ambient_p_value: f64,
    pub q_value: f64,
    pub called: bool,
}

#[derive(Debug, Clone)]
pub struct GuideCalls {
    /// Multi-guide by construction: every cell owns zero or more independent calls.
    pub by_cell: HashMap<u64, Vec<GuideCall>>,
    pub flat: Vec<GuideCall>,
}

impl GuideCalls {
        pub fn from_model(
        model: &FittedModel,
        cfg: &CallConfig,
    ) -> Self {
        let pvalues: Vec<f64> = model
            .observations
            .iter()
            .map(|state| {
                let obs = state.observation;

                let lambda =
                    model.lambda_by_cell[&obs.cell_id];

                let expected =
                    lambda * model.ambient.p(obs.guide_id);

                poisson_upper_tail(
                    obs.count,
                    expected,
                )
            })
            .collect();

        let qvalues =
            benjamini_hochberg(&pvalues);

        let mut flat =
            Vec::with_capacity(
                model.observations.len()
            );

        for ((state, p_value), q_value) in model
            .observations
            .iter()
            .zip(pvalues)
            .zip(qvalues)
        {
            let obs =
                state.observation;

            let lambda =
                model.lambda_by_cell[&obs.cell_id];

            let pg =
                model.ambient.p(obs.guide_id);

            let expected =
                lambda * pg;

            let posterior =
                state.evidence;

            let called =
                posterior.probability >= cfg.minimum_posterior
                    && q_value <= cfg.maximum_fdr;

            flat.push(
                GuideCall {
                    cell_id: obs.cell_id,
                    guide_id: obs.guide_id,
                    count: obs.count,

                    lambda_cell: lambda,
                    ambient_probability: pg,
                    expected_ambient: expected,

                    posterior,

                    ambient_p_value: p_value,
                    q_value,
                    called,
                }
            );
        }

        let mut by_cell:
            HashMap<u64, Vec<GuideCall>> =
            HashMap::new();

        for call in &flat {
            by_cell
                .entry(call.cell_id)
                .or_default()
                .push(call.clone());
        }

        for calls in by_cell.values_mut() {
            calls.sort_unstable_by_key(
                |call| call.guide_id
            );
        }

        Self {
            by_cell,
            flat,
        }
    }

    pub fn called_for_cell(&self, cell_id: u64) -> impl Iterator<Item = &GuideCall> {
        self.by_cell
            .get(&cell_id)
            .into_iter()
            .flat_map(|x| x.iter())
            .filter(|x| x.called)
    }
}
