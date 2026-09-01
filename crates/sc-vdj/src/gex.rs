use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub trait ExpressionMatrix {
    fn expression(&self, cell: &str, gene: &str) -> f64;
}

/// Minimal adapter useful for testing/integration: tab-separated cell, gene, value.
/// Nelrune can implement ExpressionMatrix directly over its native expression object
/// without converting the matrix.
#[derive(Debug, Clone, Default)]
pub struct LongTsvExpression {
    values: HashMap<(String, String), f64>,
}

impl LongTsvExpression {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let r = BufReader::new(
            File::open(path)
                .with_context(|| format!("opening expression TSV {}", path.display()))?,
        );
        let mut out = Self::default();
        for (i, line) in r.lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 3 {
                continue;
            }
            let Ok(v) = f[2].parse::<f64>() else {
                if i == 0 {
                    continue;
                } else {
                    continue;
                }
            };
            out.values
                .insert((f[0].to_string(), f[1].to_ascii_uppercase()), v);
        }
        Ok(out)
    }
}
impl ExpressionMatrix for LongTsvExpression {
    fn expression(&self, cell: &str, gene: &str) -> f64 {
        *self
            .values
            .get(&(cell.to_string(), gene.to_ascii_uppercase()))
            .unwrap_or(&0.0)
    }
}
