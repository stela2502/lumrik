use std::collections::HashSet;
use std::path::Path;

use sc_beacon::{KneeCountFit, fit_knee_counts, write_count_qc};
use scdata::QuantData;

#[derive(Debug, Clone)]
pub struct BeaconCellCalling {
    pub fit: KneeCountFit,
    pub cells: Vec<(u64, u32, bool)>,
    pub retained: HashSet<u64>,
}

impl BeaconCellCalling {
    pub fn write_qc<P: AsRef<Path>>(&self, out_dir: P) -> Result<(), String> {
        let counts: Vec<u32> = self.cells.iter().map(|(_, count, _)| *count).collect();
        write_count_qc(&counts, &self.fit, out_dir).map_err(|e| e.to_string())
    }

    pub fn write_tsv<P: AsRef<Path>>(&self, path: P) -> Result<(), String> {
        use std::io::Write;

        let path = path.as_ref();
        let mut rows = self.cells.clone();
        rows.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        let mut out = std::fs::File::create(path)
            .map_err(|e| format!("creating {}: {e}", path.display()))?;
        writeln!(out, "rank\tcell_id\tumi_count\tcalled")
            .map_err(|e| format!("writing {}: {e}", path.display()))?;
        for (rank, (cell_id, count, called)) in rows.into_iter().enumerate() {
            writeln!(out, "{}\t{}\t{}\t{}", rank + 1, cell_id, count, called)
                .map_err(|e| format!("writing {}: {e}", path.display()))?;
        }
        Ok(())
    }
}

pub fn beacon_cell_calling(data: &QuantData) -> Result<BeaconCellCalling, String> {
    let cell_counts = data.cell_umi_counts(QuantData::EXONIC);
    let counts: Vec<u32> = cell_counts.iter().map(|(_, count)| *count).collect();
    let fit = fit_knee_counts(&counts).map_err(|e| e.to_string())?;
    let cells: Vec<(u64, u32, bool)> = cell_counts
        .into_iter()
        .map(|(cell_id, count)| (cell_id, count, count >= fit.umi_cutoff))
        .collect();
    let retained = cells
        .iter()
        .filter_map(|(cell_id, _, called)| called.then_some(*cell_id))
        .collect();
    Ok(BeaconCellCalling { fit, cells, retained })
}
