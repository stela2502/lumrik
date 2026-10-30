//! Cell-cycle classification using expression-matched module scores.
use crate::{module_score, SingleCellData};
use anyhow::{ensure, Result};
use std::collections::BTreeMap;

// Regev/Tirosh cell-cycle genes, in Scanpy's standard order.
const S: &str = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8";
const G2M: &str = "HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA";

pub struct CellCycleScores {
    pub s: Vec<f64>,
    pub g2m: Vec<f64>,
    pub phase: Vec<String>,
    pub s_matched: usize,
    pub g2m_matched: usize,
}

pub fn score(data: &SingleCellData) -> Result<CellCycleScores> {
    let lookup: std::collections::HashSet<&str> = data.features.iter().map(String::as_str).collect();
    let s: Vec<_> = S.split_whitespace().filter(|g| lookup.contains(*g)).collect();
    let g2m: Vec<_> = G2M.split_whitespace().filter(|g| lookup.contains(*g)).collect();
    ensure!(!s.is_empty() && !g2m.is_empty(), "insufficient cell-cycle genes: S={}/{}, G2M={}/{}", s.len(), 43, g2m.len(), 54);
    let mut modules = BTreeMap::new();
    for gene in &s { modules.insert((*gene).to_string(), 1); }
    for gene in &g2m { modules.insert((*gene).to_string(), 2); }
    let result = module_score::score(data, &modules)?;
    let s_index = result.ids.iter().position(|id| *id == 1).expect("S module");
    let g_index = result.ids.iter().position(|id| *id == 2).expect("G2M module");
    let mut s_scores = Vec::with_capacity(result.scores.len());
    let mut g_scores = Vec::with_capacity(result.scores.len());
    let mut phase = Vec::with_capacity(result.scores.len());
    for row in &result.scores {
        let (a,b) = (row[s_index], row[g_index]);
        s_scores.push(a); g_scores.push(b);
        phase.push(if a < 0.0 && b < 0.0 { "G1" } else if a > b { "S" } else { "G2M" }.to_string());
    }
    Ok(CellCycleScores { s: s_scores, g2m: g_scores, phase, s_matched: s.len(), g2m_matched: g2m.len() })
}

#[cfg(test)]
mod tests {
    #[test]
    fn marker_counts() {
        assert_eq!(super::S.split_whitespace().count(), 43);
        assert_eq!(super::G2M.split_whitespace().count(), 54);
    }
}
