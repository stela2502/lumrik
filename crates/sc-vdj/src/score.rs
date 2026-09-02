use crate::gex::ExpressionMatrix;

#[derive(Debug, Clone)]
pub struct MarkerContribution {
    pub gene: &'static str,
    pub expression: f64,
    pub activity: f64,
}

/// Deliberately narrow evidence for ongoing antigen-receptor recombination.
///
/// This does not attempt to assign a cell type or developmental stage. RAG1 and
/// RAG2 are the direct signal; DNTT is retained as supporting junctional-diversity
/// evidence. `rag_pair_detected` is intentionally literal: both RAG genes have
/// non-zero expression in the supplied expression matrix.
#[derive(Debug, Clone)]
pub struct RecombinationActivityEvidence {
    pub rag1_expression: f64,
    pub rag2_expression: f64,
    pub dntt_expression: f64,
    pub rag_activity: f64,
    pub dntt_activity: f64,
    pub rag_pair_detected: bool,
    pub contributions: Vec<MarkerContribution>,
}

pub fn score_recombination_activity<E: ExpressionMatrix>(
    gex: &E,
    cell: &str,
) -> RecombinationActivityEvidence {
    let rag1_expression = gex.expression(cell, "RAG1");
    let rag2_expression = gex.expression(cell, "RAG2");
    let dntt_expression = gex.expression(cell, "DNTT");
    let rag1_activity = squash(rag1_expression);
    let rag2_activity = squash(rag2_expression);
    let dntt_activity = squash(dntt_expression);
    let rag_activity = (rag1_activity * rag2_activity).sqrt();
    let contributions = [
        ("RAG1", rag1_expression, rag1_activity),
        ("RAG2", rag2_expression, rag2_activity),
        ("DNTT", dntt_expression, dntt_activity),
    ]
    .into_iter()
    .map(|(gene, expression, activity)| MarkerContribution {
        gene,
        expression,
        activity,
    })
    .collect();

    RecombinationActivityEvidence {
        rag1_expression,
        rag2_expression,
        dntt_expression,
        rag_activity,
        dntt_activity,
        rag_pair_detected: rag1_expression > 0.0 && rag2_expression > 0.0,
        contributions,
    }
}

fn squash(x: f64) -> f64 {
    if x <= 0.0 {
        0.0
    } else {
        1.0 - (-x / 2.0).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    struct TestExpression(HashMap<(String, String), f64>);
    impl ExpressionMatrix for TestExpression {
        fn expression(&self, cell: &str, gene: &str) -> f64 {
            *self.0.get(&(cell.to_string(), gene.to_string())).unwrap_or(&0.0)
        }
    }

    #[test]
    fn rag_activity_requires_both_rag_genes_for_pair_detection() {
        let gex = TestExpression(HashMap::from([
            (("c".into(), "RAG1".into()), 4.0),
            (("c".into(), "RAG2".into()), 0.0),
            (("c".into(), "DNTT".into()), 3.0),
        ]));
        let evidence = score_recombination_activity(&gex, "c");
        assert!(!evidence.rag_pair_detected);
        assert_eq!(evidence.rag_activity, 0.0);
        assert!(evidence.dntt_activity > 0.0);
    }

    #[test]
    fn paired_rag_expression_produces_activity() {
        let gex = TestExpression(HashMap::from([
            (("c".into(), "RAG1".into()), 2.0),
            (("c".into(), "RAG2".into()), 2.0),
        ]));
        let evidence = score_recombination_activity(&gex, "c");
        assert!(evidence.rag_pair_detected);
        assert!(evidence.rag_activity > 0.0);
    }
}
