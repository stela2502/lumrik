//! Framework-independent regulatory graph and single-cell observation contracts.
//!
//! TF expression is not required: regulatory activities will be latent variables.
//! Training backends (e.g. Burn) can be introduced without changing graph IDs.

use std::collections::HashSet;

/// Stable index into `RegulatoryGraph::nodes`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct NodeId(pub usize);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NodeKind {
    TranscriptionFactor,
    ProteinComplex,
    BindingSite,
    Enhancer,
    Promoter,
    Transcript,
    Gene,
}

#[derive(Clone, Debug)]
pub struct Node {
    pub id: NodeId,
    pub kind: NodeKind,
    /// Source identifier (e.g. Ensembl, ENCODE, FANTOM); not assumed unique
    /// across distinct node kinds or sources.
    pub external_id: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Relation {
    InteractsWith,
    BindsAt,
    AssociatedWith,
    Regulates,
    Initiates,
    Encodes,
    BelongsToGene,
}

#[derive(Clone, Debug)]
pub struct Edge {
    pub from: NodeId,
    pub to: NodeId,
    pub relation: Relation,
    /// Evidence weight is nonnegative; it is NOT a learned effect size.
    pub evidence_weight: f32,
    pub evidence_source: String,
}

#[derive(Clone, Debug, Default)]
pub struct RegulatoryGraph {
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
}

#[derive(Debug, PartialEq, Eq)]
pub enum ValidationError {
    InvalidNodeId,
    DuplicateNodeId,
    InvalidEdgeEndpoint,
    InvalidEvidenceWeight,
}

impl RegulatoryGraph {
    pub fn validate(&self) -> Result<(), ValidationError> {
        let mut ids = HashSet::new();
        for (position, node) in self.nodes.iter().enumerate() {
            if node.id.0 != position {
                return Err(ValidationError::InvalidNodeId);
            }
            if !ids.insert(node.id) {
                return Err(ValidationError::DuplicateNodeId);
            }
        }
        for edge in &self.edges {
            if edge.from.0 >= self.nodes.len() || edge.to.0 >= self.nodes.len() {
                return Err(ValidationError::InvalidEdgeEndpoint);
            }
            if !edge.evidence_weight.is_finite() || edge.evidence_weight < 0.0 {
                return Err(ValidationError::InvalidEvidenceWeight);
            }
        }
        Ok(())
    }
}

/// Sparse observations use explicit missingness: absence from the list is NOT
/// interpreted as a biological zero. A recorded zero is an observed zero.
#[derive(Clone, Debug)]
pub struct Observation {
    pub cell_index: usize,
    pub gene: NodeId,
    pub count: u32,
}

#[derive(Clone, Debug, Default)]
pub struct CellObservations {
    pub cell_ids: Vec<String>,
    pub observations: Vec<Observation>,
}

impl CellObservations {
    pub fn validate(&self, graph: &RegulatoryGraph) -> Result<(), ValidationError> {
        for observation in &self.observations {
            if observation.cell_index >= self.cell_ids.len()
                || graph.nodes.get(observation.gene.0).map(|n| n.kind) != Some(NodeKind::Gene)
            {
                return Err(ValidationError::InvalidEdgeEndpoint);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accepts_empty_graph() {
        assert_eq!(RegulatoryGraph::default().validate(), Ok(()));
    }

    #[test]
    fn rejects_dangling_edges() {
        let graph = RegulatoryGraph {
            nodes: vec![Node {
                id: NodeId(0),
                kind: NodeKind::Gene,
                external_id: "G".into(),
            }],
            edges: vec![Edge {
                from: NodeId(1),
                to: NodeId(0),
                relation: Relation::Regulates,
                evidence_weight: 1.0,
                evidence_source: "example".into(),
            }],
        };
        assert_eq!(graph.validate(), Err(ValidationError::InvalidEdgeEndpoint));
    }

    #[test]
    fn distinguishes_missing_from_observed_zero() {
        let graph = RegulatoryGraph {
            nodes: vec![Node {
                id: NodeId(0),
                kind: NodeKind::Gene,
                external_id: "G".into(),
            }],
            edges: vec![],
        };
        let data = CellObservations {
            cell_ids: vec!["cell-1".into()],
            observations: vec![Observation {
                cell_index: 0,
                gene: NodeId(0),
                count: 0,
            }],
        };
        assert_eq!(data.validate(&graph), Ok(()));
        assert_eq!(data.observations.len(), 1);
    }
}
