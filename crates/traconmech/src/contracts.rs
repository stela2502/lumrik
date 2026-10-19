//! Declarative contract registry. No Ollama response is trusted or executed here.
use anyhow::{Result, ensure};
use serde_json::{Value, json};
use std::collections::HashSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Availability {
    Available,
    NotConfigured,
    Planned,
}
impl Availability {
    fn as_str(self) -> &'static str {
        match self {
            Self::Available => "available",
            Self::NotConfigured => "not_configured",
            Self::Planned => "planned",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Entry {
    pub id: &'static str,
    pub title: &'static str,
    pub description: &'static str,
    pub status: Availability,
}
impl Entry {
    fn json(self) -> Value {
        json!({"id":self.id,"title":self.title,"description":self.description,"status":self.status.as_str()})
    }
}

pub const RESULTS: &[Entry] = &[
    Entry {
        id: "column_order",
        title: "Ordered categories",
        description: "Propose an order for existing observation labels without changing their identifiers.",
        status: Availability::Planned,
    },
    Entry {
        id: "annotation_column",
        title: "Observation annotation",
        description: "Propose cell-type or stage metadata, subject to validation and approval.",
        status: Availability::Planned,
    },
    Entry {
        id: "analysis_target",
        title: "Analysis target",
        description: "Propose a population, gene or comparison for a registered analysis.",
        status: Availability::Planned,
    },
    Entry {
        id: "gene_set",
        title: "Gene set",
        description: "Propose a set of genes with provenance and resolvable identifiers.",
        status: Availability::Planned,
    },
    Entry {
        id: "developmental_graph",
        title: "Developmental graph",
        description: "Propose directed and possibly branching population relationships.",
        status: Availability::Planned,
    },
];
pub const CALLBACKS: &[Entry] = &[
    Entry {
        id: "propose_column_order",
        title: "Review column order",
        description: "Validate existing levels and request approval before reordering.",
        status: Availability::Planned,
    },
    Entry {
        id: "propose_annotation",
        title: "Review annotation",
        description: "Validate cell or population mappings and request approval before writing.",
        status: Availability::Planned,
    },
    Entry {
        id: "propose_analysis",
        title: "Review analysis",
        description: "Validate requested target and require approval before running.",
        status: Availability::Planned,
    },
    Entry {
        id: "propose_gene_set",
        title: "Review gene set",
        description: "Resolve genes and require approval before storing.",
        status: Availability::Planned,
    },
    Entry {
        id: "propose_developmental_graph",
        title: "Review developmental graph",
        description: "Validate edges and request approval before accepting.",
        status: Availability::Planned,
    },
];
pub const PROVIDERS: &[Entry] = &[
    Entry {
        id: "cluster_signatures",
        title: "Cluster signatures",
        description: "Read measured sc-analysis marker summary and top-20 genes.",
        status: Availability::Available,
    },
    Entry {
        id: "marker_expression",
        title: "Marker expression",
        description: "Compute mean, SD and detection fraction per observed population.",
        status: Availability::Available,
    },
    Entry {
        id: "ommverse",
        title: "Ommverse reference",
        description: "Genome and protein reference; provider API integration pending.",
        status: Availability::Planned,
    },
    Entry {
        id: "string",
        title: "STRING associations",
        description: "Protein association evidence with source scores and provenance.",
        status: Availability::Planned,
    },
    Entry {
        id: "gsea",
        title: "GSEA / MSigDB",
        description: "Gene-set and enrichment evidence.",
        status: Availability::Planned,
    },
    Entry {
        id: "encode",
        title: "ENCODE",
        description: "Regulatory element and TF-binding evidence.",
        status: Availability::Planned,
    },
    Entry {
        id: "interpro",
        title: "InterPro",
        description: "Protein domains via Ommverse reference queries.",
        status: Availability::Planned,
    },
    Entry {
        id: "variants",
        title: "Variant evidence",
        description: "Variant, transcript and protein consequence evidence.",
        status: Availability::Planned,
    },
    Entry {
        id: "images",
        title: "Image evidence",
        description: "Image metadata, references and derived visual evidence.",
        status: Availability::Planned,
    },
];

pub struct Template {
    pub entry: Entry,
    pub inputs: &'static [&'static str],
    pub outputs: &'static [&'static str],
    pub callbacks: &'static [&'static str],
    pub prompt: &'static str,
}
pub const TEMPLATES: &[Template] = &[
    Template {
        entry: Entry {
            id: "identify_cell_types",
            title: "Identify cell types",
            description: "Assign tentative identities from measured cluster signatures.",
            status: Availability::Planned,
        },
        inputs: &["cluster_signatures"],
        outputs: &["annotation_column"],
        callbacks: &["propose_annotation"],
        prompt: "Use attached cluster signatures to propose cell identities. Cite supporting measured genes; distinguish hypotheses from observations.",
    },
    Template {
        entry: Entry {
            id: "developmental_order",
            title: "Developmental order",
            description: "Propose ordered stages and branching relationships between real population IDs.",
            status: Availability::Planned,
        },
        inputs: &["cluster_signatures"],
        outputs: &["column_order", "developmental_graph"],
        callbacks: &["propose_column_order", "propose_developmental_graph"],
        prompt: "Propose developmental ordering and branching using actual population IDs and measured marker evidence; label uncertain edges.",
    },
    Template {
        entry: Entry {
            id: "investigate_dba_tp53",
            title: "Investigate DBA and TP53",
            description: "Prioritize erythroid populations for a measured DBA/TP53 investigation.",
            status: Availability::Planned,
        },
        inputs: &["cluster_signatures"],
        outputs: &["analysis_target"],
        callbacks: &["propose_analysis"],
        prompt: "Given approved developmental annotations, propose a population for DBA and TP53 investigation. Explain required case/control and sample-level evidence; do not invent measurements.",
    },
];

pub fn validate() -> Result<()> {
    let mut ids = HashSet::new();
    for (kind, entries) in [
        ("result", RESULTS),
        ("callback", CALLBACKS),
        ("provider", PROVIDERS),
    ] {
        ids.clear();
        for e in entries {
            ensure!(
                !e.id.is_empty() && ids.insert(e.id),
                "duplicate/empty {kind} id: {}",
                e.id
            );
        }
    }
    ids.clear();
    for t in TEMPLATES {
        ensure!(ids.insert(t.entry.id), "duplicate template: {}", t.entry.id);
        ensure!(
            !t.outputs.is_empty() && !t.prompt.is_empty(),
            "incomplete template: {}",
            t.entry.id
        );
        for id in t.inputs {
            ensure!(
                PROVIDERS.iter().any(|p| p.id == *id),
                "unknown provider {id}"
            );
        }
        for id in t.outputs {
            ensure!(RESULTS.iter().any(|p| p.id == *id), "unknown result {id}");
        }
        for id in t.callbacks {
            ensure!(
                CALLBACKS.iter().any(|p| p.id == *id),
                "unknown callback {id}"
            );
        }
    }
    Ok(())
}

pub fn registry_json() -> Result<Value> {
    validate()?;
    Ok(
        json!({"schema_version":1,"result_types":RESULTS.iter().map(|x|x.json()).collect::<Vec<_>>(),
        "callbacks":CALLBACKS.iter().map(|x|x.json()).collect::<Vec<_>>(),
        "evidence_providers":PROVIDERS.iter().map(|x|x.json()).collect::<Vec<_>>(),
        "templates":TEMPLATES.iter().map(|t|json!({"id":t.entry.id,"title":t.entry.title,
            "description":t.entry.description,"status":t.entry.status.as_str(),
            "inputs":t.inputs,"returns":t.outputs,"callbacks":t.callbacks,"prompt":t.prompt})).collect::<Vec<_>>() }),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn contracts_are_valid() {
        validate().unwrap();
    }
    #[test]
    fn registry_has_templates_and_providers() {
        let v = registry_json().unwrap();
        assert!(v["templates"].as_array().unwrap().len() >= 3);
        assert_eq!(v["schema_version"], 1);
    }
}
