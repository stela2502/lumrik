impl Ommverse {
    pub fn protein_id(&self, name_or_accession: &str) -> Option<usize> {
        self.protein_by_accession
            .get(name_or_accession)
            .or_else(|| self.protein_by_accession.get(strip_version(name_or_accession)))
            .copied()
    }

    pub fn protein(&self, name_or_accession: &str) -> Option<&Protein> {
        self.protein_id(name_or_accession).map(|i| &self.proteins[i])
    }

    /// Resolve a gene by any retained symbol, alias or stable accession.
    pub fn gene(&self, name_or_accession: &str) -> Option<&Gene> {
        self.gene_by_name
            .get(&name_or_accession.trim().to_ascii_lowercase())
            .and_then(|&id| self.splice.genes.get(id))
    }

    /// Resolve a transcript using the splice index's retained transcript names/accessions.
    pub fn transcript(&self, name_or_accession: &str) -> Option<&Transcript> {
        self.splice.transcript_by_name(name_or_accession).ok()
    }

    /// Follow a typed protein -> transcript -> gene relationship.
    pub fn gene_for_protein(&self, accession: &str) -> Option<&Gene> {
        let protein = self.protein(accession)?;
        let tx_id = *protein.transcript_ids.first()?;
        let tx = self.splice.transcripts.get(tx_id)?;
        self.splice.genes.get(tx.gene_id)
    }

    pub fn protein_ids_for_gene(&self, name_or_accession: &str) -> Vec<usize> {
        let mut protein_ids = self
            .gene(name_or_accession)
            .into_iter()
            .flat_map(|gene| gene.transcript_ids().iter().copied())
            .flat_map(|tx_id| self.proteins_by_transcript.get(&tx_id).into_iter().flatten().copied())
            .collect::<Vec<_>>();
        protein_ids.sort_unstable();
        protein_ids.dedup();
        protein_ids
    }

    pub fn proteins_for_gene(&self, name_or_accession: &str) -> impl Iterator<Item = &Protein> {
        self.protein_ids_for_gene(name_or_accession)
            .into_iter()
            .map(|protein_id| &self.proteins[protein_id])
    }

    pub fn protein_ids_for_transcript(&self, name_or_accession: &str) -> Vec<usize> {
        self.transcript(name_or_accession)
            .and_then(|tx| self.proteins_by_transcript.get(&tx.id))
            .cloned()
            .unwrap_or_default()
    }

    /// Proteins explicitly linked to a retained transcript name/accession.
    /// This is intentionally transcript-specific: callers must not project a
    /// gene-level protein feature onto every isoform.
    pub fn proteins_for_transcript(&self, name_or_accession: &str) -> impl Iterator<Item = &Protein> {
        self.protein_ids_for_transcript(name_or_accession)
            .into_iter()
            .map(|protein_id| &self.proteins[protein_id])
    }

    /// Resolve every transcript explicitly linked to a protein identity.
    pub fn transcripts_for_protein(&self, accession: &str) -> impl Iterator<Item = &Transcript> {
        let transcript_ids = self
            .protein(accession)
            .map(|protein| protein.transcript_ids.clone())
            .unwrap_or_default();
        transcript_ids.into_iter().filter_map(|tx_id| self.splice.transcripts.get(tx_id))
    }

    pub fn proteins_with_feature(
        &self,
        kind: ProteinFeatureKind,
    ) -> impl Iterator<Item = &Protein> {
        self.proteins
            .iter()
            .filter(move |p| p.features.iter().any(|f| f.kind == kind))
    }

}
