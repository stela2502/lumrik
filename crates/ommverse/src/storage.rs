impl Ommverse {
    pub fn save(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        let mut file = File::create(path)?;
        file.write_all(MAGIC)?;
        file.write_all(&OMMVERSE_FORMAT_VERSION.to_le_bytes())?;
        bincode::serialize_into(file, self)?;
        crate::schema::OmmverseSchema::for_ommverse(self).save_for_index(path)?;
        Ok(())
    }

    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let mut file = File::open(path.as_ref())?;
        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;
        if &magic != MAGIC {
            bail!("not an Ommverse index");
        }
        let mut version = [0u8; 4];
        file.read_exact(&mut version)?;
        let version = u32::from_le_bytes(version);
        let mut out = match version {
            OMMVERSE_FORMAT_VERSION => bincode::deserialize_from(file)?,
            9 => {
                let old: OmmverseV9 = bincode::deserialize_from(file)?;
                Self {
                    assembly: old.assembly, source_root: old.source_root, genome_twobit: old.genome_twobit,
                    splice: old.splice, proteins: old.proteins, protein_features: ProteinFeatureIndex::default(), report: old.report,
                    interpro_entries: old.interpro_entries, chromatin: old.chromatin,
                    protein_binding: old.protein_binding, ctcf: old.ctcf, experimental_loops: old.experimental_loops,
                    sources: old.sources, protein_by_accession: HashMap::new(), proteins_by_gene: HashMap::new(), proteins_by_transcript: HashMap::new(), gene_by_name: HashMap::new(),
                }
            }
            7 => {
                let old: OmmverseV7 = bincode::deserialize_from(file)?;
                Self {
                    assembly: old.assembly, source_root: old.source_root, genome_twobit: old.genome_twobit,
                    splice: old.splice, proteins: old.proteins, protein_features: ProteinFeatureIndex::default(), report: old.report,
                    interpro_entries: old.interpro_entries, chromatin: old.chromatin,
                    protein_binding: old.protein_binding, ctcf: old.ctcf, experimental_loops: old.experimental_loops,
                    sources: None, protein_by_accession: HashMap::new(), proteins_by_gene: HashMap::new(), proteins_by_transcript: HashMap::new(), gene_by_name: HashMap::new(),
                }
            }
            6 => {
                let old: OmmverseV6 = bincode::deserialize_from(file)?;
                Self {
                    assembly: old.assembly,
                    source_root: old.source_root,
                    genome_twobit: old.genome_twobit,
                    splice: old.splice,
                    proteins: old.proteins,
                    protein_features: ProteinFeatureIndex::default(),
                    report: old.report,
                    interpro_entries: old.interpro_entries,
                    chromatin: old.chromatin,
                    protein_binding: old.protein_binding,
                    ctcf: old.ctcf,
                    experimental_loops: ExperimentalLoopArchitecture::default(),
                    sources: None,
                    protein_by_accession: HashMap::new(),
                    proteins_by_gene: HashMap::new(),
                    proteins_by_transcript: HashMap::new(),
                    gene_by_name: HashMap::new(),
                }
            }
            5 => {
                let old: OmmverseV5 = bincode::deserialize_from(file)?;
                Self {
                    assembly: old.assembly,
                    source_root: old.source_root,
                    genome_twobit: old.genome_twobit,
                    splice: old.splice,
                    proteins: old.proteins,
                    protein_features: ProteinFeatureIndex::default(),
                    report: old.report,
                    interpro_entries: old.interpro_entries,
                    chromatin: old.chromatin,
                    protein_binding: old.protein_binding,
                    ctcf: CtcfArchitecture::default(),
                    experimental_loops: ExperimentalLoopArchitecture::default(),
                    sources: None,
                    protein_by_accession: HashMap::new(),
                    proteins_by_gene: HashMap::new(),
                    proteins_by_transcript: HashMap::new(),
                    gene_by_name: HashMap::new(),
                }
            }
            4 => {
                let old: OmmverseV4 = bincode::deserialize_from(file)?;
                // v4 embedded all rPeaks.  Drop that experiment-level payload on
                // load; rebuilding v5 reconstructs the compact union from cache.
                Self {
                    assembly: old.assembly,
                    source_root: old.source_root,
                    genome_twobit: old.genome_twobit,
                    splice: old.splice,
                    proteins: old.proteins,
                    protein_features: ProteinFeatureIndex::default(),
                    report: old.report,
                    interpro_entries: old.interpro_entries,
                    chromatin: old.chromatin,
                    protein_binding: ProteinBindingUnion::default(),
                    ctcf: CtcfArchitecture::default(),
                    experimental_loops: ExperimentalLoopArchitecture::default(),
                    sources: None,
                    protein_by_accession: HashMap::new(),
                    proteins_by_gene: HashMap::new(),
                    proteins_by_transcript: HashMap::new(),
                    gene_by_name: HashMap::new(),
                }
            }
            3 => {
                let old: OmmverseV3 = bincode::deserialize_from(file)?;
                Self {
                    assembly: old.assembly,
                    source_root: old.source_root,
                    genome_twobit: old.genome_twobit,
                    splice: old.splice,
                    proteins: old.proteins,
                    protein_features: ProteinFeatureIndex::default(),
                    report: old.report,
                    interpro_entries: old.interpro_entries,
                    chromatin: old.chromatin,
                    protein_binding: ProteinBindingUnion::default(),
                    ctcf: CtcfArchitecture::default(),
                    experimental_loops: ExperimentalLoopArchitecture::default(),
                    sources: None,
                    protein_by_accession: HashMap::new(),
                    proteins_by_gene: HashMap::new(),
                    proteins_by_transcript: HashMap::new(),
                    gene_by_name: HashMap::new(),
                }
            }
            2 => {
                let old: OmmverseV2 = bincode::deserialize_from(file)?;
                Self {
                    assembly: old.assembly,
                    source_root: old.source_root,
                    genome_twobit: old.genome_twobit,
                    splice: old.splice,
                    proteins: old.proteins,
                    protein_features: ProteinFeatureIndex::default(),
                    report: old.report,
                    interpro_entries: old.interpro_entries,
                    chromatin: Vec::new(),
                    protein_binding: ProteinBindingUnion::default(),
                    ctcf: CtcfArchitecture::default(),
                    experimental_loops: ExperimentalLoopArchitecture::default(),
                    sources: None,
                    protein_by_accession: HashMap::new(),
                    proteins_by_gene: HashMap::new(),
                    proteins_by_transcript: HashMap::new(),
                    gene_by_name: HashMap::new(),
                }
            }
            1 => {
                let old: OmmverseV1 = bincode::deserialize_from(file)?;
                Self {
                    assembly: old.assembly,
                    source_root: old.source_root,
                    genome_twobit: old.genome_twobit,
                    splice: old.splice,
                    proteins: old.proteins,
                    protein_features: ProteinFeatureIndex::default(),
                    report: old.report,
                    interpro_entries: HashMap::new(),
                    chromatin: Vec::new(),
                    protein_binding: ProteinBindingUnion::default(),
                    ctcf: CtcfArchitecture::default(),
                    experimental_loops: ExperimentalLoopArchitecture::default(),
                    sources: None,
                    protein_by_accession: HashMap::new(),
                    proteins_by_gene: HashMap::new(),
                    proteins_by_transcript: HashMap::new(),
                    gene_by_name: HashMap::new(),
                }
            }
            _ => bail!("unsupported Ommverse format version {version}"),
        };
        out.reindex();
        Ok(out)
    }

    fn reindex(&mut self) {
        self.protein_by_accession.clear();
        self.proteins_by_gene.clear();
        self.proteins_by_transcript.clear();
        self.gene_by_name.clear();
        for gene in &self.splice.genes {
            for name in &gene.names {
                let key = name.trim().to_ascii_lowercase();
                if !key.is_empty() {
                    self.gene_by_name.entry(key).or_insert(gene.id);
                }
            }
        }
        for (i, p) in self.proteins.iter().enumerate() {
            for identifier in &p.identifiers {
                self.protein_by_accession.entry(identifier.clone()).or_insert(i);
                self.protein_by_accession.entry(strip_version(identifier).to_owned()).or_insert(i);
            }
            self.protein_by_accession.entry(p.accession.clone()).or_insert(i);
            for &tx_id in &p.transcript_ids {
                self.proteins_by_transcript.entry(tx_id).or_default().push(i);
            }
            if !p.gene_symbol.is_empty() {
                self.proteins_by_gene
                    .entry(p.gene_symbol.to_ascii_lowercase())
                    .or_default()
                    .push(i);
            }
        }
    }
}

