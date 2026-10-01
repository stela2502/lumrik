impl Ommverse {
    /// Fetch an arbitrary genomic interval from the assembly twoBit sequence truth.
    pub fn genomic_sequence(&self, chromosome: &str, start: u32, end: u32) -> Result<IntToDna> {
        let mut genome = TwoBitReader::open(&self.genome_twobit)?;
        Ok(genome.sequence(chromosome, start, end)?)
    }

    /// Materialize a gene's genomic span in transcription orientation.
    pub fn gene_sequence(&self, name_or_accession: &str) -> Result<IntToDna> {
        let gene = self
            .gene(name_or_accession)
            .with_context(|| format!("unknown gene {name_or_accession}"))?;
        let mut start = u32::MAX;
        let mut end = 0u32;
        let mut chr_id = None;
        let mut strand = Strand::Unknown;
        for &tx_id in gene.transcript_ids() {
            let tx = &self.splice.transcripts[tx_id];
            chr_id.get_or_insert(tx.chr_id);
            start = start.min(self.splice.tx_span_start[tx_id]);
            end = end.max(self.splice.tx_span_end[tx_id]);
            if strand == Strand::Unknown {
                strand = tx.strand;
            }
        }
        let chr_id = chr_id.context("gene has no transcripts")?;
        let chromosome = self.splice.chr_names.get(chr_id).context("gene chromosome missing")?;
        let mut genome = TwoBitReader::open(&self.genome_twobit)?;
        let dna = genome.sequence(chromosome, start, end)?;
        if strand != Strand::Minus {
            return Ok(dna);
        }
        let text = dna.to_string((end - start) as usize);
        let reverse = text.bytes().rev().map(complement).collect::<Vec<_>>();
        IntToDna::try_new(reverse).map_err(anyhow::Error::msg)
    }

    /// Materialize a spliced transcript directly from the assembly twoBit source.
    pub fn transcript_sequence(&self, name_or_accession: &str) -> Result<IntToDna> {
        let tx = self
            .transcript(name_or_accession)
            .with_context(|| format!("unknown transcript {name_or_accession}"))?;
        let mut genome = TwoBitReader::open(&self.genome_twobit)?;
        self.transcript_sequence_with_reader(tx, &mut genome)
    }

    fn transcript_sequence_with_reader(
        &self,
        tx: &Transcript,
        genome: &mut TwoBitReader,
    ) -> Result<IntToDna> {
        let chr = self.splice.chr_names.get(tx.chr_id).context("transcript chromosome missing")?;
        let mut cdna = Vec::<u8>::with_capacity(tx.transcript_len());
        match tx.strand {
            Strand::Plus | Strand::Unknown => {
                for exon in tx.exons() {
                    let dna = genome.sequence(chr, exon.start, exon.end)?;
                    cdna.extend_from_slice(dna.to_string(exon.len() as usize).as_bytes());
                }
            }
            Strand::Minus => {
                for exon in tx.exons().iter().rev() {
                    let dna = genome.sequence(chr, exon.start, exon.end)?;
                    let text = dna.to_string(exon.len() as usize);
                    cdna.extend(text.bytes().rev().map(complement));
                }
            }
        }
        IntToDna::try_new(cdna).map_err(anyhow::Error::msg)
    }

    /// Materialize only the annotated CDS of a transcript from the global twoBit source.
    pub fn transcript_coding_sequence(&self, name_or_accession: &str) -> Result<Option<IntToDna>> {
        let tx = self
            .transcript(name_or_accession)
            .with_context(|| format!("unknown transcript {name_or_accession}"))?;
        let Some((start, end)) = tx.cds_transcript_span() else {
            return Ok(None);
        };
        let mut genome = TwoBitReader::open(&self.genome_twobit)?;
        let cdna = self.transcript_sequence_with_reader(tx, &mut genome)?;
        let text = cdna.to_string(cdna.size);
        Ok(Some(IntToDna::try_new(&text.as_bytes()[start..end]).map_err(anyhow::Error::msg)?))
    }

    /// Reconstruct the protein sequence from the UCSC twoBit genome and the
    /// linked GTF transcript. Protein sequence is deliberately not persisted.
    pub fn protein_sequence(&self, accession: &str) -> Result<IntToProt> {
        let protein = self
            .protein(accession)
            .with_context(|| format!("unknown protein {accession}"))?;
        let mut genome = TwoBitReader::open(&self.genome_twobit)?;
        self.protein_sequence_with_reader(protein, &mut genome)
    }

    fn protein_sequence_with_reader(
        &self,
        protein: &Protein,
        genome: &mut TwoBitReader,
    ) -> Result<IntToProt> {
        let &tx_id = protein
            .transcript_ids
            .first()
            .context("protein has no linked transcript")?;
        let tx = &self.splice.transcripts[tx_id];
        let cdna = self.transcript_sequence_with_reader(tx, genome)?;
        let (start, end) = tx
            .cds_transcript_span()
            .context("transcript has no usable CDS")?;
        let text = cdna.to_string(cdna.size);
        if end > text.len() {
            bail!(
                "CDS {}..{} exceeds reconstructed cDNA length {}",
                start,
                end,
                text.len()
            );
        }
        Ok(IntToDna::try_new(&text.as_bytes()[start..end])
            .map_err(anyhow::Error::msg)?
            .translate())
    }

}
