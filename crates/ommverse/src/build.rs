fn find_gtf(genes_dir: &Path) -> Result<PathBuf> {
    if !genes_dir.is_dir() {
        bail!(
            "required UCSC Genes directory missing: {}",
            genes_dir.display()
        );
    }

    let mut gtfs = std::fs::read_dir(genes_dir)
        .with_context(|| format!("reading Genes directory {}", genes_dir.display()))?
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .filter(|path| {
            let Some(name) = path.file_name().and_then(|name| name.to_str()) else {
                return false;
            };
            path.is_file() && (name.ends_with(".gtf") || name.ends_with(".gtf.gz"))
        })
        .collect::<Vec<_>>();
    // knownGene is transcript-centric: in its UCSC GTF export the `gene_id`
    // can effectively identify individual knownGene transcripts.  That is useful
    // for alignment, but disastrous for Ommverse gene identity (hg38 produced
    // ~236k one-name "genes").  Prefer annotations that preserve a real gene
    // hierarchy and human-readable gene names.
    gtfs.sort_by_key(|path| {
        let name = path.file_name().and_then(|x| x.to_str()).unwrap_or("");
        if name.contains("gencode") {
            0
        } else if name.contains("ncbiRefSeq") {
            1
        } else if name == "refGene.gtf.gz" || name == "refGene.gtf" {
            2
        } else if name.contains("knownGene") {
            3
        } else {
            4
        }
    });

    gtfs.into_iter().next().with_context(|| {
        format!(
            "no GTF annotation (*.gtf or *.gtf.gz) found in {}",
            genes_dir.display()
        )
    })
}

#[derive(Debug, Clone, Default)]
pub struct InterProImportReport {
    pub entries_loaded: usize,
    pub records_streamed: usize,
    pub matched_records: usize,
    pub features_added: usize,
    pub projected_features_added: usize,
    pub duplicate_features: usize,
    pub malformed_records: usize,
    pub unknown_entries: usize,
}

fn interpro_feature_kind(value: &str) -> ProteinFeatureKind {
    match value
        .to_ascii_lowercase()
        .replace(['-', '_', ' '], "")
        .as_str()
    {
        "domain" => ProteinFeatureKind::Domain,
        "family" => ProteinFeatureKind::ProteinFamily,
        "homologoussuperfamily" => ProteinFeatureKind::HomologousSuperfamily,
        "repeat" => ProteinFeatureKind::Repeat,
        "activesite" => ProteinFeatureKind::ActiveSite,
        "bindingsite" => ProteinFeatureKind::BindingSite,
        "conservedsite" => ProteinFeatureKind::ConservedSite,
        "ptm" | "ptmsite" => ProteinFeatureKind::PtmSite,
        _ => ProteinFeatureKind::Other,
    }
}

pub fn ingest_interpro_many<F>(
    ommverses: &mut [Ommverse],
    protein2ipr: &Path,
    entry_list: &Path,
    parent_child_tree: Option<&Path>,
    mut progress: F,
) -> Result<Vec<InterProImportReport>>
where
    F: FnMut(usize, &[InterProImportReport]),
{
    use flate2::read::MultiGzDecoder;
    use std::io::{BufRead, BufReader};

    let mut entries = HashMap::<String, (ProteinFeatureKind, String)>::new();
    let reader = BufReader::new(
        File::open(entry_list).with_context(|| format!("opening {}", entry_list.display()))?,
    );
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let mut f = line.split('\t');
        let Some(ipr) = f.next().map(str::trim).filter(|x| !x.is_empty()) else {
            continue;
        };
        let kind = interpro_feature_kind(f.next().unwrap_or("").trim());
        let name = f.next().unwrap_or("").trim().to_owned();
        entries.insert(ipr.to_owned(), (kind, name));
    }

    for omm in ommverses.iter_mut() {
        for (ipr, (kind, name)) in &entries {
            omm.interpro_entries
                .entry(ipr.clone())
                .or_insert_with(|| InterProEntry {
                    kind: *kind,
                    name: name.clone(),
                    parents: Vec::new(),
                    children: Vec::new(),
                });
        }
        if let Some(tree) = parent_child_tree {
            let reader = BufReader::new(
                File::open(tree).with_context(|| format!("opening {}", tree.display()))?,
            );
            let mut stack: Vec<String> = Vec::new();
            for line in reader.lines() {
                let line = line?;
                if line.trim().is_empty() || line.starts_with('#') {
                    continue;
                }
                let mut depth = 0usize;
                let bytes = line.as_bytes();
                while bytes.get(depth * 2..depth * 2 + 2) == Some(b"--") {
                    depth += 1;
                }
                let body = &line[depth * 2..];
                let Some(ipr) = body
                    .split("::")
                    .next()
                    .map(str::trim)
                    .filter(|x| x.starts_with("IPR"))
                else {
                    continue;
                };
                stack.truncate(depth);
                if depth > 0 {
                    if let Some(parent) = stack.get(depth - 1).cloned() {
                        if let Some(entry) = omm.interpro_entries.get_mut(ipr) {
                            if !entry.parents.contains(&parent) {
                                entry.parents.push(parent.clone());
                            }
                        }
                        if let Some(entry) = omm.interpro_entries.get_mut(&parent) {
                            if !entry.children.iter().any(|x| x == ipr) {
                                entry.children.push(ipr.to_owned());
                            }
                        }
                    }
                }
                stack.push(ipr.to_owned());
            }
        }
    }

    // InterPro is keyed by source protein accessions (normally UniProt), while
    // Ommverse proteins may have been born from GENCODE/ENSP and only acquired
    // that accession as an identifier during UCSC protein mapping.  Build this
    // import bridge from every protein identity, not just Protein::accession.
    // Keep it local to the importer: InterPro does not get to pollute or redefine
    // the central protein identity registry.
    let mut locations = HashMap::<String, Vec<(usize, usize)>>::new();
    for (oi, omm) in ommverses.iter().enumerate() {
        for (pi, protein) in omm.proteins.iter().enumerate() {
            let mut add = |identifier: &str| {
                if identifier.is_empty() {
                    return;
                }
                let targets = locations.entry(identifier.to_owned()).or_default();
                if !targets.contains(&(oi, pi)) {
                    targets.push((oi, pi));
                }
                let stable = strip_version(identifier);
                if stable != identifier {
                    let targets = locations.entry(stable.to_owned()).or_default();
                    if !targets.contains(&(oi, pi)) {
                        targets.push((oi, pi));
                    }
                }
            };
            add(&protein.accession);
            for identifier in &protein.identifiers {
                add(identifier);
            }
        }
    }
    // STRING's alias table is used only as an import-time identity bridge.
    // In particular, UniProt_AC rows connect protein2ipr accessions such as
    // O14770 to the GENCODE/ENSP proteins already present in Ommverse.  The
    // aliases are not inserted into the central identity registry: once the
    // InterPro interval has been projected to the genome, genomic placement is
    // the durable identity.
    for (oi, omm) in ommverses.iter().enumerate() {
        let interactions = omm.source_root.join("Interactions");
        let aliases = std::fs::read_dir(&interactions).ok().and_then(|entries| {
            let mut paths = entries
                .filter_map(|e| e.ok().map(|e| e.path()))
                .filter(|p| p.file_name().and_then(|n| n.to_str()).is_some_and(|n| {
                    n.starts_with("protein.aliases") && (n.ends_with(".txt") || n.ends_with(".txt.gz"))
                }))
                .collect::<Vec<_>>();
            paths.sort();
            paths.pop()
        });
        let Some(path) = aliases else { continue; };
        let file = File::open(&path).with_context(|| format!("opening STRING aliases {}", path.display()))?;
        let gz = path.extension().is_some_and(|x| x == "gz");
        let input: Box<dyn Read> = if gz {
            Box::new(MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };
        let reader = BufReader::with_capacity(1024 * 1024, input);
        for line in reader.lines() {
            let line = line?;
            let mut f = line.split('\t');
            let Some(node) = f.next() else { continue; };
            let Some(alias) = f.next() else { continue; };
            let source = f.next().unwrap_or("");
            if source != "UniProt_AC" {
                continue;
            }
            let ensp = node.split_once('.').map(|(_, id)| id).unwrap_or(node);
            let stable = strip_version(ensp);
            let Some(targets) = locations.get(stable).cloned().or_else(|| locations.get(ensp).cloned()) else {
                continue;
            };
            for target in targets.into_iter().filter(|(target_oi, _)| *target_oi == oi) {
                let out = locations.entry(alias.to_owned()).or_default();
                if !out.contains(&target) {
                    out.push(target);
                }
            }
        }
    }

    // Rebuild the projected index from this import.  Base Ommverse files have
    // no InterPro projection; enriched files should contain exactly the
    // geometry produced by the current source stream.
    for omm in ommverses.iter_mut() {
        omm.protein_features = ProteinFeatureIndex::new(omm.splice.bin_width, omm.splice.chr_names.len());
    }

    let mut reports = vec![InterProImportReport::default(); ommverses.len()];
    for report in &mut reports {
        report.entries_loaded = entries.len();
    }

    let file =
        File::open(protein2ipr).with_context(|| format!("opening {}", protein2ipr.display()))?;
    let gz = protein2ipr.extension().is_some_and(|x| x == "gz");
    let input: Box<dyn Read> = if gz {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let reader = BufReader::with_capacity(1024 * 1024, input);
    let mut streamed = 0usize;
    for line in reader.lines() {
        let line = line?;
        streamed += 1;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 6 {
            for r in &mut reports {
                r.malformed_records += 1;
            }
            continue;
        }
        let accession = f[0].trim();
        let Some(targets) = locations.get(accession) else {
            if streamed % 1_000_000 == 0 {
                progress(streamed, &reports);
            }
            continue;
        };
        let ipr = f[1].trim();
        let signature = f[3].trim();
        let start_1: u32 = match f[4].trim().parse() {
            Ok(v) if v > 0 => v,
            _ => {
                for &(oi, _) in targets {
                    reports[oi].malformed_records += 1;
                }
                continue;
            }
        };
        let end_1: u32 = match f[5].trim().parse() {
            Ok(v) if v >= start_1 => v,
            _ => {
                for &(oi, _) in targets {
                    reports[oi].malformed_records += 1;
                }
                continue;
            }
        };
        let known = entries.get(ipr);
        let (kind, entry_name) = known
            .cloned()
            .unwrap_or_else(|| (ProteinFeatureKind::Other, f[2].trim().to_owned()));
        let description = if signature.is_empty() {
            entry_name
        } else if entry_name.is_empty() {
            format!("member signature {signature}")
        } else {
            format!("{entry_name} [{signature}]")
        };
        for &(oi, pi) in targets {
            let report = &mut reports[oi];
            report.records_streamed = streamed;
            report.matched_records += 1;
            if known.is_none() {
                report.unknown_entries += 1;
            }
            let feature = ProteinFeature {
                kind,
                protein_range: Some((start_1 - 1, end_1)),
                label: ipr.to_owned(),
                description: description.clone(),
                chromosome: String::new(),
                genomic_start: 0,
                genomic_end: 0,
                source_db: "InterPro".to_owned(),
                review_status: ReviewStatus::Other,
            };
            if !ommverses[oi].proteins[pi].features.contains(&feature) {
                ommverses[oi].proteins[pi].features.push(feature);
                report.features_added += 1;
            } else {
                report.duplicate_features += 1;
            }

            // Bring the annotation home: project the source protein interval
            // through each GENCODE transcript linked to this protein and store
            // the resulting spliced genomic model.
            let transcript_ids = ommverses[oi].proteins[pi].transcript_ids.clone();
            for tx_id in transcript_ids {
                let Some(tx) = ommverses[oi].splice.transcripts.get(tx_id) else { continue; };
                if ommverses[oi].protein_features.add_from_transcript(
                    tx,
                    accession,
                    ipr,
                    kind,
                    description.clone(),
                    signature,
                    start_1 - 1,
                    end_1,
                ).is_some() {
                    report.projected_features_added += 1;
                }
            }
        }
        if streamed % 1_000_000 == 0 {
            progress(streamed, &reports);
        }
    }
    for (omm, report) in ommverses.iter_mut().zip(reports.iter_mut()) {
        report.records_streamed = streamed;
        omm.protein_features.finalize();
        omm.report.feature_records += report.features_added;
    }
    progress(streamed, &reports);
    Ok(reports)
}


impl Ommverse {
    /// Fetch all currently supported reference sources for an assembly into a cache,
    /// then build the integrated Ommverse from those source files.
    pub fn fetch_and_build(assembly: &str, cache_root: impl AsRef<Path>) -> Result<Self> {
        Self::fetch_and_build_with_progress(assembly, cache_root, |_, _| {})
    }

    /// Fetch supported sources and build while reporting coarse build stages and
    /// record counters. The callback is intentionally cheap so callers can feed
    /// the shared Lumrik status server without putting it on hot per-record paths.
    pub fn fetch_and_build_with_progress(
        assembly: &str,
        cache_root: impl AsRef<Path>,
        progress: impl FnMut(&str, usize),
    ) -> Result<Self> {
        Self::fetch_and_build_with_annotation_and_progress(
            assembly, cache_root, sources::AnnotationOptions::default(), progress,
        )
    }

    pub fn fetch_and_build_with_annotation_and_progress(
        assembly: &str,
        cache_root: impl AsRef<Path>,
        annotation: sources::AnnotationOptions,
        mut progress: impl FnMut(&str, usize),
    ) -> Result<Self> {
        progress("resolving reference sources", 0);
        let root = sources::fetch_assembly_with_annotation(assembly, cache_root.as_ref(), &annotation)?;
        Self::build_ucsc_with_gtf_and_progress(root, annotation.explicit_gtf.as_deref(), false, progress)
    }

    fn ingest_fantom5_if_present(&mut self) -> Result<()> {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let Some(url) = sources::fantom5::enhancer_url(&self.assembly) else {
            return Ok(());
        };
        let path = self
            .source_root
            .join("Chromatin")
            .join("FANTOM5")
            .join("F5.hg38.enhancers.bed.gz");
        if !path.is_file() {
            return Ok(());
        }
        let file = File::open(&path).with_context(|| format!("opening {}", path.display()))?;
        let reader = BufReader::new(MultiGzDecoder::new(file));
        for (line_no, line) in reader.lines().enumerate() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 12 {
                bail!(
                    "{}:{}: expected BED12, got {} fields",
                    path.display(),
                    line_no + 1,
                    f.len()
                );
            }
            let start: u32 = f[1].parse()?;
            let end: u32 = f[2].parse()?;
            if start >= end {
                bail!(
                    "{}:{}: invalid interval {start}-{end}",
                    path.display(),
                    line_no + 1
                );
            }
            let block_count: usize = f[9].parse()?;
            let sizes: Vec<u32> = f[10]
                .trim_end_matches(',')
                .split(',')
                .map(str::parse)
                .collect::<std::result::Result<_, _>>()?;
            let starts: Vec<u32> = f[11]
                .trim_end_matches(',')
                .split(',')
                .map(str::parse)
                .collect::<std::result::Result<_, _>>()?;
            if sizes.len() != block_count || starts.len() != block_count {
                bail!(
                    "{}:{}: BED block count mismatch",
                    path.display(),
                    line_no + 1
                );
            }
            let mut blocks = Vec::with_capacity(block_count);
            for (&size, &offset) in sizes.iter().zip(&starts) {
                let block_start = start
                    .checked_add(offset)
                    .context("FANTOM5 block start overflow")?;
                let block_end = block_start
                    .checked_add(size)
                    .context("FANTOM5 block end overflow")?;
                if block_start < start || block_end > end || block_start >= block_end {
                    bail!("{}:{}: invalid BED block", path.display(), line_no + 1);
                }
                blocks.push(RefBlock::new(block_start, block_end));
            }
            self.chromatin.push(ChromatinElement {
                chromosome: f[0].to_owned(),
                region: RefBlock::new(start, end),
                name: f[3].to_owned(),
                score: f[4].parse().unwrap_or(0),
                blocks,
                source: "FANTOM5 enhancer".to_owned(),
                source_url: url.to_owned(),
            });
        }
        Ok(())
    }

    fn ingest_encode4_binding_if_present(
        &mut self,
        progress: &mut dyn FnMut(&str, usize),
    ) -> Result<()> {
        let Some(url) = sources::encode4::tf_rpeaks_url(&self.assembly) else {
            return Ok(());
        };
        let path = self
            .source_root
            .join("Chromatin")
            .join("ENCODE4")
            .join("TFrPeakClusters.bb");
        if !path.is_file() {
            return Ok(());
        }

        // The bigBed is coordinate sorted.  Collapse all overlapping/touching
        // rPeaks while streaming, so memory scales with the union rather than
        // with the ~22 million source observations.
        let mut chromosome_ids = HashMap::<String, u16>::new();
        let mut index = ProteinBindingUnion {
            source: "ENCODE4 TF rPeak union".to_owned(),
            source_url: url.to_owned(),
            ..ProteinBindingUnion::default()
        };
        let mut current: Option<ProteinBindingRegion> = None;
        let mut ctcf_current: Option<CtcfAnchor> = None;
        let mut ctcf_anchors = Vec::<CtcfAnchor>::new();
        let mut ctcf_source_peaks = 0usize;
        let mut imported = 0usize;

        read_all_bigbed(&path, |chrom, start, end, rest| {
            let chromosome_id = if let Some(&id) = chromosome_ids.get(chrom) {
                id
            } else {
                let id = u16::try_from(index.chromosomes.len())
                    .context("too many ENCODE chromosomes")?;
                index.chromosomes.push(chrom.to_owned());
                chromosome_ids.insert(chrom.to_owned(), id);
                id
            };
            let next = RefBlock::new(start, end);
            match current.as_mut() {
                Some(region)
                    if region.chromosome_id == chromosome_id && next.start <= region.region.end =>
                {
                    region.region.end = region.region.end.max(next.end);
                    region.source_peak_count = region.source_peak_count.saturating_add(1);
                }
                _ => {
                    if let Some(region) = current.take() {
                        index.regions.push(region);
                    }
                    current = Some(ProteinBindingRegion {
                        chromosome_id,
                        region: next,
                        source_peak_count: 1,
                    });
                }
            }
            // Do not depend on a hard-coded rPeak column number: the source
            // record carries CTCF as its own tab-delimited factor value. If the
            // upstream schema changes so this no longer works, the zero-anchor
            // guard below makes the build fail loudly.
            if rest
                .split('\t')
                .any(|field| field.trim().eq_ignore_ascii_case("CTCF"))
            {
                ctcf_source_peaks += 1;
                match ctcf_current.as_mut() {
                    Some(anchor)
                        if anchor.chromosome_id == chromosome_id
                            && next.start <= anchor.region.end =>
                    {
                        anchor.region.end = anchor.region.end.max(next.end);
                        anchor.source_peak_count = anchor.source_peak_count.saturating_add(1);
                    }
                    _ => {
                        if let Some(anchor) = ctcf_current.take() {
                            ctcf_anchors.push(anchor);
                        }
                        ctcf_current = Some(CtcfAnchor {
                            chromosome_id,
                            region: next,
                            source_peak_count: 1,
                        });
                    }
                }
            }

            imported += 1;
            if imported % 1_000_000 == 0 {
                eprintln!(
                    "[ommverse] ENCODE4 TF rPeaks: {imported} scanned, {} union regions",
                    index.regions.len()
                );
                progress("building ENCODE4 binding union", imported);
            }
            Ok(())
        })?;
        if let Some(region) = current.take() {
            index.regions.push(region);
        }
        if let Some(anchor) = ctcf_current.take() {
            ctcf_anchors.push(anchor);
        }
        index.source_peak_count = imported;
        self.protein_binding = index;

        if imported > 0 && ctcf_anchors.is_empty() {
            bail!(
                "ENCODE4 rPeaks were present but no CTCF factor records were recognized; refusing to build an empty CTCF architecture"
            );
        }
        let mut domains = Vec::<CtcfDomain>::new();
        for (left_idx, pair) in ctcf_anchors.windows(2).enumerate() {
            let left = pair[0];
            let right = pair[1];
            if left.chromosome_id != right.chromosome_id {
                continue;
            }
            let start = left.region.start + (left.region.end - left.region.start) / 2;
            let end = right.region.start + (right.region.end - right.region.start) / 2;
            if start >= end {
                continue;
            }
            domains.push(CtcfDomain {
                chromosome_id: left.chromosome_id,
                region: RefBlock::new(start, end),
                left_anchor: left_idx as u32,
                right_anchor: (left_idx + 1) as u32,
            });
        }
        eprintln!(
            "[ommverse] CTCF architecture: {ctcf_source_peaks} source rPeaks -> {} anchors -> {} adjacent domains",
            ctcf_anchors.len(),
            domains.len()
        );
        self.ctcf = CtcfArchitecture {
            chromosomes: self.protein_binding.chromosomes.clone(),
            anchors: ctcf_anchors,
            domains,
            source_peak_count: ctcf_source_peaks,
            source: "ENCODE4 CTCF rPeaks; adjacent-anchor candidate domains".to_owned(),
            source_url: url.to_owned(),
        };
        progress("building ENCODE4 CTCF architecture", ctcf_source_peaks);
        progress("building ENCODE4 binding union", imported);
        Ok(())
    }

    /// Import compact recurrence from raw Loop Catalog loop files.
    ///
    /// Cache layout:
    ///   Chromatin/LoopCatalog/experiments.tsv
    /// Each non-comment row is: EXPERIMENT_ID<TAB>RELATIVE_LOOP_FILE.
    /// Repeating EXPERIMENT_ID across callers/resolutions is intentional: all
    /// calls are unioned first, so that experiment contributes at most one vote
    /// to an anchor or loop. Loop files may be plain text or .gz and only need
    /// BEDPE-like chromosome/start/end/chromosome/start/end as their first six
    /// columns.
    fn ingest_loop_catalog_if_present(
        &mut self,
        progress: &mut dyn FnMut(&str, usize),
    ) -> Result<()> {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        const BIN: u32 = 25_000;
        let root = self.source_root.join("Chromatin").join("LoopCatalog");
        let manifest = root.join("experiments.tsv");
        if !manifest.is_file() {
            return Ok(());
        }

        let mut experiment_files = HashMap::<String, Vec<PathBuf>>::new();
        for line in BufReader::new(File::open(&manifest)?).lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut fields = line.split('\t');
            let experiment = fields.next().unwrap_or("").trim();
            let relative = fields.next().unwrap_or("").trim();
            if experiment.is_empty() || relative.is_empty() {
                bail!(
                    "invalid Loop Catalog manifest row {line:?}; expected experiment_id<TAB>relative_file"
                );
            }
            experiment_files
                .entry(experiment.to_owned())
                .or_default()
                .push(root.join(relative));
        }

        type Endpoint = (String, u32);
        type LoopKey = (Endpoint, Endpoint);
        let mut anchor_support = HashMap::<Endpoint, u32>::new();
        let mut loop_support = HashMap::<LoopKey, u32>::new();
        let mut source_loop_count = 0usize;

        for (experiment_idx, (_experiment, files)) in experiment_files.iter().enumerate() {
            let mut experiment_loops = HashSet::<LoopKey>::new();
            for path in files {
                let file = File::open(path)
                    .with_context(|| format!("opening Loop Catalog calls {}", path.display()))?;
                let input: Box<dyn Read> = if path.extension().is_some_and(|x| x == "gz") {
                    Box::new(MultiGzDecoder::new(file))
                } else {
                    Box::new(file)
                };
                for line in BufReader::with_capacity(1024 * 1024, input).lines() {
                    let line = line?;
                    if line.is_empty() || line.starts_with('#') {
                        continue;
                    }
                    let f = line.split_whitespace().collect::<Vec<_>>();
                    if f.len() < 6 {
                        continue;
                    }
                    let (Ok(s1), Ok(e1), Ok(s2), Ok(e2)) = (
                        f[1].parse::<u32>(),
                        f[2].parse::<u32>(),
                        f[4].parse::<u32>(),
                        f[5].parse::<u32>(),
                    ) else {
                        continue;
                    };
                    let midpoint1 = s1.saturating_add(e1.saturating_sub(s1) / 2);
                    let midpoint2 = s2.saturating_add(e2.saturating_sub(s2) / 2);
                    let a = (f[0].to_owned(), midpoint1 / BIN);
                    let b = (f[3].to_owned(), midpoint2 / BIN);
                    let key = if a <= b { (a, b) } else { (b, a) };
                    experiment_loops.insert(key);
                    source_loop_count += 1;
                }
            }
            let mut experiment_anchors = HashSet::<Endpoint>::new();
            for (a, b) in &experiment_loops {
                experiment_anchors.insert(a.clone());
                experiment_anchors.insert(b.clone());
            }
            for anchor in experiment_anchors {
                *anchor_support.entry(anchor).or_default() += 1;
            }
            for loop_key in experiment_loops {
                *loop_support.entry(loop_key).or_default() += 1;
            }
            if (experiment_idx + 1) % 25 == 0 {
                progress("importing experimental loop recurrence", experiment_idx + 1);
            }
        }

        let mut chromosomes = anchor_support
            .keys()
            .map(|x| x.0.clone())
            .collect::<Vec<_>>();
        chromosomes.sort();
        chromosomes.dedup();
        let chromosome_ids = chromosomes
            .iter()
            .enumerate()
            .map(|(i, c)| (c.clone(), i as u16))
            .collect::<HashMap<_, _>>();
        let mut endpoint_keys = anchor_support.keys().cloned().collect::<Vec<_>>();
        endpoint_keys.sort();
        let mut endpoint_to_anchor = HashMap::<Endpoint, u32>::new();
        let mut anchors = Vec::with_capacity(endpoint_keys.len());
        for endpoint in endpoint_keys {
            let id = anchors.len() as u32;
            let start = endpoint.1.saturating_mul(BIN);
            anchors.push(ExperimentalLoopAnchor {
                chromosome_id: chromosome_ids[&endpoint.0],
                region: RefBlock::new(start, start.saturating_add(BIN)),
                experiment_count: anchor_support[&endpoint],
            });
            endpoint_to_anchor.insert(endpoint, id);
        }
        let mut loops = loop_support
            .into_iter()
            .map(|((a, b), experiment_count)| ExperimentalLoop {
                left_anchor: endpoint_to_anchor[&a],
                right_anchor: endpoint_to_anchor[&b],
                experiment_count,
            })
            .collect::<Vec<_>>();
        loops.sort_by_key(|x| (x.left_anchor, x.right_anchor));
        eprintln!(
            "[ommverse] experimental loops: {} experiments, {source_loop_count} source calls -> {} recurrent anchor bins -> {} distinct loops",
            experiment_files.len(),
            anchors.len(),
            loops.len()
        );
        self.experimental_loops = ExperimentalLoopArchitecture {
            chromosomes,
            anchors,
            loops,
            experiment_count: experiment_files.len(),
            source_loop_count,
            anchor_bin_size: BIN,
            source: "Loop Catalog experimental loop recurrence".to_owned(),
            source_url: "https://loopcatalog.lji.org/".to_owned(),
        };
        progress(
            "importing experimental loop recurrence",
            experiment_files.len(),
        );
        Ok(())
    }

    /// Build an Ommverse directly from the directory layout produced by
    /// scripts/download_ucsc_reference.sh.
    pub fn build_ucsc(root: impl AsRef<Path>) -> Result<Self> {
        Self::build_ucsc_with_debug_and_progress(root, false, |_, _| {})
    }

    /// Build from UCSC and optionally report a bounded sample of feature records
    /// that could not be attached to a mapped protein.
    pub fn build_ucsc_with_debug(
        root: impl AsRef<Path>,
        debug_failed_mappings: bool,
    ) -> Result<Self> {
        Self::build_ucsc_with_debug_and_progress(root, debug_failed_mappings, |_, _| {})
    }

    pub fn build_ucsc_with_debug_and_progress(
        root: impl AsRef<Path>,
        debug_failed_mappings: bool,
        progress: impl FnMut(&str, usize),
    ) -> Result<Self> {
        Self::build_ucsc_with_gtf_and_progress(root, None, debug_failed_mappings, progress)
    }

    pub fn build_ucsc_with_gtf_and_progress(
        root: impl AsRef<Path>,
        explicit_gtf: Option<&Path>,
        debug_failed_mappings: bool,
        mut progress: impl FnMut(&str, usize),
    ) -> Result<Self> {
        let root = root
            .as_ref()
            .canonicalize()
            .with_context(|| format!("reference root {}", root.as_ref().display()))?;
        let assembly = root
            .file_name()
            .and_then(|x| x.to_str())
            .context("reference root has no assembly name")?
            .to_owned();
        let genes_dir = root.join("Genes");
        let gtf = match explicit_gtf {
            Some(path) => path.canonicalize().with_context(|| format!("GTF {}", path.display()))?,
            None => find_gtf(&genes_dir)?,
        };
        let twobit = root.join("Genome").join(format!("{assembly}.2bit"));
        let protein_dir = root.join("Protein");
        for required in [&twobit, &protein_dir] {
            if !required.exists() {
                bail!(
                    "required UCSC reference component missing: {}",
                    required.display()
                );
            }
        }

        let gtf_name = gtf
            .file_name()
            .and_then(|x| x.to_str())
            .unwrap_or("annotation.gtf");
        progress(&format!("building transcript index ({gtf_name})"), 0);
        let splice = SpliceIndex::from_path(&gtf, 100_000, IdNameKeys::default())
            .with_context(|| format!("building splice index from {}", gtf.display()))?;
        let mut tx_by_stable = HashMap::<String, TranscriptId>::new();
        for tx in &splice.transcripts {
            for name in &tx.names {
                tx_by_stable
                    .entry(strip_version(name).to_owned())
                    .or_insert(tx.id);
            }
        }

        progress("mapping proteins", 0);
        let (mapping, mapping_schema) = find_swissprot_mapping_bigbed(&protein_dir)?;
        let mut proteins = Vec::<Protein>::new();
        let mut protein_by_accession = HashMap::<String, usize>::new();
        let mut linked = HashSet::new();
        let mut unlinked_mapping_records = 0usize;

        // GENCODE/GTF is authoritative for transcript -> protein identity.
        // Seed the protein universe from every explicit transcript protein_id
        // before UniProt/UCSC enrichment.  Versioned and unversioned aliases
        // resolve to the same Ommverse protein object, and multiple transcripts
        // may deliberately converge on that object.
        for tx in &splice.transcripts {
            if tx.protein_aliases().is_empty() {
                continue;
            }
            let existing = tx
                .protein_aliases()
                .iter()
                .find_map(|name| protein_by_accession.get(name).copied())
                .or_else(|| {
                    tx.protein_aliases()
                        .iter()
                        .find_map(|name| protein_by_accession.get(strip_version(name)).copied())
                });
            let idx = existing.unwrap_or_else(|| {
                let idx = proteins.len();
                proteins.push(Protein {
                    accession: tx.protein_aliases()[0].clone(),
                    entry_name: String::new(),
                    name: String::new(),
                    gene_symbol: String::new(),
                    identifiers: Vec::new(),
                    aliases: Vec::new(),
                    ensembl_gene: None,
                    ensembl_protein: tx
                        .protein_aliases()
                        .iter()
                        .find(|name| name.starts_with("ENSP"))
                        .cloned(),
                    transcript_ids: Vec::new(),
                    review_status: ReviewStatus::Other,
                    features: Vec::new(),
                });
                idx
            });
            if !proteins[idx].transcript_ids.contains(&tx.id) {
                proteins[idx].transcript_ids.push(tx.id);
            }
            linked.insert(tx.id);
            for name in tx.protein_aliases() {
                add_identifier(&mut proteins[idx].identifiers, name);
                protein_by_accession.entry(name.clone()).or_insert(idx);
                protein_by_accession
                    .entry(strip_version(name).to_owned())
                    .or_insert(idx);
            }
        }

        read_all_bigbed(&mapping, |_, _, _, rest| {
            let f: Vec<&str> = rest.split('\t').collect();
            // UCSC mapping schemas vary between assemblies and annotation
            // namespaces.  Do not assume the transcript identifier lives in a
            // fixed column: search from right to left because the most specific
            // cross-reference fields normally live at the end of these records.
            // Only exact GTF-known identifiers (with an optional numeric version
            // stripped) are accepted.
            if f.is_empty() {
                return Ok(());
            }
            // ensGene_*.swissprot.bb is a transcript mapping table, while
            // unipAliSwissprot.bb is a bigPsl protein-to-genome alignment.  In
            // bigPsl the BED name (rest field 0) is the UniProt accession; the
            // later fields are alignment statistics, not UniProt metadata.
            let accession = f[0].trim();
            if accession.is_empty() {
                return Ok(());
            }

            let tx_id = if mapping_schema == MappingSchema::EnsGene {
                let Some(tx_id) = f.iter().rev().find_map(|value| {
                    value
                        .split(|c: char| {
                            !(c.is_ascii_alphanumeric() || c == '_' || c == '.' || c == '-')
                        })
                        .rev()
                        .filter(|candidate| !candidate.is_empty())
                        .find_map(|candidate| tx_by_stable.get(strip_version(candidate)).copied())
                }) else {
                    unlinked_mapping_records += 1;
                    return Ok(());
                };
                linked.insert(tx_id);
                Some(tx_id)
            } else {
                // The canonical UCSC track already proves that this UniProt
                // protein aligns to this assembly.  Transcript linkage is added
                // below from protMapInfo.tsv, whose job is to bridge the protein
                // and assembly transcript namespaces.
                None
            };

            let seeded_idx = tx_id.and_then(|tx_id| {
                splice.transcripts.get(tx_id).and_then(|tx| {
                    tx.protein_aliases().iter().find_map(|name| {
                        protein_by_accession
                            .get(name)
                            .or_else(|| protein_by_accession.get(strip_version(name)))
                            .copied()
                    })
                })
            });
            let idx = if let Some(idx) = seeded_idx {
                protein_by_accession.entry(accession.to_owned()).or_insert(idx);
                idx
            } else {
                *protein_by_accession
                    .entry(accession.to_owned())
                    .or_insert_with(|| {
                        let idx = proteins.len();
                        proteins.push(Protein {
                            accession: accession.to_owned(),
                        entry_name: if mapping_schema == MappingSchema::EnsGene {
                            f.get(22).copied().unwrap_or("").to_owned()
                        } else {
                            String::new()
                        },
                        review_status: if mapping_schema == MappingSchema::EnsGene {
                            review_status(f.get(23).copied().unwrap_or(""))
                        } else {
                            ReviewStatus::SwissProt
                        },
                        name: if mapping_schema == MappingSchema::EnsGene {
                            f.get(26).copied().unwrap_or("").to_owned()
                        } else {
                            String::new()
                        },
                        gene_symbol: if mapping_schema == MappingSchema::EnsGene {
                            f.get(27).copied().unwrap_or("").to_owned()
                        } else {
                            String::new()
                        },
                        identifiers: Vec::new(),
                        aliases: if mapping_schema == MappingSchema::EnsGene {
                            parse_aliases(
                                f.get(30).copied().unwrap_or(""),
                                f.get(31).copied().unwrap_or(""),
                            )
                        } else {
                            Vec::new()
                        },
                        ensembl_gene: if mapping_schema == MappingSchema::EnsGene {
                            nonempty(f.get(f.len().saturating_sub(3)).copied().unwrap_or(""))
                        } else {
                            None
                        },
                        ensembl_protein: if mapping_schema == MappingSchema::EnsGene {
                            nonempty(f.get(f.len().saturating_sub(2)).copied().unwrap_or(""))
                        } else {
                            None
                        },
                        transcript_ids: Vec::new(),
                        features: Vec::new(),
                        });
                        idx
                    })
            };
            add_identifier(&mut proteins[idx].identifiers, accession);
            if proteins[idx].entry_name.is_empty() && mapping_schema == MappingSchema::EnsGene {
                proteins[idx].entry_name = f.get(22).copied().unwrap_or("").to_owned();
            }
            if proteins[idx].name.is_empty() && mapping_schema == MappingSchema::EnsGene {
                proteins[idx].name = f.get(26).copied().unwrap_or("").to_owned();
            }
            if proteins[idx].gene_symbol.is_empty() && mapping_schema == MappingSchema::EnsGene {
                proteins[idx].gene_symbol = f.get(27).copied().unwrap_or("").to_owned();
            }
            if proteins[idx].ensembl_gene.is_none() && mapping_schema == MappingSchema::EnsGene {
                proteins[idx].ensembl_gene = nonempty(f.get(f.len().saturating_sub(3)).copied().unwrap_or(""));
            }
            if proteins[idx].ensembl_protein.is_none() && mapping_schema == MappingSchema::EnsGene {
                proteins[idx].ensembl_protein = nonempty(f.get(f.len().saturating_sub(2)).copied().unwrap_or(""));
            }
            if let Some(tx_id) = tx_id {
                if !proteins[idx].transcript_ids.contains(&tx_id) {
                    proteins[idx].transcript_ids.push(tx_id);
                }
            }
            Ok(())
        })?;

        // UCSC's protMapInfo.tsv is a source-agnostic rescue bridge between
        // UniProt accessions and the transcript namespace used for this
        // assembly.  Some assemblies (notably RefSeq-backed ones) do not have
        // an ensGene Swiss-Prot mapping track, and unipAli*.bb may use isoform
        // accessions as their BED names.  Prefer the mapping BigBed above, but
        // augment it here with canonical UniProt accessions whenever an exact
        // GTF-known transcript identifier is present in protMapInfo.tsv.
        //
        // Example UCSC row:
        // Q9BL78  trembl  ...  NM_058267.6  NM_058267.6  chrI:55339-64021
        //
        // This is intentionally still an exact identifier rescue: no fuzzy
        // gene-name or coordinate inference is performed here.
        let prot_map_info = protein_dir.join("protMapInfo.tsv");
        if prot_map_info.is_file() {
            let text = fs::read_to_string(&prot_map_info)
                .with_context(|| format!("reading {}", prot_map_info.display()))?;
            for line in text.lines() {
                if line.trim().is_empty() || line.starts_with('#') {
                    continue;
                }
                let f: Vec<&str> = line.split('\t').collect();
                let accession = f.first().copied().unwrap_or("").trim();
                if accession.is_empty() {
                    continue;
                }

                let Some(tx_id) = f.iter().rev().find_map(|value| {
                    value
                        .split(|c: char| {
                            !(c.is_ascii_alphanumeric() || c == '_' || c == '.' || c == '-')
                        })
                        .rev()
                        .filter(|candidate| !candidate.is_empty())
                        .find_map(|candidate| tx_by_stable.get(strip_version(candidate)).copied())
                }) else {
                    continue;
                };

                linked.insert(tx_id);
                let seeded_idx = splice.transcripts.get(tx_id).and_then(|tx| {
                    tx.protein_aliases().iter().find_map(|name| {
                        protein_by_accession
                            .get(name)
                            .or_else(|| protein_by_accession.get(strip_version(name)))
                            .copied()
                    })
                });
                let idx = if let Some(idx) = seeded_idx {
                    protein_by_accession.entry(accession.to_owned()).or_insert(idx);
                    idx
                } else {
                    *protein_by_accession
                        .entry(accession.to_owned())
                        .or_insert_with(|| {
                            let idx = proteins.len();
                            proteins.push(Protein {
                                accession: accession.to_owned(),
                            entry_name: String::new(),
                            name: String::new(),
                            gene_symbol: String::new(),
                            identifiers: Vec::new(),
                            aliases: Vec::new(),
                            ensembl_gene: None,
                            ensembl_protein: None,
                            transcript_ids: Vec::new(),
                            review_status: review_status(f.get(1).copied().unwrap_or("")),
                            features: Vec::new(),
                            });
                            idx
                        })
                };
                add_identifier(&mut proteins[idx].identifiers, accession);
                if !proteins[idx].transcript_ids.contains(&tx_id) {
                    proteins[idx].transcript_ids.push(tx_id);
                }
            }
        }

        let feature_tracks = [
            ("unipLocTransMemb.bb", ProteinFeatureKind::Transmembrane),
            ("unipDomain.bb", ProteinFeatureKind::Domain),
            ("unipLocSignal.bb", ProteinFeatureKind::SignalPeptide),
            ("unipLocCytopl.bb", ProteinFeatureKind::Cytoplasmic),
            ("unipLocExtra.bb", ProteinFeatureKind::Extracellular),
            ("unipModif.bb", ProteinFeatureKind::ModifiedResidue),
            ("unipDisulfBond.bb", ProteinFeatureKind::Disulfide),
            ("unipRepeat.bb", ProteinFeatureKind::Repeat),
            ("unipChain.bb", ProteinFeatureKind::Chain),
            ("unipConflict.bb", ProteinFeatureKind::Conflict),
            ("unipInterest.bb", ProteinFeatureKind::Interest),
            ("unipMut.bb", ProteinFeatureKind::Mutagenesis),
            ("unipSplice.bb", ProteinFeatureKind::SpliceVariant),
            ("unipStruct.bb", ProteinFeatureKind::Structure),
        ];
        let mut feature_records = 0usize;
        let mut orphan_features = 0usize;
        let mut warnings = Vec::new();
        for (file, kind) in feature_tracks {
            let path = protein_dir.join(file);
            if !path.exists() {
                warnings.push(format!("optional protein track missing: {file}"));
                continue;
            }
            let mut failed_debug_printed = 0usize;
            read_all_bigbed(&path, |chrom, start, end, rest| {
                let f: Vec<&str> = rest.split('\t').collect();
                if f.is_empty() {
                    return Ok(());
                }

                // UCSC UniProt feature bigBeds are BED12+14 after BigBed's BED3
                // prefix.  The explicit UniProt accession is therefore f[24].
                // f[19] independently states "... on protein ACCESSION" and is
                // used as a schema/identity cross-check rather than as a rescue.
                let accession = f.get(24).copied().unwrap_or("").trim();
                let aa_text = f.get(19).copied().unwrap_or("");
                let described_accession = protein_accession_from_feature_text(aa_text);
                let identity_consistent = !accession.is_empty()
                    && described_accession.map_or(true, |described| described == accession);

                let Some(&idx) = identity_consistent
                    .then(|| protein_by_accession.get(accession))
                    .flatten()
                else {
                    orphan_features += 1;
                    if debug_failed_mappings && failed_debug_printed < 10 {
                        eprintln!(
                            "[ommverse] unmapped feature {file} {chrom}:{start}-{end} accession={accession:?} described={described_accession:?} fields={} label={:?}",
                            f.len(),
                            f.get(18).copied().unwrap_or("")
                        );
                        failed_debug_printed += 1;
                    }
                    return Ok(());
                };

                let source_db = f.get(13).copied().unwrap_or("").to_owned();
                let status_text = f.get(17).copied().unwrap_or("");
                let label = f.get(18).copied().unwrap_or("").to_owned();
                let description = aa_text.to_owned();
                // Feature tracks carry the UniProt/gene metadata that the
                // canonical bigPsl alignment deliberately does not.  Use it to
                // enrich canonical protein records without changing identity.
                let protein = &mut proteins[idx];
                if protein.gene_symbol.is_empty() {
                    protein.gene_symbol = f.get(14).copied().unwrap_or("").to_owned();
                }
                if protein.name.is_empty() {
                    protein.name = f.get(20).copied().unwrap_or("").to_owned();
                }
                if protein.aliases.is_empty() {
                    protein.aliases = parse_aliases(f.get(21).copied().unwrap_or(""), "");
                }
                protein.features.push(ProteinFeature {
                    kind,
                    protein_range: parse_amino_acid_range(aa_text),
                    label,
                    description,
                    chromosome: chrom.to_owned(),
                    genomic_start: start,
                    genomic_end: end,
                    source_db,
                    review_status: review_status(status_text),
                });
                feature_records += 1;
                Ok(())
            })?;
        }

        // Greedily retain every protein identity supplied by protein sources and
        // by the GTF transcript->protein relationship.
        for protein in &mut proteins {
            let accession = protein.accession.clone();
            let entry_name = protein.entry_name.clone();
            let ensembl_protein = protein.ensembl_protein.clone();
            add_identifier(&mut protein.identifiers, &accession);
            add_identifier(&mut protein.identifiers, &entry_name);
            if let Some(ensembl) = ensembl_protein.as_deref() {
                add_identifier(&mut protein.identifiers, ensembl);
            }
            for &tx_id in &protein.transcript_ids {
                if let Some(tx) = splice.transcripts.get(tx_id) {
                    for protein_name in tx.protein_aliases() {
                        add_identifier(&mut protein.identifiers, protein_name);
                    }
                }
            }
        }

        // Protein sources carry additional gene namespaces (HGNC-style symbol,
        // aliases and, where available, Ensembl gene accession).  Attach those
        // names to the gene already established by the exact transcript link.
        // Protein accessions themselves remain protein identifiers; the
        // transcript relation provides the typed bridge back to the gene.
        let mut splice = splice;
        for protein in &proteins {
            for &tx_id in &protein.transcript_ids {
                let Some(tx) = splice.transcripts.get(tx_id) else {
                    continue;
                };
                let Some(gene) = splice.genes.get_mut(tx.gene_id) else {
                    continue;
                };
                gene.add_name(&protein.gene_symbol);
                for alias in &protein.aliases {
                    gene.add_name(alias);
                }
                if let Some(accession) = &protein.ensembl_gene {
                    gene.add_name(accession);
                }
            }
        }

        let report = BuildReport {
            transcripts: splice.transcripts.len(),
            mapped_proteins: proteins.len(),
            linked_transcripts: linked.len(),
            feature_records,
            unlinked_mapping_records,
            feature_records_without_protein: orphan_features,
            warnings,
        };
        let sources_path = root.join("sources.yaml");
        let sources = sources::SourceManifest::load_optional(&sources_path)?;
        let protein_features = ProteinFeatureIndex::new(splice.bin_width, splice.chr_names.len());
        let mut out = Self {
            assembly,
            source_root: root,
            genome_twobit: twobit,
            splice,
            proteins,
            protein_features,
            report,
            interpro_entries: HashMap::new(),
            chromatin: Vec::new(),
            protein_binding: ProteinBindingUnion::default(),
            ctcf: CtcfArchitecture::default(),
            experimental_loops: ExperimentalLoopArchitecture::default(),
            sources,
            protein_by_accession: HashMap::new(),
            proteins_by_gene: HashMap::new(),
            proteins_by_transcript: HashMap::new(),
            gene_by_name: HashMap::new(),
        };
        progress("building lookup indexes", 0);
        out.reindex();
        progress("importing FANTOM5 chromatin", 0);
        out.ingest_fantom5_if_present()?;
        progress("importing ENCODE4 TF rPeaks", 0);
        out.ingest_encode4_binding_if_present(&mut progress)?;
        progress("importing experimental loop recurrence", 0);
        out.ingest_loop_catalog_if_present(&mut progress)?;
        progress("build complete", out.protein_binding.regions.len());
        Ok(out)
    }

    /// Add assembly-relevant InterPro annotations to an existing Ommverse.
    ///
    /// `protein2ipr.dat.gz` is streamed and only records whose UniProt accession
    /// already exists in this Ommverse are retained. InterPro coordinates are
    /// 1-based inclusive and are converted to Ommverse's 0-based half-open
    /// protein coordinates. The stable IPR accession is stored in `label`, the
    /// human-readable entry name plus contributing member signature in
    /// `description`, and `source_db` is `InterPro`.
    pub fn ingest_interpro(
        &mut self,
        protein2ipr: impl AsRef<Path>,
        entry_list: impl AsRef<Path>,
        parent_child_tree: Option<&Path>,
    ) -> Result<InterProImportReport> {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};

        let entry_list = entry_list.as_ref();
        let mut entries = HashMap::<String, (ProteinFeatureKind, String)>::new();
        let reader = BufReader::new(
            File::open(entry_list).with_context(|| format!("opening {}", entry_list.display()))?,
        );
        for line in reader.lines() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let mut f = line.split('\t');
            let Some(ipr) = f.next().map(str::trim).filter(|x| !x.is_empty()) else {
                continue;
            };
            let kind_text = f.next().unwrap_or("").trim();
            let name = f.next().unwrap_or("").trim();
            let kind = interpro_feature_kind(kind_text);
            entries.insert(ipr.to_owned(), (kind, name.to_owned()));
            self.interpro_entries
                .entry(ipr.to_owned())
                .or_insert_with(|| InterProEntry {
                    kind,
                    name: name.to_owned(),
                    parents: Vec::new(),
                    children: Vec::new(),
                });
        }
        if let Some(tree) = parent_child_tree {
            let reader = BufReader::new(
                File::open(tree).with_context(|| format!("opening {}", tree.display()))?,
            );
            let mut stack: Vec<String> = Vec::new();
            for line in reader.lines() {
                let line = line?;
                if line.trim().is_empty() || line.starts_with('#') {
                    continue;
                }
                let mut depth = 0usize;
                let bytes = line.as_bytes();
                while bytes.get(depth * 2..depth * 2 + 2) == Some(b"--") {
                    depth += 1;
                }
                let body = &line[depth * 2..];
                let Some(ipr) = body
                    .split("::")
                    .next()
                    .map(str::trim)
                    .filter(|x| x.starts_with("IPR"))
                else {
                    continue;
                };
                stack.truncate(depth);
                if depth > 0 {
                    if let Some(parent) = stack.get(depth - 1).cloned() {
                        if let Some(entry) = self.interpro_entries.get_mut(ipr) {
                            if !entry.parents.contains(&parent) {
                                entry.parents.push(parent.clone());
                            }
                        }
                        if let Some(entry) = self.interpro_entries.get_mut(&parent) {
                            if !entry.children.iter().any(|x| x == ipr) {
                                entry.children.push(ipr.to_owned());
                            }
                        }
                    }
                }
                stack.push(ipr.to_owned());
            }
        }

        let mut report = InterProImportReport {
            entries_loaded: entries.len(),
            ..Default::default()
        };
        let path = protein2ipr.as_ref();
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let gz = path.extension().is_some_and(|x| x == "gz");
        let input: Box<dyn Read> = if gz {
            Box::new(MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };
        let reader = BufReader::with_capacity(1024 * 1024, input);

        for line in reader.lines() {
            let line = line?;
            report.records_streamed += 1;
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 6 {
                report.malformed_records += 1;
                continue;
            }
            let accession = f[0].trim();
            let Some(&idx) = self.protein_by_accession.get(accession) else {
                continue;
            };
            report.matched_records += 1;

            let ipr = f[1].trim();
            let signature = f[3].trim();
            let start_1: u32 = match f[4].trim().parse() {
                Ok(v) if v > 0 => v,
                _ => {
                    report.malformed_records += 1;
                    continue;
                }
            };
            let end_1: u32 = match f[5].trim().parse() {
                Ok(v) if v >= start_1 => v,
                _ => {
                    report.malformed_records += 1;
                    continue;
                }
            };
            let (kind, entry_name) = entries
                .get(ipr)
                .cloned()
                .unwrap_or_else(|| (ProteinFeatureKind::Other, f[2].trim().to_owned()));
            if !entries.contains_key(ipr) {
                report.unknown_entries += 1;
            }

            let description = if signature.is_empty() {
                entry_name
            } else if entry_name.is_empty() {
                format!("member signature {signature}")
            } else {
                format!("{entry_name} [{signature}]")
            };
            let feature = ProteinFeature {
                kind,
                protein_range: Some((start_1 - 1, end_1)),
                label: ipr.to_owned(),
                description,
                chromosome: String::new(),
                genomic_start: 0,
                genomic_end: 0,
                source_db: "InterPro".to_owned(),
                review_status: ReviewStatus::Other,
            };
            if !self.proteins[idx].features.contains(&feature) {
                self.proteins[idx].features.push(feature);
                report.features_added += 1;
            } else {
                report.duplicate_features += 1;
            }
        }
        self.report.feature_records += report.features_added;
        Ok(report)
    }
}
