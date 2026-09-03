# sc-te

`sc-te` is the transposable-element analysis component of **Lumrik**.

It is designed for single-cell RNA-seq data where reads originating from transposable elements (TEs) frequently map to multiple genomic locations.

Rather than assuming that every short read can be assigned to one exact TE insertion, `sc-te` represents TE expression at a configurable spatial resolution:

```text
chromosome + genomic bin + TE subfamily
```

For example, using the current default 1 Mb genomic bins:

```text
chr1:100000000-101000000:L1M2
```

This provides spatial information while avoiding the assumption that short-read sequencing can reliably distinguish millions of highly similar individual TE loci.

The TE implementation builds on existing Lumrik components:

* `gtf-splice-index` for genomic annotation and spatial indexing
* `scdata` for sparse single-cell count storage
* `sc-mapper` for selective remapping
* STAR as the external mapper
* `sc-te` for TE-specific candidate collection and multimapper resolution

## Current status

The basic analysis pipeline is implemented.

It can:

* convert the UCSC RepeatMasker table into a GTF suitable for Lumrik
* index millions of RepeatMasker loci
* represent TE features as chromosome + genomic bin + TE subfamily
* identify unambiguous TE assignments
* reuse low-order multimapping information already present in a name-sorted BAM
* selectively remap difficult reads with STAR
* collect ambiguous TE candidates
* resolve multimappers using EM
* produce sparse single-cell TE count matrices

The current EM implementation operates on the multimapper pool itself.

Unique TE anchors are collected and tracked, but **they are not currently used to influence the EM probabilities**. Anchor-informed multimapper resolution is therefore an experimental direction rather than part of the current algorithm.

---

# Why spatial TE features?

RepeatMasker annotations contain millions of individual TE loci.

For example, the current mouse mm10/GRCm38 RepeatMasker index contains approximately:

```text
4,261 TE genes/subfamilies
5.3 million RepeatMasker loci
```

Trying to quantify every locus independently creates two related problems.

First, short reads often do not contain enough sequence information to distinguish closely related TE copies.

Second, treating every possible locus as an independent single-cell feature produces an extremely sparse and highly ambiguous feature space.

`sc-te` therefore uses the genomic bins already provided by `gtf-splice-index`.

The current feature identity is:

```text
(chr_id, bin_id, gene_id)
```

where `gene_id` corresponds to the RepeatMasker TE name/subfamily.

This means that multiple copies of the same TE subfamily within the same genomic bin contribute to the same spatial TE feature.

The current bin width is:

```text
1,000,000 bp
```

Importantly, **1 Mb is currently an experimental resolution choice rather than a claim that this is the biologically optimal scale**.

Because binning belongs to the genomic index rather than the TE algorithm itself, alternative resolutions can be tested without redesigning `sc-te`.

---

# Building the TE reference

## UCSC RepeatMasker input

`sc-te` can convert the official UCSC `rmsk` table into a GTF suitable for `gtf-splice-index`.

For mouse mm10/GRCm38, download:

```text
https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz
```

Then convert it:

```bash
sc-te ucsc-rmsk-to-gtf \
    --input rmsk.txt.gz \
    --output mm10_rmsk_TE.gtf
```

The converter maps RepeatMasker information approximately as:

```text
repName   -> gene_id / gene_name
repFamily -> family_id
repClass  -> class_id
```

Each RepeatMasker annotation is represented as a GTF feature with a unique transcript identifier.

UCSC genomic coordinates are converted from 0-based half-open coordinates to the 1-based inclusive coordinates expected by GTF.

## Build the spatial index

The resulting GTF can be indexed with `gtf-splice-index`:

```bash
gtf-splice-index build \
    --annotation mm10_rmsk_TE.gtf \
    --index mm10_rmsk_TE.gtf.dat
```

The resulting `.dat` file is the TE annotation used by `sc-te` and `nelrune-te`.

---

# Input BAM

The production `nelrune-te` analysis uses a **query-name-sorted BAM**.

This is important because normal alignment BAMs can contain primary and secondary mappings for the same read at widely separated genomic coordinates.

Grouping by QNAME allows `nelrune-te` to inspect the complete mapping information already present for a read before deciding whether remapping is necessary.

A Cell Ranger BAM can for example be prepared with:

```bash
samtools sort -n \
    -@ 16 \
    -T "$SNIC_TMP/nelrune_te_sort" \
    -o "$SNIC_TMP/input.name_sorted.bam" \
    input.bam
```

The BAM is expected to contain the normal single-cell barcode and UMI tags, typically:

```text
CB
UB
```

and mapping multiplicity information through `NH`.

---

# Reusing existing alignments

A major design goal of `nelrune-te` is to avoid unnecessarily remapping reads whose useful alignment information is already present in the original BAM.

For each QNAME group, the existing primary and secondary mappings are inspected.

Low-order multimappers can therefore be converted directly into spatial TE candidates.

The current default strategy is conceptually:

```text
complete low-NH group
        |
        +--> reuse existing BAM alignments

high-NH / incomplete / unmapped group
        |
        +--> submit one representative read to STAR
```

The remapping threshold is configurable with:

```text
--remap-nh-above
```

with the current experimental default:

```text
NH > 5
```

A low-NH group is only reused when the BAM appears to contain the complete set of mappings reported by `NH`.

Otherwise the read is treated as requiring remapping.

This avoids interpreting an incomplete collection of secondary alignments as the complete candidate space.

---

# Spatial TE candidates

Each mapped alignment is intersected with the TE splice index.

Matching CIGAR operations contribute genomic overlap, while deletions and skipped regions only advance the reference position.

The resulting TE candidates are collapsed to unique spatial features:

```text
chromosome + bin + TE subfamily
```

Consequently, several individual RepeatMasker loci can collapse to the same spatial TE candidate.

This distinction is important:

```text
number of genomic alignments != number of spatial TE candidates
```

A read may therefore have several genomic mappings while still resolving to a single spatial TE feature.

---

# Anchor evidence

When a read or molecule resolves unambiguously to one spatial TE feature, it contributes anchor evidence.

Anchor counts are stored separately in `Scdata`.

Conceptually:

```text
read
 |
 +--> one spatial TE candidate
         |
         +--> anchor
```

These anchors provide direct evidence that a TE subfamily is expressed within a particular genomic region.

They also provide potentially useful information for resolving ambiguous molecules.

At present, however, the EM algorithm does **not** use anchor abundance as a prior.

---

# Selective STAR remapping

Reads that cannot safely reuse their original BAM mappings are submitted to STAR.

This includes groups such as:

* unmapped reads
* high-order multimappers
* groups whose available secondary alignments do not appear to represent the complete `NH` candidate set

Only one representative sequence is submitted for each QNAME group.

The representative preferentially comes from the primary alignment.

When necessary, reverse-strand BAM sequence and quality information is transformed back into the original read orientation before being sent to STAR.

The current STAR multimapping ceiling is deliberately finite rather than attempting to enumerate arbitrarily many mappings.

The current default is:

```text
--max-multimap 100
```

The intention is to capture useful candidate information without spending excessive computational effort enumerating hundreds or thousands of nearly equivalent genomic placements.

---

# Multimapper resolution

Reads or molecules with more than one spatial TE candidate enter the multimapper pool.

`TeCollector` currently performs an EM-based assignment of this pool.

Conceptually:

```text
candidate sets
     |
     v
multimapper EM
     |
     v
fractional spatial TE counts
```

The resulting counts are stored separately from direct anchor evidence.

Current outputs distinguish:

* direct anchor signal
* uniquely rescued signal
* multimapper EM signal
* EM signal associated with features that have anchor support
* EM signal associated with features without anchor support

A combined matrix can then be generated from the individual components.

## Important current limitation

The distinction between anchored and unanchored features is currently useful for reporting and downstream inspection, but **anchor abundance does not currently determine the EM probabilities**.

The current implementation should therefore not be described as anchor-informed EM.

---

# Diagnostic analysis

`nelrune-te-tests` can be used to inspect the structure of TE ambiguity before running the complete analysis.

It reports information including:

* mapped and unmapped records
* availability of cell barcode and UMI tags
* NH distribution
* TE overlap
* spatial TE anchors
* TE features observed only through multimappers
* number of spatial candidates per multimapping read
* whether ambiguity occurs within a bin, between bins, or between chromosomes
* processing throughput

For example:

```bash
nelrune-te-tests \
    --bam input.bam \
    --te-index mm10_rmsk_TE.gtf.dat \
    --max-reads 5000000 \
    --out te_diagnostic
```

A current 5-million-record mouse test produced approximately:

```text
5.0 M BAM records
4.26 k TE genes/subfamilies
5.33 M RepeatMasker loci

~859 k unambiguous TE anchor reads
~98 k spatial features with anchor support
~15 k features observed only through multimappers
```

Among TE-associated multimappers in this test:

```text
30.7%  -> 1 spatial candidate
56.9%  -> <= 2
68.2%  -> <= 3
74.4%  -> <= 4
78.6%  -> <= 5
```

Most remaining ambiguity was between genomic locations on different chromosomes.

These numbers are dataset-specific and should be treated as diagnostic observations rather than general properties of TE expression.

---

# Output interpretation

The central goal of `sc-te` is not necessarily to claim expression from one exact RepeatMasker insertion.

Instead, the default interpretation is:

> evidence for expression of a TE subfamily within a genomic region.

This distinction becomes particularly important for young and highly repetitive TE families where the short-read sequence may fundamentally not contain enough information to identify the originating locus.

The resulting feature matrix can be treated similarly to other single-cell feature matrices and used for downstream analyses such as:

* cell-type-specific TE expression
* differential TE activity
* regional TE activation
* clustering based on TE activity
* association between gene-expression state and TE activity
* comparison of TE programs between experimental conditions

---

# Outlook

Several aspects of the current approach are deliberately experimental.

## Anchor-informed multimapper resolution

The most immediate extension is to incorporate anchor evidence directly into multimapper assignment.

For example, consider a molecule with candidates:

```text
chr1:100-101Mb:L1M2
chr7:42-43Mb:L1M2
chr12:88-89Mb:L1M2
```

If the surrounding dataset contains strong unambiguous evidence for the first feature and essentially none for the other two, that information could contribute to the probability of assigning the ambiguous molecule.

A future EM model could therefore initialize or constrain candidate probabilities using uniquely assigned molecules.

This would turn the current:

```text
multimappers -> EM
```

into something closer to:

```text
unique spatial evidence
          |
          v
multimappers -> anchor-informed EM
```

Whether this improves biological accuracy needs to be tested rather than assumed.

## Determine the useful spatial resolution

The current 1 Mb bin size is intentionally simple.

Potential resolutions include:

```text
1 Mb
500 kb
100 kb
individual RepeatMasker loci
```

Smaller bins increase spatial resolution but also increase ambiguity and sparsity.

Larger bins sacrifice locus precision but allow more genomic mappings to collapse into interpretable regional signals.

The optimal resolution may also differ between TE families.

A useful future benchmark is therefore to quantify how candidate multiplicity, anchor density and biological reproducibility change with bin size.

## Determine whether extreme multimappers are useful

The most repetitive reads may have tens or hundreds of plausible genomic origins.

Such reads are expensive to map and may contain little useful spatial information.

The selective-remapping architecture makes it possible to test this directly.

For example:

```text
A: anchors / unique spatial assignments only

B: A + low-order multimappers from the original BAM

C: B + high-order and unmapped reads rescued with STAR
```

If B and C produce essentially the same biological structure while C requires substantially more computation, aggressive remapping of the highly repetitive tail may not be justified.

Conversely, if particular young TE families gain substantial reproducible signal only in C, those reads may need a different treatment.

## Separate spatial and non-spatial TE information

A read can be impossible to assign spatially while still providing useful evidence for TE-family or TE-subfamily expression.

Future versions could therefore maintain two related representations:

```text
spatial:
chromosome + bin + TE subfamily

global:
TE subfamily
```

Highly repetitive reads that cannot be resolved spatially could still contribute to global TE abundance instead of being discarded.

This may be particularly useful for young TE families.

## Family-specific resolution

Different TE families have very different sequence uniqueness.

An eventual model could therefore allow the spatial resolution to depend on how much information the reads actually provide.

Older/diverged elements might support relatively fine spatial localization, while young/high-copy elements may only support regional or global subfamily abundance.

This should only be introduced if empirical benchmarking demonstrates a clear advantage over the simpler fixed-bin model.

---

# Design principle

The central principle of `sc-te` is:

> Do not claim more spatial resolution than the sequencing data supports.

The method therefore attempts to preserve useful genomic information while explicitly representing ambiguity.

Rather than forcing every TE-derived short read onto one exact RepeatMasker locus, `sc-te` asks what level of genomic localization can actually be supported by the observed alignments.

The current 1 Mb spatial model, selective reuse of existing BAM mappings, selective STAR remapping and multimapper EM are the first implementation of that idea.

