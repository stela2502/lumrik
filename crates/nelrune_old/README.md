# Nelrune

**Nelrune** is the high-level single-cell sequencing pipeline of the Lumrik ecosystem.

Its goal is to process raw sequencing data into a clean, analysis-ready `scdata` dataset without tying the analysis to a particular sequencing platform, library chemistry, mapper, or vendor pipeline.

Nelrune itself should contain as little sequencing-specific processing logic as possible. It is primarily an **orchestration layer** that connects the specialized Lumrik crates.

The long-term goal is to provide a platform- and chemistry-independent alternative to vendor-specific processing pipelines such as Cell Ranger and the BD Rhapsody pipeline.

> **Raw sequencing data should be interpreted according to the experiment, not according to the vendor that produced it.**

## Status

Nelrune is currently at the initial design stage.

The architecture described here is provisional and records the design decisions made before implementation.

---

# Core idea

A Nelrune analysis is primarily defined by two things:

1. **Sequencing system / chemistry**
2. **Collection strategy**

The sequencing system describes **how the reads are structured**.

The collection strategy describes **what information should be extracted from those reads and how it should be processed**.

These definitions are combined with run-specific information such as:

- input FASTQ/BAM files
- genomic reference/index
- annotations
- output location
- mapper configuration
- additional reference sequences
- processing policies

Together they define the processing path through the Lumrik tools.

Conceptually:

    sequencing system
           |
           v
    primer/read structure
           |
           +-------------------+
                               |
    collection strategy        |
           |                   |
           +---------+---------+
                     |
                     v
              resolved pipeline
                     |
          +----------+----------+
          |          |          |
          v          v          v
      sc_primer  fast_mapper  sc_mapper
          |          |          |
          +----------+----------+
                     |
                     v
                   scdata
                     |
                     v
                 sc-beacon
                     |
                     v
              analysis-ready data

---

# Lumrik components

Nelrune should reuse existing Lumrik crates rather than reimplement their functionality.

## `sc_primer`

`sc_primer` interprets the physical structure of sequencing reads.

It is responsible for identifying elements such as:

- cell barcodes
- UMIs
- fixed primer/adaptor sequences
- genomic inserts
- feature sequences
- sample-tag sequences
- other chemistry-defined read segments

The same abstraction should support both short-read and long-read sequencing.

For standard Illumina experiments, the outer input may be:

    R1 + optional R2

For Oxford Nanopore and other long-read data, a single physical read may contain structures that `sc_primer` can identify and potentially split into multiple logical pieces.

Nelrune should therefore not encode Illumina-specific assumptions into the downstream pipeline.

---

# Sequencing systems

Known sequencing chemistries should be available as presets.

A provisional structure could resemble:

    pub enum SequencingSystem {
        Tenx(TenxChemistry),
        Bd(BdChemistry),
        Ont(OntChemistry),
        Custom(CustomChemistry),
    }

These presets should ultimately resolve into generic `sc_primer` descriptions rather than requiring separate processing implementations.

The purpose of a sequencing-system definition is to answer:

> **How do we interpret the sequences that were produced by this experiment?**

Vendor names should therefore exist mainly at the configuration/preset level, not throughout the processing code.

A custom description layer should allow unsupported or experimental chemistries to be defined without modifying Nelrune itself.

---

# Collection strategy

The collection layer answers a different question:

> **What information do we want to extract from the interpreted reads?**

This should be kept independent of sequencing chemistry.

For example, the same conceptual gene-expression collection may be applicable to:

    10x 3' + gene expression
    BD Rhapsody + gene expression
    ONT + gene expression

The read interpretation differs, while much of the downstream biological collection does not.

Collection definitions determine things such as:

- which `sc_primer` segment is sent to genomic mapping
- which external mapper is used
- what genomic information is collected
- which segments are searched against small reference collections
- how resulting information is stored in `scdata`
- which cleanup/inference steps are performed

---

# Genomic mapping

Large genomic or transcriptomic references are handled through `sc_mapper`.

Nelrune should orchestrate the mapper rather than implement mapping itself.

`sc_mapper` currently provides a common interface around external mappers such as:

- minimap2
- STAR
- BWA

The selected sequencing system and collection strategy may define a sensible default mapper, while the user should ultimately be able to override mapper configuration where appropriate.

Genomic collection may eventually include information such as:

- gene assignment
- exon/intron evidence
- transcript assignment
- variants
- multimapping information

VDJ, translocations, fusion detection and other more specialized analyses are intentionally left for later design.

---

# Collectables

A major design decision is to treat small sets of known, non-genomic sequences generically.

Examples include:

- sample tags
- HTOs
- CRISPR guides
- feature barcodes
- antibody-derived tags
- spike-ins
- custom experimental sequences

Nelrune should **not** require separate fundamental implementations for each of these.

They share the same underlying problem:

> A sequence is compared against a relatively small collection of known sequences, where identifying the matching reference carries biological or experimental meaning.

These are provisionally called **collectables**.

A minimal representation could be:

    pub struct Collectable {
        pub name: String,
        pub reference: PathBuf,
        pub source: SegmentName,
    }

where:

- `name` identifies the resulting dataset/modality
- `reference` describes the sequences to recognize
- `source` identifies the sequence segment produced by `sc_primer` that should be searched

For example:

    [[collectables]]
    name = "Sample Tag"
    reference = "sample_tags.fa"
    source = "SAMPLE_TAG"

    [[collectables]]
    name = "CRISPR Guide"
    reference = "guides.fa"
    source = "FEATURE"

The exact TOML schema remains to be designed.

---

# `fast_mapper`

Collectables should be identified using `fast_mapper`.

This avoids repeatedly constructing or loading heavyweight genomic mapper indices for tiny sets of highly conserved known sequences.

From the mapper's perspective:

    Sample Tag
    CRISPR Guide
    HTO
    Antibody Tag
    Custom Feature

are not fundamentally different concepts.

They are all:

    sequence
       |
       v
    fast_mapper
       |
       v
    reference identity

The semantic meaning comes from the **collectable definition**, not from special-purpose mapper logic.

---

# Named data in `scdata`

The `Collectable.name` should be preserved in the resulting `scdata` representation.

For example:

    [[collectables]]
    name = "CRISPR Guide"
    reference = "guides.fa"
    source = "FEATURE"

with:

    >TP53_g1
    ACGT...

    >TP53_g2
    TGCA...

could naturally produce a named data collection conceptually containing:

    CRISPR Guide / TP53_g1
    CRISPR Guide / TP53_g2

Likewise:

    [[collectables]]
    name = "Sample Tag"
    reference = "sample_tags.fa"
    source = "SAMPLE_TAG"

could produce:

    Sample Tag / Donor_A
    Sample Tag / Donor_B

This allows the name provided by the experimental configuration to propagate naturally into MMX/scdata output.

Importantly, `name` should remain a `String`, not an enum.

Nelrune may provide conventional presets such as `"Sample Tag"` or `"CRISPR Guide"`, but a user should be free to define:

    name = "Spatial Barcode"

or:

    name = "My Experimental Feature"

without requiring a new Nelrune release.

---

# `sc-beacon`

Raw collectable matches should generally not be exposed directly as final biological assignments.

Small feature-reference experiments frequently contain:

- ambient signal
- sequencing errors
- low-level false matches
- ambiguous assignments
- genuine multiple assignments
- multiplets
- experiment-specific background

This applies to CRISPR guides, sample tags, HTOs and many other collectable types.

Therefore collectable data should be suitable for processing by `sc-beacon`.

Conceptually:

    fast_mapper
         |
         v
    raw collectable counts
         |
         v
      scdata
         |
         v
     sc-beacon
         |
         v
    resolved assignments

This may be particularly useful for problematic sample-tag assignments, including cases where vendor pipelines produce substantial ambiguous or unassigned populations.

The intention is for Beacon to solve the **general statistical assignment problem**, rather than having separate HTO-, CRISPR-, or sample-tag-specific cleanup implementations whenever the underlying statistical problem is equivalent.

---

# Input and output

The outer Nelrune CLI should deal with physical run-level information.

At minimum this is expected to include:

- R1
- optional R2
- BAM where applicable
- sequencing-system/chemistry definition
- collection definition
- genomic index/reference
- annotation
- additional collectable references
- output location
- threads/resources
- mapper overrides

A provisional invocation might eventually resemble:

    nelrune \
        --system tenx-3pv3 \
        --collection experiment.toml \
        --r1 sample_R1.fastq.gz \
        --r2 sample_R2.fastq.gz \
        --index genome-index \
        --out sample/

The actual CLI has not yet been designed.

---

# Configuration philosophy

Nelrune should distinguish between **presets** and the generic descriptions they resolve into.

For example:

    Tenx3PrimeV3

is convenient for users, but internally it should resolve into a generic description of the relevant primer/read structure.

Likewise, a known assay may provide a predefined collection configuration, but that configuration should use the same primitives available to custom experiments.

The objective is:

> **Known technologies are presets, not hard-coded pipeline architectures.**

This is essential for supporting new chemistries, custom assays and future sequencing technologies without continually adding special cases.

---

# Provisional run configuration

The resolved internal configuration may eventually resemble:

    pub struct NelruneRun {
        pub input: InputConfig,

        pub chemistry: ChemistryConfig,
        pub collection: CollectionConfig,

        pub mapping: MappingConfig,
        pub references: References,

        pub policy: CollectionPolicy,

        pub output: OutputConfig,
    }

This is deliberately provisional.

Users should normally not need to manually construct the entire resolved configuration. Sequencing-system presets, collection descriptions and CLI arguments should be combined to produce it.

---

# Design principles

### Nelrune is an orchestrator

Processing functionality belongs in reusable Lumrik crates whenever possible.

Nelrune primarily determines which tools run, how they are connected, and where their results go.

### Chemistry and biological collection are independent

The sequencing system describes **how to interpret the reads**.

The collection strategy describes **what information to extract from them**.

These should not become one giant assay enum.

### Vendor technologies are presets

10x, BD and future platform definitions should resolve into generic structures.

Vendor-specific assumptions should not leak unnecessarily into the core pipeline.

### Collectables are generic

HTOs, CRISPR guides, sample tags and similar known sequence sets are all instances of the same underlying abstraction unless a genuine biological/technical difference requires otherwise.

### Preserve semantic names

User- or preset-defined collectable names should propagate into `scdata`/MMX rather than being translated into hard-coded Nelrune categories.

### Clean before presenting

Raw assignment ambiguity and background should be retained long enough for `sc-beacon` to make informed decisions rather than prematurely collapsing information.

### Avoid premature information loss

Mapping ambiguity, barcode evidence and collectable evidence should be retained where practical until the component responsible for resolving that ambiguity can make the decision.

### Do not design around Illumina

Paired FASTQ is an input format, not the fundamental abstraction of the pipeline.

The architecture should remain capable of processing long-read data where `sc_primer` may interpret or split a single physical read differently.

---

# Current scope

The initial implementation should concentrate on the path:

    FASTQ
      |
      v
    sc_primer
      |
      +--------> fast_mapper ----> collectables
      |
      +--------> sc_mapper ------> genomic evidence
                                   |
                                   v
                                 scdata
                                   |
                                   v
                               sc-beacon
                                   |
                                   v
                                 output

The immediate goal is not to reproduce every feature of existing vendor pipelines.

The goal is to establish a generic architecture capable of expressing their common processing requirements without being constrained by their assumptions.

## Later

Potential future extensions include:

- VDJ processing
- translocation detection
- fusion detection
- richer variant-aware processing
- long-read transcript/isoform analysis
- hybrid short-/long-read experiments

These should be considered when designing interfaces, but should **not** complicate the first implementation unnecessarily.

---

# Long-term goal

Nelrune should make the sequencing platform and library chemistry configuration details rather than architectural boundaries.

Ideally, adding a new experiment eventually consists primarily of describing:

1. **what the reads look like**
2. **which pieces contain useful information**
3. **where those pieces should be mapped**
4. **what evidence should be collected**
5. **how that evidence should be resolved**

rather than implementing another pipeline.

If that succeeds, 10x, BD Rhapsody, Oxford Nanopore and custom single-cell experiments can all travel through the same underlying system while retaining the flexibility to treat their data appropriately.