# reference_curator

Persistent, evidence-backed reference curation for Lumrik.

`reference_curator` remembers reconstructed sequences that are not adequately
represented by the reference used for an analysis. A candidate may remain
unresolved indefinitely. Later runs can add genomic or external alignment
evidence and, when justified, attach a resolved reference representation
without replacing the sequence that was originally observed.

The canonical genome/annotation is never modified by the store.

## Data model

Each exact observed sequence has one stable `ReferenceCandidate`:

- `sequence` — immutable reconstructed/observed sequence.
- `resolved` — optional curated reference representation.
- `observations` — who observed the sequence, in which sample/run, and with
  what support.
- `annotations` — `HashMap<ReferenceId, Vec<Annotation>>`. Evidence is grouped
  by the exact reference/assembly against which it was generated.
- every `Annotation` records the producer (STAR, BWA, minimap2, BLAST, ...),
  producer version/parameters, target coordinates, strand and alignment
  details.

Exact sequence identity is intentionally strict. Biological/fuzzy grouping is
the producer's job. For example, ClonoMap may decide that several receptor
observations represent one reconstructed V model before submitting that
sequence to the curator.

## Typical workflow

The intended flow is:

```text
ClonoMap / another producer
        |
        | observe(sequence, Observation)
        v
reference_curator store
        |
        | export unresolved candidates as FASTA
        v
sc-mapper / STAR / BWA / minimap2
        |
        | BAM -> collect-evidence
        v
Annotation objects
        |
        | add_annotation(candidate, reference_id, annotation)
        v
reference_curator store
        |
        +-- convincing evidence -> resolve(...)
        |
        `-- insufficient evidence -> keep unresolved
```

A missing genomic match is not an error. Such a candidate simply remains in
the persistent store. It can later be tested against a newer genome assembly
or against other evidence sources such as BLAST/RNA databases.

`sc-mapper` owns mapper execution and BAM parsing. `reference_curator` owns
candidate identity, persistence, accumulated evidence and curation decisions.

## Using it from Rust

Open an existing store, or create a new one automatically if the path does not
exist:

```rust
use reference_curator::{Observation, ReferenceCurator};

let path = "reference_curator.bin";
let mut curator = ReferenceCurator::open(path)?;

let candidate_id = curator.observe(
    b"ACGTACGTACGT",
    Observation {
        source: "clonomap_family".into(),
        sample: Some("ZD-4631-Undetermined".into()),
        run: None,
        kind: "reconstructed_v".into(),
        support: 42,
    },
)?;

curator.save(path)?;
```

Calling `observe()` again with the exact same normalized sequence returns the
same candidate ID and appends the new observation. The curator does not
fuzzy-merge similar sequences.

### Adding mapping evidence

Evidence is attached to both the candidate and the exact reference used:

```rust
use reference_curator::{Annotation, Producer, Strand};

let mut annotation = Annotation::new(
    Producer {
        name: "STAR".into(),
        version: Some("2.7.11b".into()),
        parameters: Some("--outSAMtype BAM Unsorted".into()),
    },
    "chr16",
    12_345_678,
    12_345_897,
    Strand::Forward,
    220,
);

annotation.cigar = Some("89M76N119M".into());
annotation.mapq = Some(60);
annotation.edit_distance = Some(0);

let evidence_id = curator.add_annotation(
    &candidate_id,
    "GRCm39_M39",
    annotation,
)?;

curator.save(path)?;
```

The reference ID should identify the actual reference/assembly used, not merely
the organism. Multiple annotations are retained for one reference, so
secondary/supplementary alignments and evidence from different mappers do not
overwrite each other.

### Resolving a candidate

Resolution is an explicit curation decision. It never changes
`ReferenceCandidate.sequence`:

```rust
use reference_curator::ResolvedReference;

curator.resolve(
    &candidate_id,
    ResolvedReference {
        sequence: b"ACGTACGTACGT".to_vec(),
        source: "genomic placement on GRCm39_M39".into(),
        evidence_ids: vec![evidence_id],
    },
)?;

curator.save(path)?;
```

Evidence IDs used for resolution must belong to that candidate.

## Command-line use

Build the binary with the crate:

```bash
cargo build -r -p reference_curator
```

Inspect a store:

```bash
target/release/reference-curator inspect \
    --store reference_curator.bin
```

This reports the number of candidates and how many are resolved/unresolved.

Export all candidates to FASTA. Resolved candidates use their curated
`resolved.sequence`; unresolved candidates use their original observed
`sequence`:

```bash
target/release/reference-curator export-fasta \
    --store reference_curator.bin \
    --out candidates.fa
```

Export only resolved candidates:

```bash
target/release/reference-curator export-fasta \
    --store reference_curator.bin \
    --out resolved.fa \
    --resolved-only true
```

The last flag is optional:

```
no flag            export all candidates
--unresolved-only  export candidates still needing evidence
--resolved-only    export curated/resolved candidates
```

The CLI is intentionally small at present. Candidate discovery currently comes
through library users such as ClonoMap. Mapper/BAM evidence collection belongs
in `sc-mapper`; as that interface is added, its parsed annotations are fed back
through `add_annotation()`.

## Persistence contract

The store is a versioned binary file (`LUMREFC1`) and is written atomically via
a temporary file followed by rename.

The store is the authoritative persistent state. FASTA, BAM, TSV, AIRR and
future supplemental GTF files are interchange/report products, not competing
databases.

The important invariants are:

1. The observed `sequence` is never silently replaced.
2. The same exact sequence reuses its stable candidate identity.
3. Mapping evidence remembers both *who mapped it* and *what reference was
   mapped against*.
4. Multiple pieces of evidence are retained rather than overwritten.
5. `resolved` is first-class but optional.
6. Unresolved candidates are valid and may remain unresolved indefinitely.
