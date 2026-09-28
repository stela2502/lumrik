# How to query Ommverse

Ommverse is Lumrik's integrated biological reference model. There are currently two ways to ask it questions:

1. the general SQL-like `ommverse query` interface, currently implemented for genomic variants; and
2. existing typed Ommverse commands for genes, proteins, model training, and model inspection.

The long-term direction is to move biological lookup behind the general query engine while keeping specialized computational workflows as explicit commands.

## Build or update a reference first

For a supported assembly, Ommverse owns the assembly-specific directory below the cache root:

```bash
target/release/ommverse build --assembly hg38 --cache /data1/UCSC
```

With no explicit `--out`, this writes the model to:

```text
/data1/UCSC/hg38/hg38.ommverse
```

The build reuses source files already present in the reference tree and fetches missing configured resources. The resolved source manifest is retained with the Ommverse model so query execution can locate source-backed resources such as the indexed variant VCF.

An existing downloaded reference tree can instead be built explicitly:

```bash
target/release/ommverse build \
    --reference /data1/UCSC/hg38 \
    --out /data1/UCSC/hg38/hg38.ommverse
```

## General query language

The current grammar is deliberately small:

```text
SELECT <entity> WHERE <field> = <value> [AND <field> = <value> ...]
```

Keyword matching is case-insensitive. The first supported entity is `variant` (`variants` is accepted as an alias).

### Exact genomic position

```bash
target/release/ommverse query \
    --index /data1/UCSC/hg38/hg38.ommverse \
    'SELECT variant WHERE chromosome = "chr17" AND position = 7674894'
```

`chromosome` may also be written as `chrom` or `chr`; `position` may be written as `pos`.

### Genomic range

A colon denotes an inclusive biological coordinate range:

```bash
target/release/ommverse query \
    --index /data1/UCSC/hg38/hg38.ommverse \
    'SELECT variant WHERE chromosome = "chr17" AND position = 7674894:7679000'
```

`7674894:7679000` means positions 7,674,894 through 7,679,000 inclusive. User-facing positions are 1-based. The query backend converts them to the coordinate convention required by the indexed source.

Current variant output is tab-separated:

```text
chromosome    position    identifier    reference    alternate(s)
```

For example:

```text
chr17    7678991    rs1160313639    C    T
```

### What the variant query does internally

The query is not a linear scan of dbSNP and it does not require all dbSNP records to be stored inside the `.ommverse` file. Ommverse:

1. parses the query into its internal `Query`/`Predicate` representation;
2. resolves the configured `variants` resource from the source manifest stored with the model;
3. translates Ommverse chromosome naming to the source representation when required (for hg38, for example, `chr17` is resolved to `NC_000017.11` for the current dbSNP source);
4. uses the indexed VCF/TBI backend for the requested interval; and
5. translates the result back to Ommverse-facing coordinates/names.

This is intentional: the query language describes biology while the query engine chooses the appropriate backing resource/index.

### Query grammar currently supported

Supported now:

```text
SELECT variant WHERE chromosome = "chr17" AND position = 7674894
SELECT variants WHERE chr = "chr17" AND pos = 7674894:7679000
SELECT variant WHERE chromosome = "chr17" AND position = 7674894:7679000 AND clinical_effect IS NOT NULL
SELECT variant WHERE chromosome = "chr17" AND position = 7674894:7679000 AND clinical_effect IS NULL
```

Values currently supported are:

```text
"text"
12345
12345:67890
```

The numeric range is inclusive. Multiple predicates are combined with `AND`. `clinical_effect IS NOT NULL` keeps variants that have a matching ClinVar allele in the same interval; `clinical_effect IS NULL` keeps variants without a matching ClinVar allele. Ommverse queries the indexed dbSNP and ClinVar resources separately and intersects them by position, REF and ALT.

Not yet implemented in the query grammar include `OR`, comparison operators, `ORDER BY`, `GROUP BY`, generic joins, gene-based variant lookup, consequence lookup, effect lookup, or arbitrary VCF INFO-field predicates. Do not assume SQL features that are not listed above.

## Existing typed lookup commands

The general query language is new. Ommverse already contains useful typed lookup commands that remain available.

### Gene lookup

```bash
target/release/ommverse gene \
    --index /data1/UCSC/hg38/hg38.ommverse \
    TP53
```

`gene` accepts a literal gene symbol. Existing Ommverse gene-input support can also autodetect supported file input forms where present in the current build, including one-symbol-per-line input and 10x-style `features.tsv[.gz]` input.

### Protein lookup

```bash
target/release/ommverse protein \
    --index /data1/UCSC/hg38/hg38.ommverse \
    <ACCESSION>
```

To also reconstruct/show sequence information:

```bash
target/release/ommverse protein \
    --index /data1/UCSC/hg38/hg38.ommverse \
    <ACCESSION> \
    --sequence
```

Protein sequence is reconstructed from the biological reference rather than treated as an unrelated duplicate identity.

## InterPro enrichment

InterPro annotations can be streamed into an existing Ommverse index:

```bash
target/release/ommverse ingest-interpro \
    --index /data1/UCSC/hg38/hg38.ommverse \
    --protein2ipr <protein2ipr.dat.gz> \
    --entry-list <entry.list> \
    --out /data1/UCSC/hg38/hg38.interpro.ommverse
```

An optional `--parent-child-tree` can retain the InterPro hierarchy. `--index-dir` supports recursively enriching base `.ommverse` indexes in a directory instead of specifying one `--index`.

## Protein-feature model commands

These commands are not database lookups, but they are part of the current Ommverse binary and operate on information contained in an Ommverse index.

### Inspect a supervised feature corpus

```bash
target/release/ommverse training-corpus \
    --index <reference.ommverse> \
    --feature <FEATURE>
```

Optional `--flank-width` controls the context width used when reporting collisions between adjacent feature flanks.

### Train an exact-amino-acid feature HMM

```bash
target/release/ommverse train-exact-aa \
    --index <reference.ommverse> \
    --feature <FEATURE> \
    --out <model>
```

Useful options include `--flank-width` and `--scan-unannotated`.

### Train a model vault

```bash
target/release/ommverse train-model-vault \
    --index <reference.ommverse> \
    --out <vault>
```

### Inspect a model vault

```bash
target/release/ommverse model-vault-info --vault <vault>
```

### Evaluate a model vault on another Ommverse reference

```bash
target/release/ommverse evaluate-model-vault \
    --index <target.ommverse> \
    --vault <vault>
```

### Scan unannotated proteins with a frozen exact-AA model

```bash
target/release/ommverse scan-exact-aa \
    --index <reference.ommverse> \
    --model <model>
```

## Where the query language should go next

The useful next step is not to reproduce a generic SQL database. It is to add biological access paths and let the query planner choose the appropriate index.

Examples of intended future questions are:

```sql
SELECT variant WHERE gene = "TP53"
SELECT variant WHERE gene = "TP53" AND clinical_effect IS NOT NULL
SELECT variant WHERE gene = "TP53" AND protein_feature = "DNA-binding domain"
SELECT effect WHERE variant = "rs..."
```

The corresponding indexes should represent biological relationships such as:

```text
gene -> genomic intervals
gene -> transcripts
transcript -> protein
protein -> functional features
variant -> transcript consequences
variant -> effects
external variant identifier -> genomic allele
```

Genomic interval queries should continue to use specialized interval indexes such as tabix where they are already the correct data structure. The Ommverse query planner should compose those indexes with the cross-resource biological relationships rather than replacing them with one monolithic storage engine.

## Discover the exact CLI available in a build

The binary remains authoritative for the exact command surface:

```bash
target/release/ommverse --help
target/release/ommverse query --help
target/release/ommverse gene --help
target/release/ommverse protein --help
```

Use the corresponding `<subcommand> --help` for the training/model commands as their options evolve.
