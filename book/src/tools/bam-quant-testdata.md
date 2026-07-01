# bam-quant-testdata

`bam-quant-testdata` generates synthetic BAM data together with ground-truth quantification results.

It is used to test and validate `bam-quant` end-to-end.

## Status

Internal / developer-focused tool.

Not intended for general users, but highly useful for testing, debugging, and method development.

## What it does

`bam-quant-testdata` creates:

- synthetic BAM files
- matching annotation (implicit or explicit)
- ground-truth count matrices

The generated data is fully controlled, allowing exact validation of quantification results.

## Why it exists

Testing single-cell quantification is difficult because:

- real data has no ground truth
- annotation may be incomplete
- alignments are noisy
- expected counts are not known exactly

This tool solves that by generating:

👉 reads with known origin  
👉 exact expected counts  
👉 controlled edge cases  

## Typical use

1. Generate test data:

```bash
bam-quant-testdata \
  --out testdata/
```

2. Run quantification:

```bash
bam-quant \
  --bam testdata/reads.bam \
  --index testdata/index.splice.idx \
  --outpath quant
```

3. Compare results:

```bash
# using internal comparison tools / tests
```

## What is encoded

The generated dataset can include:

- known gene assignments
- known transcript structures
- controlled intronic reads
- known SNP positions (optional)
- defined CB / UB combinations

This allows validation of:

- gene-level counts
- transcript-level counts
- intronic classification
- SNP ref/alt assignment

## Role in development

`bam-quant-testdata` is used to:

- build regression tests
- validate refactoring
- debug matching logic
- ensure correctness of edge cases

It enables true end-to-end testing:

```
synthetic truth → BAM → bam-quant → output → comparison
```

## Design philosophy

Instead of testing isolated components, this tool enables:

- full pipeline validation
- biologically meaningful test scenarios
- reproducible debugging

## Limitations

- synthetic data cannot capture all real-world complexity
- alignment artifacts may differ from real data
- test scenarios must be explicitly designed

## Relationship to bam-quant

This tool exists solely to validate `bam-quant`.

It ensures that changes in:

- matching logic
- SNP handling
- intronic classification

do not silently break expected behaviour.

## See also

- [`bam-quant`](./bam-quant.md)