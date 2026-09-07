# Nelrune posterior integration

Recommended position in the pipeline:

```text
normalization -> STAR -> BAM quantification -> expression matrix
                       |                    |
                       +------ retained ----+
                                  |
                                  v
                               sc-vdj
```

The posterior analyzer should receive the *same* cell/UMI identity used by Nelrune quantification. Do not introduce a second barcode parser.

Minimal integration shape:

```rust
struct NelruneIdentity;
impl sc_vdj::BamIdentityResolver for NelruneIdentity {
    fn resolve(&self, record: &rust_htslib::bam::Record) -> Option<(String, String)> {
        // Call the existing Nelrune QNAME decoder here.
        todo!()
    }
}

struct NelruneExpression<'a>(&'a Scdata);
impl sc_vdj::ExpressionMatrix for NelruneExpression<'_> {
    fn expression(&self, cell: &str, gene: &str) -> f64 {
        // Thin adapter onto existing expression storage.
        todo!()
    }
}

let reference = sc_vdj::VdjReferenceBuilder::default().build(gtf, genome)?;
let mapper = sc_vdj::VdjMapper::new(reference.clone(), Default::default());
let reads = sc_vdj::read_bam(mapper_bam, &NelruneIdentity)?;
let analyzer = sc_vdj::PosteriorAnalyzer::new(&reference, &mapper, Default::default());
let cells = analyzer.analyze(reads, &NelruneExpression(&scdata));
sc_vdj::output::write_reports(out_dir.join("vdj"), &cells)?;
```

The two `todo!()` blocks are intentionally the only Nelrune-specific glue: use the code that already exists in the workspace rather than inventing competing cell IDs or matrix formats.
