use anyhow::{Context, Result};

use fast_tag_mapper::{
    BuiltinTagSet,
    FastTagMapper,
};

use mapping_info::MappingInfo;

use scdata::{
    GeneUmiHash,
    MatrixValueType,
    Scdata,
};

use crate::cli::{
    Cli,
    FastFeatureSource,
};

use crate::fastq::ParsedPair;


pub struct FastFeatures {
    pub mapper: FastTagMapper,
    pub data: Scdata,
    pub report: MappingInfo,
}

impl FastFeatures {
    pub fn new(
        mapper: FastTagMapper,
        threads: usize,
    ) -> Self {
        Self {
            mapper,
            data: Scdata::new(
                threads.max(1),
                MatrixValueType::Real,
            ),
            report: MappingInfo::new(
                None,
                0.0,
                usize::MAX,
            ),
        }
    }
}

pub fn build_fast_features(
    args: &Cli,
) -> Result<Option<FastFeatures>> {
    if args.fast_features.is_empty() {
        return Ok(None);
    }

    let mut mapper = FastTagMapper::new();

    for source in
        &args.fast_features
    {
        match source {
            FastFeatureSource::BdSampleHuman => {
                let n =
                    mapper.add_builtin(
                        BuiltinTagSet::Human,
                    );

                println!(
                    "Loaded {n} built-in BD human sample tags"
                );
            }

            FastFeatureSource::BdSampleMouse => {
                let n =
                    mapper.add_builtin(
                        BuiltinTagSet::Mouse,
                    );

                println!(
                    "Loaded {n} built-in BD mouse sample tags"
                );
            }

            FastFeatureSource::Fasta(path) => {
                let n =
                    mapper
                        .load_fasta(path)
                        .with_context(|| {
                            format!(
                                "loading fast feature FASTA {}",
                                path.display()
                            )
                        })?;

                println!(
                    "Loaded {n} fast features from {}",
                    path.display()
                );
            }
        }
    }

    mapper =
        mapper.with_min_hits(
            args.fast_feature_min_hits,
        );

    println!(
        "Fast feature mapper contains {} features",
        mapper.feature_count()
    );

    Ok(Some(
        FastFeatures::new(
            mapper,
            args.threads,
        )
    ))
}

pub fn process_fast_feature_read(
    read: &ParsedPair,
    features: Option<&mut FastFeatures>,
) {
    let Some(features) =
        features
    else {
        return;
    };

    /*
     * Initial Nelrune model:
     *
     * R2 is searched against every requested short-feature
     * reference.
     *
     * No distinction exists here between guide / HTO /
     * sample tag / feature barcode / extra synthetic gene.
     */

    let Some(feature_id) =
        features
            .mapper
            .map_feature_id(
                &read.r2.seq,
                &mut features.report,
            )
    else {
        return;
    };

    features.data.try_insert(
        &read.cell_id,
        GeneUmiHash(
            feature_id,
            read.umi_id,
        ),
        1.0,
        &mut features.report,
    );
}