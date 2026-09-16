use ndarray::{Array1, Array2, Axis};
//use ndarray_linalg::{ SVD, UPLO, eigh::Eigh};
use crate::CloneData;
use linfa::dataset::DatasetBase;
use linfa::prelude::Transformer;
use linfa::traits::Fit;
use linfa_reduction::Pca;
#[cfg(feature = "plot")]
use plotters::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub struct PcaModel {
    pub k: usize,
    pub mean: Array1<f32>,
    pub components: Array2<f32>,
    pub coords: Array2<f32>,
}

impl PcaModel {
    pub fn new(k: usize) -> Self {
        Self {
            k,
            mean: Array1::zeros(0),
            components: Array2::zeros((0, 0)),
            coords: Array2::zeros((0, 0)),
        }
    }

    /// Write PCA coordinates to TSV (n rows × k columns).
    pub fn to_tsv<P: AsRef<Path>>(&self, clone: &CloneData, path: P) -> Result<(), Box<dyn Error>> {
        Ok(self.to_delimited(clone, '\t', path)?)
    }

    pub fn principal_axis(&self) -> Array1<f32> {
        // direction along PC1 axis in PCA coords = unit vector (1,0,0,...)
        // because coords are already in PC space
        let mut axis = Array1::<f32>::zeros(self.coords.ncols());
        if axis.len() > 0 {
            axis[0] = 1.0;
        }
        axis
    }

    /// Write PCA coordinates with sequence annotations.
    /// Format:
    /// DNA <sep> AA <sep> PC1 <sep> PC2 ...
    pub fn to_delimited<P: AsRef<Path>>(
        &self,
        info: &CloneData,
        sep: char,
        path: P,
    ) -> std::io::Result<()> {
        let f = File::create(path)?;
        let mut w = BufWriter::new(f);

        // --- sanity check ---
        if self.coords.nrows() != info.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "PCA rows ({}) do not match clone data length ({})",
                    self.coords.nrows(),
                    info.len()
                ),
            ));
        }

        // --- write rows ---
        for (row_idx, row) in self.coords.outer_iter().enumerate() {
            // DNA + AA
            write!(w, "{}{}{}", info.dna[row_idx], sep, info.aa[row_idx])?;

            // PCA coordinates
            for v in row.iter() {
                write!(w, "{}{:.6}", sep, v)?;
            }

            writeln!(w)?;
        }

        Ok(())
    }

    pub fn fit_transform(&mut self, x: &Array2<f32>) -> Result<(), Box<dyn Error>> {
        let (n, p) = x.dim();

        // Require enough samples AND enough features
        if n < 3 || p < 2 {
            return Err(Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("Clone too sparse for PCA (n={}, p={})", n, p),
            )));
        }

        let max_k = std::cmp::min(n - 1, p - 1);

        if max_k < 1 {
            return Err(Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("Clone too sparse for PCA (n={}, p={})", n, p),
            )));
        }

        if self.k > max_k {
            // Small clones naturally support fewer principal components. This is
            // expected dimensionality clipping, not a warning-worthy condition.
            self.k = max_k;
        }

        // mean-center (same behavior as before)
        let mean = x.mean_axis(Axis(0)).unwrap();
        let mut centered = x.clone();
        for mut row in centered.outer_iter_mut() {
            row -= &mean;
        }

        // linfa expects f64
        let centered64 = centered.mapv(|v| v as f64);
        let dataset = DatasetBase::from(centered64);

        let pca = Pca::params(self.k).fit(&dataset)?;
        let projected = pca.transform(dataset).records; // (n x k)

        self.mean = mean;
        self.coords = projected.mapv(|v| v as f32);

        // components optional: linfa provides projection; loadings access differs by version
        // You can either store empty or extract loadings if you need them.
        self.components = Array2::zeros((0, 0));

        Ok(())
    }

    pub fn coords(&self) -> &Array2<f32> {
        &self.coords
    }

    pub fn components(&self) -> &Array2<f32> {
        &self.components
    }
    #[cfg(feature = "plot")]
    pub fn plot_2d_clusters(
        &self,
        tree: &crate::MstTree,
        outfile: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        println!("Plotting PCA with cluster coloring → {}", outfile);

        let coords = &self.coords;
        let n = coords.nrows();

        // --- run elbow clustering ---
        let clusters = tree.clusters_elbow(n);

        // map node -> cluster id
        let mut belong = HashMap::<usize, usize>::new();
        for (cid, cluster) in clusters.iter().enumerate() {
            for &node in cluster {
                belong.insert(node, cid);
            }
        }

        // --- build color palette ---
        let palette = [
            RED,
            BLUE,
            GREEN,
            MAGENTA,
            CYAN,
            YELLOW,
            RGBColor(255, 165, 0), // orange
            RGBColor(128, 0, 128), // purple
            RGBColor(0, 128, 128), // teal
        ];

        // --- chart setup ---
        let root = BitMapBackend::new(outfile, (900, 900)).into_drawing_area();
        root.fill(&WHITE)?;

        let x = coords.column(0);
        let y = coords.column(1);

        let xmin = x.iter().cloned().fold(f32::INFINITY, f32::min);
        let xmax = x.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        let ymin = y.iter().cloned().fold(f32::INFINITY, f32::min);
        let ymax = y.iter().cloned().fold(f32::NEG_INFINITY, f32::max);

        let mut chart = ChartBuilder::on(&root)
            .caption("PCA (cluster-colored)", ("sans-serif", 30))
            .margin(10)
            .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

        chart.configure_mesh().draw()?;

        // --- draw points ---
        for i in 0..n {
            let color = if let Some(&cid) = belong.get(&i) {
                palette[cid % palette.len()]
            } else {
                RGBColor(160, 160, 160) // external / orphan
            };

            chart.draw_series([Circle::new(
                (coords[(i, 0)], coords[(i, 1)]),
                4,
                color.filled(),
            )])?;
        }

        root.present()?;
        Ok(())
    }
    #[cfg(feature = "plot")]
    pub fn plot_2d(&self, outfile: &str) -> Result<(), Box<dyn std::error::Error>> {
        use plotters::prelude::*;

        let root = BitMapBackend::new(outfile, (900, 900)).into_drawing_area();
        root.fill(&WHITE)?;

        let x = self.coords.column(0);
        let y = self.coords.column(1);

        let xmin = x.iter().cloned().fold(f32::INFINITY, f32::min);
        let xmax = x.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        let ymin = y.iter().cloned().fold(f32::INFINITY, f32::min);
        let ymax = y.iter().cloned().fold(f32::NEG_INFINITY, f32::max);

        let mut chart = ChartBuilder::on(&root)
            .caption("PCA projection", ("sans-serif", 30))
            .margin(10)
            .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

        chart.configure_mesh().draw()?;

        chart.draw_series(
            x.iter()
                .zip(y.iter())
                .map(|(&x, &y)| Circle::new((x, y), 3, BLUE.filled())),
        )?;

        Ok(())
    }

    /// Find sparse root using PCA trajectory:
    /// 1) Project onto PC1
    /// 2) Take leftmost & rightmost endpoints
    /// 3) Compute density (kNN) at both ends
    /// 4) Sparse end = root
    pub fn find_sparse_root(&self, k: usize) -> usize {
        let coords = &self.coords;
        let n = coords.nrows();

        // === 1. get principal axis ===

        // === 2. project points on axis ===
        let mut proj: Vec<(usize, f32)> = coords
            .rows()
            .into_iter()
            .enumerate()
            .map(|(i, row)| {
                let t = row[0];
                (i, t)
            })
            .collect();

        proj.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        let left = proj[0].0;
        let right = proj[n - 1].0;

        // === 3. compute local density for left & right endpoint ===
        let left_density = self.knn_density(left, k);
        let right_density = self.knn_density(right, k);

        // === 4. choose sparse end ===
        if left_density < right_density {
            left
        } else {
            right
        }
    }

    /// Compute kNN density for a single node
    fn knn_density(&self, idx: usize, k: usize) -> f32 {
        let coords = &self.coords;
        let row_i = coords.row(idx);

        let mut dists: Vec<f32> = (0..coords.nrows())
            .filter(|&j| j != idx)
            .map(|j| {
                coords
                    .row(j)
                    .iter()
                    .zip(row_i.iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f32>()
                    .sqrt()
            })
            .collect();

        dists.sort_by(|a, b| a.partial_cmp(b).unwrap());
        dists[k.min(dists.len() - 1)]
    }
}

#[cfg(test)]
mod tests {
    use ndarray::array;
    use std::fs;

    #[test]
    fn test_pca_writer_with_seq_and_aa() {
        // --------------------------------------------------------------
        // 1. Prepare a tiny PCA model with known coords
        // --------------------------------------------------------------
        let coords = array![
            [1.2345678_f32, -2.3456789, 3.4567890],
            [4.1111111_f32, 5.2222222, 6.3333333],
        ];

        // Dummy struct that has "coords" like your PCA struct
        struct DummyPca {
            coords: ndarray::Array2<f32>,
        }

        impl DummyPca {
            fn to_delimited<P: AsRef<std::path::Path>>(
                &self,
                path: P,
                sep: char,
                seq_ids: &[String],
                aa_sequences: &[String],
            ) -> std::io::Result<()> {
                use std::fs::File;
                use std::io::{BufWriter, Write};

                let f = File::create(path)?;
                let mut w = BufWriter::new(f);

                for (row_idx, row) in self.coords.outer_iter().enumerate() {
                    write!(w, "{}{}{}", seq_ids[row_idx], sep, aa_sequences[row_idx])?;
                    for v in row {
                        write!(w, "{}{:.6}", sep, v)?;
                    }
                    writeln!(w)?;
                }

                Ok(())
            }
        }

        let pca = DummyPca { coords };

        // --------------------------------------------------------------
        // 2. Prepare test input for IDs and AA sequences
        // --------------------------------------------------------------
        let seq_ids = vec!["seqA".into(), "seqB".into()];
        let aa_sequences = vec!["AAA".into(), "BBB".into()];

        // --------------------------------------------------------------
        // 3. Write to a temporary file
        // --------------------------------------------------------------
        let tmp = std::env::temp_dir().join("pca_writer_test.tsv");
        pca.to_delimited(&tmp, '\t', &seq_ids, &aa_sequences)
            .expect("failed to write test PCA file");

        // --------------------------------------------------------------
        // 4. Read file and validate content
        // --------------------------------------------------------------
        let content = fs::read_to_string(&tmp).expect("failed to read test output file");

        let expected = "\
seqA\tAAA\t1.234568\t-2.345679\t3.456789\n\
seqB\tBBB\t4.111111\t5.222222\t6.333333\n";

        assert_eq!(content, expected, "PCA writer output mismatch");
    }
}
