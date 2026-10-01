#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn strips_only_numeric_versions() {
        assert_eq!(strip_version("ENSMUST1.5"), "ENSMUST1");
        assert_eq!(strip_version("Q9EST3-1"), "Q9EST3-1");
    }
    #[test]
    fn parses_feature_protein_accession() {
        assert_eq!(
            protein_accession_from_feature_text("amino acids 485-507 on protein Q5GH67"),
            Some("Q5GH67")
        );
        assert_eq!(
            protein_accession_from_feature_text("not a protein range"),
            None
        );
    }
    #[test]
    fn parses_ucsc_uniprot_coordinates() {
        assert_eq!(
            parse_amino_acid_range("amino acids 485-507 on protein Q5GH67"),
            Some((484, 507))
        );
        assert_eq!(
            parse_amino_acid_range("amino acids 42 on protein X"),
            Some((41, 42))
        );
        assert_eq!(
            parse_amino_acid_range("amino acid 197 on protein X"),
            Some((196, 197))
        );
    }
    #[test]
    fn summarizes_geometry_distribution() {
        let s = distribution_summary(&[1, 2, 3, 4, 5]);
        assert_eq!(s.count, 5);
        assert_eq!(s.mean, 3.0);
        assert_eq!(s.median, 3.0);
        assert_eq!(s.q1, 2.0);
        assert_eq!(s.q3, 4.0);
        assert_eq!(s.min, 1);
        assert_eq!(s.max, 5);
        assert!((s.sd - 2.0f64.sqrt()).abs() < 1e-12);
    }
    #[test]
    fn builds_pre_and_post_feature_context() {
        let truth = [false, false, false, true, true, false, false, false];
        assert_eq!(
            feature_context_states(&truth, 2),
            vec![0, 1, 1, 2, 2, 3, 3, 0]
        );
    }
    #[test]
    fn splits_short_inter_feature_gap_by_nearest_boundary() {
        let truth = [false, true, true, false, false, false, true, true, false];
        assert_eq!(
            feature_context_states(&truth, 3),
            vec![1, 2, 2, 3, 1, 1, 2, 2, 3]
        );
    }
}
