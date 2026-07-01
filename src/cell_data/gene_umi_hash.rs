use core::fmt;

/// GeneUmiHash (gene_id, umi)
#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct GeneUmiHash(pub u64, pub u64);

impl fmt::Display for GeneUmiHash {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "GeneUmiHash(gene_id: {}, umi: {})", self.0, self.1)
    }
}

#[test]
fn equality() {
    let a = GeneUmiHash(5, 10);
    let b = GeneUmiHash(5, 10);
    let c = GeneUmiHash(5, 11);

    assert_eq!(a, b);
    assert_ne!(a, c);
}

#[test]
fn ordering() {
    let a = GeneUmiHash(1, 5);
    let b = GeneUmiHash(2, 1);
    let c = GeneUmiHash(2, 10);

    assert!(a < b);
    assert!(b < c);
}
