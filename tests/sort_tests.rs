use bitvec::vec::BitVec;
use kira_spatial_io::determinism::sort::sort_bins;
use kira_spatial_io::{CoordSystem, SpatialDomain};

#[test]
fn sort_bins_permutation_is_consistent() {
    let mut domain = SpatialDomain::new(
        vec![30.0, 10.0, 20.0],
        vec![3.0, 1.0, 2.0],
        Some(vec![3, 1, 2]),
        Some(vec![30, 10, 20]),
        vec![100, 101, 102],
        BitVec::from_iter([true, false, true]),
        CoordSystem::Pixel,
        0,
    )
    .expect("domain");

    let result = sort_bins(&mut domain).expect("sort");

    assert_eq!(result.old_to_new, vec![2, 0, 1]);
    assert_eq!(domain.grid_row.as_deref(), Some(&[1, 2, 3][..]));
    assert_eq!(domain.grid_col.as_deref(), Some(&[10, 20, 30][..]));
    assert_eq!(domain.x, vec![10.0, 20.0, 30.0]);
    assert_eq!(domain.y, vec![1.0, 2.0, 3.0]);
    assert!(!domain.tissue_mask[0]);
    assert!(domain.tissue_mask[1]);
    assert!(domain.tissue_mask[2]);
    assert_eq!(domain.bin_id, vec![0, 1, 2]);
}
