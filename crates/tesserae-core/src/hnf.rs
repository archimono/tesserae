//! Hermite Normal Form enumeration for derivative structure generation.
//!
//! Hart, G. L. W.; Forcade, R. W. (2008). Algorithm for generating
//! derivative structures. *Phys. Rev. B* 77, 224115.
//! <https://doi.org/10.1103/PhysRevB.77.224115>

/// Enumerate all 2×2 Hermite Normal Form matrices with the given determinant.
///
/// Returns lower-triangular matrices `[[a, 0], [b, c]]` where
/// `a * c = det`, `0 ≤ b < c`, following the Hart–Forcade convention.
/// The count equals σ(det), the sum of divisors of det.
pub fn hnf_2d(det: usize) -> Vec<Vec<Vec<i64>>> {
    let mut results = Vec::new();
    for a in 1..=det {
        if !det.is_multiple_of(a) {
            continue;
        }
        let c = (det / a) as i64;
        let a = a as i64;
        for b in 0..c {
            results.push(vec![vec![a, 0], vec![b, c]]);
        }
    }
    results
}

/// Enumerate all 3×3 Hermite Normal Form matrices with the given determinant.
///
/// Returns lower-triangular matrices `[[a,0,0],[b,c,0],[d,e,f]]` where
/// `a * c * f = det`, `0 ≤ b < c`, `0 ≤ d < f`, `0 ≤ e < f`,
/// following the Hart–Forcade convention.
pub fn hnf_3d(det: usize) -> Vec<Vec<Vec<i64>>> {
    let mut results = Vec::new();
    for a in 1..=det {
        if !det.is_multiple_of(a) {
            continue;
        }
        let cf = det / a;
        for c in 1..=cf {
            if !cf.is_multiple_of(c) {
                continue;
            }
            let f = cf / c;
            for b in 0..c {
                for d in 0..f {
                    for e in 0..f {
                        results.push(vec![
                            vec![a as i64, 0, 0],
                            vec![b as i64, c as i64, 0],
                            vec![d as i64, e as i64, f as i64],
                        ]);
                    }
                }
            }
        }
    }
    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    fn sigma(n: usize) -> usize {
        (1..=n).filter(|d| n.is_multiple_of(*d)).sum()
    }

    fn expected_hnf_3d_count(n: usize) -> usize {
        (1..=n)
            .filter(|d| n.is_multiple_of(*d))
            .map(|d| d * sigma(d))
            .sum()
    }

    fn det_2x2(m: &[Vec<i64>]) -> i64 {
        m[0][0] * m[1][1] - m[0][1] * m[1][0]
    }

    fn det_3x3(m: &[Vec<i64>]) -> i64 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    #[test]
    fn count_2d_matches_sigma() {
        for det in 1..=20 {
            assert_eq!(
                hnf_2d(det).len(),
                sigma(det),
                "2D HNF count mismatch at det={det}"
            );
        }
    }

    #[test]
    fn count_3d_matches_formula() {
        for det in 1..=12 {
            assert_eq!(
                hnf_3d(det).len(),
                expected_hnf_3d_count(det),
                "3D HNF count mismatch at det={det}"
            );
        }
    }

    #[test]
    fn count_3d_known_values() {
        let known = [(1, 1), (2, 7), (3, 13), (4, 35), (5, 31), (6, 91)];
        for (det, expected) in known {
            assert_eq!(hnf_3d(det).len(), expected, "3D count for det={det}");
        }
    }

    #[test]
    fn determinants_2d() {
        for det in 1..=10 {
            for m in hnf_2d(det) {
                assert_eq!(det_2x2(&m), det as i64);
            }
        }
    }

    #[test]
    fn determinants_3d() {
        for det in 1..=8 {
            for m in hnf_3d(det) {
                assert_eq!(det_3x3(&m), det as i64);
            }
        }
    }

    #[test]
    fn constraints_2d() {
        for det in 1..=10 {
            for m in hnf_2d(det) {
                assert!(m[0][0] >= 1);
                assert_eq!(m[0][1], 0);
                assert!(m[1][0] >= 0 && m[1][0] < m[1][1]);
                assert!(m[1][1] >= 1);
            }
        }
    }

    #[test]
    fn constraints_3d() {
        for det in 1..=6 {
            for m in hnf_3d(det) {
                assert!(m[0][0] >= 1);
                assert!(m[1][1] >= 1);
                assert!(m[2][2] >= 1);
                assert_eq!(m[0][1], 0);
                assert_eq!(m[0][2], 0);
                assert_eq!(m[1][2], 0);
                assert!(m[1][0] >= 0 && m[1][0] < m[1][1]);
                assert!(m[2][0] >= 0 && m[2][0] < m[2][2]);
                assert!(m[2][1] >= 0 && m[2][1] < m[2][2]);
            }
        }
    }

    #[test]
    fn no_duplicates_2d() {
        for det in 1..=10 {
            let matrices = hnf_2d(det);
            let unique: HashSet<Vec<Vec<i64>>> = matrices.iter().cloned().collect();
            assert_eq!(matrices.len(), unique.len(), "duplicate at det={det}");
        }
    }

    #[test]
    fn no_duplicates_3d() {
        for det in 1..=6 {
            let matrices = hnf_3d(det);
            let unique: HashSet<Vec<Vec<i64>>> = matrices.iter().cloned().collect();
            assert_eq!(matrices.len(), unique.len(), "duplicate at det={det}");
        }
    }
}
