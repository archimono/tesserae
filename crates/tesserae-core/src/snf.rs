/// Result of a Smith Normal Form decomposition: D = L · M · R,
/// where D = diag(d0, d1, ...) with d0 | d1 | d2 | ...
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmithNormalForm {
    /// Diagonal entries (invariant factors), with d[i] | d[i+1].
    pub diagonal: Vec<i64>,
    /// Left transformation matrix (row operations).
    pub left: Vec<Vec<i64>>,
    /// Right transformation matrix (column operations).
    pub right: Vec<Vec<i64>>,
}

/// Compute the Smith Normal Form of an n×n integer matrix.
///
/// Returns `SmithNormalForm { diagonal, left, right }` such that
/// `left * matrix * right = diag(diagonal)` and `diagonal[i]` divides
/// `diagonal[i+1]`. All diagonal entries are non-negative.
///
/// # Panics
/// Panics if `matrix` is not square.
pub fn smith_normal_form(matrix: &[Vec<i64>]) -> SmithNormalForm {
    let n = matrix.len();
    assert!(n > 0);
    for row in matrix {
        assert_eq!(row.len(), n, "matrix must be square");
    }

    let mut d: Vec<Vec<i64>> = matrix.to_vec();
    let mut left = identity_matrix(n);
    let mut right = identity_matrix(n);

    for pivot in 0..n {
        if !eliminate_pivot(&mut d, &mut left, &mut right, pivot, n) {
            continue;
        }
    }

    // Ensure diagonal entries are non-negative.
    for i in 0..n {
        if d[i][i] < 0 {
            d[i][i] = -d[i][i];
            for val in &mut left[i] {
                *val = -*val;
            }
        }
    }

    // Enforce divisibility chain: d[i] | d[i+1].
    enforce_divisibility(&mut d, &mut left, &mut right, n);

    let diagonal = (0..n).map(|i| d[i][i]).collect();
    SmithNormalForm {
        diagonal,
        left,
        right,
    }
}

fn eliminate_pivot(
    d: &mut [Vec<i64>],
    left: &mut [Vec<i64>],
    right: &mut [Vec<i64>],
    pivot: usize,
    n: usize,
) -> bool {
    // Find a nonzero entry in the submatrix d[pivot..n][pivot..n].
    let mut found = false;
    'search: for i in pivot..n {
        for j in pivot..n {
            if d[i][j] != 0 {
                if i != pivot {
                    swap_rows(d, left, i, pivot);
                }
                if j != pivot {
                    swap_cols(d, right, j, pivot);
                }
                found = true;
                break 'search;
            }
        }
    }
    if !found {
        return false;
    }

    loop {
        let mut changed = false;

        // Eliminate column entries below pivot.
        for i in (pivot + 1)..n {
            if d[i][pivot] != 0 {
                let q = d[i][pivot] / d[pivot][pivot];
                add_row_multiple(d, left, i, pivot, -q, n);
                if d[i][pivot] != 0 {
                    swap_rows(d, left, i, pivot);
                }
                changed = true;
            }
        }

        // Eliminate row entries right of pivot.
        for j in (pivot + 1)..n {
            if d[pivot][j] != 0 {
                let q = d[pivot][j] / d[pivot][pivot];
                add_col_multiple(d, right, j, pivot, -q, n);
                if d[pivot][j] != 0 {
                    swap_cols(d, right, j, pivot);
                }
                changed = true;
            }
        }

        if !changed {
            break;
        }
    }

    true
}

fn enforce_divisibility(
    d: &mut [Vec<i64>],
    left: &mut [Vec<i64>],
    right: &mut [Vec<i64>],
    n: usize,
) {
    for i in 0..n.saturating_sub(1) {
        if d[i][i] == 0 {
            continue;
        }
        if d[i + 1][i + 1] == 0 {
            continue;
        }
        if d[i + 1][i + 1] % d[i][i] != 0 {
            // Add row i+1 to row i, then re-eliminate.
            add_row_multiple(d, left, i, i + 1, 1, n);
            eliminate_pivot(d, left, right, i, n);

            if d[i][i] < 0 {
                d[i][i] = -d[i][i];
                for val in &mut left[i] {
                    *val = -*val;
                }
            }
            if d[i + 1][i + 1] < 0 {
                d[i + 1][i + 1] = -d[i + 1][i + 1];
                for val in &mut left[i + 1] {
                    *val = -*val;
                }
            }

            // Restart divisibility check from beginning.
            enforce_divisibility(d, left, right, n);
            return;
        }
    }
}

fn identity_matrix(n: usize) -> Vec<Vec<i64>> {
    let mut m = vec![vec![0i64; n]; n];
    for (i, row) in m.iter_mut().enumerate() {
        row[i] = 1;
    }
    m
}

fn swap_rows(d: &mut [Vec<i64>], left: &mut [Vec<i64>], i: usize, j: usize) {
    d.swap(i, j);
    left.swap(i, j);
}

fn swap_cols(d: &mut [Vec<i64>], right: &mut [Vec<i64>], i: usize, j: usize) {
    for row in d.iter_mut() {
        row.swap(i, j);
    }
    for row in right.iter_mut() {
        row.swap(i, j);
    }
}

fn add_row_multiple(
    d: &mut [Vec<i64>],
    left: &mut [Vec<i64>],
    target: usize,
    source: usize,
    factor: i64,
    n: usize,
) {
    for j in 0..n {
        d[target][j] += factor * d[source][j];
        left[target][j] += factor * left[source][j];
    }
}

fn add_col_multiple(
    d: &mut [Vec<i64>],
    right: &mut [Vec<i64>],
    target: usize,
    source: usize,
    factor: i64,
    _n: usize,
) {
    for row in d.iter_mut() {
        row[target] += factor * row[source];
    }
    for row in right.iter_mut() {
        row[target] += factor * row[source];
    }
}

#[cfg(test)]
fn mat_mul(a: &[Vec<i64>], b: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let n = a.len();
    let mut c = vec![vec![0i64; n]; n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    fn diag_matrix(diagonal: &[i64]) -> Vec<Vec<i64>> {
        let n = diagonal.len();
        let mut m = vec![vec![0i64; n]; n];
        for i in 0..n {
            m[i][i] = diagonal[i];
        }
        m
    }

    fn verify_snf(matrix: &[Vec<i64>], snf: &SmithNormalForm) {
        let n = matrix.len();

        // L * M * R == D
        let lm = mat_mul(&snf.left, matrix);
        let lmr = mat_mul(&lm, &snf.right);
        let expected = diag_matrix(&snf.diagonal);
        assert_eq!(lmr, expected, "L * M * R != diag(D)");

        // All diagonal entries non-negative.
        for &d in &snf.diagonal {
            assert!(d >= 0, "negative diagonal entry: {d}");
        }

        // Divisibility chain.
        for i in 0..n.saturating_sub(1) {
            if snf.diagonal[i] == 0 {
                continue;
            }
            assert_eq!(
                snf.diagonal[i + 1] % snf.diagonal[i],
                0,
                "divisibility violated: {} does not divide {}",
                snf.diagonal[i],
                snf.diagonal[i + 1]
            );
        }
    }

    #[test]
    fn identity_2x2() {
        let m = vec![vec![1, 0], vec![0, 1]];
        let snf = smith_normal_form(&m);
        assert_eq!(snf.diagonal, vec![1, 1]);
        verify_snf(&m, &snf);
    }

    #[test]
    fn identity_3x3() {
        let m = vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
        let snf = smith_normal_form(&m);
        assert_eq!(snf.diagonal, vec![1, 1, 1]);
        verify_snf(&m, &snf);
    }

    #[test]
    fn diagonal_2x2() {
        let m = vec![vec![6, 0], vec![0, 4]];
        let snf = smith_normal_form(&m);
        assert_eq!(snf.diagonal, vec![2, 12]);
        verify_snf(&m, &snf);
    }

    #[test]
    fn supercell_2x2() {
        // A 2×2 supercell matrix with det=4
        let m = vec![vec![2, 0], vec![0, 2]];
        let snf = smith_normal_form(&m);
        assert_eq!(snf.diagonal, vec![2, 2]);
        verify_snf(&m, &snf);
    }

    #[test]
    fn off_diagonal_2x2() {
        let m = vec![vec![2, 4], vec![1, 3]];
        let snf = smith_normal_form(&m);
        // det = 6-4 = 2, so diagonal product = 2
        assert_eq!(snf.diagonal.iter().product::<i64>(), 2);
        verify_snf(&m, &snf);
    }

    #[test]
    fn supercell_3x3_diagonal() {
        let m = vec![vec![2, 0, 0], vec![0, 3, 0], vec![0, 0, 5]];
        let snf = smith_normal_form(&m);
        assert_eq!(snf.diagonal, vec![1, 1, 30]);
        verify_snf(&m, &snf);
    }

    #[test]
    fn supercell_3x3_index_8() {
        // FCC conventional cell: index 8 supercell
        let m = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 2]];
        let snf = smith_normal_form(&m);
        assert_eq!(snf.diagonal, vec![2, 2, 2]);
        verify_snf(&m, &snf);
    }

    #[test]
    fn non_diagonal_3x3() {
        let m = vec![vec![2, 4, 4], vec![-6, 6, 12], vec![10, -4, -16]];
        let snf = smith_normal_form(&m);
        verify_snf(&m, &snf);
        // det = 2*(6*(-16)-12*(-4)) - 4*((-6)*(-16)-12*10) + 4*((-6)*(-4)-6*10)
        //     = 2*(-96+48) - 4*(96-120) + 4*(24-60)
        //     = 2*(-48) - 4*(-24) + 4*(-36)
        //     = -96 + 96 - 144 = -144
        // |det| = 144, product of diagonal must be 144
        assert_eq!(snf.diagonal.iter().product::<i64>(), 144);
    }

    #[test]
    fn zero_matrix_2x2() {
        let m = vec![vec![0, 0], vec![0, 0]];
        let snf = smith_normal_form(&m);
        assert_eq!(snf.diagonal, vec![0, 0]);
        verify_snf(&m, &snf);
    }

    #[test]
    fn singular_3x3() {
        let m = vec![vec![1, 2, 3], vec![4, 5, 6], vec![7, 8, 9]];
        let snf = smith_normal_form(&m);
        verify_snf(&m, &snf);
        // Rank 2, so last diagonal is 0
        assert_eq!(snf.diagonal[2], 0);
    }
}
