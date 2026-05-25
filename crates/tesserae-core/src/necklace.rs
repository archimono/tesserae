//! Booth's O(n) algorithm for lexicographically least circular substrings.
//!
//! Booth, K. S. (1980). Lexicographically least circular substrings.
//! *Inf. Process. Lett.* 10(4-5), 240–242.
//! <https://doi.org/10.1016/0020-0190(80)90149-0>

/// Find the starting index of the lexicographically smallest rotation
/// of `s` in O(n) time.
///
/// Returns `k` such that `s[k..] ++ s[..k]` is the lex-min rotation.
/// For an empty slice, returns 0.
///
/// Uses a two-pointer technique: maintain candidates `i` and `j`,
/// compare them character-by-character, and advance the loser past
/// the mismatch. Total work is O(n) since `i + j + k` increases
/// by at least 1 each iteration.
pub fn lex_min_rotation<T: Ord>(s: &[T]) -> usize {
    let n = s.len();
    if n == 0 {
        return 0;
    }
    let mut i = 0;
    let mut j = 1;
    let mut k = 0;
    while i < n && j < n && k < n {
        let a = &s[(i + k) % n];
        let b = &s[(j + k) % n];
        match a.cmp(b) {
            std::cmp::Ordering::Equal => k += 1,
            std::cmp::Ordering::Less => {
                j = (j + k + 1).max(i + 1);
                k = 0;
            }
            std::cmp::Ordering::Greater => {
                i = (i + k + 1).max(j + 1);
                k = 0;
            }
        }
    }
    i.min(j)
}

/// Like [`lex_min_rotation`] but with a custom comparison function.
///
/// `cmp(i, j)` compares the element at position `i` with the element
/// at position `j` in the circular sequence of length `n`.
pub fn lex_min_rotation_by(n: usize, cmp: impl Fn(usize, usize) -> std::cmp::Ordering) -> usize {
    if n == 0 {
        return 0;
    }
    let mut i = 0;
    let mut j = 1;
    let mut k = 0;
    while i < n && j < n && k < n {
        match cmp((i + k) % n, (j + k) % n) {
            std::cmp::Ordering::Equal => k += 1,
            std::cmp::Ordering::Less => {
                j = (j + k + 1).max(i + 1);
                k = 0;
            }
            std::cmp::Ordering::Greater => {
                i = (i + k + 1).max(j + 1);
                k = 0;
            }
        }
    }
    i.min(j)
}

/// Return the lex-min rotation of `s` as a new `Vec`.
pub fn canonical_rotation<T: Ord + Clone>(s: &[T]) -> Vec<T> {
    let k = lex_min_rotation(s);
    let mut result = Vec::with_capacity(s.len());
    result.extend_from_slice(&s[k..]);
    result.extend_from_slice(&s[..k]);
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty() {
        let s: &[u8] = &[];
        assert_eq!(lex_min_rotation(s), 0);
    }

    #[test]
    fn single_element() {
        assert_eq!(lex_min_rotation(&[5]), 0);
    }

    #[test]
    fn already_minimal() {
        assert_eq!(lex_min_rotation(&[0, 1, 2, 3]), 0);
    }

    #[test]
    fn simple_rotation() {
        // [2, 0, 1] → lex min is [0, 1, 2] starting at index 1
        assert_eq!(lex_min_rotation(&[2, 0, 1]), 1);
    }

    #[test]
    fn all_same() {
        assert_eq!(lex_min_rotation(&[1, 1, 1, 1]), 0);
    }

    #[test]
    fn binary_necklace() {
        // [1, 0, 0, 1, 0] → lex min [0, 0, 1, 0, 1] at index 1
        assert_eq!(lex_min_rotation(&[1, 0, 0, 1, 0]), 1);
        assert_eq!(canonical_rotation(&[1, 0, 0, 1, 0]), vec![0, 0, 1, 0, 1]);
    }

    #[test]
    fn idempotent() {
        let cases: Vec<Vec<u8>> = vec![
            vec![3, 1, 2],
            vec![0, 1, 0, 1],
            vec![2, 2, 1, 2],
            vec![5, 3, 1, 4, 2],
        ];
        for s in &cases {
            let canon = canonical_rotation(s);
            let canon2 = canonical_rotation(&canon);
            assert_eq!(canon, canon2, "not idempotent for {s:?}");
        }
    }

    #[test]
    fn all_rotations_agree() {
        let s = vec![3u8, 1, 4, 1, 5];
        let canon = canonical_rotation(&s);
        let n = s.len();
        for start in 0..n {
            let rotated: Vec<u8> = s[start..]
                .iter()
                .chain(s[..start].iter())
                .copied()
                .collect();
            assert_eq!(
                canonical_rotation(&rotated),
                canon,
                "rotation at {start} gives different canonical form"
            );
        }
    }

    #[test]
    fn brute_force_cross_check() {
        let cases: Vec<Vec<u8>> = vec![
            vec![2, 1, 3],
            vec![1, 3, 2, 1],
            vec![0, 0, 1, 0, 0, 1],
            vec![4, 3, 2, 1, 0],
        ];
        for s in &cases {
            let n = s.len();
            let mut min_rotation = s.clone();
            for start in 1..n {
                let rotated: Vec<u8> = s[start..]
                    .iter()
                    .chain(s[..start].iter())
                    .copied()
                    .collect();
                if rotated < min_rotation {
                    min_rotation = rotated;
                }
            }
            assert_eq!(
                canonical_rotation(s),
                min_rotation,
                "Booth disagrees with brute force for {s:?}"
            );
        }
    }
}
