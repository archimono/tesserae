//! T-canonicalization: lex-min representative under the translation subgroup.
//!
//! 1D uses Booth's O(n) algorithm (see [`crate::necklace`]).
//! 2D/3D use an axis-by-axis approach: Booth on the innermost axis prunes
//! the outer-axis search space. 4D+ falls back to brute force over T.

use crate::necklace::{canonical_rotation, lex_min_rotation};
use crate::perm::Permutation;

/// Return the T-canonical form for cyclic T = ℤ_n using O(n) lex-min rotation.
///
/// Equivalent to `t_canonical(coloring, &cyclic_translations(n))` but runs
/// in O(n) instead of O(n²).
pub fn t_canonical_cyclic(coloring: &[u32]) -> Vec<u32> {
    canonical_rotation(coloring)
}

/// Return the T-canonical form for a product translation group
/// ℤ_{dims[0]} × ℤ_{dims[1]} × ... using axis-by-axis lex-min rotation.
///
/// For 1D, this is O(n) via lex-min rotation.
/// For 2D (dims = [d0, d1]), this is O(d0² × d1) amortized — saves factor d1
/// vs brute force O(d0² × d1²).
/// For 3D+, falls back to brute force.
pub fn t_canonical_product(coloring: &[u32], dims: &[usize]) -> Vec<u32> {
    debug_assert_eq!(coloring.len(), dims.iter().product::<usize>());
    match dims.len() {
        0 => coloring.to_vec(),
        1 => canonical_rotation(coloring),
        2 => t_canonical_product_2d(coloring, dims[0], dims[1]),
        3 => t_canonical_product_3d(coloring, dims[0], dims[1], dims[2]),
        _ => {
            let translations = product_translations(dims);
            t_canonical(coloring, &translations)
        }
    }
}

fn t_canonical_product_2d(coloring: &[u32], d0: usize, d1: usize) -> Vec<u32> {
    let n = d0 * d1;
    let mut best = coloring.to_vec();

    for s0 in 0..d0 {
        let row0 = &coloring[s0 * d1..s0 * d1 + d1];
        let primary_s1 = lex_min_rotation(row0);
        let period = string_period(row0);

        let mut s1 = primary_s1;
        loop {
            let mut candidate = Vec::with_capacity(n);
            for a in 0..d0 {
                let row = (a + s0) % d0;
                for b in 0..d1 {
                    candidate.push(coloring[row * d1 + (b + s1) % d1]);
                }
            }

            if candidate < best {
                best = candidate;
            }

            s1 = (s1 + period) % d1;
            if s1 == primary_s1 {
                break;
            }
        }
    }

    best
}

fn t_canonical_product_3d(coloring: &[u32], d0: usize, d1: usize, d2: usize) -> Vec<u32> {
    let n = d0 * d1 * d2;
    let mut best = coloring.to_vec();

    for s0 in 0..d0 {
        for s1 in 0..d1 {
            let fiber_base = s0 * d1 * d2 + s1 * d2;
            let fiber = &coloring[fiber_base..fiber_base + d2];
            let primary_s2 = lex_min_rotation(fiber);
            let period = string_period(fiber);

            let mut s2 = primary_s2;
            loop {
                let mut candidate = Vec::with_capacity(n);
                for a in 0..d0 {
                    for b in 0..d1 {
                        let src = ((a + s0) % d0) * d1 * d2 + ((b + s1) % d1) * d2;
                        for c in 0..d2 {
                            candidate.push(coloring[src + (c + s2) % d2]);
                        }
                    }
                }

                if candidate < best {
                    best = candidate;
                }

                s2 = (s2 + period) % d2;
                if s2 == primary_s2 {
                    break;
                }
            }
        }
    }

    best
}

fn string_period(s: &[u32]) -> usize {
    let n = s.len();
    for p in 1..n {
        if !n.is_multiple_of(p) {
            continue;
        }
        if (0..n).all(|i| s[i] == s[i % p]) {
            return p;
        }
    }
    n
}

/// Generate multi-site translation permutations for
/// ℤ_{dims[0]} × ℤ_{dims[1]} × ... acting on `n_orbits * product(dims)` sites.
///
/// Each translation shifts ALL orbit blocks by the same amount.
/// Flat indexing: site `(j, y)` → `j * V + row_major_index(y)` where `V = product(dims)`.
pub fn product_translations_multisite(dims: &[usize], n_orbits: usize) -> Vec<Permutation> {
    let single = product_translations(dims);
    let v: usize = dims.iter().product();
    let n_total = n_orbits * v;

    single
        .into_iter()
        .map(|t| {
            let image: Vec<u32> = (0..n_total)
                .map(|flat| {
                    let orbit = flat / v;
                    let site = flat % v;
                    (orbit * v + t.apply(site as u32) as usize) as u32
                })
                .collect();
            Permutation::new(image)
        })
        .collect()
}

/// Generate all translation permutations for a product group
/// ℤ_{dims[0]} × ℤ_{dims[1]} × ... acting on `product(dims)` sites.
///
/// Sites are indexed in row-major order.
pub fn product_translations(dims: &[usize]) -> Vec<Permutation> {
    let n: usize = dims.iter().product();
    if n == 0 {
        return vec![Permutation::identity(0)];
    }

    let strides: Vec<usize> = {
        let mut s = vec![1usize; dims.len()];
        for k in (0..dims.len() - 1).rev() {
            s[k] = s[k + 1] * dims[k + 1];
        }
        s
    };

    let num_translations: usize = dims.iter().product();
    let mut translations = Vec::with_capacity(num_translations);

    let mut shift = vec![0usize; dims.len()];
    loop {
        let image: Vec<u32> = (0..n)
            .map(|idx| {
                let mut new_idx = 0;
                let mut remaining = idx;
                for k in 0..dims.len() {
                    let coord = remaining / strides[k];
                    remaining %= strides[k];
                    new_idx += ((coord + shift[k]) % dims[k]) * strides[k];
                }
                new_idx as u32
            })
            .collect();
        translations.push(Permutation::new(image));

        let mut carry = true;
        for k in (0..dims.len()).rev() {
            if carry {
                shift[k] += 1;
                if shift[k] >= dims[k] {
                    shift[k] = 0;
                } else {
                    carry = false;
                }
            }
        }
        if carry {
            break;
        }
    }

    translations
}

/// Return the lexicographically smallest coloring reachable from `coloring`
/// under any translation in `translations`.
///
/// `translations` should contain all elements of the translation subgroup T
/// (including the identity). The result is the T-canonical representative.
pub fn t_canonical(coloring: &[u32], translations: &[Permutation]) -> Vec<u32> {
    debug_assert!(
        !translations.is_empty(),
        "translation group must be non-empty"
    );

    let n = coloring.len();
    let mut best = coloring.to_vec();

    for t in translations {
        debug_assert_eq!(t.len(), n);
        let mut is_better = false;
        for i in 0..n {
            let val = coloring[t.apply(i as u32) as usize];
            match val.cmp(&best[i]) {
                std::cmp::Ordering::Less => {
                    is_better = true;
                    break;
                }
                std::cmp::Ordering::Greater => break,
                std::cmp::Ordering::Equal => {}
            }
        }
        if is_better {
            for i in 0..n {
                best[i] = coloring[t.apply(i as u32) as usize];
            }
        }
    }

    best
}

/// Generate all translation permutations for a 1D cyclic group ℤ_n.
pub fn cyclic_translations(n: usize) -> Vec<Permutation> {
    (0..n)
        .map(|shift| {
            let image: Vec<u32> = (0..n).map(|i| ((i + shift) % n) as u32).collect();
            Permutation::new(image)
        })
        .collect()
}

/// Generate all translation permutations for a product group
/// ℤ_{d0} × ℤ_{d1} acting on d0*d1 sites.
///
/// Sites are indexed in row-major order: site `(a, b)` maps to index `a * d1 + b`.
pub fn product_translations_2d(d0: usize, d1: usize) -> Vec<Permutation> {
    let n = d0 * d1;
    let mut translations = Vec::with_capacity(n);
    for s0 in 0..d0 {
        for s1 in 0..d1 {
            let image: Vec<u32> = (0..d0)
                .flat_map(|a| {
                    (0..d1).map(move |b| {
                        let new_a = (a + s0) % d0;
                        let new_b = (b + s1) % d1;
                        (new_a * d1 + new_b) as u32
                    })
                })
                .collect();
            debug_assert_eq!(image.len(), n);
            translations.push(Permutation::new(image));
        }
    }
    translations
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polya::polya_count_binary;
    use std::collections::HashSet;

    #[test]
    fn identity_only() {
        let t = vec![Permutation::identity(4)];
        let coloring = vec![1, 0, 1, 0];
        assert_eq!(t_canonical(&coloring, &t), coloring);
    }

    #[test]
    fn idempotent_cyclic() {
        let translations = cyclic_translations(6);
        let colorings: Vec<Vec<u32>> = vec![
            vec![1, 0, 0, 1, 0, 0],
            vec![0, 1, 0, 1, 0, 1],
            vec![1, 1, 0, 0, 0, 0],
        ];
        for c in &colorings {
            let canon = t_canonical(c, &translations);
            let canon2 = t_canonical(&canon, &translations);
            assert_eq!(canon, canon2, "not idempotent for {c:?}");
        }
    }

    #[test]
    fn cyclic_orbit_count_matches_polya() {
        for n in [4, 6, 8, 10] {
            let translations = cyclic_translations(n);
            for n0 in 0..=n {
                let expected = polya_count_binary(&translations, n0, n - n0);
                let actual = count_distinct_canonicals(n, n0, &translations);
                assert_eq!(
                    actual,
                    expected as usize,
                    "C_{n}: orbit mismatch at composition ({n0},{})",
                    n - n0
                );
            }
        }
    }

    #[test]
    fn trivial_group_orbit_count() {
        let n = 5;
        let translations = vec![Permutation::identity(n)];
        for n0 in 0..=n {
            let expected = polya_count_binary(&translations, n0, n - n0);
            let actual = count_distinct_canonicals(n, n0, &translations);
            assert_eq!(actual, expected as usize);
        }
    }

    #[test]
    fn product_2d_orbit_count_matches_polya() {
        // ℤ_2 × ℤ_3 on 6 sites
        let translations = product_translations_2d(2, 3);
        assert_eq!(translations.len(), 6);
        for n0 in 0..=6 {
            let expected = polya_count_binary(&translations, n0, 6 - n0);
            let actual = count_distinct_canonicals(6, n0, &translations);
            assert_eq!(
                actual,
                expected as usize,
                "ℤ_2×ℤ_3: orbit mismatch at composition ({n0},{})",
                6 - n0
            );
        }
    }

    #[test]
    fn product_3x3_orbit_count_matches_polya() {
        // ℤ_3 × ℤ_3 on 9 sites
        let translations = product_translations_2d(3, 3);
        assert_eq!(translations.len(), 9);
        for n0 in 0..=9 {
            let expected = polya_count_binary(&translations, n0, 9 - n0);
            let actual = count_distinct_canonicals(9, n0, &translations);
            assert_eq!(
                actual,
                expected as usize,
                "ℤ_3×ℤ_3: orbit mismatch at composition ({n0},{})",
                9 - n0
            );
        }
    }

    /// Enumerate all binary colorings with composition (n0, n-n0),
    /// canonicalize each, and count distinct results.
    fn count_distinct_canonicals(n: usize, n0: usize, translations: &[Permutation]) -> usize {
        let mut seen = HashSet::new();
        for_each_binary_coloring(n, n0, &mut |c| {
            seen.insert(t_canonical(c, translations));
        });
        seen.len()
    }

    #[test]
    fn booth_agrees_with_brute_force() {
        for n in [4, 6, 8, 10, 12] {
            let translations = cyclic_translations(n);
            for n0 in 0..=n {
                for_each_binary_coloring(n, n0, &mut |c| {
                    let brute = t_canonical(c, &translations);
                    let booth = t_canonical_cyclic(c);
                    assert_eq!(brute, booth, "mismatch at n={n}, coloring={c:?}");
                });
            }
        }
    }

    #[test]
    fn product_1d_agrees_with_cyclic() {
        for n in [4, 6, 8] {
            for n0 in 0..=n {
                for_each_binary_coloring(n, n0, &mut |c| {
                    let cyclic = t_canonical_cyclic(c);
                    let product = t_canonical_product(c, &[n]);
                    assert_eq!(
                        cyclic, product,
                        "1D product mismatch at n={n}, coloring={c:?}"
                    );
                });
            }
        }
    }

    #[test]
    fn product_2d_agrees_with_brute_force() {
        for (d0, d1) in [(2, 3), (3, 3), (2, 4), (3, 4)] {
            let n = d0 * d1;
            let translations = product_translations_2d(d0, d1);
            for n0 in 0..=n {
                for_each_binary_coloring(n, n0, &mut |c| {
                    let brute = t_canonical(c, &translations);
                    let product = t_canonical_product(c, &[d0, d1]);
                    assert_eq!(
                        brute, product,
                        "2D product mismatch at {d0}x{d1}, coloring={c:?}"
                    );
                });
            }
        }
    }

    #[test]
    fn product_3d_agrees_with_brute_force() {
        for (d0, d1, d2) in [(2, 2, 2), (1, 2, 3), (2, 2, 3)] {
            let n = d0 * d1 * d2;
            let translations = product_translations(&[d0, d1, d2]);
            for n0 in 0..=n {
                for_each_binary_coloring(n, n0, &mut |c| {
                    let brute = t_canonical(c, &translations);
                    let product = t_canonical_product(c, &[d0, d1, d2]);
                    assert_eq!(
                        brute, product,
                        "3D product mismatch at {d0}x{d1}x{d2}, coloring={c:?}"
                    );
                });
            }
        }
    }

    #[test]
    fn product_translations_agrees_with_2d() {
        for (d0, d1) in [(2, 3), (3, 3)] {
            let specific = product_translations_2d(d0, d1);
            let general = product_translations(&[d0, d1]);
            assert_eq!(specific.len(), general.len());
            let specific_set: HashSet<_> = specific.into_iter().collect();
            let general_set: HashSet<_> = general.into_iter().collect();
            assert_eq!(specific_set, general_set, "mismatch at {d0}x{d1}");
        }
    }

    /// Iterate over all binary colorings of length `n` with exactly `n0` zeros.
    fn for_each_binary_coloring(n: usize, n0: usize, f: &mut impl FnMut(&[u32])) {
        let mut coloring = vec![0u32; n];
        for c in coloring.iter_mut().skip(n0) {
            *c = 1;
        }
        loop {
            f(&coloring);
            if !next_permutation(&mut coloring) {
                break;
            }
        }
    }

    fn next_permutation(arr: &mut [u32]) -> bool {
        let n = arr.len();
        if n <= 1 {
            return false;
        }
        let mut i = n - 1;
        while i > 0 && arr[i - 1] >= arr[i] {
            i -= 1;
        }
        if i == 0 {
            return false;
        }
        let mut j = n - 1;
        while arr[j] <= arr[i - 1] {
            j -= 1;
        }
        arr.swap(i - 1, j);
        arr[i..].reverse();
        true
    }
}
