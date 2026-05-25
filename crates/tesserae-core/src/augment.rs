//! McKay-style canonical augmentation for isomorph-free exhaustive generation.
//!
//! McKay, B. D. (1998). Isomorph-free exhaustive generation.
//! *J. Algorithms* 26(2), 306–324.
//! <https://doi.org/10.1006/jagm.1997.0898>

use crate::perm::Permutation;

/// Enumerate all canonical colorings under `group` with the given `composition`.
///
/// `composition[k]` is the number of sites to color with color `k`.
/// `group` must contain the identity. Calls `callback` with each canonical
/// coloring exactly once.
///
/// Uses McKay's canonical augmentation: colorings are built site-by-site,
/// maintaining an active set of group elements that could still make the
/// partial coloring lex-smaller. A branch is pruned as soon as any active
/// element produces a smaller image.
///
/// # Panics
///
/// Panics (debug) if `group` is empty or its elements have length
/// different from `composition.iter().sum()`.
pub fn enumerate_canonical(
    composition: &[usize],
    group: &[Permutation],
    callback: &mut impl FnMut(&[u32]),
) {
    let n: usize = composition.iter().sum();
    debug_assert!(!group.is_empty());
    debug_assert!(group.iter().all(|g| g.len() == n));

    let active: Vec<usize> = (0..group.len()).collect();
    let mut partial = Vec::with_capacity(n);
    let mut remaining = composition.to_vec();

    augment(&mut partial, &mut remaining, &active, group, n, callback);
}

fn augment(
    partial: &mut Vec<u32>,
    remaining: &mut [usize],
    active: &[usize],
    group: &[Permutation],
    n: usize,
    callback: &mut impl FnMut(&[u32]),
) {
    if partial.len() == n {
        callback(partial);
        return;
    }

    for color in 0..remaining.len() {
        if remaining[color] == 0 {
            continue;
        }

        partial.push(color as u32);
        let level = partial.len() - 1;

        let mut new_active = Vec::new();
        let mut canonical = true;

        for &g_idx in active {
            match partial_compare(partial, &group[g_idx], level) {
                PartialCmp::Smaller => {
                    canonical = false;
                    break;
                }
                PartialCmp::Larger => {}
                PartialCmp::Undecided => {
                    new_active.push(g_idx);
                }
            }
        }

        if canonical {
            remaining[color] -= 1;
            augment(partial, remaining, &new_active, group, n, callback);
            remaining[color] += 1;
        }

        partial.pop();
    }
}

/// Enumerate colorings that are canonical under `point_group` and pruned
/// against `translations`, using dual active sets.
///
/// At each tree level, branches are pruned if any P-element or any T-element
/// produces a lex-smaller partial image. This eliminates most non-G-canonical
/// subtrees without ever building the full group G = T ⋊ P.
///
/// Complete colorings that survive pruning are passed to `callback`.
/// Not all survivors are G-canonical — the caller must apply a final
/// `is_g_canonical` check.
pub fn enumerate_t_pruned(
    composition: &[usize],
    point_group: &[Permutation],
    translations: &[Permutation],
    callback: &mut impl FnMut(&[u32]),
) {
    let n: usize = composition.iter().sum();
    debug_assert!(!point_group.is_empty());
    debug_assert!(!translations.is_empty());
    debug_assert!(point_group.iter().all(|g| g.len() == n));
    debug_assert!(translations.iter().all(|g| g.len() == n));

    let p_active: Vec<usize> = (0..point_group.len()).collect();
    let t_active: Vec<usize> = (0..translations.len()).collect();
    let mut partial = Vec::with_capacity(n);
    let mut remaining = composition.to_vec();

    augment_t_pruned(
        &mut partial,
        &mut remaining,
        &p_active,
        &t_active,
        point_group,
        translations,
        n,
        callback,
    );
}

#[allow(clippy::too_many_arguments)]
fn augment_t_pruned(
    partial: &mut Vec<u32>,
    remaining: &mut [usize],
    p_active: &[usize],
    t_active: &[usize],
    point_group: &[Permutation],
    translations: &[Permutation],
    n: usize,
    callback: &mut impl FnMut(&[u32]),
) {
    if partial.len() == n {
        callback(partial);
        return;
    }

    for color in 0..remaining.len() {
        if remaining[color] == 0 {
            continue;
        }

        partial.push(color as u32);
        let level = partial.len() - 1;

        let mut new_p_active = Vec::new();
        let mut canonical = true;

        for &g_idx in p_active {
            match partial_compare(partial, &point_group[g_idx], level) {
                PartialCmp::Smaller => {
                    canonical = false;
                    break;
                }
                PartialCmp::Larger => {}
                PartialCmp::Undecided => {
                    new_p_active.push(g_idx);
                }
            }
        }

        if canonical {
            let mut new_t_active = Vec::new();

            for &t_idx in t_active {
                match partial_compare(partial, &translations[t_idx], level) {
                    PartialCmp::Smaller => {
                        canonical = false;
                        break;
                    }
                    PartialCmp::Larger => {}
                    PartialCmp::Undecided => {
                        new_t_active.push(t_idx);
                    }
                }
            }

            if canonical {
                remaining[color] -= 1;
                augment_t_pruned(
                    partial,
                    remaining,
                    &new_p_active,
                    &new_t_active,
                    point_group,
                    translations,
                    n,
                    callback,
                );
                remaining[color] += 1;
            }
        }

        partial.pop();
    }
}

/// Like [`enumerate_canonical`] but only places colors allowed at each site.
///
/// `allowed[i][c]` must be `true` iff color `c` may be placed on site `i`.
/// The allowed sets must be G-equivariant (see [`polya_count_constrained`]).
pub fn enumerate_canonical_constrained(
    composition: &[usize],
    group: &[Permutation],
    allowed: &[Vec<bool>],
    callback: &mut impl FnMut(&[u32]),
) {
    let n: usize = composition.iter().sum();
    debug_assert!(!group.is_empty());
    debug_assert!(group.iter().all(|g| g.len() == n));
    debug_assert_eq!(allowed.len(), n);

    let active: Vec<usize> = (0..group.len()).collect();
    let mut partial = Vec::with_capacity(n);
    let mut remaining = composition.to_vec();

    augment_constrained(
        &mut partial,
        &mut remaining,
        &active,
        group,
        allowed,
        n,
        callback,
    );
}

#[allow(clippy::too_many_arguments)]
fn augment_constrained(
    partial: &mut Vec<u32>,
    remaining: &mut [usize],
    active: &[usize],
    group: &[Permutation],
    allowed: &[Vec<bool>],
    n: usize,
    callback: &mut impl FnMut(&[u32]),
) {
    if partial.len() == n {
        callback(partial);
        return;
    }

    let level = partial.len();

    for color in 0..remaining.len() {
        if remaining[color] == 0 || !allowed[level][color] {
            continue;
        }

        partial.push(color as u32);

        let mut new_active = Vec::new();
        let mut canonical = true;

        for &g_idx in active {
            match partial_compare(partial, &group[g_idx], level) {
                PartialCmp::Smaller => {
                    canonical = false;
                    break;
                }
                PartialCmp::Larger => {}
                PartialCmp::Undecided => {
                    new_active.push(g_idx);
                }
            }
        }

        if canonical {
            remaining[color] -= 1;
            augment_constrained(partial, remaining, &new_active, group, allowed, n, callback);
            remaining[color] += 1;
        }

        partial.pop();
    }
}

/// Like [`enumerate_t_pruned`] but only places colors allowed at each site.
pub fn enumerate_t_pruned_constrained(
    composition: &[usize],
    point_group: &[Permutation],
    translations: &[Permutation],
    allowed: &[Vec<bool>],
    callback: &mut impl FnMut(&[u32]),
) {
    let n: usize = composition.iter().sum();
    debug_assert!(!point_group.is_empty());
    debug_assert!(!translations.is_empty());
    debug_assert!(point_group.iter().all(|g| g.len() == n));
    debug_assert!(translations.iter().all(|g| g.len() == n));
    debug_assert_eq!(allowed.len(), n);

    let p_active: Vec<usize> = (0..point_group.len()).collect();
    let t_active: Vec<usize> = (0..translations.len()).collect();
    let mut partial = Vec::with_capacity(n);
    let mut remaining = composition.to_vec();

    augment_t_pruned_constrained(
        &mut partial,
        &mut remaining,
        &p_active,
        &t_active,
        point_group,
        translations,
        allowed,
        n,
        callback,
    );
}

#[allow(clippy::too_many_arguments)]
fn augment_t_pruned_constrained(
    partial: &mut Vec<u32>,
    remaining: &mut [usize],
    p_active: &[usize],
    t_active: &[usize],
    point_group: &[Permutation],
    translations: &[Permutation],
    allowed: &[Vec<bool>],
    n: usize,
    callback: &mut impl FnMut(&[u32]),
) {
    if partial.len() == n {
        callback(partial);
        return;
    }

    let level = partial.len();

    for color in 0..remaining.len() {
        if remaining[color] == 0 || !allowed[level][color] {
            continue;
        }

        partial.push(color as u32);

        let mut new_p_active = Vec::new();
        let mut canonical = true;

        for &g_idx in p_active {
            match partial_compare(partial, &point_group[g_idx], level) {
                PartialCmp::Smaller => {
                    canonical = false;
                    break;
                }
                PartialCmp::Larger => {}
                PartialCmp::Undecided => {
                    new_p_active.push(g_idx);
                }
            }
        }

        if canonical {
            let mut new_t_active = Vec::new();

            for &t_idx in t_active {
                match partial_compare(partial, &translations[t_idx], level) {
                    PartialCmp::Smaller => {
                        canonical = false;
                        break;
                    }
                    PartialCmp::Larger => {}
                    PartialCmp::Undecided => {
                        new_t_active.push(t_idx);
                    }
                }
            }

            if canonical {
                remaining[color] -= 1;
                augment_t_pruned_constrained(
                    partial,
                    remaining,
                    &new_p_active,
                    &new_t_active,
                    point_group,
                    translations,
                    allowed,
                    n,
                    callback,
                );
                remaining[color] += 1;
            }
        }

        partial.pop();
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PartialCmp {
    Smaller,
    Larger,
    Undecided,
}

/// Compare partial coloring `c` with its image `g·c` where `(g·c)[i] = c[g(i)]`.
///
/// Scans positions `0..=level`. If `g(i)` is beyond the assigned range, the
/// comparison is unresolvable and we return `Undecided`. Otherwise, the first
/// position where `c[g(i)] ≠ c[i]` determines the result.
fn partial_compare(partial: &[u32], g: &Permutation, level: usize) -> PartialCmp {
    for i in 0..=level {
        let j = g.apply(i as u32) as usize;
        if j > level {
            return PartialCmp::Undecided;
        }
        match partial[j].cmp(&partial[i]) {
            std::cmp::Ordering::Less => return PartialCmp::Smaller,
            std::cmp::Ordering::Greater => return PartialCmp::Larger,
            std::cmp::Ordering::Equal => {}
        }
    }
    PartialCmp::Undecided
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::perm::generate_group;
    use crate::polya::polya_count_binary;
    use std::collections::HashSet;

    fn count_canonical(composition: &[usize], group: &[Permutation]) -> usize {
        let mut count = 0;
        enumerate_canonical(composition, group, &mut |_| count += 1);
        count
    }

    #[test]
    fn trivial_group() {
        let group = vec![Permutation::identity(4)];
        assert_eq!(count_canonical(&[2, 2], &group), 6);
    }

    #[test]
    fn cyclic_matches_polya() {
        for n in [3, 4, 5, 6, 8] {
            let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
            let group = generate_group(n, &[rot]);
            for n0 in 0..=n {
                let expected = polya_count_binary(&group, n0, n - n0) as usize;
                let actual = count_canonical(&[n0, n - n0], &group);
                assert_eq!(actual, expected, "C_{n}: mismatch at ({n0}, {})", n - n0);
            }
        }
    }

    #[test]
    fn dihedral_matches_polya() {
        for n in [3, 4, 5, 6] {
            let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
            let refl = Permutation::new((0..n).map(|i| ((n - i) % n) as u32).collect());
            let group = generate_group(n, &[rot, refl]);
            for n0 in 0..=n {
                let expected = polya_count_binary(&group, n0, n - n0) as usize;
                let actual = count_canonical(&[n0, n - n0], &group);
                assert_eq!(actual, expected, "D_{n}: mismatch at ({n0}, {})", n - n0);
            }
        }
    }

    #[test]
    fn symmetric_group_matches_polya() {
        let s01 = Permutation::new(vec![1, 0, 2]);
        let s12 = Permutation::new(vec![0, 2, 1]);
        let group = generate_group(3, &[s01, s12]);
        assert_eq!(group.len(), 6);
        for n0 in 0..=3 {
            let expected = polya_count_binary(&group, n0, 3 - n0) as usize;
            let actual = count_canonical(&[n0, 3 - n0], &group);
            assert_eq!(actual, expected, "S_3: ({n0}, {})", 3 - n0);
        }
    }

    #[test]
    fn outputs_are_lex_minimal() {
        let n = 5;
        let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
        let group = generate_group(n, &[rot]);

        enumerate_canonical(&[2, 3], &group, &mut |coloring| {
            for g in &group {
                let image: Vec<u32> = (0..n)
                    .map(|i| coloring[g.apply(i as u32) as usize])
                    .collect();
                assert!(
                    coloring <= &image[..],
                    "{coloring:?} not lex-min: g gives {image:?}"
                );
            }
        });
    }

    #[test]
    fn no_duplicates() {
        let n = 6;
        let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
        let refl = Permutation::new((0..n).map(|i| ((n - i) % n) as u32).collect());
        let group = generate_group(n, &[rot, refl]);

        let mut results = Vec::new();
        enumerate_canonical(&[3, 3], &group, &mut |c| results.push(c.to_vec()));

        let unique: HashSet<Vec<u32>> = results.iter().cloned().collect();
        assert_eq!(results.len(), unique.len());
    }

    #[test]
    fn ternary_cyclic() {
        let n = 3;
        let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
        let group = generate_group(n, &[rot]);
        // 3!/3 = 2 orbits: [0,1,2] and [0,2,1]
        assert_eq!(count_canonical(&[1, 1, 1], &group), 2);
    }

    #[test]
    fn ternary_symmetric() {
        let s01 = Permutation::new(vec![1, 0, 2]);
        let s12 = Permutation::new(vec![0, 2, 1]);
        let group = generate_group(3, &[s01, s12]);
        // S_3 on (1,1,1): all 6 colorings in one orbit
        assert_eq!(count_canonical(&[1, 1, 1], &group), 1);
    }
}
