use crate::augment::{
    enumerate_canonical, enumerate_canonical_constrained, enumerate_t_pruned,
    enumerate_t_pruned_constrained,
};
use crate::perm::Permutation;
use crate::polya::{polya_count, polya_count_constrained};

/// Enumerate all symmetry-inequivalent colorings with the given composition
/// under a permutation group.
///
/// Calls `callback` with each inequivalent coloring exactly once.
/// The enumeration count is cross-checked against the Pólya cycle-index
/// value; a mismatch indicates a bug and triggers a panic.
///
/// # Panics
///
/// Panics if `group` is empty, if the composition sum differs from
/// the permutation degree, or if the enumeration count disagrees
/// with the Pólya prediction.
pub fn enumerate(composition: &[usize], group: &[Permutation], callback: &mut impl FnMut(&[u32])) {
    let expected = polya_count(group, composition);
    let mut count = 0u128;
    enumerate_canonical(composition, group, &mut |coloring| {
        callback(coloring);
        count += 1;
    });
    assert_eq!(
        count, expected,
        "enumeration produced {count} colorings but Pólya predicts {expected}"
    );
}

/// Enumerate using an explicit T ⋊ P decomposition.
///
/// `translations` must be the full translation subgroup T (including
/// the identity). `point_group` must be the full point group P
/// (including the identity). The full group G = T ⋊ P is constructed
/// as {t ∘ p : t ∈ T, p ∈ P} and passed to [`enumerate`].
///
/// # Panics
///
/// Panics (debug) if the product has duplicates (meaning T ∩ P ≠ {e}).
pub fn enumerate_with_decomposition(
    composition: &[usize],
    translations: &[Permutation],
    point_group: &[Permutation],
    callback: &mut impl FnMut(&[u32]),
) {
    let full_group = build_full_group(translations, point_group);
    enumerate(composition, &full_group, callback);
}

/// Enumerate using the T ⋊ P decomposition with T-pruned augmentation.
///
/// Uses dual active sets (P-active + T-active) during McKay tree traversal
/// to prune branches that are not canonical under either P or T. Surviving
/// leaves are checked against the full G-canonical form via a short-circuit
/// predicate and emitted only if canonical.
///
/// No HashSet. No unbounded memory. Streaming-compatible.
///
/// # Panics
///
/// Panics (debug) if the product T × P has duplicates, or if the final
/// count disagrees with the Pólya prediction for the full group G = T ⋊ P.
pub fn enumerate_decomposed(
    composition: &[usize],
    translations: &[Permutation],
    point_group: &[Permutation],
    callback: &mut impl FnMut(&[u32]),
) {
    let n: usize = composition.iter().sum();
    let full_group = build_full_group(translations, point_group);
    let expected = polya_count(&full_group, composition);

    let mut count = 0u128;
    enumerate_t_pruned(composition, point_group, translations, &mut |coloring| {
        if is_g_canonical(coloring, translations, point_group, n) {
            callback(coloring);
            count += 1;
        }
    });

    assert_eq!(
        count, expected,
        "decomposed enumeration produced {count} colorings but Pólya predicts {expected}"
    );
}

/// Like [`enumerate_decomposed`] but infers the cyclic translation group
/// ℤ_n from `n = sum(composition)`.
///
/// # Panics
///
/// Same as [`enumerate_decomposed`].
pub fn enumerate_decomposed_cyclic(
    composition: &[usize],
    point_group: &[Permutation],
    callback: &mut impl FnMut(&[u32]),
) {
    let n: usize = composition.iter().sum();
    let translations = crate::t_canon::cyclic_translations(n);
    enumerate_decomposed(composition, &translations, point_group, callback);
}

fn stabilizer_size(coloring: &[u32], group: &[Permutation]) -> u128 {
    let n = coloring.len();
    group
        .iter()
        .filter(|g| (0..n).all(|i| coloring[g.apply(i as u32) as usize] == coloring[i]))
        .count() as u128
}

fn multinomial(composition: &[usize]) -> u128 {
    let n: usize = composition.iter().sum();
    let mut result: u128 = 1;
    let mut remaining = n;
    for &k in composition {
        for j in 1..=k {
            result = result
                .checked_mul(remaining as u128)
                .expect("multinomial overflow")
                / j as u128;
            remaining -= 1;
        }
    }
    result
}

/// Enumerate all symmetry-inequivalent colorings with multiplicity.
///
/// Like [`enumerate`] but calls `callback` with `(coloring, orbit_size)`.
/// The orbit size is `|G| / |Stab(c)|` — the number of distinct images of
/// the coloring under the group. Asserts that the sum of all orbit sizes
/// equals the multinomial coefficient n! / ∏(comp[i]!).
pub fn enumerate_with_multiplicity(
    composition: &[usize],
    group: &[Permutation],
    callback: &mut impl FnMut(&[u32], u128),
) {
    let expected = polya_count(group, composition);
    let group_size = group.len() as u128;
    let multi = multinomial(composition);
    let mut count = 0u128;
    let mut weight_sum = 0u128;
    enumerate_canonical(composition, group, &mut |coloring| {
        let stab = stabilizer_size(coloring, group);
        let orbit_size = group_size / stab;
        callback(coloring, orbit_size);
        count += 1;
        weight_sum += orbit_size;
    });
    assert_eq!(
        count, expected,
        "enumeration produced {count} colorings but Pólya predicts {expected}"
    );
    assert_eq!(
        weight_sum, multi,
        "orbit size sum {weight_sum} != multinomial {multi}"
    );
}

/// Enumerate using T ⋊ P decomposition with multiplicity.
///
/// Like [`enumerate_decomposed`] but calls `callback` with `(coloring, orbit_size)`.
/// Asserts both the Pólya orbit count and the multinomial sum.
pub fn enumerate_decomposed_with_multiplicity(
    composition: &[usize],
    translations: &[Permutation],
    point_group: &[Permutation],
    callback: &mut impl FnMut(&[u32], u128),
) {
    let n: usize = composition.iter().sum();
    let full_group = build_full_group(translations, point_group);
    let expected = polya_count(&full_group, composition);
    let group_size = full_group.len() as u128;
    let multi = multinomial(composition);

    let mut count = 0u128;
    let mut weight_sum = 0u128;
    enumerate_t_pruned(composition, point_group, translations, &mut |coloring| {
        if is_g_canonical(coloring, translations, point_group, n) {
            let stab = stabilizer_size(coloring, &full_group);
            let orbit_size = group_size / stab;
            callback(coloring, orbit_size);
            count += 1;
            weight_sum += orbit_size;
        }
    });

    assert_eq!(
        count, expected,
        "decomposed enumeration produced {count} colorings but Pólya predicts {expected}"
    );
    assert_eq!(
        weight_sum, multi,
        "orbit size sum {weight_sum} != multinomial {multi}"
    );
}

/// Enumerate primitive (non-superperiodic) colorings with multiplicity.
///
/// Like [`enumerate_decomposed_primitive`] but calls `callback` with
/// `(coloring, orbit_size)`. The Pólya cross-check validates the full
/// G-canonical count (before primitive filtering). The multinomial sum
/// is NOT asserted because superperiodic orbits are excluded.
pub fn enumerate_decomposed_primitive_with_multiplicity(
    composition: &[usize],
    translations: &[Permutation],
    point_group: &[Permutation],
    callback: &mut impl FnMut(&[u32], u128),
) {
    let n: usize = composition.iter().sum();
    let full_group = build_full_group(translations, point_group);
    let expected = polya_count(&full_group, composition);
    let group_size = full_group.len() as u128;

    let mut count = 0u128;
    enumerate_t_pruned(composition, point_group, translations, &mut |coloring| {
        if is_g_canonical(coloring, translations, point_group, n) {
            count += 1;
            if !is_superperiodic(coloring, translations) {
                let stab = stabilizer_size(coloring, &full_group);
                let orbit_size = group_size / stab;
                callback(coloring, orbit_size);
            }
        }
    });

    assert_eq!(
        count, expected,
        "primitive enumeration: {count} G-canonical colorings but Pólya predicts {expected}"
    );
}

/// Enumerate symmetry-inequivalent colorings with per-site species constraints.
///
/// Like [`enumerate`] but only places colors allowed at each site. The Pólya
/// cross-check uses the constrained cycle-index DP.
///
/// `allowed[i][c]` must be `true` iff color `c` may be placed on site `i`.
/// The allowed sets must be G-equivariant.
///
/// # Panics
///
/// Panics if the enumeration count disagrees with the constrained Pólya
/// prediction.
pub fn enumerate_constrained(
    composition: &[usize],
    group: &[Permutation],
    allowed: &[Vec<bool>],
    callback: &mut impl FnMut(&[u32]),
) {
    let expected = polya_count_constrained(group, composition, allowed);
    let mut count = 0u128;
    enumerate_canonical_constrained(composition, group, allowed, &mut |coloring| {
        callback(coloring);
        count += 1;
    });
    assert_eq!(
        count, expected,
        "constrained enumeration produced {count} colorings but Pólya predicts {expected}"
    );
}

/// Enumerate constrained colorings using T ⋊ P decomposition.
///
/// Like [`enumerate_decomposed`] but with per-site species constraints.
///
/// # Panics
///
/// Same as [`enumerate_decomposed`] plus constraint-related mismatches.
pub fn enumerate_decomposed_constrained(
    composition: &[usize],
    translations: &[Permutation],
    point_group: &[Permutation],
    allowed: &[Vec<bool>],
    callback: &mut impl FnMut(&[u32]),
) {
    let n: usize = composition.iter().sum();
    let full_group = build_full_group(translations, point_group);
    let expected = polya_count_constrained(&full_group, composition, allowed);

    let mut count = 0u128;
    enumerate_t_pruned_constrained(
        composition,
        point_group,
        translations,
        allowed,
        &mut |coloring| {
            if is_g_canonical(coloring, translations, point_group, n) {
                callback(coloring);
                count += 1;
            }
        },
    );

    assert_eq!(
        count, expected,
        "constrained decomposed enumeration produced {count} colorings but Pólya predicts {expected}"
    );
}

/// Reservoir sampling adapter for streaming enumeration.
///
/// Wraps a callback to collect at most `k` samples uniformly at random
/// from the stream, using Vitter's Algorithm R. After enumeration completes,
/// call `into_samples()` to get the collected samples.
pub struct ReservoirSampler {
    reservoir: Vec<Vec<u32>>,
    k: usize,
    seen: u64,
    rng_state: u64,
}

impl ReservoirSampler {
    pub fn new(k: usize, seed: u64) -> Self {
        Self {
            reservoir: Vec::with_capacity(k),
            k,
            seen: 0,
            rng_state: if seed == 0 { 1 } else { seed },
        }
    }

    pub fn observe(&mut self, coloring: &[u32]) {
        self.seen += 1;
        if self.reservoir.len() < self.k {
            self.reservoir.push(coloring.to_vec());
        } else {
            let j = self.bounded_random(self.seen);
            if (j as usize) < self.k {
                self.reservoir[j as usize] = coloring.to_vec();
            }
        }
    }

    pub fn into_samples(self) -> Vec<Vec<u32>> {
        self.reservoir
    }

    pub fn samples(&self) -> &[Vec<u32>] {
        &self.reservoir
    }

    fn next_u64(&mut self) -> u64 {
        self.rng_state ^= self.rng_state << 13;
        self.rng_state ^= self.rng_state >> 7;
        self.rng_state ^= self.rng_state << 17;
        self.rng_state
    }

    fn bounded_random(&mut self, bound: u64) -> u64 {
        let threshold = bound.wrapping_neg() % bound;
        loop {
            let r = self.next_u64();
            if r >= threshold {
                return r % bound;
            }
        }
    }
}

/// Generate all compositions of `n` into `k` parts where each part
/// satisfies `min[i] <= part[i] <= max[i]`.
///
/// Returns an empty vec if no valid compositions exist.
///
/// # Panics
///
/// Panics if `min.len() != k` or `max.len() != k`.
pub fn compositions_in_range(n: usize, min: &[usize], max: &[usize]) -> Vec<Vec<usize>> {
    let k = min.len();
    assert_eq!(max.len(), k);

    let mut result = Vec::new();
    let mut current = vec![0usize; k];
    fn recurse(
        k: usize,
        pos: usize,
        remaining: usize,
        min: &[usize],
        max: &[usize],
        current: &mut Vec<usize>,
        result: &mut Vec<Vec<usize>>,
    ) {
        if pos == k - 1 {
            if remaining >= min[pos] && remaining <= max[pos] {
                current[pos] = remaining;
                result.push(current.clone());
            }
            return;
        }
        let min_rest: usize = min[pos + 1..].iter().sum();
        let max_rest: usize = max[pos + 1..].iter().sum();
        let lo = min[pos].max(remaining.saturating_sub(max_rest));
        let hi = max[pos].min(remaining.saturating_sub(min_rest));
        for i in lo..=hi {
            current[pos] = i;
            recurse(k, pos + 1, remaining - i, min, max, current, result);
        }
    }
    if k == 0 {
        if n == 0 {
            result.push(Vec::new());
        }
        return result;
    }
    recurse(k, 0, n, min, max, &mut current, &mut result);
    result
}

/// Check if `coloring` is its own G-canonical form, short-circuiting
/// as soon as any (p, t) pair produces a lex-smaller image.
fn is_g_canonical(
    coloring: &[u32],
    translations: &[Permutation],
    point_group: &[Permutation],
    n: usize,
) -> bool {
    for p in point_group {
        let p_image: Vec<u32> = (0..n)
            .map(|i| coloring[p.apply(i as u32) as usize])
            .collect();
        for t in translations {
            let mut decided = false;
            for i in 0..n {
                let val = p_image[t.apply(i as u32) as usize];
                match val.cmp(&coloring[i]) {
                    std::cmp::Ordering::Less => return false,
                    std::cmp::Ordering::Greater => {
                        decided = true;
                        break;
                    }
                    std::cmp::Ordering::Equal => {}
                }
            }
            if !decided {
                // t(p·c) == c, which is fine (identity or stabilizer element)
            }
        }
    }
    true
}

/// Enumerate primitive (non-superperiodic) symmetry-inequivalent colorings.
///
/// Like [`enumerate_decomposed`] but skips colorings that are invariant under
/// any non-identity translation in T. Such colorings are representable by a
/// smaller supercell and are excluded by tools like enumlib. The Pólya
/// cross-check validates the full G-canonical count (before primitive
/// filtering) so the enumeration engine is still verified on every call.
///
/// # Panics
///
/// Same as [`enumerate_decomposed`].
pub fn enumerate_decomposed_primitive(
    composition: &[usize],
    translations: &[Permutation],
    point_group: &[Permutation],
    callback: &mut impl FnMut(&[u32]),
) {
    let n: usize = composition.iter().sum();
    let full_group = build_full_group(translations, point_group);
    let expected = polya_count(&full_group, composition);

    let mut count = 0u128;
    enumerate_t_pruned(composition, point_group, translations, &mut |coloring| {
        if is_g_canonical(coloring, translations, point_group, n) {
            count += 1;
            if !is_superperiodic(coloring, translations) {
                callback(coloring);
            }
        }
    });

    assert_eq!(
        count, expected,
        "primitive enumeration: {count} G-canonical colorings but Pólya predicts {expected}"
    );
}

/// Return true if `coloring` is invariant under any non-identity translation.
///
/// `translations[0]` must be the identity permutation (convention from
/// `product_translations`). Superperiodic colorings can be represented by a
/// smaller supercell and are excluded from primitive enumeration.
fn is_superperiodic(coloring: &[u32], translations: &[Permutation]) -> bool {
    translations[1..].iter().any(|tau| {
        tau.image()
            .iter()
            .enumerate()
            .all(|(i, &ti)| coloring[ti as usize] == coloring[i])
    })
}

fn build_full_group(translations: &[Permutation], point_group: &[Permutation]) -> Vec<Permutation> {
    let mut group = Vec::with_capacity(translations.len() * point_group.len());
    for t in translations {
        for p in point_group {
            group.push(p.compose(t));
        }
    }
    debug_assert_eq!(
        group.len(),
        {
            let set: std::collections::HashSet<&Permutation> = group.iter().collect();
            set.len()
        },
        "T × P product has duplicates — T ∩ P ≠ {{e}}"
    );
    group
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::perm::generate_group;
    use crate::t_canon::{cyclic_translations, product_translations_2d};

    fn collect_colorings(composition: &[usize], group: &[Permutation]) -> Vec<Vec<u32>> {
        let mut results = Vec::new();
        enumerate(composition, group, &mut |c| results.push(c.to_vec()));
        results
    }

    // --- Full-group enumeration tests ---

    #[test]
    fn cyclic_binary() {
        let n = 6;
        let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
        let group = generate_group(n, &[rot]);
        let results = collect_colorings(&[3, 3], &group);
        assert_eq!(results.len(), 4);
    }

    #[test]
    fn dihedral_ternary() {
        let n = 4;
        let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
        let refl = Permutation::new((0..n).map(|i| ((n - i) % n) as u32).collect());
        let group = generate_group(n, &[rot, refl]);
        let results = collect_colorings(&[2, 1, 1], &group);
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn single_color() {
        let group = vec![Permutation::identity(5)];
        let results = collect_colorings(&[5], &group);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0], vec![0, 0, 0, 0, 0]);
    }

    // --- T ⋊ P decomposition tests ---

    #[test]
    fn decomposition_cyclic_trivial_p() {
        let n = 6;
        let translations = cyclic_translations(n);
        let point_group = vec![Permutation::identity(n)];
        let mut count = 0;
        enumerate_with_decomposition(&[3, 3], &translations, &point_group, &mut |_| count += 1);
        let expected = polya_count(&translations, &[3, 3]);
        assert_eq!(count as u128, expected);
    }

    #[test]
    fn decomposition_dihedral_as_cyclic_times_reflection() {
        let n = 6;
        let translations = cyclic_translations(n);
        let refl = Permutation::new((0..n).map(|i| ((n - i) % n) as u32).collect());
        let point_group = vec![Permutation::identity(n), refl.clone()];

        // Verify the full group is D_6
        let full = build_full_group(&translations, &point_group);
        assert_eq!(full.len(), 2 * n);

        let mut decomp_results = Vec::new();
        enumerate_with_decomposition(&[3, 3], &translations, &point_group, &mut |c| {
            decomp_results.push(c.to_vec())
        });

        // Cross-check: build D_6 directly and enumerate
        let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
        let direct_group = generate_group(n, &[rot, refl]);
        assert_eq!(direct_group.len(), 2 * n);

        let direct_results = collect_colorings(&[3, 3], &direct_group);
        assert_eq!(decomp_results.len(), direct_results.len());
    }

    #[test]
    fn decomposition_2d_product_translations() {
        let translations = product_translations_2d(2, 3);
        assert_eq!(translations.len(), 6);
        let point_group = vec![Permutation::identity(6)];

        let mut count = 0;
        enumerate_with_decomposition(&[3, 3], &translations, &point_group, &mut |_| count += 1);
        let expected = polya_count(&translations, &[3, 3]);
        assert_eq!(count as u128, expected);
    }

    #[test]
    fn decomposition_2d_with_point_group() {
        // ℤ_2 × ℤ_2 on 4 sites with a point-group rotation
        // that swaps the two axes: (a,b) → (b,a).
        // Site index = a * 2 + b for a ∈ {0,1}, b ∈ {0,1}.
        // Swap axes: (a,b) → (b,a) maps index a*2+b → b*2+a.
        let translations = product_translations_2d(2, 2);
        assert_eq!(translations.len(), 4);

        // Axis-swap permutation on 4 sites: 0→0, 1→2, 2→1, 3→3
        let swap = Permutation::new(vec![0, 2, 1, 3]);
        let point_group = vec![Permutation::identity(4), swap];

        let mut count = 0;
        enumerate_with_decomposition(&[2, 2], &translations, &point_group, &mut |_| count += 1);
        // Pólya cross-check is built into enumerate, so reaching here = correct
        assert!(count > 0);
    }

    // --- Decomposed (two-phase) enumeration tests ---

    fn collect_decomposed(
        composition: &[usize],
        translations: &[Permutation],
        point_group: &[Permutation],
    ) -> Vec<Vec<u32>> {
        let mut results = Vec::new();
        enumerate_decomposed(composition, translations, point_group, &mut |c| {
            results.push(c.to_vec())
        });
        results
    }

    #[test]
    fn decomposed_cyclic_trivial_p() {
        let n = 6;
        let translations = cyclic_translations(n);
        let point_group = vec![Permutation::identity(n)];
        let results = collect_decomposed(&[3, 3], &translations, &point_group);
        let full = build_full_group(&translations, &point_group);
        let expected = collect_colorings(&[3, 3], &full);
        assert_eq!(results.len(), expected.len());
    }

    #[test]
    fn decomposed_dihedral() {
        let n = 6;
        let translations = cyclic_translations(n);
        let refl = Permutation::new((0..n).map(|i| ((n - i) % n) as u32).collect());
        let point_group = vec![Permutation::identity(n), refl.clone()];
        let results = collect_decomposed(&[3, 3], &translations, &point_group);

        let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
        let direct = generate_group(n, &[rot, refl]);
        let expected = collect_colorings(&[3, 3], &direct);
        assert_eq!(results.len(), expected.len());
    }

    #[test]
    fn decomposed_2d_product() {
        let translations = product_translations_2d(3, 3);
        let point_group = vec![Permutation::identity(9)];
        let results = collect_decomposed(&[3, 3, 3], &translations, &point_group);
        let expected = polya_count(&translations, &[3, 3, 3]);
        assert_eq!(results.len() as u128, expected);
    }

    #[test]
    fn decomposed_agrees_with_full_group() {
        for n in [4, 6, 8, 10] {
            let translations = cyclic_translations(n);
            let refl = Permutation::new((0..n).map(|i| ((n - i) % n) as u32).collect());
            let point_group = vec![Permutation::identity(n), refl.clone()];

            let decomposed = collect_decomposed(&[n / 2, n - n / 2], &translations, &point_group);

            let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
            let full = generate_group(n, &[rot, refl]);
            let direct = collect_colorings(&[n / 2, n - n / 2], &full);

            assert_eq!(
                decomposed.len(),
                direct.len(),
                "decomposed vs full mismatch at n={n}"
            );
        }
    }

    // --- Decomposed cyclic (Booth fast path) tests ---

    fn collect_decomposed_cyclic(
        composition: &[usize],
        point_group: &[Permutation],
    ) -> Vec<Vec<u32>> {
        let mut results = Vec::new();
        enumerate_decomposed_cyclic(composition, point_group, &mut |c| results.push(c.to_vec()));
        results
    }

    #[test]
    fn decomposed_cyclic_booth_trivial_p() {
        let n = 6;
        let point_group = vec![Permutation::identity(n)];
        let results = collect_decomposed_cyclic(&[3, 3], &point_group);
        let translations = cyclic_translations(n);
        let expected = collect_decomposed(&[3, 3], &translations, &point_group);
        assert_eq!(results.len(), expected.len());
    }

    #[test]
    fn decomposed_cyclic_booth_dihedral() {
        for n in [4, 6, 8, 12] {
            let refl = Permutation::new((0..n).map(|i| ((n - i) % n) as u32).collect());
            let point_group = vec![Permutation::identity(n), refl.clone()];
            let results = collect_decomposed_cyclic(&[n / 2, n - n / 2], &point_group);

            let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
            let full = generate_group(n, &[rot, refl]);
            let direct = collect_colorings(&[n / 2, n - n / 2], &full);
            assert_eq!(
                results.len(),
                direct.len(),
                "cyclic Booth vs full mismatch at n={n}"
            );
        }
    }

    #[test]
    fn decomposed_cyclic_booth_ternary() {
        for n in [6, 9] {
            let point_group = vec![Permutation::identity(n)];
            let comp = vec![n / 3, n / 3, n - 2 * (n / 3)];
            let results = collect_decomposed_cyclic(&comp, &point_group);
            let translations = cyclic_translations(n);
            let expected = polya_count(&translations, &comp);
            assert_eq!(
                results.len() as u128,
                expected,
                "ternary Booth mismatch at n={n}"
            );
        }
    }

    #[test]
    fn decomposition_ternary_3x3() {
        let translations = product_translations_2d(3, 3);
        assert_eq!(translations.len(), 9);
        let point_group = vec![Permutation::identity(9)];

        let mut count = 0;
        enumerate_with_decomposition(&[3, 3, 3], &translations, &point_group, &mut |_| count += 1);
        let expected = polya_count(&translations, &[3, 3, 3]);
        assert_eq!(count as u128, expected);
    }

    #[test]
    fn decomposed_2d_with_point_group() {
        let translations = product_translations_2d(2, 2);
        let swap = Permutation::new(vec![0, 2, 1, 3]);
        let point_group = vec![Permutation::identity(4), swap];
        let results = collect_decomposed(&[2, 2], &translations, &point_group);

        let full = build_full_group(&translations, &point_group);
        let expected = collect_colorings(&[2, 2], &full);
        assert_eq!(results.len(), expected.len());
    }

    #[test]
    fn decomposed_ternary() {
        for n in [6, 9] {
            let translations = cyclic_translations(n);
            let point_group = vec![Permutation::identity(n)];
            let comp = vec![n / 3, n / 3, n - 2 * (n / 3)];
            let results = collect_decomposed(&comp, &translations, &point_group);
            let expected = polya_count(&translations, &comp);
            assert_eq!(
                results.len() as u128,
                expected,
                "decomposed ternary mismatch at n={n}"
            );
        }
    }

    // --- Multiplicity tests ---

    #[test]
    fn multiplicity_sum_equals_multinomial() {
        let n = 6;
        let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
        let refl = Permutation::new((0..n).map(|i| ((n - i) % n) as u32).collect());
        let group = generate_group(n, &[rot, refl]);

        let mut weight_sum = 0u128;
        enumerate_with_multiplicity(&[3, 3], &group, &mut |_, mult| {
            weight_sum += mult;
        });
        assert_eq!(weight_sum, multinomial(&[3, 3]));
    }

    #[test]
    fn multiplicity_trivial_stabilizer() {
        let group = vec![Permutation::identity(4)];
        let mut results = Vec::new();
        enumerate_with_multiplicity(&[2, 1, 1], &group, &mut |c, mult| {
            results.push((c.to_vec(), mult));
        });
        for (_, mult) in &results {
            assert_eq!(*mult, 1);
        }
    }

    #[test]
    fn multiplicity_single_color() {
        let n = 6;
        let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
        let group = generate_group(n, &[rot]);
        let mut results = Vec::new();
        enumerate_with_multiplicity(&[6], &group, &mut |c, mult| {
            results.push((c.to_vec(), mult));
        });
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].1, 1);
    }

    #[test]
    fn multiplicity_decomposed_sum() {
        let n = 8;
        let translations = cyclic_translations(n);
        let refl = Permutation::new((0..n).map(|i| ((n - i) % n) as u32).collect());
        let point_group = vec![Permutation::identity(n), refl];

        let mut weight_sum = 0u128;
        enumerate_decomposed_with_multiplicity(
            &[4, 4],
            &translations,
            &point_group,
            &mut |_, mult| {
                weight_sum += mult;
            },
        );
        assert_eq!(weight_sum, multinomial(&[4, 4]));
    }

    #[test]
    fn multiplicity_primitive_runs() {
        let n = 6;
        let translations = cyclic_translations(n);
        let point_group = vec![Permutation::identity(n)];
        let mut count = 0;
        enumerate_decomposed_primitive_with_multiplicity(
            &[3, 3],
            &translations,
            &point_group,
            &mut |_, mult| {
                assert!(mult > 0);
                count += 1;
            },
        );
        assert!(count > 0);
    }

    #[test]
    fn multinomial_values() {
        assert_eq!(multinomial(&[3, 3]), 20); // 6!/(3!3!)
        assert_eq!(multinomial(&[4, 4]), 70); // 8!/(4!4!)
        assert_eq!(multinomial(&[2, 1, 1]), 12); // 4!/(2!1!1!)
        assert_eq!(multinomial(&[5]), 1);
    }

    // --- Concentration range tests ---

    #[test]
    fn compositions_in_range_exact() {
        let comps = compositions_in_range(6, &[3, 3], &[3, 3]);
        assert_eq!(comps, vec![vec![3, 3]]);
    }

    #[test]
    fn compositions_in_range_binary() {
        let comps = compositions_in_range(4, &[1, 1], &[3, 3]);
        assert_eq!(comps, vec![vec![1, 3], vec![2, 2], vec![3, 1]]);
    }

    #[test]
    fn compositions_in_range_ternary() {
        let comps = compositions_in_range(4, &[0, 0, 0], &[4, 4, 4]);
        let mut expected = Vec::new();
        for a in 0..=4usize {
            for b in 0..=4 - a {
                expected.push(vec![a, b, 4 - a - b]);
            }
        }
        assert_eq!(comps.len(), expected.len());
    }

    #[test]
    fn compositions_in_range_empty() {
        let comps = compositions_in_range(4, &[3, 3], &[4, 4]);
        assert!(comps.is_empty());
    }

    #[test]
    fn compositions_in_range_full_sweep_matches_all() {
        let n = 6;
        let comps = compositions_in_range(n, &[0, 0], &[n, n]);
        assert_eq!(comps.len(), n + 1);
    }

    // --- Reservoir sampler tests ---

    #[test]
    fn sampler_collects_all_when_fewer_than_k() {
        let n = 6;
        let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
        let group = generate_group(n, &[rot]);
        let mut sampler = ReservoirSampler::new(100, 42);
        enumerate(&[3, 3], &group, &mut |c| sampler.observe(c));
        assert_eq!(sampler.into_samples().len(), 4);
    }

    #[test]
    fn sampler_limits_to_k() {
        let group = vec![Permutation::identity(8)];
        let mut sampler = ReservoirSampler::new(5, 42);
        enumerate(&[4, 4], &group, &mut |c| sampler.observe(c));
        assert_eq!(sampler.into_samples().len(), 5);
    }

    #[test]
    fn sampler_deterministic_with_same_seed() {
        let group = vec![Permutation::identity(6)];
        let mut s1 = ReservoirSampler::new(3, 123);
        enumerate(&[3, 3], &group, &mut |c| s1.observe(c));
        let mut s2 = ReservoirSampler::new(3, 123);
        enumerate(&[3, 3], &group, &mut |c| s2.observe(c));
        assert_eq!(s1.into_samples(), s2.into_samples());
    }

    // --- Constrained enumeration tests ---

    #[test]
    fn constrained_identity_disjoint() {
        // 4 sites, 4 colors. Sites 0,1 ∈ {0,1}, sites 2,3 ∈ {2,3}.
        let group = vec![Permutation::identity(4)];
        let allowed = vec![
            vec![true, true, false, false],
            vec![true, true, false, false],
            vec![false, false, true, true],
            vec![false, false, true, true],
        ];
        let mut results = Vec::new();
        enumerate_constrained(&[1, 1, 1, 1], &group, &allowed, &mut |c| {
            results.push(c.to_vec())
        });
        assert_eq!(results.len(), 4);
    }

    #[test]
    fn constrained_decomposed_two_orbits() {
        // 4 sites as 2×ℤ_2. Translations: (0↔1)(2↔3).
        // Orbit 0 (sites 0,1): {0,1}. Orbit 1 (sites 2,3): {0}.
        // P = identity only.
        // comp = [3, 1, 0]: orbit 1 forced to color 0 (2 sites),
        //   orbit 0 must have 1×color 0 + 1×color 1.
        //   Under ℤ_2 translation: 1 orbit.
        let trans = vec![Permutation::identity(4), Permutation::new(vec![1, 0, 3, 2])];
        let pg = vec![Permutation::identity(4)];
        let allowed = vec![
            vec![true, true, false],
            vec![true, true, false],
            vec![true, false, false],
            vec![true, false, false],
        ];
        let mut results = Vec::new();
        enumerate_decomposed_constrained(&[3, 1, 0], &trans, &pg, &allowed, &mut |c| {
            results.push(c.to_vec())
        });
        assert_eq!(results.len(), 1);
    }

    #[test]
    fn constrained_agrees_with_unconstrained_when_all_allowed() {
        let n = 6;
        let translations = cyclic_translations(n);
        let refl = Permutation::new((0..n).map(|i| ((n - i) % n) as u32).collect());
        let point_group = vec![Permutation::identity(n), refl];
        let all_allowed = vec![vec![true, true]; n];

        let mut uc_results = Vec::new();
        enumerate_decomposed(&[3, 3], &translations, &point_group, &mut |c| {
            uc_results.push(c.to_vec())
        });

        let mut co_results = Vec::new();
        enumerate_decomposed_constrained(
            &[3, 3],
            &translations,
            &point_group,
            &all_allowed,
            &mut |c| co_results.push(c.to_vec()),
        );

        assert_eq!(uc_results.len(), co_results.len());
    }
}
