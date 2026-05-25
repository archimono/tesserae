use crate::perm::Permutation;

/// Count the number of distinct colorings of `n` sites with two species
/// at composition `(n0, n1)` under the action of the permutation
/// group `group`, using the Pólya cycle-index method.
///
/// Returns the exact orbit count. Uses `u128` to avoid overflow.
///
/// # Panics
/// Panics if `n0 + n1` does not equal the permutation degree,
/// or if `group` is empty.
pub fn polya_count_binary(group: &[Permutation], n0: usize, n1: usize) -> u128 {
    let n = n0 + n1;
    assert!(!group.is_empty(), "group must be non-empty");
    assert_eq!(
        group[0].len(),
        n,
        "composition must match permutation degree"
    );

    let group_order = group.len() as u128;
    let mut total: u128 = 0;

    for g in group {
        total += fixed_colorings(g, n0);
    }

    assert_eq!(
        total % group_order,
        0,
        "Burnside sum must be divisible by |G|"
    );
    total / group_order
}

/// Count colorings with exactly `n0` species-0 sites that are fixed by
/// permutation `g`. A coloring is fixed iff every cycle is monochromatic.
/// We choose which cycles get species 0 such that their total size equals `n0`.
///
/// This is a subset-sum DP over cycle lengths.
fn fixed_colorings(g: &Permutation, n0: usize) -> u128 {
    let cycle_lengths: Vec<usize> = g.cycles().iter().map(|c| c.len()).collect();
    let n = g.len();

    // dp[w] = number of ways to assign cycles so that exactly w sites get species 0.
    let mut dp = vec![0u128; n + 1];
    dp[0] = 1;

    for &cl in &cycle_lengths {
        // Traverse in reverse to avoid reusing the same cycle.
        for w in (0..=n).rev() {
            if dp[w] == 0 {
                continue;
            }
            if w + cl <= n {
                dp[w + cl] += dp[w];
            }
        }
    }

    dp[n0]
}

/// Count the number of distinct k-ary colorings at a given composition
/// under a permutation group, using the Pólya cycle-index method.
///
/// `composition[j]` is the number of sites colored with color `j`.
/// Returns the exact orbit count via Burnside's lemma.
///
/// Internally uses a flattened `(k−1)`-dimensional knapsack DP where
/// the last color is implicit. For `k = 2` the result matches
/// [`polya_count_binary`].
///
/// # Panics
///
/// Panics if `group` is empty or if the composition sum differs from
/// the permutation degree.
pub fn polya_count(group: &[Permutation], composition: &[usize]) -> u128 {
    let k = composition.len();
    let n: usize = composition.iter().sum();

    assert!(!group.is_empty(), "group must be non-empty");
    assert_eq!(
        group[0].len(),
        n,
        "composition must match permutation degree"
    );

    if k <= 1 {
        return 1;
    }

    // Track k-1 dimensions; the last color count is implicit.
    let tracked = k - 1;
    let mut strides = vec![1usize; tracked];
    for i in 1..tracked {
        strides[i] = strides[i - 1] * (composition[i - 1] + 1);
    }
    let table_size = strides[tracked - 1] * (composition[tracked - 1] + 1);
    let target: usize = (0..tracked).map(|i| composition[i] * strides[i]).sum();

    let group_order = group.len() as u128;
    let mut total_fixed: u128 = 0;

    for g in group {
        debug_assert_eq!(g.len(), n);
        let cycle_lengths: Vec<usize> = g.cycles().iter().map(|c| c.len()).collect();

        let mut dp = vec![0u128; table_size];
        dp[0] = 1;

        for &cl in &cycle_lengths {
            for idx in (0..table_size).rev() {
                if dp[idx] == 0 {
                    continue;
                }
                let val = dp[idx];
                for j in 0..tracked {
                    let w_j = (idx / strides[j]) % (composition[j] + 1);
                    if w_j + cl <= composition[j] {
                        dp[idx + cl * strides[j]] += val;
                    }
                }
            }
        }

        total_fixed += dp[target];
    }

    assert_eq!(
        total_fixed % group_order,
        0,
        "Burnside sum must be divisible by |G|"
    );
    total_fixed / group_order
}

/// Count distinct k-ary colorings under a permutation group with per-site
/// species constraints.
///
/// `allowed[i][c]` is `true` iff color `c` may be placed on site `i`.
/// The allowed sets must be G-equivariant: for every `g` in `group` and
/// every site `i`, `allowed[g(i)] == allowed[i]`. Violating this
/// precondition will trigger a debug panic and produce nonsensical results.
///
/// The algorithm modifies the standard cycle-index DP: for each cycle of
/// a group element, only colors allowed on ALL sites in the cycle can
/// contribute a monochromatic coloring.
///
/// # Panics
///
/// Panics if `group` is empty, if the composition sum differs from the
/// permutation degree, or if `allowed` dimensions are inconsistent.
pub fn polya_count_constrained(
    group: &[Permutation],
    composition: &[usize],
    allowed: &[Vec<bool>],
) -> u128 {
    let k = composition.len();
    let n: usize = composition.iter().sum();

    assert!(!group.is_empty(), "group must be non-empty");
    assert_eq!(
        group[0].len(),
        n,
        "composition must match permutation degree"
    );
    assert_eq!(allowed.len(), n, "allowed must have one entry per site");
    assert!(
        allowed.iter().all(|a| a.len() == k),
        "each allowed entry must have length k"
    );
    debug_assert!(
        group
            .iter()
            .all(|g| (0..n).all(|i| allowed[g.apply(i as u32) as usize] == allowed[i])),
        "allowed sets must be G-equivariant"
    );

    if k == 0 {
        return if n == 0 { 1 } else { 0 };
    }
    if k == 1 {
        return if allowed.iter().all(|a| a[0]) { 1 } else { 0 };
    }

    let tracked = k - 1;
    let mut strides = vec![1usize; tracked];
    for i in 1..tracked {
        strides[i] = strides[i - 1] * (composition[i - 1] + 1);
    }
    let table_size = strides[tracked - 1] * (composition[tracked - 1] + 1);
    let target: usize = (0..tracked).map(|i| composition[i] * strides[i]).sum();

    let group_order = group.len() as u128;
    let mut total_fixed: u128 = 0;

    for g in group {
        debug_assert_eq!(g.len(), n);
        let cycles = g.cycles();

        let cycle_allowed: Vec<Vec<bool>> = cycles
            .iter()
            .map(|cycle| {
                let mut colors = vec![true; k];
                for &site in cycle {
                    for c in 0..k {
                        if !allowed[site as usize][c] {
                            colors[c] = false;
                        }
                    }
                }
                colors
            })
            .collect();
        let cycle_lengths: Vec<usize> = cycles.iter().map(|c| c.len()).collect();

        let mut dp = vec![0u128; table_size];
        dp[0] = 1;

        for (ci, &cl) in cycle_lengths.iter().enumerate() {
            let mut new_dp = vec![0u128; table_size];
            for idx in 0..table_size {
                if dp[idx] == 0 {
                    continue;
                }
                let val = dp[idx];
                if cycle_allowed[ci][tracked] {
                    new_dp[idx] += val;
                }
                for j in 0..tracked {
                    if !cycle_allowed[ci][j] {
                        continue;
                    }
                    let w_j = (idx / strides[j]) % (composition[j] + 1);
                    if w_j + cl <= composition[j] {
                        new_dp[idx + cl * strides[j]] += val;
                    }
                }
            }
            dp = new_dp;
        }

        total_fixed += dp[target];
    }

    assert_eq!(
        total_fixed % group_order,
        0,
        "Burnside sum must be divisible by |G|"
    );
    total_fixed / group_order
}

#[cfg(test)]
fn euler_totient(mut n: u64) -> u64 {
    let mut result = n;
    let mut p = 2u64;
    while p * p <= n {
        if n.is_multiple_of(p) {
            while n.is_multiple_of(p) {
                n /= p;
            }
            result -= result / p;
        }
        p += 1;
    }
    if n > 1 {
        result -= result / n;
    }
    result
}

#[cfg(test)]
fn necklace_count_closed_form(n: u64) -> u128 {
    let mut sum: u128 = 0;
    for d in 1..=n {
        if n.is_multiple_of(d) {
            let phi = euler_totient(d) as u128;
            let exp = (n / d) as u32;
            sum += phi * 2u128.pow(exp);
        }
    }
    sum / n as u128
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::perm::generate_group;

    fn cyclic_group(n: usize) -> Vec<Permutation> {
        let rot = Permutation::new((1..n as u32).chain(std::iter::once(0)).collect());
        generate_group(n, &[rot])
    }

    fn dihedral_group(n: usize) -> Vec<Permutation> {
        let rot = Permutation::new((1..n as u32).chain(std::iter::once(0)).collect());
        let refl_image: Vec<u32> = (0..n as u32).rev().collect();
        let refl = Permutation::new(refl_image);
        generate_group(n, &[rot, refl])
    }

    fn symmetric_group(n: usize) -> Vec<Permutation> {
        if n <= 1 {
            return vec![Permutation::identity(n)];
        }
        // Generated by adjacent transposition (0 1) and n-cycle (0 1 2 … n-1).
        let adj = {
            let mut img: Vec<u32> = (0..n as u32).collect();
            img[0] = 1;
            img[1] = 0;
            Permutation::new(img)
        };
        let cyc = Permutation::new((1..n as u32).chain(std::iter::once(0)).collect());
        generate_group(n, &[adj, cyc])
    }

    // --- Cyclic group tests (cross-check against closed-form necklace count) ---

    #[test]
    fn cyclic_total_colorings_match_necklace_formula() {
        for n in 2..=12usize {
            let group = cyclic_group(n);
            let expected = necklace_count_closed_form(n as u64);
            let mut computed: u128 = 0;
            for w in 0..=n {
                computed += polya_count_binary(&group, w, n - w);
            }
            assert_eq!(computed, expected, "C_{n}: total orbits mismatch (n={n})");
        }
    }

    #[test]
    fn cyclic_specific_compositions() {
        // C_6, k=2: known necklace counts by composition
        let group = cyclic_group(6);
        assert_eq!(polya_count_binary(&group, 6, 0), 1);
        assert_eq!(polya_count_binary(&group, 0, 6), 1);
        assert_eq!(polya_count_binary(&group, 1, 5), 1);
        // (3,3): Burnside gives (C(6,3) + 0 + 2·C(2,1) + 0 + 0 + 0) / 6 = 4
        assert_eq!(polya_count_binary(&group, 3, 3), 4);
    }

    // --- Trivial group tests ---

    #[test]
    fn trivial_group_gives_binomial() {
        for n in 1..=8usize {
            let group = vec![Permutation::identity(n)];
            for w in 0..=n {
                let expected = binomial(n, w);
                assert_eq!(
                    polya_count_binary(&group, w, n - w),
                    expected,
                    "trivial group: C({n},{w}) mismatch"
                );
            }
        }
    }

    fn binomial(n: usize, k: usize) -> u128 {
        if k > n {
            return 0;
        }
        let mut result: u128 = 1;
        for i in 0..k {
            result = result * (n - i) as u128 / (i + 1) as u128;
        }
        result
    }

    // --- Symmetric group tests ---

    #[test]
    fn symmetric_group_gives_one_orbit() {
        // S_n acts on n sites with composition (w, n-w): exactly 1 orbit.
        for n in 2..=6usize {
            let group = symmetric_group(n);
            assert_eq!(group.len(), (1..=n).product::<usize>()); // |S_n| = n!
            for w in 0..=n {
                assert_eq!(
                    polya_count_binary(&group, w, n - w),
                    1,
                    "S_{n} with composition ({w},{}) should be 1 orbit",
                    n - w
                );
            }
        }
    }

    // --- Dihedral group tests (bracelet counts) ---

    #[test]
    fn dihedral_total_colorings() {
        // OEIS A000029: total binary bracelets of length n
        // Start at n=3 because D_2 on 2 points collapses to C_2 (not faithful).
        // n:   3  4   5   6   7    8
        let expected = [4, 6, 8, 13, 18, 30];
        for (i, &exp) in expected.iter().enumerate() {
            let n = i + 3;
            let group = dihedral_group(n);
            assert_eq!(group.len(), 2 * n); // |D_n| = 2n
            let mut total: u128 = 0;
            for w in 0..=n {
                total += polya_count_binary(&group, w, n - w);
            }
            assert_eq!(total, exp, "D_{n}: bracelet count mismatch (n={n})");
        }
    }

    #[test]
    fn dihedral_specific_compositions() {
        // D_4 (square), 4 sites, composition (2,2): 2 bracelets
        let group = dihedral_group(4);
        assert_eq!(polya_count_binary(&group, 2, 2), 2);
    }

    // --- Edge cases ---

    #[test]
    fn single_site() {
        let group = vec![Permutation::identity(1)];
        assert_eq!(polya_count_binary(&group, 1, 0), 1);
        assert_eq!(polya_count_binary(&group, 0, 1), 1);
    }

    // --- k-ary Pólya tests ---

    #[test]
    fn k_ary_agrees_with_binary() {
        for n in 2..=8usize {
            let group = cyclic_group(n);
            for n0 in 0..=n {
                let binary = polya_count_binary(&group, n0, n - n0);
                let k_ary = polya_count(&group, &[n0, n - n0]);
                assert_eq!(binary, k_ary, "k=2 mismatch: C_{n}, ({n0}, {})", n - n0);
            }
        }
        for n in 3..=6usize {
            let group = dihedral_group(n);
            for n0 in 0..=n {
                let binary = polya_count_binary(&group, n0, n - n0);
                let k_ary = polya_count(&group, &[n0, n - n0]);
                assert_eq!(binary, k_ary, "k=2 mismatch: D_{n}, ({n0}, {})", n - n0);
            }
        }
    }

    #[test]
    fn ternary_cyclic() {
        let group = cyclic_group(3);
        // C_3 on (1,1,1): 3!/3 = 2 orbits
        assert_eq!(polya_count(&group, &[1, 1, 1]), 2);
        // C_3 on (2,1,0): 3!/3 = 1 orbit (only one binary pattern up to rotation)
        assert_eq!(polya_count(&group, &[2, 1, 0]), 1);
    }

    #[test]
    fn ternary_symmetric() {
        let group = symmetric_group(3);
        // S_3 on (1,1,1): 1 orbit
        assert_eq!(polya_count(&group, &[1, 1, 1]), 1);
    }

    #[test]
    fn ternary_trivial_gives_multinomial() {
        let n = 4;
        let group = vec![Permutation::identity(n)];
        // Trivial group: count = multinomial(4; 2,1,1) = 4!/(2!1!1!) = 12
        assert_eq!(polya_count(&group, &[2, 1, 1]), 12);
        // multinomial(4; 1,1,1,1) = 24
        assert_eq!(polya_count(&group, &[1, 1, 1, 1]), 24);
    }

    #[test]
    fn ternary_dihedral() {
        // D_4 on 4 sites, composition (2,1,1)
        // Burnside: e→12, r→0, r²→0, r³→0, s→0, sr→2, sr²→0, sr³→2 = 16/8 = 2
        let group = dihedral_group(4);
        assert_eq!(polya_count(&group, &[2, 1, 1]), 2);
    }

    #[test]
    fn k_ary_augment_cross_check() {
        use crate::augment::enumerate_canonical;

        // Cross-check polya_count against enumerate_canonical for ternary
        for n in [3, 4, 5, 6] {
            let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
            let group = generate_group(n, &[rot]);
            // Test a few ternary compositions
            for &comp in &[[1, 1, n - 2], [2, 1, n - 3]] {
                if comp.iter().any(|&c| c > n) {
                    continue;
                }
                if comp.iter().sum::<usize>() != n {
                    continue;
                }
                let polya = polya_count(&group, &comp);
                let mut aug_count = 0u128;
                enumerate_canonical(&comp, &group, &mut |_| aug_count += 1);
                assert_eq!(
                    aug_count, polya,
                    "C_{n} ternary: polya={polya}, augment={aug_count}, comp={comp:?}"
                );
            }
        }
    }

    #[test]
    fn single_color() {
        let group = cyclic_group(5);
        assert_eq!(polya_count(&group, &[5]), 1);
    }

    // --- Constrained Pólya tests ---

    #[test]
    fn constrained_unconstrained_equiv() {
        // All sites allow all colors → constrained == unconstrained
        for n in [3, 4, 6] {
            let group = cyclic_group(n);
            let all_allowed = vec![vec![true, true]; n];
            for n0 in 0..=n {
                let uc = polya_count(&group, &[n0, n - n0]);
                let co = polya_count_constrained(&group, &[n0, n - n0], &all_allowed);
                assert_eq!(
                    uc,
                    co,
                    "C_{n}, ({n0}, {}): constrained != unconstrained",
                    n - n0
                );
            }
        }
    }

    #[test]
    fn constrained_disjoint_identity() {
        // 4 sites, 4 colors, identity group.
        // Sites 0,1 allow {0,1}, sites 2,3 allow {2,3}.
        // Composition [1,1,1,1]: site 0 and 1 each get one of {0,1},
        // sites 2,3 each get one of {2,3} → 2 * 2 = 4.
        let group = vec![Permutation::identity(4)];
        let allowed = vec![
            vec![true, true, false, false],
            vec![true, true, false, false],
            vec![false, false, true, true],
            vec![false, false, true, true],
        ];
        assert_eq!(polya_count_constrained(&group, &[1, 1, 1, 1], &allowed), 4);
    }

    #[test]
    fn constrained_forced_zero() {
        // 2 sites, 2 colors. Site 0 allows {0}, site 1 allows {1}.
        // Composition [1,1] → exactly one coloring: [0, 1].
        let group = vec![Permutation::identity(2)];
        let allowed = vec![vec![true, false], vec![false, true]];
        assert_eq!(polya_count_constrained(&group, &[1, 1], &allowed), 1);

        // Composition [2,0] → impossible (site 1 can't take color 0).
        assert_eq!(polya_count_constrained(&group, &[2, 0], &allowed), 0);
    }

    #[test]
    fn constrained_cyclic_all_same_allowed() {
        // C_4, all sites allow {0, 1} of 3 colors, composition [2,2,0].
        // Since no site allows color 2 and composition[2]=0, this should
        // equal the C_4 binary count for (2,2).
        let group = cyclic_group(4);
        let allowed = vec![vec![true, true, false]; 4];
        let constrained = polya_count_constrained(&group, &[2, 2, 0], &allowed);
        let binary = polya_count_binary(&group, 2, 2);
        assert_eq!(constrained, binary);
        assert_eq!(constrained, 2);
    }

    #[test]
    fn constrained_two_orbits_block_diagonal() {
        // 4 sites as 2 orbits of 2. Group swaps within each orbit but
        // never mixes orbits: perm = (0↔1)(2↔3).
        // Orbit 0 (sites 0,1): allowed = {0, 1}
        // Orbit 1 (sites 2,3): allowed = {0}
        // 3 colors, composition [3, 1, 0].
        // Orbit 1 is forced: both sites = color 0 (contributes 2 to color 0).
        // Orbit 0 must contribute 1 to color 0 and 1 to color 1.
        // Under swap (0↔1): (0,1) and (1,0) are equivalent → 1 orbit.
        let swap = Permutation::new(vec![1, 0, 3, 2]);
        let group = vec![Permutation::identity(4), swap];
        let allowed = vec![
            vec![true, true, false],
            vec![true, true, false],
            vec![true, false, false],
            vec![true, false, false],
        ];
        assert_eq!(polya_count_constrained(&group, &[3, 1, 0], &allowed), 1);
    }
}
