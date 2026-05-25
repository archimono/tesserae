use crate::enumerate::{
    compositions_in_range, enumerate_decomposed, enumerate_decomposed_constrained,
    enumerate_decomposed_primitive, enumerate_decomposed_primitive_with_multiplicity,
    enumerate_decomposed_with_multiplicity,
};
use crate::hnf::{hnf_2d, hnf_3d};
use crate::perm::Permutation;
use crate::snf::{SmithNormalForm, smith_normal_form};
use crate::t_canon::{product_translations, product_translations_multisite};

/// Enumerate symmetry-inequivalent colorings of a supercell.
///
/// Given a d×d integer supercell matrix, a composition summing to
/// |det(supercell_matrix)|, and the parent point group as d×d integer
/// matrices, enumerates all colorings up to the full space group
/// G = T ⋊ P acting on the supercell sites.
///
/// The translation subgroup T is inferred from the Smith Normal Form
/// of the supercell matrix. Parent point group elements that do not
/// preserve the supercell quotient are silently filtered.
///
/// # Panics
///
/// Panics if the supercell matrix has non-positive determinant, if the
/// composition sum differs from the number of sites, or if the Pólya
/// cross-check fails.
pub fn enumerate_supercell(
    supercell_matrix: &[Vec<i64>],
    composition: &[usize],
    parent_point_group: &[Vec<Vec<i64>>],
    callback: &mut impl FnMut(&[u32]),
) {
    let snf = smith_normal_form(supercell_matrix);
    let dims: Vec<usize> = snf
        .diagonal
        .iter()
        .map(|&d| {
            assert!(d > 0, "supercell matrix must be non-singular");
            d as usize
        })
        .collect();
    let n: usize = dims.iter().product();
    assert_eq!(
        composition.iter().sum::<usize>(),
        n,
        "composition sum {} != site count {n}",
        composition.iter().sum::<usize>()
    );

    let translations = product_translations(&dims);
    let point_group = transform_point_group(parent_point_group, &snf);

    enumerate_decomposed(composition, &translations, &point_group, callback);
}

/// Enumerate symmetry-inequivalent colorings across all inequivalent
/// supercell shapes at a given volume index.
///
/// Generates all HNF matrices at the given determinant, filters to
/// symmetry-inequivalent representatives under the parent point group
/// (Hart-Forcade 2008), then enumerates colorings for each.
///
/// # Panics
///
/// Panics if `parent_point_group` contains matrices not of dimension 2 or 3,
/// or if any `enumerate_supercell` call fails.
pub fn enumerate_at_index(
    parent_point_group: &[Vec<Vec<i64>>],
    index: usize,
    composition: &[usize],
    callback: &mut impl FnMut(&[u32]),
) {
    let dim = parent_point_group[0].len();
    let hnfs = match dim {
        2 => hnf_2d(index),
        3 => hnf_3d(index),
        _ => panic!("only 2D and 3D supercells are supported"),
    };

    let inequiv = inequivalent_hnfs(&hnfs, parent_point_group);

    for s in &inequiv {
        enumerate_supercell(s, composition, parent_point_group, callback);
    }
}

/// Enumerate primitive (non-superperiodic) symmetry-inequivalent colorings
/// of a supercell.
///
/// Like [`enumerate_supercell`] but excludes colorings that are invariant
/// under any non-identity translation in T (i.e., colorings representable
/// by a smaller supercell). This matches the enumlib convention.
///
/// # Panics
///
/// Same as [`enumerate_supercell`].
pub fn enumerate_supercell_primitive(
    supercell_matrix: &[Vec<i64>],
    composition: &[usize],
    parent_point_group: &[Vec<Vec<i64>>],
    callback: &mut impl FnMut(&[u32]),
) {
    let snf = smith_normal_form(supercell_matrix);
    let dims: Vec<usize> = snf
        .diagonal
        .iter()
        .map(|&d| {
            assert!(d > 0, "supercell matrix must be non-singular");
            d as usize
        })
        .collect();
    let n: usize = dims.iter().product();
    assert_eq!(
        composition.iter().sum::<usize>(),
        n,
        "composition sum {} != site count {n}",
        composition.iter().sum::<usize>()
    );

    let translations = product_translations(&dims);
    let point_group = transform_point_group(parent_point_group, &snf);

    enumerate_decomposed_primitive(composition, &translations, &point_group, callback);
}

/// Enumerate primitive (non-superperiodic) symmetry-inequivalent colorings
/// across all inequivalent supercell shapes at a given volume index.
///
/// Like [`enumerate_at_index`] but calls [`enumerate_supercell_primitive`]
/// for each inequivalent HNF, excluding superperiodic structures.
///
/// # Panics
///
/// Same as [`enumerate_at_index`].
pub fn enumerate_at_index_primitive(
    parent_point_group: &[Vec<Vec<i64>>],
    index: usize,
    composition: &[usize],
    callback: &mut impl FnMut(&[u32]),
) {
    let dim = parent_point_group[0].len();
    let hnfs = match dim {
        2 => hnf_2d(index),
        3 => hnf_3d(index),
        _ => panic!("only 2D and 3D supercells are supported"),
    };

    let inequiv = inequivalent_hnfs(&hnfs, parent_point_group);

    for s in &inequiv {
        enumerate_supercell_primitive(s, composition, parent_point_group, callback);
    }
}

/// Like [`enumerate_supercell`] but calls `callback` with `(coloring, orbit_size)`.
pub fn enumerate_supercell_with_multiplicity(
    supercell_matrix: &[Vec<i64>],
    composition: &[usize],
    parent_point_group: &[Vec<Vec<i64>>],
    callback: &mut impl FnMut(&[u32], u128),
) {
    let snf = smith_normal_form(supercell_matrix);
    let dims: Vec<usize> = snf
        .diagonal
        .iter()
        .map(|&d| {
            assert!(d > 0, "supercell matrix must be non-singular");
            d as usize
        })
        .collect();
    let n: usize = dims.iter().product();
    assert_eq!(
        composition.iter().sum::<usize>(),
        n,
        "composition sum {} != site count {n}",
        composition.iter().sum::<usize>()
    );

    let translations = product_translations(&dims);
    let point_group = transform_point_group(parent_point_group, &snf);

    enumerate_decomposed_with_multiplicity(composition, &translations, &point_group, callback);
}

/// Like [`enumerate_at_index`] but calls `callback` with `(coloring, orbit_size)`.
pub fn enumerate_at_index_with_multiplicity(
    parent_point_group: &[Vec<Vec<i64>>],
    index: usize,
    composition: &[usize],
    callback: &mut impl FnMut(&[u32], u128),
) {
    let dim = parent_point_group[0].len();
    let hnfs = match dim {
        2 => hnf_2d(index),
        3 => hnf_3d(index),
        _ => panic!("only 2D and 3D supercells are supported"),
    };

    let inequiv = inequivalent_hnfs(&hnfs, parent_point_group);

    for s in &inequiv {
        enumerate_supercell_with_multiplicity(s, composition, parent_point_group, callback);
    }
}

/// Like [`enumerate_supercell_primitive`] but calls `callback` with `(coloring, orbit_size)`.
pub fn enumerate_supercell_primitive_with_multiplicity(
    supercell_matrix: &[Vec<i64>],
    composition: &[usize],
    parent_point_group: &[Vec<Vec<i64>>],
    callback: &mut impl FnMut(&[u32], u128),
) {
    let snf = smith_normal_form(supercell_matrix);
    let dims: Vec<usize> = snf
        .diagonal
        .iter()
        .map(|&d| {
            assert!(d > 0, "supercell matrix must be non-singular");
            d as usize
        })
        .collect();
    let n: usize = dims.iter().product();
    assert_eq!(
        composition.iter().sum::<usize>(),
        n,
        "composition sum {} != site count {n}",
        composition.iter().sum::<usize>()
    );

    let translations = product_translations(&dims);
    let point_group = transform_point_group(parent_point_group, &snf);

    enumerate_decomposed_primitive_with_multiplicity(
        composition,
        &translations,
        &point_group,
        callback,
    );
}

/// Like [`enumerate_at_index_primitive`] but calls `callback` with `(coloring, orbit_size)`.
pub fn enumerate_at_index_primitive_with_multiplicity(
    parent_point_group: &[Vec<Vec<i64>>],
    index: usize,
    composition: &[usize],
    callback: &mut impl FnMut(&[u32], u128),
) {
    let dim = parent_point_group[0].len();
    let hnfs = match dim {
        2 => hnf_2d(index),
        3 => hnf_3d(index),
        _ => panic!("only 2D and 3D supercells are supported"),
    };

    let inequiv = inequivalent_hnfs(&hnfs, parent_point_group);

    for s in &inequiv {
        enumerate_supercell_primitive_with_multiplicity(
            s,
            composition,
            parent_point_group,
            callback,
        );
    }
}

/// Enumerate colorings of a supercell for all compositions in a concentration range.
///
/// `min[i]` and `max[i]` bound the count of species `i`. Enumerates all valid
/// compositions and calls [`enumerate_supercell`] for each.
pub fn enumerate_supercell_range(
    supercell_matrix: &[Vec<i64>],
    min: &[usize],
    max: &[usize],
    parent_point_group: &[Vec<Vec<i64>>],
    callback: &mut impl FnMut(&[u32]),
) {
    let n: usize = determinant(supercell_matrix).unsigned_abs() as usize;
    for comp in compositions_in_range(n, min, max) {
        enumerate_supercell(supercell_matrix, &comp, parent_point_group, callback);
    }
}

/// Enumerate colorings across all supercell shapes at a volume index
/// for all compositions in a concentration range.
pub fn enumerate_at_index_range(
    parent_point_group: &[Vec<Vec<i64>>],
    index: usize,
    min: &[usize],
    max: &[usize],
    callback: &mut impl FnMut(&[u32]),
) {
    let dim = parent_point_group[0].len();
    let hnfs = match dim {
        2 => hnf_2d(index),
        3 => hnf_3d(index),
        _ => panic!("only 2D and 3D supercells are supported"),
    };
    let inequiv = inequivalent_hnfs(&hnfs, parent_point_group);
    for comp in compositions_in_range(index, min, max) {
        for s in &inequiv {
            enumerate_supercell(s, &comp, parent_point_group, callback);
        }
    }
}

/// Like [`enumerate_at_index_range`] but excludes superperiodic structures.
pub fn enumerate_at_index_range_primitive(
    parent_point_group: &[Vec<Vec<i64>>],
    index: usize,
    min: &[usize],
    max: &[usize],
    callback: &mut impl FnMut(&[u32]),
) {
    let dim = parent_point_group[0].len();
    let hnfs = match dim {
        2 => hnf_2d(index),
        3 => hnf_3d(index),
        _ => panic!("only 2D and 3D supercells are supported"),
    };
    let inequiv = inequivalent_hnfs(&hnfs, parent_point_group);
    for comp in compositions_in_range(index, min, max) {
        for s in &inequiv {
            enumerate_supercell_primitive(s, &comp, parent_point_group, callback);
        }
    }
}

/// Filter HNF matrices to symmetry-inequivalent representatives.
///
/// Two HNF matrices S₁ and S₂ define equivalent superlattices under the
/// parent point group if there exists R in the point group such that
/// S₁⁻¹·R·S₂ is an integer matrix (Hart-Forcade 2008). Since both have
/// the same determinant and det(R) = ±1, the product is automatically
/// unimodular.
///
/// Uses integer arithmetic via adj(S₁)·R·S₂ / det(S₁) to avoid
/// floating-point issues.
pub fn inequivalent_hnfs(
    hnfs: &[Vec<Vec<i64>>],
    parent_pg: &[Vec<Vec<i64>>],
) -> Vec<Vec<Vec<i64>>> {
    let dim = if hnfs.is_empty() { 0 } else { hnfs[0].len() };

    match dim {
        2 => inequivalent_hnfs_2d(hnfs, parent_pg),
        3 => inequivalent_hnfs_3d(hnfs, parent_pg),
        _ => inequivalent_hnfs_generic(hnfs, parent_pg),
    }
}

/// 2D fast path: uses stack-allocated `[[i64;2];2]` matrices.
fn inequivalent_hnfs_2d(hnfs: &[Vec<Vec<i64>>], parent_pg: &[Vec<Vec<i64>>]) -> Vec<Vec<Vec<i64>>> {
    type M2 = [[i64; 2]; 2];

    #[inline]
    fn from_vecs(m: &[Vec<i64>]) -> M2 {
        [[m[0][0], m[0][1]], [m[1][0], m[1][1]]]
    }

    #[inline]
    fn mul2(a: &M2, b: &M2) -> M2 {
        [
            [
                a[0][0] * b[0][0] + a[0][1] * b[1][0],
                a[0][0] * b[0][1] + a[0][1] * b[1][1],
            ],
            [
                a[1][0] * b[0][0] + a[1][1] * b[1][0],
                a[1][0] * b[0][1] + a[1][1] * b[1][1],
            ],
        ]
    }

    #[inline]
    fn det2(m: &M2) -> i64 {
        m[0][0] * m[1][1] - m[0][1] * m[1][0]
    }

    // adj([[a,b],[c,d]]) = [[d,-b],[-c,a]]
    #[inline]
    fn adj2(m: &M2) -> M2 {
        [[m[1][1], -m[0][1]], [-m[1][0], m[0][0]]]
    }

    #[inline]
    fn divisible(product: &M2, det: i64) -> bool {
        product.iter().all(|row| row.iter().all(|&v| v % det == 0))
    }

    let pg: Vec<M2> = parent_pg.iter().map(|r| from_vecs(r)).collect();

    let mut rep_mats: Vec<M2> = Vec::new();
    let mut rep_vecs: Vec<Vec<Vec<i64>>> = Vec::new();

    for s_vec in hnfs {
        let s = from_vecs(s_vec);
        let rs_vals: Vec<M2> = pg.iter().map(|r| mul2(r, &s)).collect();

        let is_equivalent = rep_mats.iter().any(|rep| {
            let det = det2(rep);
            let adj = adj2(rep);
            rs_vals.iter().any(|rs| divisible(&mul2(&adj, rs), det))
        });

        if !is_equivalent {
            rep_mats.push(s);
            rep_vecs.push(s_vec.clone());
        }
    }

    rep_vecs
}

/// 3D fast path: uses stack-allocated `[[i64;3];3]` matrices to eliminate
/// heap allocation in the hot loop (adjugate + two mat-muls per check).
fn inequivalent_hnfs_3d(hnfs: &[Vec<Vec<i64>>], parent_pg: &[Vec<Vec<i64>>]) -> Vec<Vec<Vec<i64>>> {
    type M3 = [[i64; 3]; 3];

    #[inline]
    fn from_vecs(m: &[Vec<i64>]) -> M3 {
        [
            [m[0][0], m[0][1], m[0][2]],
            [m[1][0], m[1][1], m[1][2]],
            [m[2][0], m[2][1], m[2][2]],
        ]
    }

    #[inline]
    fn mul3(a: &M3, b: &M3) -> M3 {
        let mut c = [[0i64; 3]; 3];
        for i in 0..3 {
            for k in 0..3 {
                if a[i][k] == 0 {
                    continue;
                }
                for j in 0..3 {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        c
    }

    #[inline]
    fn det3(m: &M3) -> i64 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    // Adjugate (classical adjoint) = transpose of cofactor matrix.
    #[inline]
    fn adj3(m: &M3) -> M3 {
        let [a, b, c] = [m[0][0], m[0][1], m[0][2]];
        let [d, e, f] = [m[1][0], m[1][1], m[1][2]];
        let [g, h, k] = [m[2][0], m[2][1], m[2][2]];
        [
            [e * k - f * h, c * h - b * k, b * f - c * e],
            [f * g - d * k, a * k - c * g, c * d - a * f],
            [d * h - e * g, b * g - a * h, a * e - b * d],
        ]
    }

    #[inline]
    fn divisible(product: &M3, det: i64) -> bool {
        product.iter().all(|row| row.iter().all(|&v| v % det == 0))
    }

    let pg: Vec<M3> = parent_pg.iter().map(|r| from_vecs(r)).collect();

    // For each HNF s, precompute r*s for all r so the rep loop doesn't redo it.
    let mut rep_mats: Vec<M3> = Vec::new();
    let mut rep_vecs: Vec<Vec<Vec<i64>>> = Vec::new();

    for s_vec in hnfs {
        let s = from_vecs(s_vec);
        // Precompute r*s for every r — this hoists the (r,s) mul out of the rep loop.
        let rs_vals: Vec<M3> = pg.iter().map(|r| mul3(r, &s)).collect();

        let is_equivalent = rep_mats.iter().any(|rep| {
            let det = det3(rep);
            let adj = adj3(rep);
            rs_vals.iter().any(|rs| divisible(&mul3(&adj, rs), det))
        });

        if !is_equivalent {
            rep_mats.push(s);
            rep_vecs.push(s_vec.clone());
        }
    }

    rep_vecs
}

/// Generic fallback for 2D (or other dimensions) using Vec<Vec<i64>>.
fn inequivalent_hnfs_generic(
    hnfs: &[Vec<Vec<i64>>],
    parent_pg: &[Vec<Vec<i64>>],
) -> Vec<Vec<Vec<i64>>> {
    let mut representatives: Vec<Vec<Vec<i64>>> = Vec::new();

    for s in hnfs {
        let is_equivalent = representatives.iter().any(|rep| {
            let det = determinant(rep);
            let adj = adjugate(rep);
            parent_pg.iter().any(|r| {
                let rs = mat_mul(r, s);
                let product = mat_mul(&adj, &rs);
                product.iter().all(|row| row.iter().all(|&v| v % det == 0))
            })
        });
        if !is_equivalent {
            representatives.push(s.clone());
        }
    }

    representatives
}

#[cfg(test)]
fn hnf_equivalent(rep: &[Vec<i64>], r: &[Vec<i64>], s: &[Vec<i64>]) -> bool {
    let det = determinant(rep);
    debug_assert!(det != 0);
    let adj = adjugate(rep);
    let rs = mat_mul(r, s);
    let product = mat_mul(&adj, &rs);
    product.iter().all(|row| row.iter().all(|&v| v % det == 0))
}

/// Transform parent point group matrices into permutations on SNF-basis
/// site indices.
///
/// Given d×d integer matrices representing the parent point group and the
/// Smith Normal Form of the supercell matrix (D = L·S·R), computes
/// M = L·p·L⁻¹ for each parent element p, filters elements that don't
/// preserve the quotient group ℤ_{d₀} × ℤ_{d₁} × ..., and converts
/// survivors into permutations on row-major site indices.
///
/// # Preconditions
///
/// `parent_pg` must contain the identity matrix and be closed under
/// matrix multiplication (i.e., form a valid group). Violating these
/// will produce an incomplete or invalid subgroup.
///
/// The result always includes the identity and forms a valid subgroup.
pub fn transform_point_group(
    parent_pg: &[Vec<Vec<i64>>],
    snf: &SmithNormalForm,
) -> Vec<Permutation> {
    let dims: Vec<usize> = snf.diagonal.iter().map(|&x| x as usize).collect();
    let n: usize = dims.iter().product();
    let l_inv = unimodular_inverse(&snf.left);

    let mut perms = Vec::new();
    for p in parent_pg {
        let lp = mat_mul(&snf.left, p);
        let m = mat_mul(&lp, &l_inv);

        if !preserves_quotient(&m, &snf.diagonal) {
            continue;
        }

        let image: Vec<u32> = (0..n)
            .map(|idx| {
                let coords = index_to_coords(idx, &dims);
                let new_coords = mat_vec_mod(&m, &coords, &dims);
                coords_to_index(&new_coords, &dims) as u32
            })
            .collect();

        perms.push(Permutation::new(image));
    }

    perms.sort_by(|a, b| a.image().cmp(b.image()));
    perms.dedup();

    debug_assert!(
        perms.contains(&Permutation::identity(n)),
        "transformed point group must contain the identity"
    );

    perms
}

/// Transform space group operations into permutations on multi-site supercells.
///
/// For a parent structure with `n_orbits` Wyckoff orbits (atoms in the unit cell),
/// the supercell has `n_total = n_orbits * V` sites where `V = product(snf.diagonal)`.
/// Sites use flat indexing: site `(j, y)` maps to index `j * V + row_major(y)`.
///
/// Each space group operation `(R, t)` induces a permutation:
///   `(j, y) → (orbit_maps[op][j], M_snf · y + offset[op][j] mod D)`
/// where `M_snf = L · R · L⁻¹` in the SNF basis, and `offset[op][j]` is the
/// SNF-basis integer offset from the translation part.
///
/// # Arguments
///
/// * `rotations` — Rotation matrices R for each space group operation (d×d).
/// * `orbit_maps` — For each operation, `orbit_maps[op][j]` = index of orbit
///   that orbit `j` maps to under this operation.
/// * `offsets` — For each operation, `offsets[op][j]` = integer translation offset
///   `n_jk = round(R·x_j + t - x_k)` where `k = orbit_maps[op][j]`.
///   These are in the parent lattice basis.
/// * `snf` — Smith Normal Form of the supercell matrix.
/// * `n_orbits` — Number of atoms (Wyckoff orbits) in the parent unit cell.
pub fn transform_space_group_multisite(
    rotations: &[Vec<Vec<i64>>],
    orbit_maps: &[Vec<usize>],
    offsets: &[Vec<Vec<i64>>],
    snf: &SmithNormalForm,
    n_orbits: usize,
) -> Vec<Permutation> {
    let dims: Vec<usize> = snf.diagonal.iter().map(|&x| x as usize).collect();
    let v: usize = dims.iter().product();
    let n_total = n_orbits * v;
    let d = dims.len();
    let l_inv = unimodular_inverse(&snf.left);

    let mut perms = Vec::new();
    for (op_idx, rot) in rotations.iter().enumerate() {
        let lp = mat_mul(&snf.left, rot);
        let m_snf = mat_mul(&lp, &l_inv);

        if !preserves_quotient(&m_snf, &snf.diagonal) {
            continue;
        }

        let image: Vec<u32> = (0..n_total)
            .map(|flat| {
                let orbit_j = flat / v;
                let site_y = flat % v;
                let coords_y = index_to_coords(site_y, &dims);

                let orbit_k = orbit_maps[op_idx][orbit_j];

                // SNF-basis offset: L · n_jk
                let n_jk = &offsets[op_idx][orbit_j];
                let l_offset: Vec<i64> = (0..d)
                    .map(|i| (0..d).map(|j| snf.left[i][j] * n_jk[j]).sum())
                    .collect();

                // new_coords = M_snf · y + L · n_jk  (mod D)
                let mut new_coords = mat_vec_mod(&m_snf, &coords_y, &dims);
                for i in 0..d {
                    new_coords[i] = (new_coords[i] + l_offset[i]).rem_euclid(dims[i] as i64);
                }

                let new_site = coords_to_index(&new_coords, &dims);
                (orbit_k * v + new_site) as u32
            })
            .collect();

        perms.push(Permutation::new(image));
    }

    perms.sort_by(|a, b| a.image().cmp(b.image()));
    perms.dedup();

    debug_assert!(
        perms.contains(&Permutation::identity(n_total)),
        "transformed space group must contain the identity"
    );

    perms
}

/// Enumerate symmetry-inequivalent colorings of a multi-site supercell.
///
/// For parent structures with multiple atoms per unit cell, the space group
/// `G = T ⋊ P` acts on all Wyckoff orbits simultaneously. The total number
/// of sites is `n_orbits * |det(supercell_matrix)|`.
///
/// # Arguments
///
/// * `supercell_matrix` — d×d integer matrix defining the supercell.
/// * `composition` — Number of sites per species, summing to `n_orbits * |det(M)|`.
/// * `rotations` — Rotation matrices R for each space group operation.
/// * `orbit_maps` — For each operation, the permutation of orbits.
/// * `offsets` — For each operation, integer offsets per orbit (parent basis).
/// * `n_orbits` — Number of atoms in the parent unit cell.
/// * `callback` — Called once per inequivalent coloring.
pub fn enumerate_supercell_multisite(
    supercell_matrix: &[Vec<i64>],
    composition: &[usize],
    rotations: &[Vec<Vec<i64>>],
    orbit_maps: &[Vec<usize>],
    offsets: &[Vec<Vec<i64>>],
    n_orbits: usize,
    callback: &mut impl FnMut(&[u32]),
) {
    let snf = smith_normal_form(supercell_matrix);
    let dims: Vec<usize> = snf
        .diagonal
        .iter()
        .map(|&d| {
            assert!(d > 0, "supercell matrix must be non-singular");
            d as usize
        })
        .collect();
    let v: usize = dims.iter().product();
    let n_total = n_orbits * v;
    assert_eq!(
        composition.iter().sum::<usize>(),
        n_total,
        "composition sum {} != total site count {n_total}",
        composition.iter().sum::<usize>()
    );

    let translations = product_translations_multisite(&dims, n_orbits);
    let point_group =
        transform_space_group_multisite(rotations, orbit_maps, offsets, &snf, n_orbits);

    enumerate_decomposed(composition, &translations, &point_group, callback);
}

/// Enumerate primitive (non-superperiodic) colorings of a multi-site supercell.
///
/// Like [`enumerate_supercell_multisite`] but excludes colorings that are
/// invariant under any non-identity supercell translation. Such colorings
/// can be represented by a smaller supercell and would be counted at a
/// smaller volume.
///
/// # Panics
///
/// Panics if dimensions are inconsistent or if the Pólya cross-check fails.
pub fn enumerate_supercell_multisite_primitive(
    supercell_matrix: &[Vec<i64>],
    composition: &[usize],
    rotations: &[Vec<Vec<i64>>],
    orbit_maps: &[Vec<usize>],
    offsets: &[Vec<Vec<i64>>],
    n_orbits: usize,
    callback: &mut impl FnMut(&[u32]),
) {
    let snf = smith_normal_form(supercell_matrix);
    let dims: Vec<usize> = snf
        .diagonal
        .iter()
        .map(|&d| {
            assert!(d > 0, "supercell matrix must be non-singular");
            d as usize
        })
        .collect();
    let v: usize = dims.iter().product();
    let n_total = n_orbits * v;
    assert_eq!(
        composition.iter().sum::<usize>(),
        n_total,
        "composition sum {} != total site count {n_total}",
        composition.iter().sum::<usize>()
    );

    let translations = product_translations_multisite(&dims, n_orbits);
    let point_group =
        transform_space_group_multisite(rotations, orbit_maps, offsets, &snf, n_orbits);

    enumerate_decomposed_primitive(composition, &translations, &point_group, callback);
}

/// Enumerate constrained colorings of a multi-site supercell.
///
/// Like [`enumerate_supercell_multisite`] but each Wyckoff orbit has an
/// independent set of allowed species. `orbit_allowed[j][c]` is `true`
/// iff species `c` may occupy orbit `j`. All `V` supercell images of
/// orbit `j` share the same allowed set.
///
/// # Panics
///
/// Panics if dimensions are inconsistent or if the Pólya cross-check fails.
#[allow(clippy::too_many_arguments)]
pub fn enumerate_supercell_multisite_constrained(
    supercell_matrix: &[Vec<i64>],
    composition: &[usize],
    rotations: &[Vec<Vec<i64>>],
    orbit_maps: &[Vec<usize>],
    offsets: &[Vec<Vec<i64>>],
    n_orbits: usize,
    orbit_allowed: &[Vec<bool>],
    callback: &mut impl FnMut(&[u32]),
) {
    assert_eq!(
        orbit_allowed.len(),
        n_orbits,
        "orbit_allowed must have one entry per orbit"
    );
    let snf = smith_normal_form(supercell_matrix);
    let dims: Vec<usize> = snf
        .diagonal
        .iter()
        .map(|&d| {
            assert!(d > 0, "supercell matrix must be non-singular");
            d as usize
        })
        .collect();
    let v: usize = dims.iter().product();
    let n_total = n_orbits * v;
    assert_eq!(
        composition.iter().sum::<usize>(),
        n_total,
        "composition sum {} != total site count {n_total}",
        composition.iter().sum::<usize>()
    );

    let allowed: Vec<Vec<bool>> = (0..n_total)
        .map(|flat| orbit_allowed[flat / v].clone())
        .collect();

    let translations = product_translations_multisite(&dims, n_orbits);
    let point_group =
        transform_space_group_multisite(rotations, orbit_maps, offsets, &snf, n_orbits);

    enumerate_decomposed_constrained(composition, &translations, &point_group, &allowed, callback);
}

/// Compute fractional positions for all sites in a multi-site supercell.
///
/// Returns `n_orbits * V` positions in fractional coordinates of the supercell
/// lattice. Site ordering matches the flat indexing: site `(j, y)` maps to
/// index `j * V + row_major(y)`, and its position is `(x_j + grid_point_y) · M⁻¹`.
///
/// # Arguments
///
/// * `supercell_matrix` — d×d integer matrix defining the supercell.
/// * `parent_positions` — Fractional positions of the `n_orbits` parent atoms
///   in the parent lattice basis.
pub fn supercell_fractional_positions_multisite(
    supercell_matrix: &[Vec<i64>],
    parent_positions: &[[f64; 3]],
) -> Vec<Vec<f64>> {
    let snf = smith_normal_form(supercell_matrix);
    let dims: Vec<usize> = snf.diagonal.iter().map(|&d| d as usize).collect();
    let v: usize = dims.iter().product();
    let d = dims.len();
    let r_inv = unimodular_inverse(&snf.right);
    let det = determinant(supercell_matrix);
    let adj_m = adjugate(supercell_matrix);

    let n_orbits = parent_positions.len();
    let mut positions = Vec::with_capacity(n_orbits * v);

    for parent_pos in parent_positions {
        for idx in 0..v {
            let coords = index_to_coords(idx, &dims);
            let mut parent_int = vec![0f64; d];
            for i in 0..d {
                parent_int[i] = parent_pos[i];
                for j in 0..d {
                    parent_int[i] += r_inv[j][i] as f64 * coords[j] as f64;
                }
            }
            let mut frac = vec![0.0f64; d];
            for i in 0..d {
                let mut num = 0.0f64;
                for j in 0..d {
                    num += adj_m[j][i] as f64 * parent_int[j];
                }
                frac[i] = (num / det as f64).rem_euclid(1.0);
            }
            positions.push(frac);
        }
    }

    positions
}

/// Generate a matrix group from generators by BFS closure.
///
/// Starting from the identity and the given generators, multiplies
/// all pairs until no new elements are found. Returns the full group
/// as a list of d×d integer matrices.
pub fn generate_matrix_group(generators: &[Vec<Vec<i64>>], dim: usize) -> Vec<Vec<Vec<i64>>> {
    let identity: Vec<Vec<i64>> = (0..dim)
        .map(|i| {
            let mut row = vec![0i64; dim];
            row[i] = 1;
            row
        })
        .collect();

    crate::perm::generate_group_elements(identity, generators, |a, b| mat_mul(a, b))
}

pub fn mat_mul(a: &[Vec<i64>], b: &[Vec<i64>]) -> Vec<Vec<i64>> {
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

pub(crate) fn determinant(m: &[Vec<i64>]) -> i64 {
    match m.len() {
        0 => 1,
        1 => m[0][0],
        2 => m[0][0] * m[1][1] - m[0][1] * m[1][0],
        3 => {
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        }
        n => {
            let mut det = 0i64;
            for j in 0..n {
                let minor = minor_matrix(m, 0, j);
                let sign = if j % 2 == 0 { 1 } else { -1 };
                det += sign * m[0][j] * determinant(&minor);
            }
            det
        }
    }
}

pub(crate) fn minor_matrix(m: &[Vec<i64>], row: usize, col: usize) -> Vec<Vec<i64>> {
    m.iter()
        .enumerate()
        .filter(|&(i, _)| i != row)
        .map(|(_, r)| {
            r.iter()
                .enumerate()
                .filter(|&(j, _)| j != col)
                .map(|(_, &v)| v)
                .collect()
        })
        .collect()
}

pub(crate) fn adjugate(m: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let n = m.len();
    (0..n)
        .map(|j| {
            (0..n)
                .map(|i| {
                    let sign = if (i + j) % 2 == 0 { 1 } else { -1 };
                    sign * determinant(&minor_matrix(m, i, j))
                })
                .collect()
        })
        .collect()
}

fn unimodular_inverse(m: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let det = determinant(m);
    assert!(det == 1 || det == -1, "matrix is not unimodular: det={det}");
    let adj = adjugate(m);
    if det == -1 {
        adj.into_iter()
            .map(|row| row.into_iter().map(|v| -v).collect())
            .collect()
    } else {
        adj
    }
}

/// Check if M preserves the quotient ℤ_{d₀} × ... × ℤ_{d_{k-1}}.
///
/// The condition is that D⁻¹·M·D has integer entries, where
/// D = diag(diagonal). Entry (i,j) of D⁻¹·M·D is M[i][j]·d[j]/d[i],
/// so we check d[i] | (M[i][j]·d[j]).
fn preserves_quotient(m: &[Vec<i64>], diagonal: &[i64]) -> bool {
    let n = m.len();
    for i in 0..n {
        for j in 0..n {
            if (m[i][j] * diagonal[j]) % diagonal[i] != 0 {
                return false;
            }
        }
    }
    true
}

/// Compute the fractional positions of all supercell sites.
///
/// Returns `n` positions in fractional coordinates of the **supercell**
/// lattice (L_super = M * L_parent). Site ordering matches the coloring
/// index used by the enumeration functions.
pub fn supercell_fractional_positions(supercell_matrix: &[Vec<i64>]) -> Vec<Vec<f64>> {
    let snf = smith_normal_form(supercell_matrix);
    let dims: Vec<usize> = snf.diagonal.iter().map(|&d| d as usize).collect();
    let n: usize = dims.iter().product();
    let d = dims.len();
    let r_inv = unimodular_inverse(&snf.right);
    let det = determinant(supercell_matrix);
    let adj_m = adjugate(supercell_matrix);

    (0..n)
        .map(|idx| {
            let coords = index_to_coords(idx, &dims);
            // Row-convention quotient: x · R ≡ y mod D, so x = y · R⁻¹.
            // In column indexing: parent_int[i] = Σ_j R⁻¹[j][i] · coords[j]
            let mut parent_int = vec![0i64; d];
            for i in 0..d {
                for j in 0..d {
                    parent_int[i] += r_inv[j][i] * coords[j];
                }
            }
            // Supercell fractional: r = x · M⁻¹ (row convention)
            // r[i] = Σ_j x[j] · adj(M)[j][i] / det(M)
            let mut frac = vec![0.0f64; d];
            for i in 0..d {
                let mut num = 0i64;
                for j in 0..d {
                    num += adj_m[j][i] * parent_int[j];
                }
                frac[i] = (num as f64 / det as f64).rem_euclid(1.0);
            }
            frac
        })
        .collect()
}

fn index_to_coords(mut idx: usize, dims: &[usize]) -> Vec<i64> {
    let d = dims.len();
    let mut coords = vec![0i64; d];
    for k in (0..d).rev() {
        coords[k] = (idx % dims[k]) as i64;
        idx /= dims[k];
    }
    coords
}

fn coords_to_index(coords: &[i64], dims: &[usize]) -> usize {
    let mut idx = 0;
    for k in 0..dims.len() {
        idx = idx * dims[k] + coords[k] as usize;
    }
    idx
}

fn mat_vec_mod(m: &[Vec<i64>], v: &[i64], dims: &[usize]) -> Vec<i64> {
    let n = m.len();
    let mut result = vec![0i64; n];
    for i in 0..n {
        for j in 0..n {
            result[i] += m[i][j] * v[j];
        }
        result[i] = result[i].rem_euclid(dims[i] as i64);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enumerate::enumerate;
    use crate::t_canon::product_translations_2d;

    fn identity_matrix(n: usize) -> Vec<Vec<i64>> {
        (0..n)
            .map(|i| {
                let mut row = vec![0i64; n];
                row[i] = 1;
                row
            })
            .collect()
    }

    // --- Unimodular inverse ---

    #[test]
    fn inverse_identity_2d() {
        let id = identity_matrix(2);
        let inv = unimodular_inverse(&id);
        assert_eq!(inv, id);
    }

    #[test]
    fn inverse_identity_3d() {
        let id = identity_matrix(3);
        let inv = unimodular_inverse(&id);
        assert_eq!(inv, id);
    }

    #[test]
    fn inverse_from_snf_2d() {
        for m in [
            vec![vec![2, 0], vec![0, 3]],
            vec![vec![2, 1], vec![0, 3]],
            vec![vec![3, 1], vec![1, 2]],
        ] {
            let snf = smith_normal_form(&m);
            let l_inv = unimodular_inverse(&snf.left);
            let product = mat_mul(&snf.left, &l_inv);
            assert_eq!(product, identity_matrix(2), "L·L⁻¹ ≠ I for {m:?}");

            let r_inv = unimodular_inverse(&snf.right);
            let product = mat_mul(&snf.right, &r_inv);
            assert_eq!(product, identity_matrix(2), "R·R⁻¹ ≠ I for {m:?}");
        }
    }

    #[test]
    fn inverse_from_snf_3d() {
        for m in [
            vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 2]],
            vec![vec![2, 1, 0], vec![0, 3, 0], vec![0, 0, 1]],
            vec![vec![2, 4, 4], vec![-6, 6, 12], vec![10, -4, -16]],
        ] {
            let snf = smith_normal_form(&m);
            let l_inv = unimodular_inverse(&snf.left);
            let product = mat_mul(&snf.left, &l_inv);
            assert_eq!(product, identity_matrix(3), "L·L⁻¹ ≠ I for {m:?}");
        }
    }

    // --- Point group transformation ---

    #[test]
    fn transform_trivial_pg() {
        let snf = smith_normal_form(&[vec![2, 0], vec![0, 3]]);
        let pg = vec![identity_matrix(2)];
        let perms = transform_point_group(&pg, &snf);
        assert_eq!(perms.len(), 1);
        assert_eq!(perms[0], Permutation::identity(6));
    }

    #[test]
    fn transform_inversion_2d() {
        // S = diag(2,3), SNF = diag(1,6) or diag(2,3) depending on matrix
        // Use S = diag(2,3) → SNF diagonal = [1,6] with L that diagonalizes
        // Actually for diagonal S, SNF is diag(gcd structure)...
        // Let's use S = diag(2,3): det=6, SNF should be diag(1,6) since gcd(2,3)=1.
        // Actually no: for S=diag(2,3), SNF diagonal = [1, 6] since d0|d1 and d0*d1=6.
        // Hmm let me just compute it.
        let s = vec![vec![2, 0], vec![0, 3]];
        let snf = smith_normal_form(&s);
        // SNF of diag(2,3): since gcd(2,3)=1, SNF = diag(1,6)
        assert_eq!(snf.diagonal, vec![1, 6]);

        // Inversion: -I
        let neg_i = vec![vec![-1, 0], vec![0, -1]];
        let pg = vec![identity_matrix(2), neg_i];
        let perms = transform_point_group(&pg, &snf);

        // Both should preserve quotient; -I on ℤ_1×ℤ_6 maps (a,b)→(-a,-b)=(0, 6-b mod 6)
        assert_eq!(perms.len(), 2);
        // The non-identity permutation should map index k → (6-k) mod 6
        // Since dims=[1,6], sites are just 0..5 (first coord always 0)
        // -I maps (0,b) → (0, -b mod 6) = (0, (6-b) mod 6)
        let expected_image: Vec<u32> = (0..6).map(|b| ((6 - b) % 6) as u32).collect();
        let inv_perm = perms
            .iter()
            .find(|p| p != &&Permutation::identity(6))
            .unwrap();
        assert_eq!(inv_perm.image(), &expected_image);
    }

    #[test]
    fn transform_c4_3x3() {
        let s = vec![vec![3, 0], vec![0, 3]];
        let snf = smith_normal_form(&s);
        assert_eq!(snf.diagonal, vec![3, 3]);

        // C4 rotation: [[0,-1],[1,0]]
        let c4 = vec![vec![0, -1], vec![1, 0]];
        let pg = vec![identity_matrix(2), c4];
        let perms = transform_point_group(&pg, &snf);

        // C4 should preserve ℤ_3×ℤ_3, giving 2 distinct permutations
        assert_eq!(perms.len(), 2);

        // Verify the C4 permutation: (a,b) → (-b mod 3, a mod 3)
        // For L=I (diagonal S), M = L·C4·L⁻¹ = C4
        let c4_perm = perms
            .iter()
            .find(|p| p != &&Permutation::identity(9))
            .unwrap();
        for a in 0..3i64 {
            for b in 0..3i64 {
                let idx = (a * 3 + b) as usize;
                let new_a = (-b).rem_euclid(3);
                let new_b = a;
                let expected = (new_a * 3 + new_b) as u32;
                assert_eq!(
                    c4_perm.apply(idx as u32),
                    expected,
                    "C4 mismatch at ({a},{b})"
                );
            }
        }
    }

    #[test]
    fn transform_filters_incompatible() {
        // S = diag(2,3): SNF = diag(1,6)
        // C4 = [[0,-1],[1,0]]. Check: does C4 preserve ℤ_1×ℤ_6?
        // M = L·C4·L⁻¹. Need to check D⁻¹·M·D integer where D=diag(1,6).
        // If L=I (not necessarily true), M=C4=[[0,-1],[1,0]].
        // D⁻¹·M·D: entry (0,1) = (-1)*6/1 = -6 (int), entry (1,0) = 1*1/6 = 1/6 (NOT int)
        // So C4 should be filtered out for the 2×3 supercell.
        let s = vec![vec![2, 0], vec![0, 3]];
        let snf = smith_normal_form(&s);

        let c4 = vec![vec![0, -1], vec![1, 0]];
        let pg = vec![identity_matrix(2), c4];
        let perms = transform_point_group(&pg, &snf);

        // Only identity should survive
        assert_eq!(perms.len(), 1);
        assert_eq!(perms[0], Permutation::identity(6));
    }

    // --- enumerate_supercell ---

    #[test]
    fn supercell_diagonal_trivial_pg() {
        let s = vec![vec![2, 0], vec![0, 3]];
        let pg = vec![identity_matrix(2)];
        let mut count = 0;
        enumerate_supercell(&s, &[3, 3], &pg, &mut |_| count += 1);
        // Pólya cross-check is built into enumerate_decomposed
        assert!(count > 0);
    }

    #[test]
    fn supercell_with_inversion_2d() {
        let s = vec![vec![3, 0], vec![0, 3]];
        let neg_i = vec![vec![-1, 0], vec![0, -1]];
        let pg = vec![identity_matrix(2), neg_i];

        let mut supercell_count = 0;
        enumerate_supercell(&s, &[4, 5], &pg, &mut |_| supercell_count += 1);

        // Cross-check: build full group manually
        let translations = product_translations_2d(3, 3);
        let snf = smith_normal_form(&s);
        let point_perms = transform_point_group(&pg, &snf);
        let mut full_group = Vec::new();
        for t in &translations {
            for p in &point_perms {
                full_group.push(p.compose(t));
            }
        }
        // Deduplicate
        full_group.sort_by(|a, b| a.image().cmp(b.image()));
        full_group.dedup();

        let mut direct_count = 0;
        enumerate(&[4, 5], &full_group, &mut |_| direct_count += 1);
        assert_eq!(supercell_count, direct_count);
    }

    #[test]
    fn supercell_c4_3x3() {
        let s = vec![vec![3, 0], vec![0, 3]];
        // C4 group: generated by [[0,-1],[1,0]]
        let c4 = vec![vec![0, -1], vec![1, 0]];
        let c2 = vec![vec![-1, 0], vec![0, -1]];
        let c4_inv = vec![vec![0, 1], vec![-1, 0]];
        let pg = vec![identity_matrix(2), c4, c2, c4_inv];

        let mut supercell_count = 0;
        enumerate_supercell(&s, &[4, 5], &pg, &mut |_| supercell_count += 1);

        // Cross-check via full group
        let snf = smith_normal_form(&s);
        let translations = product_translations_2d(3, 3);
        let point_perms = transform_point_group(&pg, &snf);
        let mut full_group = Vec::new();
        for t in &translations {
            for p in &point_perms {
                full_group.push(p.compose(t));
            }
        }
        full_group.sort_by(|a, b| a.image().cmp(b.image()));
        full_group.dedup();

        let mut direct_count = 0;
        enumerate(&[4, 5], &full_group, &mut |_| direct_count += 1);
        assert_eq!(supercell_count, direct_count);
    }

    #[test]
    fn supercell_off_diagonal() {
        // Non-diagonal supercell: S = [[2,1],[0,2]], det=4
        let s = vec![vec![2, 1], vec![0, 2]];
        let pg = vec![identity_matrix(2)];
        let mut count = 0;
        enumerate_supercell(&s, &[2, 2], &pg, &mut |_| count += 1);
        assert!(count > 0);
    }

    #[test]
    fn supercell_3d_diagonal() {
        let s = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 2]];
        let pg = vec![identity_matrix(3)];
        let mut count = 0;
        enumerate_supercell(&s, &[4, 4], &pg, &mut |_| count += 1);

        // Cross-check: this is ℤ_2×ℤ_2×ℤ_2 translation group
        let translations = product_translations(&[2, 2, 2]);
        let point_group = vec![Permutation::identity(8)];
        let mut direct_count = 0;
        enumerate_decomposed(&[4, 4], &translations, &point_group, &mut |_| {
            direct_count += 1;
        });
        assert_eq!(count, direct_count);
    }

    #[test]
    fn supercell_3d_with_inversion() {
        let s = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 2]];
        let neg_i = vec![vec![-1, 0, 0], vec![0, -1, 0], vec![0, 0, -1]];
        let pg = vec![identity_matrix(3), neg_i];

        let mut count = 0;
        enumerate_supercell(&s, &[4, 4], &pg, &mut |_| count += 1);

        // For ℤ_2³, -I acts trivially, so point group contributes nothing new
        let translations = product_translations(&[2, 2, 2]);
        let mut direct_count = 0;
        enumerate_decomposed(
            &[4, 4],
            &translations,
            &[Permutation::identity(8)],
            &mut |_| direct_count += 1,
        );
        assert_eq!(count, direct_count);
    }

    #[test]
    fn supercell_1d_with_inversion() {
        // 1D supercell: S = [[6]], parent PG = {I, -I}
        // This should reproduce dihedral D_6 enumeration on 6 sites.
        let s = vec![vec![6]];
        let neg_i = vec![vec![-1]];
        let pg = vec![vec![vec![1]], neg_i];

        let snf = smith_normal_form(&s);
        assert_eq!(snf.diagonal, vec![6]);

        let point_perms = transform_point_group(&pg, &snf);
        assert_eq!(point_perms.len(), 2);

        // Cross-check: build full group = T ⋊ P manually
        let translations = product_translations(&[6]);
        let mut full_group = Vec::new();
        for t in &translations {
            for p in &point_perms {
                full_group.push(p.compose(t));
            }
        }
        full_group.sort_by(|a, b| a.image().cmp(b.image()));
        full_group.dedup();
        assert_eq!(full_group.len(), 12); // |D_6| = 12

        for comp in [&[3, 3][..], &[2, 4], &[1, 5], &[0, 6]] {
            let mut supercell_count = 0;
            enumerate_supercell(&s, comp, &pg, &mut |_| supercell_count += 1);

            let mut direct_count = 0;
            enumerate(comp, &full_group, &mut |_| direct_count += 1);

            assert_eq!(
                supercell_count, direct_count,
                "1D supercell mismatch for composition {comp:?}"
            );
        }
    }

    // --- HNF equivalence ---

    fn generate_oh() -> Vec<Vec<Vec<i64>>> {
        let perms: [[usize; 3]; 6] = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];
        let mut ops = Vec::new();
        for perm in &perms {
            for signs in 0..8u8 {
                let mut m = vec![vec![0i64; 3]; 3];
                for i in 0..3 {
                    let sign = if signs & (1 << i) != 0 { -1 } else { 1 };
                    m[i][perm[i]] = sign;
                }
                ops.push(m);
            }
        }
        ops
    }

    fn brute_force_orbit_count(hnfs: &[Vec<Vec<i64>>], pg: &[Vec<Vec<i64>>]) -> usize {
        let n = hnfs.len();
        let mut parent = (0..n).collect::<Vec<_>>();

        fn find(parent: &mut [usize], x: usize) -> usize {
            if parent[x] != x {
                parent[x] = find(parent, parent[x]);
            }
            parent[x]
        }

        fn union(parent: &mut [usize], a: usize, b: usize) {
            let ra = find(parent, a);
            let rb = find(parent, b);
            if ra != rb {
                parent[ra] = rb;
            }
        }

        for i in 0..n {
            for j in (i + 1)..n {
                for r in pg {
                    if hnf_equivalent(&hnfs[i], r, &hnfs[j]) {
                        union(&mut parent, i, j);
                        break;
                    }
                }
            }
        }

        (0..n).filter(|&i| find(&mut parent, i) == i).count()
    }

    #[test]
    fn hnf_identity_pg_is_noop_2d() {
        let pg = vec![identity_matrix(2)];
        for det in 1..=8 {
            let hnfs = crate::hnf::hnf_2d(det);
            let inequiv = inequivalent_hnfs(&hnfs, &pg);
            assert_eq!(
                inequiv.len(),
                hnfs.len(),
                "identity PG should keep all HNFs at det={det}"
            );
        }
    }

    #[test]
    fn hnf_identity_pg_is_noop_3d() {
        let pg = vec![identity_matrix(3)];
        for det in 1..=6 {
            let hnfs = crate::hnf::hnf_3d(det);
            let inequiv = inequivalent_hnfs(&hnfs, &pg);
            assert_eq!(
                inequiv.len(),
                hnfs.len(),
                "identity PG should keep all HNFs at det={det}"
            );
        }
    }

    #[test]
    fn hnf_inversion_same_as_identity_3d() {
        let pg_id = vec![identity_matrix(3)];
        let neg_i = vec![vec![-1, 0, 0], vec![0, -1, 0], vec![0, 0, -1]];
        let pg_inv = vec![identity_matrix(3), neg_i];
        for det in 1..=6 {
            let hnfs = crate::hnf::hnf_3d(det);
            let count_id = inequivalent_hnfs(&hnfs, &pg_id).len();
            let count_inv = inequivalent_hnfs(&hnfs, &pg_inv).len();
            assert_eq!(
                count_id, count_inv,
                "inversion should not reduce HNF count at det={det}"
            );
        }
    }

    #[test]
    fn hnf_filter_matches_brute_force_2d() {
        let test_pgs: Vec<Vec<Vec<Vec<i64>>>> = vec![
            vec![identity_matrix(2)],
            vec![identity_matrix(2), vec![vec![-1, 0], vec![0, -1]]],
            vec![identity_matrix(2), vec![vec![0, 1], vec![1, 0]]],
            {
                let c4 = vec![vec![0, -1], vec![1, 0]];
                let c2 = vec![vec![-1, 0], vec![0, -1]];
                let c4i = vec![vec![0, 1], vec![-1, 0]];
                vec![identity_matrix(2), c4, c2, c4i]
            },
        ];
        for pg in &test_pgs {
            for det in 1..=8 {
                let hnfs = crate::hnf::hnf_2d(det);
                let filter_count = inequivalent_hnfs(&hnfs, pg).len();
                let orbit_count = brute_force_orbit_count(&hnfs, pg);
                assert_eq!(
                    filter_count,
                    orbit_count,
                    "filter vs brute-force mismatch at 2D det={det}, PG size={}",
                    pg.len()
                );
            }
        }
    }

    #[test]
    fn hnf_filter_matches_brute_force_3d() {
        let oh = generate_oh();
        let test_pgs: Vec<Vec<Vec<Vec<i64>>>> = vec![
            vec![identity_matrix(3)],
            vec![
                identity_matrix(3),
                vec![vec![-1, 0, 0], vec![0, -1, 0], vec![0, 0, -1]],
            ],
            oh,
        ];
        for pg in &test_pgs {
            for det in 1..=4 {
                let hnfs = crate::hnf::hnf_3d(det);
                let filter_count = inequivalent_hnfs(&hnfs, pg).len();
                let orbit_count = brute_force_orbit_count(&hnfs, pg);
                assert_eq!(
                    filter_count,
                    orbit_count,
                    "filter vs brute-force mismatch at 3D det={det}, PG size={}",
                    pg.len()
                );
            }
        }
    }

    #[test]
    fn hnf_oh_inequivalent_counts_3d() {
        let oh = generate_oh();
        let mut counts = Vec::new();
        for det in 1..=6 {
            let hnfs = crate::hnf::hnf_3d(det);
            counts.push(inequivalent_hnfs(&hnfs, &oh).len());
        }
        // Cross-check first value: det=1 always has exactly 1 HNF
        assert_eq!(counts[0], 1);
        // Verify via brute-force orbit counting
        for det in 1..=6 {
            let hnfs = crate::hnf::hnf_3d(det);
            let orbit_count = brute_force_orbit_count(&hnfs, &oh);
            assert_eq!(
                counts[det - 1],
                orbit_count,
                "Oh filter vs brute-force at det={det}"
            );
        }
    }

    // --- Primitive (non-superperiodic) enumeration ---

    #[test]
    fn primitive_prime_n_equals_all() {
        // n=5 (prime): ℤ_5 has no proper subgroups, so no non-trivial
        // translation can stabilize any coloring with non-constant composition.
        // All G-canonical colorings are therefore primitive.
        let s = vec![vec![5i64]];
        let pg = vec![vec![vec![1i64]]];
        for comp in [&[1usize, 4][..], &[2, 3]] {
            let mut all = 0usize;
            let mut prim = 0usize;
            enumerate_supercell(&s, comp, &pg, &mut |_| all += 1);
            enumerate_supercell_primitive(&s, comp, &pg, &mut |_| prim += 1);
            assert_eq!(
                all, prim,
                "all structures should be primitive for n=5 (prime), comp={comp:?}"
            );
        }
    }

    #[test]
    fn primitive_filters_superperiodic() {
        // ℤ_4 with composition [2,2]: exactly 2 G-canonical colorings.
        // [0,0,1,1] is primitive; [0,1,0,1] is superperiodic (period 2).
        let s = vec![vec![4i64]];
        let pg = vec![vec![vec![1i64]]];
        let mut all = 0usize;
        let mut prim = 0usize;
        enumerate_supercell(&s, &[2, 2], &pg, &mut |_| all += 1);
        enumerate_supercell_primitive(&s, &[2, 2], &pg, &mut |_| prim += 1);
        assert_eq!(all, 2, "ℤ_4 [2,2]: 2 G-canonical colorings total");
        assert_eq!(prim, 1, "ℤ_4 [2,2]: 1 primitive structure");
    }

    #[test]
    fn primitive_subset_of_all() {
        // Primitive count ≤ total count for varied 1D and 2D cases.
        let pg1d = vec![vec![vec![1i64]]];
        let cases_1d: &[(&[usize], Vec<Vec<i64>>)] =
            &[(&[3, 3], vec![vec![6i64]]), (&[4, 4], vec![vec![8i64]])];
        for (comp, s) in cases_1d {
            let mut all = 0usize;
            let mut prim = 0usize;
            enumerate_supercell(s, comp, &pg1d, &mut |_| all += 1);
            enumerate_supercell_primitive(s, comp, &pg1d, &mut |_| prim += 1);
            assert!(prim <= all, "primitive ≤ total for comp={comp:?}");
        }

        let pg2d = vec![identity_matrix(2)];
        let s2d = vec![vec![2i64, 0], vec![0, 2i64]];
        let mut all = 0usize;
        let mut prim = 0usize;
        enumerate_supercell(&s2d, &[2, 2], &pg2d, &mut |_| all += 1);
        enumerate_supercell_primitive(&s2d, &[2, 2], &pg2d, &mut |_| prim += 1);
        assert!(prim <= all, "primitive ≤ total for 2D 2×2 [2,2]");
    }

    #[test]
    fn enumerate_at_index_trivial() {
        let pg = vec![identity_matrix(2)];
        let mut total = 0;
        enumerate_at_index(&pg, 2, &[1, 1], &mut |_| total += 1);
        // det=2 has 3 HNFs in 2D (σ(2)=3), no equivalence reduction with trivial PG
        // Each HNF with [1,1] gives colorings enumerated under T only
        assert_eq!(total, 3);
    }

    #[test]
    fn supercell_agrees_with_full_group_sweep() {
        // Sweep over several 2D supercells and point groups
        let test_cases = vec![
            // diag(2,2) with C2 rotation
            (
                vec![vec![2, 0], vec![0, 2]],
                vec![identity_matrix(2), vec![vec![-1, 0], vec![0, -1]]],
            ),
            // diag(3,3) with C4
            (vec![vec![3, 0], vec![0, 3]], {
                let c4 = vec![vec![0, -1], vec![1, 0]];
                let c2 = vec![vec![-1, 0], vec![0, -1]];
                let c4i = vec![vec![0, 1], vec![-1, 0]];
                vec![identity_matrix(2), c4, c2, c4i]
            }),
            // diag(4,4) with reflection
            (
                vec![vec![4, 0], vec![0, 4]],
                vec![identity_matrix(2), vec![vec![0, 1], vec![1, 0]]],
            ),
        ];

        for (s, pg) in &test_cases {
            let snf = smith_normal_form(s);
            let dims: Vec<usize> = snf.diagonal.iter().map(|&d| d as usize).collect();
            let n: usize = dims.iter().product();

            let translations = product_translations(&dims);
            let point_perms = transform_point_group(pg, &snf);

            let mut full_group = Vec::new();
            for t in &translations {
                for p in &point_perms {
                    full_group.push(p.compose(t));
                }
            }
            full_group.sort_by(|a, b| a.image().cmp(b.image()));
            full_group.dedup();

            let comp = vec![n / 2, n - n / 2];

            let mut supercell_count = 0;
            enumerate_supercell(s, &comp, pg, &mut |_| supercell_count += 1);

            let mut direct_count = 0;
            enumerate(&comp, &full_group, &mut |_| direct_count += 1);

            assert_eq!(
                supercell_count,
                direct_count,
                "mismatch for supercell {s:?} with PG size {}",
                pg.len()
            );
        }
    }

    #[test]
    fn supercell_positions_diagonal_2x2x2() {
        let m = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 2]];
        let positions = supercell_fractional_positions(&m);
        assert_eq!(positions.len(), 8);
        for pos in &positions {
            for &f in pos {
                assert!(
                    (0.0..1.0).contains(&f),
                    "fractional coord out of range: {f}"
                );
            }
        }
        // All positions should be distinct
        let mut sorted: Vec<_> = positions
            .iter()
            .map(|p| {
                (
                    (p[0] * 1000.0) as i64,
                    (p[1] * 1000.0) as i64,
                    (p[2] * 1000.0) as i64,
                )
            })
            .collect();
        sorted.sort();
        sorted.dedup();
        assert_eq!(sorted.len(), 8);
    }

    #[test]
    fn supercell_positions_expected_values() {
        let m = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 2]];
        let positions = supercell_fractional_positions(&m);
        // For diagonal M, each coord should be 0.0 or 0.5
        for pos in &positions {
            for &f in pos {
                assert!(
                    (f - 0.0).abs() < 1e-10 || (f - 0.5).abs() < 1e-10,
                    "expected 0.0 or 0.5, got {f}"
                );
            }
        }
    }

    #[test]
    fn supercell_positions_nondiagonal() {
        let m = vec![vec![2, 1, 0], vec![0, 2, 0], vec![0, 0, 1]];
        let positions = supercell_fractional_positions(&m);
        assert_eq!(positions.len(), 4);

        // Round-trip: r_super · M must have all-integer components
        for pos in &positions {
            #[allow(clippy::needless_range_loop)]
            for i in 0..3 {
                let mut val = 0.0f64;
                for j in 0..3 {
                    val += pos[j] * m[j][i] as f64;
                }
                let rounded = val.round();
                assert!(
                    (val - rounded).abs() < 1e-10,
                    "r_super · M not integer: pos={pos:?}, component {i} = {val}"
                );
            }
        }

        // Expected set (unordered): {(0,0,0), (0.5,0.75,0), (0,0.5,0), (0.5,0.25,0)}
        let mut actual: Vec<(i64, i64, i64)> = positions
            .iter()
            .map(|p| {
                (
                    (p[0] * 1000.0).round() as i64,
                    (p[1] * 1000.0).round() as i64,
                    (p[2] * 1000.0).round() as i64,
                )
            })
            .collect();
        actual.sort();
        let mut expected = vec![(0, 0, 0), (500, 750, 0), (0, 500, 0), (500, 250, 0)];
        expected.sort();
        assert_eq!(actual, expected, "positions don't match expected set");
    }

    // --- Multi-site enumeration ---

    #[test]
    fn multisite_1d_two_atoms_inversion() {
        // 1D chain with 2 atoms at x=0 and x=0.5.
        // Supercell matrix: [[2]]. SNF diagonal: [2]. V = 2.
        // n_total = 2 orbits * 2 = 4 sites.
        //
        // Space group: identity + inversion.
        // Identity: orbit_map = [0, 1], offsets = [[0], [0]]
        // Inversion (x → -x): maps atom 0 (x=0) → atom 0 (x=0), offset [0]
        //                      maps atom 1 (x=0.5) → atom 1 (x=-0.5 ≡ 0.5), offset [-1]
        //
        // With composition [2, 2], Pólya predicts 3 orbits.

        let supercell_matrix = vec![vec![2]];
        let rotations = vec![
            vec![vec![1]],  // identity
            vec![vec![-1]], // inversion
        ];
        let orbit_maps = vec![
            vec![0, 1], // identity
            vec![0, 1], // inversion: both atoms map to themselves
        ];
        let offsets = vec![
            vec![vec![0], vec![0]],  // identity
            vec![vec![0], vec![-1]], // inversion: atom 1 gets offset -1
        ];

        let mut count = 0;
        enumerate_supercell_multisite(
            &supercell_matrix,
            &[2, 2],
            &rotations,
            &orbit_maps,
            &offsets,
            2,
            &mut |_| count += 1,
        );
        assert_eq!(
            count, 3,
            "1D 2-atom inversion chain should give 3 colorings"
        );
    }

    #[test]
    fn multisite_identity_only_matches_single_site() {
        // With only the identity operation and 1 orbit, multi-site
        // enumeration should match single-site enumeration exactly.
        let supercell_matrix = vec![vec![2, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
        let rotations = vec![vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]]];
        let orbit_maps = vec![vec![0]];
        let offsets = vec![vec![vec![0, 0, 0]]];

        let mut multi_count = 0;
        enumerate_supercell_multisite(
            &supercell_matrix,
            &[1, 1],
            &rotations,
            &orbit_maps,
            &offsets,
            1,
            &mut |_| multi_count += 1,
        );

        let pg = vec![vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]]];
        let mut single_count = 0;
        enumerate_supercell(&supercell_matrix, &[1, 1], &pg, &mut |_| {
            single_count += 1;
        });

        assert_eq!(
            multi_count, single_count,
            "single-orbit multisite should match single-site"
        );
    }

    #[test]
    fn multisite_fractional_positions_two_atoms() {
        // 1D supercell [[2]] with parent atoms at 0.0 and 0.5.
        // Sites: (orbit 0, y=0), (orbit 0, y=1), (orbit 1, y=0), (orbit 1, y=1)
        // Supercell fracs: 0.0, 0.5, 0.25, 0.75  (after M⁻¹ transform)
        let supercell_matrix = vec![vec![2, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
        let parent_positions = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let positions =
            supercell_fractional_positions_multisite(&supercell_matrix, &parent_positions);
        assert_eq!(positions.len(), 4);

        let mut xs: Vec<i64> = positions
            .iter()
            .map(|p| (p[0] * 1000.0).round() as i64)
            .collect();
        xs.sort();
        assert_eq!(xs, vec![0, 250, 500, 750]);
    }

    #[test]
    fn multisite_constrained_disjoint_species() {
        // 1D, 2 atoms, supercell ×2, inversion symmetry.
        // 4 total sites: orbit 0 (sites 0,1), orbit 1 (sites 2,3).
        // 3 species. Orbit 0 allows {0,1}, orbit 1 allows {0,2}.
        // Composition [2, 1, 1]: 2×species0, 1×species1, 1×species2.
        // species1 must go to orbit 0, species2 must go to orbit 1.
        // Remaining: 1×species0 in orbit 0, 1×species0 in orbit 1.
        // Under inversion (swap within each orbit): 1 way.
        let supercell_matrix = vec![vec![2]];
        let rotations = vec![vec![vec![1]], vec![vec![-1]]];
        let orbit_maps = vec![vec![0, 1], vec![0, 1]];
        let offsets = vec![vec![vec![0], vec![0]], vec![vec![0], vec![0]]];
        let orbit_allowed = vec![vec![true, true, false], vec![true, false, true]];
        let mut results = Vec::new();
        enumerate_supercell_multisite_constrained(
            &supercell_matrix,
            &[2, 1, 1],
            &rotations,
            &orbit_maps,
            &offsets,
            2,
            &orbit_allowed,
            &mut |c| results.push(c.to_vec()),
        );
        // 4 valid colorings, 2 orbits under translation (0↔1)(2↔3):
        //   {[0,1,0,2], [1,0,2,0]} and {[0,1,2,0], [1,0,0,2]}
        assert_eq!(results.len(), 2);
    }
}
