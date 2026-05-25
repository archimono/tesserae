use moyo::MoyoDataset;
use moyo::base::{Cell, Lattice};
use nalgebra::{Matrix3, Vector3};

/// A symmetry operation: integer rotation matrix + fractional translation.
#[derive(Debug, Clone)]
pub struct SymmetryOperation {
    pub rotation: [[i32; 3]; 3],
    pub translation: [f64; 3],
}

/// Get all symmetry operations of a crystal structure using moyo.
///
/// `lattice` is a 3×3 matrix where each row is a lattice vector (a, b, c)
/// in Cartesian coordinates. `positions` are fractional coordinates.
/// `types` maps each atom to a species index.
///
/// Returns rotation matrices in the fractional coordinate basis (integer)
/// and fractional translation vectors.
///
/// # Panics
///
/// Panics if moyo fails to determine symmetry.
pub fn get_symmetry(
    lattice: &[[f64; 3]; 3],
    positions: &[[f64; 3]],
    types: &[i32],
    symprec: f64,
) -> Vec<SymmetryOperation> {
    assert_eq!(positions.len(), types.len());

    let lat = Lattice::new(Matrix3::new(
        lattice[0][0],
        lattice[0][1],
        lattice[0][2],
        lattice[1][0],
        lattice[1][1],
        lattice[1][2],
        lattice[2][0],
        lattice[2][1],
        lattice[2][2],
    ));

    let pos: Vec<Vector3<f64>> = positions
        .iter()
        .map(|p| Vector3::new(p[0], p[1], p[2]))
        .collect();

    let cell = Cell::new(lat, pos, types.to_vec());

    let dataset =
        MoyoDataset::with_default(&cell, symprec).expect("moyo failed to determine symmetry");

    dataset
        .operations
        .iter()
        .map(|op| {
            let r = op.rotation;
            SymmetryOperation {
                rotation: [
                    [r[(0, 0)], r[(0, 1)], r[(0, 2)]],
                    [r[(1, 0)], r[(1, 1)], r[(1, 2)]],
                    [r[(2, 0)], r[(2, 1)], r[(2, 2)]],
                ],
                translation: [op.translation[0], op.translation[1], op.translation[2]],
            }
        })
        .collect()
}

/// Extract unique rotation matrices from symmetry operations (single-site enumeration only).
///
/// Returns the crystallographic point group P by discarding translations and
/// deduplicating. Valid for Bravais lattice enumeration: every rotation R in
/// {R | τ} is a lattice symmetry regardless of τ. For multi-atom parents, use
/// [`extract_rotations`] + [`compute_orbit_maps`] with the multisite path.
pub fn extract_point_group(ops: &[SymmetryOperation]) -> Vec<Vec<Vec<i64>>> {
    let mut pg: Vec<Vec<Vec<i64>>> = Vec::new();
    for op in ops {
        let mat: Vec<Vec<i64>> = op
            .rotation
            .iter()
            .map(|row| row.iter().map(|&v| v as i64).collect())
            .collect();
        if !pg.contains(&mat) {
            pg.push(mat);
        }
    }
    pg
}

/// Extract rotation matrices as `Vec<Vec<Vec<i64>>>` from symmetry operations.
///
/// Unlike `extract_point_group`, this keeps duplicates: the returned list
/// is indexed in lockstep with `orbit_maps` and `offsets` from
/// `compute_orbit_maps`.
pub fn extract_rotations(ops: &[SymmetryOperation]) -> Vec<Vec<Vec<i64>>> {
    ops.iter()
        .map(|op| {
            op.rotation
                .iter()
                .map(|row| row.iter().map(|&v| v as i64).collect())
                .collect()
        })
        .collect()
}

/// Compute orbit permutations and integer offsets for multi-site enumeration.
///
/// For each space group operation `(R, t)`, determines which parent atom `k`
/// atom `j` maps to: `R·x_j + t ≈ x_k + n_jk` (modulo lattice), where
/// `n_jk` is an integer vector.
///
/// Returns `(orbit_maps, offsets)`:
/// - `orbit_maps[op][j]` = index `k` of the target atom
/// - `offsets[op][j]` = integer vector `n_jk` (parent lattice basis)
///
/// # Arguments
///
/// * `ops` — Space group operations from `get_symmetry`.
/// * `positions` — Fractional positions of parent atoms.
/// * `symprec` — Tolerance for matching positions modulo lattice.
pub fn compute_orbit_maps(
    ops: &[SymmetryOperation],
    positions: &[[f64; 3]],
    symprec: f64,
) -> (Vec<Vec<usize>>, Vec<Vec<Vec<i64>>>) {
    let n_atoms = positions.len();
    let mut orbit_maps = Vec::with_capacity(ops.len());
    let mut offsets = Vec::with_capacity(ops.len());

    for op in ops {
        let mut op_map = Vec::with_capacity(n_atoms);
        let mut op_offsets = Vec::with_capacity(n_atoms);

        for (j, pos_j) in positions.iter().enumerate() {
            let mapped: [f64; 3] = std::array::from_fn(|i| {
                op.translation[i]
                    + (0..3)
                        .map(|k| op.rotation[i][k] as f64 * pos_j[k])
                        .sum::<f64>()
            });

            let mut found = false;
            for (k, pos_k) in positions.iter().enumerate() {
                let diff: [f64; 3] = std::array::from_fn(|i| mapped[i] - pos_k[i]);

                let n_jk: Vec<i64> = diff.iter().map(|&d| d.round() as i64).collect();
                let is_match = (0..3).all(|i| (diff[i] - n_jk[i] as f64).abs() < symprec);

                if is_match {
                    op_map.push(k);
                    op_offsets.push(n_jk);
                    found = true;
                    break;
                }
            }
            assert!(
                found,
                "space group op did not map atom {j} to any known atom"
            );
        }

        orbit_maps.push(op_map);
        offsets.push(op_offsets);
    }

    (orbit_maps, offsets)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cubic_nacl_has_48_point_ops() {
        let a = 5.64;
        let lattice = [
            [0.0, a / 2.0, a / 2.0],
            [a / 2.0, 0.0, a / 2.0],
            [a / 2.0, a / 2.0, 0.0],
        ];
        let positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];
        let types = [1, 2];

        let ops = get_symmetry(&lattice, &positions, &types, 1e-5);
        assert_eq!(ops.len(), 48, "NaCl Fm-3m should have 48 operations");

        let pg = extract_point_group(&ops);
        assert_eq!(pg.len(), 48, "NaCl should have 48 distinct point group ops");
    }

    #[test]
    fn identity_always_present() {
        let a = 4.0;
        let lattice = [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]];
        let positions = [[0.0, 0.0, 0.0]];
        let types = [1];

        let ops = get_symmetry(&lattice, &positions, &types, 1e-5);
        let has_identity = ops
            .iter()
            .any(|op| op.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
        assert!(has_identity, "identity must be present");
    }

    #[test]
    fn hexagonal_2d() {
        let a = 3.0;
        let lattice = [
            [a, 0.0, 0.0],
            [-a / 2.0, a * 3.0_f64.sqrt() / 2.0, 0.0],
            [0.0, 0.0, 20.0],
        ];
        let positions = [[0.0, 0.0, 0.0]];
        let types = [1];

        let ops = get_symmetry(&lattice, &positions, &types, 1e-5);
        assert_eq!(ops.len(), 24);
    }

    #[test]
    fn diamond_nonsymmorphic_fd3m() {
        // Diamond cubic (Fd-3m, #227) — non-symmorphic space group.
        // FCC lattice with two atoms at (0,0,0) and (1/4,1/4,1/4).
        let a = 5.43; // Si lattice constant
        let lattice = [
            [0.0, a / 2.0, a / 2.0],
            [a / 2.0, 0.0, a / 2.0],
            [a / 2.0, a / 2.0, 0.0],
        ];
        let positions = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]];
        let types = [1, 1];

        let ops = get_symmetry(&lattice, &positions, &types, 1e-5);
        // Primitive FCC cell: 192 conventional ops / 4 centering = 48
        assert_eq!(
            ops.len(),
            48,
            "diamond primitive cell should have 48 operations"
        );

        // Non-symmorphic: some operations have non-zero fractional translations
        let has_nontrivial_translation = ops.iter().any(|op| {
            op.translation.iter().any(|&t| {
                let t_mod = (t % 1.0 + 1.0) % 1.0;
                t_mod > 1e-5 && t_mod < 1.0 - 1e-5
            })
        });
        assert!(
            has_nontrivial_translation,
            "Fd-3m must have non-trivial translations"
        );

        // extract_point_group gives unique rotations = crystallographic point group Oh
        let pg = extract_point_group(&ops);
        assert_eq!(
            pg.len(),
            48,
            "diamond point group should be Oh with 48 elements"
        );
    }
}
