use std::fmt::Write;

/// Write a POSCAR-format string for a supercell coloring.
///
/// # Arguments
///
/// * `parent_lattice` — 3×3 parent lattice vectors (rows = vectors, Å).
/// * `supercell_matrix` — 3×3 integer supercell matrix M (L_super = M * L_parent).
/// * `species_names` — Species name for each color index in the coloring.
/// * `coloring` — Species label per site (0-indexed).
/// * `comment` — POSCAR comment line.
///
/// The output lattice is M * L_parent. Fractional coordinates are in the
/// supercell basis, sorted by species.
pub fn write_poscar(
    parent_lattice: &[[f64; 3]; 3],
    supercell_matrix: &[Vec<i64>],
    species_names: &[String],
    coloring: &[u32],
    comment: &str,
) -> String {
    let n = coloring.len();
    assert!(
        supercell_matrix.len() == 3 && supercell_matrix.iter().all(|r| r.len() == 3),
        "supercell_matrix must be 3x3"
    );

    // Supercell lattice: L_super[i] = sum_j M[i][j] * L_parent[j]
    let mut super_lattice = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            for j in 0..3 {
                super_lattice[i][k] += supercell_matrix[i][j] as f64 * parent_lattice[j][k];
            }
        }
    }

    let positions = tesserae_core::supercell::supercell_fractional_positions(supercell_matrix);
    assert_eq!(positions.len(), n);

    let k = species_names.len();
    assert!(
        coloring.iter().all(|&l| (l as usize) < k),
        "coloring contains label >= species_names.len()"
    );
    let mut by_species: Vec<Vec<usize>> = vec![Vec::new(); k];
    for (i, &label) in coloring.iter().enumerate() {
        by_species[label as usize].push(i);
    }

    let mut out = String::new();
    writeln!(out, "{comment}").unwrap();
    writeln!(out, "1.0").unwrap();
    for row in &super_lattice {
        writeln!(out, "  {:.10}  {:.10}  {:.10}", row[0], row[1], row[2]).unwrap();
    }
    let names: Vec<&str> = species_names.iter().map(|s| s.as_str()).collect();
    writeln!(out, "{}", names.join("  ")).unwrap();
    let counts: Vec<String> = by_species.iter().map(|v| v.len().to_string()).collect();
    writeln!(out, "{}", counts.join("  ")).unwrap();
    writeln!(out, "Direct").unwrap();
    for group in &by_species {
        for &idx in group {
            let p = &positions[idx];
            writeln!(out, "  {:.10}  {:.10}  {:.10}", p[0], p[1], p[2]).unwrap();
        }
    }

    out
}

/// Write a CIF-format string for a supercell coloring.
///
/// Produces a minimal CIF block with lattice parameters and atom sites.
pub fn write_cif(
    parent_lattice: &[[f64; 3]; 3],
    supercell_matrix: &[Vec<i64>],
    species_names: &[String],
    coloring: &[u32],
    data_name: &str,
) -> String {
    let n = coloring.len();
    assert!(
        supercell_matrix.len() == 3 && supercell_matrix.iter().all(|r| r.len() == 3),
        "supercell_matrix must be 3x3"
    );

    // Supercell lattice
    let mut super_lattice = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            for j in 0..3 {
                super_lattice[i][k] += supercell_matrix[i][j] as f64 * parent_lattice[j][k];
            }
        }
    }

    let (a, b, c, alpha, beta, gamma) = lattice_parameters(&super_lattice);
    let positions = tesserae_core::supercell::supercell_fractional_positions(supercell_matrix);
    assert_eq!(positions.len(), n);

    let mut out = String::new();
    writeln!(out, "data_{data_name}").unwrap();
    writeln!(out, "_cell_length_a   {a:.6}").unwrap();
    writeln!(out, "_cell_length_b   {b:.6}").unwrap();
    writeln!(out, "_cell_length_c   {c:.6}").unwrap();
    writeln!(out, "_cell_angle_alpha   {alpha:.4}").unwrap();
    writeln!(out, "_cell_angle_beta    {beta:.4}").unwrap();
    writeln!(out, "_cell_angle_gamma   {gamma:.4}").unwrap();
    writeln!(out, "loop_").unwrap();
    writeln!(out, "_atom_site_label").unwrap();
    writeln!(out, "_atom_site_type_symbol").unwrap();
    writeln!(out, "_atom_site_fract_x").unwrap();
    writeln!(out, "_atom_site_fract_y").unwrap();
    writeln!(out, "_atom_site_fract_z").unwrap();

    assert!(
        coloring.iter().all(|&l| (l as usize) < species_names.len()),
        "coloring contains label >= species_names.len()"
    );
    let mut label_counts: Vec<usize> = vec![0; species_names.len()];
    for (i, &label) in coloring.iter().enumerate() {
        let sp = label as usize;
        label_counts[sp] += 1;
        let name = &species_names[sp];
        let p = &positions[i];
        writeln!(
            out,
            "{}{} {} {:.10} {:.10} {:.10}",
            name, label_counts[sp], name, p[0], p[1], p[2]
        )
        .unwrap();
    }

    out
}

fn lattice_parameters(lattice: &[[f64; 3]; 3]) -> (f64, f64, f64, f64, f64, f64) {
    let norm = |v: &[f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let dot = |a: &[f64; 3], b: &[f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let a = norm(&lattice[0]);
    let b = norm(&lattice[1]);
    let c = norm(&lattice[2]);
    let alpha = (dot(&lattice[1], &lattice[2]) / (b * c))
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees();
    let beta = (dot(&lattice[0], &lattice[2]) / (a * c))
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees();
    let gamma = (dot(&lattice[0], &lattice[1]) / (a * b))
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees();
    (a, b, c, alpha, beta, gamma)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poscar_roundtrip() {
        let lattice = [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]];
        let supercell = vec![vec![2, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
        let species = vec!["A".to_string(), "B".to_string()];
        let coloring = vec![0, 1];
        let poscar = write_poscar(&lattice, &supercell, &species, &coloring, "test");
        assert!(poscar.contains("A  B"));
        assert!(poscar.contains("1  1"));
        assert!(poscar.contains("Direct"));
        let lines: Vec<&str> = poscar.lines().collect();
        assert_eq!(lines[0], "test");
    }

    #[test]
    fn cif_output() {
        let lattice = [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]];
        let supercell = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 2]];
        let species = vec!["Na".to_string(), "Cl".to_string()];
        let coloring = vec![0, 0, 0, 0, 1, 1, 1, 1];
        let cif = write_cif(&lattice, &supercell, &species, &coloring, "test");
        assert!(cif.contains("data_test"));
        assert!(cif.contains("_cell_length_a"));
        assert!(cif.contains("_atom_site_fract_x"));
        assert!(cif.contains("Na"));
        assert!(cif.contains("Cl"));
    }

    #[test]
    fn poscar_species_grouped() {
        let lattice = [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]];
        let supercell = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 2]];
        let species = vec!["A".to_string(), "B".to_string()];
        let coloring = vec![0, 1, 0, 1, 0, 1, 0, 1];
        let poscar = write_poscar(&lattice, &supercell, &species, &coloring, "test");
        assert!(poscar.contains("4  4"));
    }
}
