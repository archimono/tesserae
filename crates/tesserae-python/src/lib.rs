use pyo3::prelude::*;

#[pyfunction]
#[pyo3(signature = (supercell_matrix, composition, parent_point_group))]
fn enumerate_supercell(
    supercell_matrix: Vec<Vec<i64>>,
    composition: Vec<usize>,
    parent_point_group: Vec<Vec<Vec<i64>>>,
) -> Vec<Vec<u32>> {
    let mut results = Vec::new();
    tesserae_core::supercell::enumerate_supercell(
        &supercell_matrix,
        &composition,
        &parent_point_group,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (parent_point_group, index, composition))]
fn enumerate_at_index(
    parent_point_group: Vec<Vec<Vec<i64>>>,
    index: usize,
    composition: Vec<usize>,
) -> Vec<Vec<u32>> {
    let mut results = Vec::new();
    tesserae_core::supercell::enumerate_at_index(
        &parent_point_group,
        index,
        &composition,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (hnfs, parent_point_group))]
fn inequivalent_hnfs(
    hnfs: Vec<Vec<Vec<i64>>>,
    parent_point_group: Vec<Vec<Vec<i64>>>,
) -> Vec<Vec<Vec<i64>>> {
    tesserae_core::supercell::inequivalent_hnfs(&hnfs, &parent_point_group)
}

#[pyfunction]
#[pyo3(signature = (poscar_path, supercell_matrix, composition, symprec=1e-5))]
fn enumerate_from_poscar(
    poscar_path: &str,
    supercell_matrix: Vec<Vec<i64>>,
    composition: Vec<usize>,
    symprec: f64,
) -> Vec<Vec<u32>> {
    let path = std::path::Path::new(poscar_path);
    let mut results = Vec::new();
    tesserae_io::enumerate_from_poscar(path, &supercell_matrix, &composition, symprec, &mut |c| {
        results.push(c.to_vec())
    });
    results
}

#[pyfunction]
#[pyo3(signature = (cif_path, supercell_matrix, composition, symprec=1e-5))]
fn enumerate_from_cif(
    cif_path: &str,
    supercell_matrix: Vec<Vec<i64>>,
    composition: Vec<usize>,
    symprec: f64,
) -> Vec<Vec<u32>> {
    let path = std::path::Path::new(cif_path);
    let mut results = Vec::new();
    tesserae_io::enumerate_from_cif(path, &supercell_matrix, &composition, symprec, &mut |c| {
        results.push(c.to_vec())
    });
    results
}

#[pyfunction]
#[pyo3(signature = (supercell_matrix, composition, parent_point_group))]
fn enumerate_supercell_primitive(
    supercell_matrix: Vec<Vec<i64>>,
    composition: Vec<usize>,
    parent_point_group: Vec<Vec<Vec<i64>>>,
) -> Vec<Vec<u32>> {
    let mut results = Vec::new();
    tesserae_core::supercell::enumerate_supercell_primitive(
        &supercell_matrix,
        &composition,
        &parent_point_group,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (parent_point_group, index, composition))]
fn enumerate_at_index_primitive(
    parent_point_group: Vec<Vec<Vec<i64>>>,
    index: usize,
    composition: Vec<usize>,
) -> Vec<Vec<u32>> {
    let mut results = Vec::new();
    tesserae_core::supercell::enumerate_at_index_primitive(
        &parent_point_group,
        index,
        &composition,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (poscar_path, supercell_matrix, composition, symprec=1e-5))]
fn enumerate_from_poscar_primitive(
    poscar_path: &str,
    supercell_matrix: Vec<Vec<i64>>,
    composition: Vec<usize>,
    symprec: f64,
) -> Vec<Vec<u32>> {
    let path = std::path::Path::new(poscar_path);
    let mut results = Vec::new();
    tesserae_io::enumerate_from_poscar_primitive(
        path,
        &supercell_matrix,
        &composition,
        symprec,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (cif_path, supercell_matrix, composition, symprec=1e-5))]
fn enumerate_from_cif_primitive(
    cif_path: &str,
    supercell_matrix: Vec<Vec<i64>>,
    composition: Vec<usize>,
    symprec: f64,
) -> Vec<Vec<u32>> {
    let path = std::path::Path::new(cif_path);
    let mut results = Vec::new();
    tesserae_io::enumerate_from_cif_primitive(
        path,
        &supercell_matrix,
        &composition,
        symprec,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (cif_path, index, composition, symprec=1e-5))]
fn enumerate_from_cif_at_index(
    cif_path: &str,
    index: usize,
    composition: Vec<usize>,
    symprec: f64,
) -> Vec<Vec<u32>> {
    let path = std::path::Path::new(cif_path);
    let mut results = Vec::new();
    tesserae_io::enumerate_from_cif_at_index(path, index, &composition, symprec, &mut |c| {
        results.push(c.to_vec())
    });
    results
}

#[pyfunction]
#[pyo3(signature = (cif_path, index, composition, symprec=1e-5))]
fn enumerate_from_cif_at_index_primitive(
    cif_path: &str,
    index: usize,
    composition: Vec<usize>,
    symprec: f64,
) -> Vec<Vec<u32>> {
    let path = std::path::Path::new(cif_path);
    let mut results = Vec::new();
    tesserae_io::enumerate_from_cif_at_index_primitive(
        path,
        index,
        &composition,
        symprec,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (poscar_path, index, composition, symprec=1e-5))]
fn enumerate_from_poscar_at_index(
    poscar_path: &str,
    index: usize,
    composition: Vec<usize>,
    symprec: f64,
) -> Vec<Vec<u32>> {
    let path = std::path::Path::new(poscar_path);
    let mut results = Vec::new();
    tesserae_io::enumerate_from_poscar_at_index(path, index, &composition, symprec, &mut |c| {
        results.push(c.to_vec())
    });
    results
}

#[pyfunction]
#[pyo3(signature = (poscar_path, index, composition, symprec=1e-5))]
fn enumerate_from_poscar_at_index_primitive(
    poscar_path: &str,
    index: usize,
    composition: Vec<usize>,
    symprec: f64,
) -> Vec<Vec<u32>> {
    let path = std::path::Path::new(poscar_path);
    let mut results = Vec::new();
    tesserae_io::enumerate_from_poscar_at_index_primitive(
        path,
        index,
        &composition,
        symprec,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (supercell_matrix, composition, parent_point_group))]
fn enumerate_supercell_with_multiplicity(
    supercell_matrix: Vec<Vec<i64>>,
    composition: Vec<usize>,
    parent_point_group: Vec<Vec<Vec<i64>>>,
) -> Vec<(Vec<u32>, u128)> {
    let mut results = Vec::new();
    tesserae_core::supercell::enumerate_supercell_with_multiplicity(
        &supercell_matrix,
        &composition,
        &parent_point_group,
        &mut |c, w| results.push((c.to_vec(), w)),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (parent_point_group, index, composition))]
fn enumerate_at_index_with_multiplicity(
    parent_point_group: Vec<Vec<Vec<i64>>>,
    index: usize,
    composition: Vec<usize>,
) -> Vec<(Vec<u32>, u128)> {
    let mut results = Vec::new();
    tesserae_core::supercell::enumerate_at_index_with_multiplicity(
        &parent_point_group,
        index,
        &composition,
        &mut |c, w| results.push((c.to_vec(), w)),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (parent_point_group, index, min_conc, max_conc))]
fn enumerate_at_index_range(
    parent_point_group: Vec<Vec<Vec<i64>>>,
    index: usize,
    min_conc: Vec<usize>,
    max_conc: Vec<usize>,
) -> Vec<Vec<u32>> {
    let mut results = Vec::new();
    tesserae_core::supercell::enumerate_at_index_range(
        &parent_point_group,
        index,
        &min_conc,
        &max_conc,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (parent_point_group, index, min_conc, max_conc))]
fn enumerate_at_index_range_primitive(
    parent_point_group: Vec<Vec<Vec<i64>>>,
    index: usize,
    min_conc: Vec<usize>,
    max_conc: Vec<usize>,
) -> Vec<Vec<u32>> {
    let mut results = Vec::new();
    tesserae_core::supercell::enumerate_at_index_range_primitive(
        &parent_point_group,
        index,
        &min_conc,
        &max_conc,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (supercell_matrix, composition, rotations, orbit_maps, offsets, n_orbits))]
fn enumerate_supercell_multisite(
    supercell_matrix: Vec<Vec<i64>>,
    composition: Vec<usize>,
    rotations: Vec<Vec<Vec<i64>>>,
    orbit_maps: Vec<Vec<usize>>,
    offsets: Vec<Vec<Vec<i64>>>,
    n_orbits: usize,
) -> Vec<Vec<u32>> {
    let mut results = Vec::new();
    tesserae_core::supercell::enumerate_supercell_multisite(
        &supercell_matrix,
        &composition,
        &rotations,
        &orbit_maps,
        &offsets,
        n_orbits,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (poscar_path, supercell_matrix, composition, symprec=1e-5))]
fn enumerate_from_poscar_multisite(
    poscar_path: &str,
    supercell_matrix: Vec<Vec<i64>>,
    composition: Vec<usize>,
    symprec: f64,
) -> Vec<Vec<u32>> {
    let path = std::path::Path::new(poscar_path);
    let mut results = Vec::new();
    tesserae_io::enumerate_from_poscar_multisite(
        path,
        &supercell_matrix,
        &composition,
        symprec,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (cif_path, supercell_matrix, composition, symprec=1e-5))]
fn enumerate_from_cif_multisite(
    cif_path: &str,
    supercell_matrix: Vec<Vec<i64>>,
    composition: Vec<usize>,
    symprec: f64,
) -> Vec<Vec<u32>> {
    let path = std::path::Path::new(cif_path);
    let mut results = Vec::new();
    tesserae_io::enumerate_from_cif_multisite(
        path,
        &supercell_matrix,
        &composition,
        symprec,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (supercell_matrix, composition, rotations, orbit_maps, offsets, n_orbits, orbit_allowed))]
fn enumerate_supercell_multisite_constrained(
    supercell_matrix: Vec<Vec<i64>>,
    composition: Vec<usize>,
    rotations: Vec<Vec<Vec<i64>>>,
    orbit_maps: Vec<Vec<usize>>,
    offsets: Vec<Vec<Vec<i64>>>,
    n_orbits: usize,
    orbit_allowed: Vec<Vec<bool>>,
) -> Vec<Vec<u32>> {
    let mut results = Vec::new();
    tesserae_core::supercell::enumerate_supercell_multisite_constrained(
        &supercell_matrix,
        &composition,
        &rotations,
        &orbit_maps,
        &offsets,
        n_orbits,
        &orbit_allowed,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (cif_path, supercell_matrix, composition, symprec=1e-5))]
fn enumerate_from_cif_constrained(
    cif_path: &str,
    supercell_matrix: Vec<Vec<i64>>,
    composition: Vec<usize>,
    symprec: f64,
) -> Vec<Vec<u32>> {
    let path = std::path::Path::new(cif_path);
    let mut results = Vec::new();
    tesserae_io::enumerate_from_cif_constrained(
        path,
        &supercell_matrix,
        &composition,
        symprec,
        &mut |c| results.push(c.to_vec()),
    );
    results
}

#[pyfunction]
#[pyo3(signature = (cif_path, symprec=1e-5))]
#[allow(clippy::type_complexity)]
fn cif_partial_occupancy_info(
    cif_path: &str,
    symprec: f64,
) -> ([[f64; 3]; 3], Vec<Vec<f64>>, Vec<String>, Vec<Vec<bool>>) {
    let path = std::path::Path::new(cif_path);
    let (lattice, positions, species, orbit_allowed) =
        tesserae_io::cif_partial_occupancy_info(path, symprec);
    let pos_vecs: Vec<Vec<f64>> = positions.iter().map(|p| p.to_vec()).collect();
    (lattice, pos_vecs, species, orbit_allowed)
}

#[pyfunction]
#[pyo3(signature = (supercell_matrix))]
fn supercell_fractional_positions(supercell_matrix: Vec<Vec<i64>>) -> Vec<Vec<f64>> {
    tesserae_core::supercell::supercell_fractional_positions(&supercell_matrix)
}

#[pyfunction]
#[pyo3(signature = (supercell_matrix, parent_positions))]
fn supercell_fractional_positions_multisite(
    supercell_matrix: Vec<Vec<i64>>,
    parent_positions: Vec<[f64; 3]>,
) -> Vec<Vec<f64>> {
    tesserae_core::supercell::supercell_fractional_positions_multisite(
        &supercell_matrix,
        &parent_positions,
    )
}

#[pyfunction]
#[pyo3(signature = (generators, dim))]
fn generate_matrix_group(generators: Vec<Vec<Vec<i64>>>, dim: usize) -> Vec<Vec<Vec<i64>>> {
    tesserae_core::supercell::generate_matrix_group(&generators, dim)
}

#[pyfunction]
#[pyo3(signature = (parent_lattice, supercell_matrix, species_names, coloring, comment))]
fn write_poscar(
    parent_lattice: [[f64; 3]; 3],
    supercell_matrix: Vec<Vec<i64>>,
    species_names: Vec<String>,
    coloring: Vec<u32>,
    comment: &str,
) -> String {
    tesserae_io::write::write_poscar(
        &parent_lattice,
        &supercell_matrix,
        &species_names,
        &coloring,
        comment,
    )
}

#[pyfunction]
#[pyo3(signature = (parent_lattice, supercell_matrix, species_names, coloring, data_name))]
fn write_cif(
    parent_lattice: [[f64; 3]; 3],
    supercell_matrix: Vec<Vec<i64>>,
    species_names: Vec<String>,
    coloring: Vec<u32>,
    data_name: &str,
) -> String {
    tesserae_io::write::write_cif(
        &parent_lattice,
        &supercell_matrix,
        &species_names,
        &coloring,
        data_name,
    )
}

#[pymodule]
fn tesserae(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(enumerate_supercell, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_at_index, m)?)?;
    m.add_function(wrap_pyfunction!(inequivalent_hnfs, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_from_poscar, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_from_cif, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_supercell_primitive, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_at_index_primitive, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_from_poscar_primitive, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_from_cif_primitive, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_from_cif_at_index, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_from_cif_at_index_primitive, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_from_poscar_at_index, m)?)?;
    m.add_function(wrap_pyfunction!(
        enumerate_from_poscar_at_index_primitive,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(enumerate_supercell_with_multiplicity, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_at_index_with_multiplicity, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_at_index_range, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_at_index_range_primitive, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_supercell_multisite, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_from_poscar_multisite, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_from_cif_multisite, m)?)?;
    m.add_function(wrap_pyfunction!(
        enumerate_supercell_multisite_constrained,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(enumerate_from_cif_constrained, m)?)?;
    m.add_function(wrap_pyfunction!(cif_partial_occupancy_info, m)?)?;
    m.add_function(wrap_pyfunction!(supercell_fractional_positions, m)?)?;
    m.add_function(wrap_pyfunction!(
        supercell_fractional_positions_multisite,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(generate_matrix_group, m)?)?;
    m.add_function(wrap_pyfunction!(write_poscar, m)?)?;
    m.add_function(wrap_pyfunction!(write_cif, m)?)?;
    Ok(())
}
