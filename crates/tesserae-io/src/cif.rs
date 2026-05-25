use std::path::Path;

/// A crystal structure parsed from a CIF file.
///
/// Handles the subset of CIF tags needed for enumeration:
/// cell parameters, atom sites (fractional coordinates), and species labels.
#[derive(Debug, Clone)]
pub struct Cif {
    pub lattice: [[f64; 3]; 3],
    pub species_names: Vec<String>,
    pub species_counts: Vec<usize>,
    pub positions: Vec<[f64; 3]>,
    pub occupancies: Vec<f64>,
    /// Per-atom species label in file order (parallel to `positions`).
    pub labels: Vec<String>,
}

impl Cif {
    pub fn parse(s: &str) -> Self {
        let mut a = 0.0_f64;
        let mut b = 0.0_f64;
        let mut c = 0.0_f64;
        let mut alpha = 90.0_f64;
        let mut beta = 90.0_f64;
        let mut gamma = 90.0_f64;

        let lines: Vec<&str> = s.lines().collect();
        let mut i = 0;

        let mut frac_x: Vec<f64> = Vec::new();
        let mut frac_y: Vec<f64> = Vec::new();
        let mut frac_z: Vec<f64> = Vec::new();
        let mut labels: Vec<String> = Vec::new();
        let mut occupancies: Vec<f64> = Vec::new();

        while i < lines.len() {
            let line = lines[i].trim();

            if line.starts_with('#') || line.is_empty() {
                i += 1;
                continue;
            }

            if line.starts_with("_cell_length_a") {
                a = parse_cif_float(tag_value(line));
            } else if line.starts_with("_cell_length_b") {
                b = parse_cif_float(tag_value(line));
            } else if line.starts_with("_cell_length_c") {
                c = parse_cif_float(tag_value(line));
            } else if line.starts_with("_cell_angle_alpha") {
                alpha = parse_cif_float(tag_value(line));
            } else if line.starts_with("_cell_angle_beta") {
                beta = parse_cif_float(tag_value(line));
            } else if line.starts_with("_cell_angle_gamma") {
                gamma = parse_cif_float(tag_value(line));
            } else if line == "loop_" {
                i += 1;
                let (cols, data_start) = parse_loop_header(&lines, i);
                if let Some(atom_loop) = AtomSiteLoop::from_columns(&cols) {
                    for line in &lines[data_start..] {
                        let row = line.trim();
                        if row.is_empty()
                            || row.starts_with('_')
                            || row.starts_with("loop_")
                            || row.starts_with("data_")
                        {
                            break;
                        }
                        let tokens: Vec<&str> = row.split_whitespace().collect();
                        if tokens.len() < cols.len() {
                            break;
                        }
                        let label = atom_loop
                            .type_symbol_idx
                            .map(|idx| tokens[idx])
                            .or(atom_loop.label_idx.map(|idx| tokens[idx]))
                            .expect("no atom label or type_symbol in CIF");
                        let species = strip_oxidation(label);
                        labels.push(species);
                        frac_x.push(parse_cif_float(tokens[atom_loop.fract_x_idx]));
                        frac_y.push(parse_cif_float(tokens[atom_loop.fract_y_idx]));
                        frac_z.push(parse_cif_float(tokens[atom_loop.fract_z_idx]));
                        let occ = atom_loop
                            .occupancy_idx
                            .map(|idx| parse_cif_float(tokens[idx]))
                            .unwrap_or(1.0);
                        occupancies.push(occ);
                    }
                }
                i = data_start;
                continue;
            }

            i += 1;
        }

        assert!(!frac_x.is_empty(), "no atom sites found in CIF");

        let lattice = lattice_from_params(a, b, c, alpha, beta, gamma);

        let mut positions = Vec::with_capacity(frac_x.len());
        for j in 0..frac_x.len() {
            positions.push([frac_x[j], frac_y[j], frac_z[j]]);
        }

        let mut species_order: Vec<String> = Vec::new();
        let mut species_counts: Vec<usize> = Vec::new();
        for label in &labels {
            if let Some(pos) = species_order.iter().position(|s| s == label) {
                species_counts[pos] += 1;
            } else {
                species_order.push(label.clone());
                species_counts.push(1);
            }
        }

        Cif {
            lattice,
            species_names: species_order,
            species_counts,
            positions,
            occupancies,
            labels,
        }
    }

    pub fn from_file(path: &Path) -> Self {
        let contents = std::fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
        Self::parse(&contents)
    }

    pub fn types(&self) -> Vec<i32> {
        self.labels
            .iter()
            .map(|label| {
                let idx = self.species_names.iter().position(|s| s == label).unwrap();
                (idx + 1) as i32
            })
            .collect()
    }

    /// Check whether this CIF has partial occupancy (any occupancy < 1.0).
    pub fn has_partial_occupancy(&self) -> bool {
        self.occupancies.iter().any(|&o| o < 1.0 - 1e-6)
    }

    /// Merge co-located atoms into unique parent sites with allowed species.
    ///
    /// Returns `(parent_positions, all_species, orbit_allowed)` where:
    /// - `parent_positions[j]` is the fractional coordinate of parent site `j`
    /// - `all_species` is the ordered list of unique species across all orbits
    /// - `orbit_allowed[j][c]` is `true` iff species `c` can occupy orbit `j`
    ///
    /// Two CIF entries are considered co-located if their fractional coordinates
    /// differ by less than `tol` modulo lattice.
    pub fn partial_occupancy_orbits(
        &self,
        tol: f64,
    ) -> (Vec<[f64; 3]>, Vec<String>, Vec<Vec<bool>>) {
        let mut parent_positions: Vec<[f64; 3]> = Vec::new();
        let mut orbit_species: Vec<Vec<String>> = Vec::new();

        for (idx, pos) in self.positions.iter().enumerate() {
            if self.occupancies[idx] <= 1e-6 {
                continue;
            }
            let mut found = None;
            for (j, parent) in parent_positions.iter().enumerate() {
                let close = (0..3).all(|d| {
                    let diff = (pos[d] - parent[d]).rem_euclid(1.0);
                    diff < tol || diff > 1.0 - tol
                });
                if close {
                    found = Some(j);
                    break;
                }
            }

            if let Some(j) = found {
                let sp = &self.labels[idx];
                if !orbit_species[j].contains(sp) {
                    orbit_species[j].push(sp.clone());
                }
            } else {
                parent_positions.push(*pos);
                orbit_species.push(vec![self.labels[idx].clone()]);
            }
        }

        let mut all_species: Vec<String> = Vec::new();
        for orbit_sp in &orbit_species {
            for sp in orbit_sp {
                if !all_species.contains(sp) {
                    all_species.push(sp.clone());
                }
            }
        }

        let k = all_species.len();
        let orbit_allowed: Vec<Vec<bool>> = orbit_species
            .iter()
            .map(|orbit_sp| {
                let mut allowed = vec![false; k];
                for sp in orbit_sp {
                    let idx = all_species.iter().position(|s| s == sp).unwrap();
                    allowed[idx] = true;
                }
                allowed
            })
            .collect();

        (parent_positions, all_species, orbit_allowed)
    }

    /// Produce a types array for the merged parent sites.
    ///
    /// Each merged parent position gets a unique type index, used by moyo
    /// for symmetry detection. Atoms at the same merged position share the
    /// same type.
    pub fn partial_occupancy_types(&self, tol: f64) -> (Vec<[f64; 3]>, Vec<i32>) {
        let (parent_positions, _, _) = self.partial_occupancy_orbits(tol);
        let n = parent_positions.len();
        let types: Vec<i32> = (1..=n as i32).collect();
        (parent_positions, types)
    }
}

fn tag_value(line: &str) -> &str {
    line.split_whitespace()
        .nth(1)
        .expect("missing value after CIF tag")
}

fn parse_cif_float(s: &str) -> f64 {
    // CIF floats may have uncertainty in parens, e.g. "5.640(1)"
    let s = if let Some(idx) = s.find('(') {
        &s[..idx]
    } else {
        s
    };
    s.parse()
        .unwrap_or_else(|_| panic!("invalid CIF float: {s}"))
}

fn strip_oxidation(s: &str) -> String {
    // Remove trailing charge like "Fe3+" or "O2-"
    let bytes = s.as_bytes();
    let mut end = bytes.len();
    if end > 0 && (bytes[end - 1] == b'+' || bytes[end - 1] == b'-') {
        end -= 1;
        while end > 0 && bytes[end - 1].is_ascii_digit() {
            end -= 1;
        }
    }
    s[..end].to_string()
}

struct AtomSiteLoop {
    label_idx: Option<usize>,
    type_symbol_idx: Option<usize>,
    fract_x_idx: usize,
    fract_y_idx: usize,
    fract_z_idx: usize,
    occupancy_idx: Option<usize>,
}

impl AtomSiteLoop {
    fn from_columns(cols: &[String]) -> Option<Self> {
        let fx = cols.iter().position(|c| c == "_atom_site_fract_x")?;
        let fy = cols.iter().position(|c| c == "_atom_site_fract_y")?;
        let fz = cols.iter().position(|c| c == "_atom_site_fract_z")?;
        let label = cols.iter().position(|c| c == "_atom_site_label");
        let type_sym = cols.iter().position(|c| c == "_atom_site_type_symbol");
        let occupancy = cols.iter().position(|c| c == "_atom_site_occupancy");
        if label.is_none() && type_sym.is_none() {
            return None;
        }
        Some(AtomSiteLoop {
            label_idx: label,
            type_symbol_idx: type_sym,
            fract_x_idx: fx,
            fract_y_idx: fy,
            fract_z_idx: fz,
            occupancy_idx: occupancy,
        })
    }
}

fn parse_loop_header(lines: &[&str], start: usize) -> (Vec<String>, usize) {
    let mut cols = Vec::new();
    let mut i = start;
    while i < lines.len() {
        let line = lines[i].trim();
        if line.starts_with('_') {
            cols.push(line.to_string());
            i += 1;
        } else {
            break;
        }
    }
    (cols, i)
}

fn lattice_from_params(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> [[f64; 3]; 3] {
    let alpha_r = alpha.to_radians();
    let beta_r = beta.to_radians();
    let gamma_r = gamma.to_radians();

    let cos_alpha = alpha_r.cos();
    let cos_beta = beta_r.cos();
    let cos_gamma = gamma_r.cos();
    let sin_gamma = gamma_r.sin();

    // Standard crystallographic convention (a along x, b in xy-plane)
    let ax = a;
    let bx = b * cos_gamma;
    let by = b * sin_gamma;
    let cx = c * cos_beta;
    let cy = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
    let cz = (c * c - cx * cx - cy * cy).sqrt();

    [[ax, 0.0, 0.0], [bx, by, 0.0], [cx, cy, cz]]
}

#[cfg(test)]
mod tests {
    use super::*;

    const FCC_CIF: &str = "\
data_fcc
_cell_length_a   4.0
_cell_length_b   4.0
_cell_length_c   4.0
_cell_angle_alpha   90.0
_cell_angle_beta    90.0
_cell_angle_gamma   90.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
A1 A 0.0 0.0 0.0
A2 A 0.5 0.5 0.0
A3 A 0.5 0.0 0.5
A4 A 0.0 0.5 0.5
";

    #[test]
    fn parse_fcc() {
        let cif = Cif::parse(FCC_CIF);
        assert_eq!(cif.positions.len(), 4);
        assert_eq!(cif.species_names, vec!["A"]);
        assert_eq!(cif.species_counts, vec![4]);
        assert!((cif.lattice[0][0] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn parse_cif_float_with_uncertainty() {
        assert!((parse_cif_float("5.640(1)") - 5.64).abs() < 1e-10);
        assert!((parse_cif_float("90.0") - 90.0).abs() < 1e-10);
    }

    #[test]
    fn strip_oxidation_state() {
        assert_eq!(strip_oxidation("Fe3+"), "Fe");
        assert_eq!(strip_oxidation("O2-"), "O");
        assert_eq!(strip_oxidation("Na"), "Na");
        assert_eq!(strip_oxidation("Fe"), "Fe");
    }

    #[test]
    fn lattice_cubic() {
        let lat = lattice_from_params(4.0, 4.0, 4.0, 90.0, 90.0, 90.0);
        assert!((lat[0][0] - 4.0).abs() < 1e-10);
        assert!((lat[1][1] - 4.0).abs() < 1e-10);
        assert!((lat[2][2] - 4.0).abs() < 1e-10);
        assert!(lat[0][1].abs() < 1e-10);
    }

    const NACL_CIF: &str = "\
data_NaCl
_cell_length_a   5.64
_cell_length_b   5.64
_cell_length_c   5.64
_cell_angle_alpha   90.0
_cell_angle_beta    90.0
_cell_angle_gamma   90.0
loop_
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na 0.0 0.0 0.0
Na 0.5 0.5 0.0
Na 0.5 0.0 0.5
Na 0.0 0.5 0.5
Cl 0.5 0.5 0.5
Cl 0.0 0.0 0.5
Cl 0.0 0.5 0.0
Cl 0.5 0.0 0.0
";

    #[test]
    fn parse_nacl_cif() {
        let cif = Cif::parse(NACL_CIF);
        assert_eq!(cif.positions.len(), 8);
        assert_eq!(cif.species_names, vec!["Na", "Cl"]);
        assert_eq!(cif.species_counts, vec![4, 4]);
        assert_eq!(cif.types(), vec![1, 1, 1, 1, 2, 2, 2, 2]);
    }

    #[test]
    fn partial_occupancy_interleaved_species() {
        // Regression: species interleaved in file order must pair correctly
        // with positions. A and B share (0,0,0); A alone at (0.5,0.5,0.5).
        let cif_str = "\
data_test
_cell_length_a   4.0
_cell_length_b   4.0
_cell_length_c   4.0
_cell_angle_alpha   90.0
_cell_angle_beta    90.0
_cell_angle_gamma   90.0
loop_
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
A 0.0 0.0 0.0 0.5
B 0.0 0.0 0.0 0.5
A 0.5 0.5 0.5 1.0
";
        let cif = Cif::parse(cif_str);
        assert_eq!(cif.labels, vec!["A", "B", "A"]);

        let (parents, all_sp, allowed) = cif.partial_occupancy_orbits(0.01);
        assert_eq!(parents.len(), 2);
        assert_eq!(all_sp, vec!["A", "B"]);
        // orbit 0 at (0,0,0): A and B allowed
        assert_eq!(allowed[0], vec![true, true]);
        // orbit 1 at (0.5,0.5,0.5): only A allowed
        assert_eq!(allowed[1], vec![true, false]);
    }
}
