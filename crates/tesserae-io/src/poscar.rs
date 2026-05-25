use std::path::Path;

/// A crystal structure parsed from a POSCAR/CONTCAR file.
#[derive(Debug, Clone)]
pub struct Poscar {
    pub comment: String,
    pub lattice: [[f64; 3]; 3],
    pub species_names: Vec<String>,
    pub species_counts: Vec<usize>,
    pub positions: Vec<[f64; 3]>,
    pub is_cartesian: bool,
}

impl Poscar {
    /// Parse a POSCAR file from a string.
    ///
    /// Supports VASP 5+ format (with species names line).
    ///
    /// # Panics
    ///
    /// Panics on malformed input: missing lines, non-numeric values,
    /// or unsupported selective dynamics.
    pub fn parse(s: &str) -> Self {
        let lines: Vec<&str> = s.lines().collect();
        assert!(lines.len() >= 7, "POSCAR too short");

        let comment = lines[0].to_string();

        let scale: f64 = lines[1].trim().parse().expect("invalid scale factor");
        assert!(scale > 0.0, "scale factor must be positive");

        let lattice = [
            parse_vec3(lines[2], scale),
            parse_vec3(lines[3], scale),
            parse_vec3(lines[4], scale),
        ];

        // VASP 5: line 5 is species names, line 6 is counts
        // VASP 4: line 5 is counts directly (all numeric)
        let (species_names, counts_line_idx) = if is_species_line(lines[5]) {
            let names: Vec<String> = lines[5].split_whitespace().map(String::from).collect();
            (names, 6)
        } else {
            panic!(
                "VASP 4 format (no species names) is not supported; add species names on line 6"
            );
        };

        let species_counts: Vec<usize> = lines[counts_line_idx]
            .split_whitespace()
            .map(|s| s.parse().expect("invalid species count"))
            .collect();
        assert_eq!(
            species_names.len(),
            species_counts.len(),
            "species names and counts must have same length"
        );

        let total_atoms: usize = species_counts.iter().sum();

        let mut coord_line = counts_line_idx + 1;
        let coord_type = lines[coord_line].trim();

        // Skip selective dynamics line if present
        if coord_type.starts_with('S') || coord_type.starts_with('s') {
            coord_line += 1;
        }

        let is_cartesian = {
            let c = lines[coord_line].trim().as_bytes()[0];
            c == b'C' || c == b'c' || c == b'K' || c == b'k'
        };
        coord_line += 1;

        let mut positions = Vec::with_capacity(total_atoms);
        for i in 0..total_atoms {
            let vals: Vec<f64> = lines[coord_line + i]
                .split_whitespace()
                .take(3)
                .map(|s| s.parse().expect("invalid coordinate"))
                .collect();
            assert_eq!(
                vals.len(),
                3,
                "expected 3 coordinates on line {}",
                coord_line + i + 1
            );
            if is_cartesian {
                positions.push([vals[0] * scale, vals[1] * scale, vals[2] * scale]);
            } else {
                positions.push([vals[0], vals[1], vals[2]]);
            }
        }

        Poscar {
            comment,
            lattice,
            species_names,
            species_counts,
            positions,
            is_cartesian,
        }
    }

    /// Parse a POSCAR file from disk.
    pub fn from_file(path: &Path) -> Self {
        let contents = std::fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
        Self::parse(&contents)
    }

    /// Convert Cartesian positions to fractional coordinates.
    /// Returns fractional positions. If already fractional, returns a clone.
    pub fn fractional_positions(&self) -> Vec<[f64; 3]> {
        if !self.is_cartesian {
            return self.positions.clone();
        }
        // VASP convention: rows of L are lattice vectors, so r = f·L (row-vector).
        // In column-vector form: r = Lᵀ·f, hence f = (Lᵀ)⁻¹·r.
        let lt = transpose_3x3(&self.lattice);
        let inv_lt = invert_3x3(&lt);
        self.positions
            .iter()
            .map(|p| mat_vec_mul(&inv_lt, p))
            .collect()
    }

    /// Species type array: maps each atom to a 1-based species index.
    pub fn types(&self) -> Vec<i32> {
        let mut types = Vec::with_capacity(self.positions.len());
        for (i, &count) in self.species_counts.iter().enumerate() {
            for _ in 0..count {
                types.push((i + 1) as i32);
            }
        }
        types
    }
}

fn parse_vec3(line: &str, scale: f64) -> [f64; 3] {
    let vals: Vec<f64> = line
        .split_whitespace()
        .map(|s| s.parse().expect("invalid lattice coordinate"))
        .collect();
    assert_eq!(vals.len(), 3, "lattice vector must have 3 components");
    [vals[0] * scale, vals[1] * scale, vals[2] * scale]
}

fn is_species_line(line: &str) -> bool {
    line.split_whitespace()
        .next()
        .is_some_and(|tok| tok.parse::<usize>().is_err())
}

fn transpose_3x3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

fn invert_3x3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    assert!(det.abs() > 1e-12, "lattice matrix is singular");
    let inv_det = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ]
}

fn mat_vec_mul(m: &[[f64; 3]; 3], v: &[f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    const NACL_POSCAR: &str = "\
NaCl rock salt
5.64
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0
Na Cl
1 1
Direct
0.0 0.0 0.0
0.5 0.5 0.5
";

    #[test]
    fn parse_nacl() {
        let p = Poscar::parse(NACL_POSCAR);
        assert_eq!(p.comment, "NaCl rock salt");
        assert_eq!(p.species_names, vec!["Na", "Cl"]);
        assert_eq!(p.species_counts, vec![1, 1]);
        assert_eq!(p.positions.len(), 2);
        assert!(!p.is_cartesian);
        assert_eq!(p.types(), vec![1, 2]);
    }

    #[test]
    fn fractional_roundtrip() {
        let p = Poscar::parse(NACL_POSCAR);
        let frac = p.fractional_positions();
        assert!((frac[0][0]).abs() < 1e-10);
        assert!((frac[1][0] - 0.5).abs() < 1e-10);
    }

    const CUBIC_POSCAR: &str = "\
Simple cubic
4.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Si
1
Direct
0.0 0.0 0.0
";

    #[test]
    fn parse_cubic() {
        let p = Poscar::parse(CUBIC_POSCAR);
        assert_eq!(p.lattice[0], [4.0, 0.0, 0.0]);
        assert_eq!(p.species_counts, vec![1]);
    }

    const SELECTIVE_POSCAR: &str = "\
With selective dynamics
1.0
4.0 0.0 0.0
0.0 4.0 0.0
0.0 0.0 4.0
Si
2
Selective dynamics
Direct
0.0 0.0 0.0 T T T
0.5 0.5 0.5 F F F
";

    #[test]
    fn parse_selective_dynamics() {
        let p = Poscar::parse(SELECTIVE_POSCAR);
        assert_eq!(p.positions.len(), 2);
        assert!((p.positions[1][0] - 0.5).abs() < 1e-10);
    }

    const CARTESIAN_POSCAR: &str = "\
Cartesian coords
1.0
4.0 0.0 0.0
0.0 4.0 0.0
0.0 0.0 4.0
Si
2
Cartesian
0.0 0.0 0.0
2.0 2.0 2.0
";

    #[test]
    fn parse_cartesian() {
        let p = Poscar::parse(CARTESIAN_POSCAR);
        assert!(p.is_cartesian);
        let frac = p.fractional_positions();
        assert!((frac[1][0] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn cartesian_hexagonal() {
        // Hexagonal lattice (non-symmetric L) with Cartesian coordinates.
        // a1 = (3, 0, 0), a2 = (-1.5, 2.598..., 0), a3 = (0, 0, 5)
        let a = 3.0_f64;
        let half_sqrt3 = a * 3.0_f64.sqrt() / 2.0;
        // Atom at Cartesian (0, 0, 0) => fractional (0, 0, 0)
        // Atom at Cartesian a2 = (-1.5, 2.598..., 0) => fractional (0, 1, 0)
        let poscar = format!(
            "\
Hexagonal
1.0
{a} 0.0 0.0
{half_a} {half_sqrt3} 0.0
0.0 0.0 5.0
A
2
Cartesian
0.0 0.0 0.0
{half_a} {half_sqrt3} 0.0
",
            half_a = -a / 2.0,
        );
        let p = Poscar::parse(&poscar);
        assert!(p.is_cartesian);
        let frac = p.fractional_positions();
        // First atom: (0,0,0)
        assert!(frac[0][0].abs() < 1e-10);
        assert!(frac[0][1].abs() < 1e-10);
        assert!(frac[0][2].abs() < 1e-10);
        // Second atom: (0, 1, 0) in fractional
        assert!(frac[1][0].abs() < 1e-10, "expected 0, got {}", frac[1][0]);
        assert!(
            (frac[1][1] - 1.0).abs() < 1e-10,
            "expected 1, got {}",
            frac[1][1]
        );
        assert!(frac[1][2].abs() < 1e-10, "expected 0, got {}", frac[1][2]);
    }

    #[test]
    fn cartesian_monoclinic_general_position() {
        // Monoclinic lattice (β ≠ 90°) with a general-position atom.
        // a=4, b=5, c=6, β=100°
        // a1 = (4, 0, 0), a2 = (0, 5, 0), a3 = (6 cos 100°, 0, 6 sin 100°)
        let beta = 100.0_f64.to_radians();
        let a3x = 6.0 * beta.cos();
        let a3z = 6.0 * beta.sin();

        // Fractional (0.3, 0.5, 0.7) => Cartesian:
        // r = 0.3*(4,0,0) + 0.5*(0,5,0) + 0.7*(a3x, 0, a3z)
        let rx = 0.3 * 4.0 + 0.7 * a3x;
        let ry = 0.5 * 5.0;
        let rz = 0.7 * a3z;

        let poscar = format!(
            "\
Monoclinic
1.0
4.0 0.0 0.0
0.0 5.0 0.0
{a3x} 0.0 {a3z}
X
1
Cartesian
{rx} {ry} {rz}
"
        );
        let p = Poscar::parse(&poscar);
        assert!(p.is_cartesian);
        let frac = p.fractional_positions();
        assert!(
            (frac[0][0] - 0.3).abs() < 1e-10,
            "expected f0=0.3, got {}",
            frac[0][0]
        );
        assert!(
            (frac[0][1] - 0.5).abs() < 1e-10,
            "expected f1=0.5, got {}",
            frac[0][1]
        );
        assert!(
            (frac[0][2] - 0.7).abs() < 1e-10,
            "expected f2=0.7, got {}",
            frac[0][2]
        );
    }
}
