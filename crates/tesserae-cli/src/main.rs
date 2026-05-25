use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser)]
#[command(
    name = "tesserae",
    about = "Symmetry-inequivalent crystal structure enumeration"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Enumerate colorings of a supercell from a CIF file
    FromCif {
        /// Path to CIF file
        cif: PathBuf,

        /// Supercell matrix as semicolon-separated rows, e.g. "2,0,0;0,2,0;0,0,2"
        #[arg(long)]
        supercell: String,

        /// Composition as comma-separated counts, e.g. "4,4"
        #[arg(long)]
        composition: String,

        /// Symmetry precision
        #[arg(long, default_value = "1e-5")]
        symprec: f64,

        /// Print each coloring (default: count only)
        #[arg(long)]
        print: bool,

        /// Exclude superperiodic structures (colorings representable by a smaller supercell).
        /// Matches the enumlib convention.
        #[arg(long)]
        filter_superperiodic: bool,

        /// Write structure files to this directory (one per coloring)
        #[arg(long)]
        output_dir: Option<PathBuf>,

        /// Output format for structure files (default: cif, matching input)
        #[arg(long, default_value = "cif")]
        output_format: OutputFormat,

        /// Species names for coloring labels, e.g. "Na,Cl" (default: A,B,C,...)
        #[arg(long)]
        species: Option<String>,

        /// Multi-site enumeration: use full space group on all Wyckoff orbits simultaneously.
        /// Required for parent structures with multiple atoms per unit cell.
        #[arg(long)]
        multisite: bool,

        /// Constrained enumeration for CIF files with partial occupancy.
        /// Each Wyckoff orbit only allows the species listed in the CIF;
        /// co-located atoms (same position, different species) define the
        /// allowed set.
        #[arg(long)]
        constrained: bool,
    },

    /// Enumerate colorings of a supercell from a POSCAR file
    FromPoscar {
        /// Path to VASP POSCAR file
        poscar: PathBuf,

        /// Supercell matrix as semicolon-separated rows, e.g. "2,0,0;0,2,0;0,0,2"
        #[arg(long)]
        supercell: String,

        /// Composition as comma-separated counts, e.g. "4,4"
        #[arg(long)]
        composition: String,

        /// Symmetry precision
        #[arg(long, default_value = "1e-5")]
        symprec: f64,

        /// Print each coloring (default: count only)
        #[arg(long)]
        print: bool,

        /// Exclude superperiodic structures (colorings representable by a smaller supercell).
        /// Matches the enumlib convention.
        #[arg(long)]
        filter_superperiodic: bool,

        /// Write structure files to this directory (one per coloring)
        #[arg(long)]
        output_dir: Option<PathBuf>,

        /// Output format for structure files (default: poscar, matching input)
        #[arg(long, default_value = "poscar")]
        output_format: OutputFormat,

        /// Species names for coloring labels, e.g. "Na,Cl" (default: A,B,C,...)
        #[arg(long)]
        species: Option<String>,

        /// Multi-site enumeration: use full space group on all Wyckoff orbits simultaneously.
        /// Required for parent structures with multiple atoms per unit cell.
        #[arg(long)]
        multisite: bool,
    },

    /// Enumerate colorings across all inequivalent supercells at a given index, reading symmetry from a CIF file
    FromCifAtIndex {
        /// Path to CIF file
        cif: PathBuf,

        /// Volume index (determinant of supercell matrix)
        index: usize,

        /// Composition as comma-separated counts, e.g. "4,4"
        #[arg(long)]
        composition: String,

        /// Symmetry precision
        #[arg(long, default_value = "1e-5")]
        symprec: f64,

        /// Print each coloring (default: count only)
        #[arg(long)]
        print: bool,

        /// Exclude superperiodic structures (colorings representable by a smaller supercell).
        /// Matches the enumlib convention.
        #[arg(long)]
        filter_superperiodic: bool,
    },

    /// Enumerate colorings across all inequivalent supercells at a given index, reading symmetry from a POSCAR file
    FromPoscarAtIndex {
        /// Path to VASP POSCAR file
        poscar: PathBuf,

        /// Volume index (determinant of supercell matrix)
        index: usize,

        /// Composition as comma-separated counts, e.g. "4,4"
        #[arg(long)]
        composition: String,

        /// Symmetry precision
        #[arg(long, default_value = "1e-5")]
        symprec: f64,

        /// Print each coloring (default: count only)
        #[arg(long)]
        print: bool,

        /// Exclude superperiodic structures (colorings representable by a smaller supercell).
        /// Matches the enumlib convention.
        #[arg(long)]
        filter_superperiodic: bool,
    },

    /// Enumerate colorings across all inequivalent supercells at a given volume index
    AtIndex {
        /// Volume index (determinant of supercell matrix)
        index: usize,

        /// Number of species (k-ary enumeration)
        #[arg(long, default_value = "2")]
        nspecies: usize,

        /// Dimension (2 or 3)
        #[arg(long, default_value = "3")]
        dim: usize,

        /// Composition as comma-separated counts, e.g. "4,4". If omitted, sweeps all compositions.
        #[arg(long)]
        composition: Option<String>,

        /// Concentration range as "min0-max0,min1-max1,...", e.g. "1-3,1-3".
        /// Overrides --composition; sweeps all valid compositions within bounds.
        #[arg(long)]
        concentration_range: Option<String>,

        /// Parent point group generators as semicolon-separated matrices.
        /// Each matrix is rows separated by '/', elements by commas.
        /// e.g. "0,-1,0/1,0,0/0,0,1;-1,0,0/0,-1,0/0,0,-1" for C4 + inversion
        #[arg(long)]
        point_group: String,

        /// Print each coloring (default: count only)
        #[arg(long)]
        print: bool,

        /// Exclude superperiodic structures (colorings representable by a smaller supercell).
        /// Matches the enumlib convention.
        #[arg(long)]
        filter_superperiodic: bool,

        /// Reservoir-sample N colorings uniformly at random instead of enumerating all.
        #[arg(long)]
        sample: Option<usize>,

        /// RNG seed for --sample (default: 42).
        #[arg(long, default_value = "42")]
        seed: u64,
    },
}

#[derive(Clone, ValueEnum)]
enum OutputFormat {
    Poscar,
    Cif,
}

fn default_species(k: usize) -> Vec<String> {
    (0..k)
        .map(|i| {
            if i < 26 {
                String::from((b'A' + i as u8) as char)
            } else {
                format!("X{i}")
            }
        })
        .collect()
}

fn parse_species(s: &str) -> Vec<String> {
    s.split(',').map(|v| v.trim().to_string()).collect()
}

fn write_structure(
    format: &OutputFormat,
    lattice: &[[f64; 3]; 3],
    supercell_matrix: &[Vec<i64>],
    species_names: &[String],
    coloring: &[u32],
    index: u64,
    dir: &std::path::Path,
) {
    let content = match format {
        OutputFormat::Poscar => {
            let comment = format!("coloring_{index}");
            tesserae_io::write::write_poscar(
                lattice,
                supercell_matrix,
                species_names,
                coloring,
                &comment,
            )
        }
        OutputFormat::Cif => {
            let data_name = format!("coloring_{index}");
            tesserae_io::write::write_cif(
                lattice,
                supercell_matrix,
                species_names,
                coloring,
                &data_name,
            )
        }
    };
    let ext = match format {
        OutputFormat::Poscar => "vasp",
        OutputFormat::Cif => "cif",
    };
    let path = dir.join(format!("coloring_{index}.{ext}"));
    std::fs::write(&path, content)
        .unwrap_or_else(|e| panic!("failed to write {}: {e}", path.display()));
}

fn parse_matrix(s: &str) -> Vec<Vec<i64>> {
    s.split(';')
        .map(|row| {
            row.split(',')
                .map(|v| v.trim().parse::<i64>().expect("invalid integer in matrix"))
                .collect()
        })
        .collect()
}

fn parse_composition(s: &str) -> Vec<usize> {
    s.split(',')
        .map(|v| {
            v.trim()
                .parse::<usize>()
                .expect("invalid integer in composition")
        })
        .collect()
}

fn parse_point_group(s: &str, dim: usize) -> Vec<Vec<Vec<i64>>> {
    let generators: Vec<Vec<Vec<i64>>> = s
        .split(';')
        .map(|mat| {
            mat.split('/')
                .map(|row| {
                    row.split(',')
                        .map(|v| {
                            v.trim()
                                .parse::<i64>()
                                .expect("invalid integer in point group")
                        })
                        .collect()
                })
                .collect()
        })
        .collect();

    tesserae_core::supercell::generate_matrix_group(&generators, dim)
}

fn parse_concentration_range(s: &str) -> (Vec<usize>, Vec<usize>) {
    let mut mins = Vec::new();
    let mut maxs = Vec::new();
    for part in s.split(',') {
        let bounds: Vec<&str> = part.trim().split('-').collect();
        match bounds.len() {
            1 => {
                let v: usize = bounds[0]
                    .parse()
                    .expect("invalid integer in concentration range");
                mins.push(v);
                maxs.push(v);
            }
            2 => {
                let lo: usize = bounds[0]
                    .parse()
                    .expect("invalid integer in concentration range");
                let hi: usize = bounds[1]
                    .parse()
                    .expect("invalid integer in concentration range");
                mins.push(lo);
                maxs.push(hi);
            }
            _ => panic!("invalid concentration range format: {part}"),
        }
    }
    (mins, maxs)
}

fn all_compositions(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut current = vec![0usize; k];
    fn recurse(
        k: usize,
        pos: usize,
        remaining: usize,
        current: &mut Vec<usize>,
        result: &mut Vec<Vec<usize>>,
    ) {
        if pos == k - 1 {
            current[pos] = remaining;
            result.push(current.clone());
            return;
        }
        for i in 0..=remaining {
            current[pos] = i;
            recurse(k, pos + 1, remaining - i, current, result);
        }
    }
    recurse(k, 0, n, &mut current, &mut result);
    result
}

#[allow(clippy::too_many_arguments)]
fn run_at_index_enumeration(
    pg: &[Vec<Vec<i64>>],
    index: usize,
    dim: usize,
    nspecies: usize,
    composition: Option<&str>,
    concentration_range: Option<&str>,
    filter_superperiodic: bool,
    cb: &mut impl FnMut(&[u32]),
) {
    if let Some(range_str) = concentration_range {
        let (mins, maxs) = parse_concentration_range(range_str);
        if mins.len() != nspecies {
            eprintln!(
                "error: --concentration-range has {} entries but --nspecies is {}",
                mins.len(),
                nspecies
            );
            std::process::exit(1);
        }
        if filter_superperiodic {
            tesserae_core::supercell::enumerate_at_index_range_primitive(
                pg, index, &mins, &maxs, cb,
            );
        } else {
            tesserae_core::supercell::enumerate_at_index_range(pg, index, &mins, &maxs, cb);
        }
    } else if let Some(comp_str) = composition {
        let comp = parse_composition(comp_str);
        if filter_superperiodic {
            tesserae_core::supercell::enumerate_at_index_primitive(pg, index, &comp, cb);
        } else {
            tesserae_core::supercell::enumerate_at_index(pg, index, &comp, cb);
        }
    } else {
        let hnfs = if dim == 3 {
            tesserae_core::hnf::hnf_3d(index)
        } else {
            tesserae_core::hnf::hnf_2d(index)
        };
        let inequiv = tesserae_core::supercell::inequivalent_hnfs(&hnfs, pg);
        for comp in all_compositions(index, nspecies) {
            for s in &inequiv {
                if filter_superperiodic {
                    tesserae_core::supercell::enumerate_supercell_primitive(s, &comp, pg, cb);
                } else {
                    tesserae_core::supercell::enumerate_supercell(s, &comp, pg, cb);
                }
            }
        }
    }
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Command::FromCif {
            cif,
            supercell,
            composition,
            symprec,
            print,
            filter_superperiodic,
            output_dir,
            output_format,
            species,
            multisite,
            constrained,
        } => {
            let supercell_matrix = parse_matrix(&supercell);
            let comp = parse_composition(&composition);
            let structure = tesserae_io::cif::Cif::from_file(&cif);
            let species_names = species
                .map(|s| parse_species(&s))
                .unwrap_or_else(|| default_species(comp.len()));
            if species_names.len() != comp.len() {
                eprintln!(
                    "error: --species has {} names but --composition has {} entries",
                    species_names.len(),
                    comp.len()
                );
                std::process::exit(1);
            }
            if output_dir.is_some() && (multisite || constrained) {
                eprintln!(
                    "error: --output-dir is not yet supported with --multisite or --constrained"
                );
                std::process::exit(1);
            }
            if let Some(ref dir) = output_dir {
                std::fs::create_dir_all(dir)
                    .unwrap_or_else(|e| panic!("failed to create {}: {e}", dir.display()));
            }
            let start = Instant::now();
            let mut count: u64 = 0;
            let cb = &mut |coloring: &[u32]| {
                count += 1;
                if print {
                    let s: Vec<String> = coloring.iter().map(|c| c.to_string()).collect();
                    println!("{}", s.join(" "));
                }
                if let Some(ref dir) = output_dir {
                    write_structure(
                        &output_format,
                        &structure.lattice,
                        &supercell_matrix,
                        &species_names,
                        coloring,
                        count,
                        dir,
                    );
                }
            };
            if constrained {
                let (_, all_species, _) = structure.partial_occupancy_orbits(symprec);
                if comp.len() != all_species.len() {
                    eprintln!(
                        "error: --composition has {} entries but CIF partial occupancy has {} species ({})",
                        comp.len(),
                        all_species.len(),
                        all_species.join(", ")
                    );
                    std::process::exit(1);
                }
                tesserae_io::enumerate_from_cif_constrained(
                    &cif,
                    &supercell_matrix,
                    &comp,
                    symprec,
                    cb,
                );
            } else if multisite && filter_superperiodic {
                tesserae_io::enumerate_from_cif_multisite_primitive(
                    &cif,
                    &supercell_matrix,
                    &comp,
                    symprec,
                    cb,
                );
            } else if multisite {
                tesserae_io::enumerate_from_cif_multisite(
                    &cif,
                    &supercell_matrix,
                    &comp,
                    symprec,
                    cb,
                );
            } else if filter_superperiodic {
                tesserae_io::enumerate_from_cif_primitive(
                    &cif,
                    &supercell_matrix,
                    &comp,
                    symprec,
                    cb,
                );
            } else {
                tesserae_io::enumerate_from_cif(&cif, &supercell_matrix, &comp, symprec, cb);
            }
            eprintln!("{count} colorings in {:.3?}", start.elapsed());
        }

        Command::FromPoscar {
            poscar,
            supercell,
            composition,
            symprec,
            print,
            filter_superperiodic,
            output_dir,
            output_format,
            species,
            multisite,
        } => {
            let supercell_matrix = parse_matrix(&supercell);
            let comp = parse_composition(&composition);
            let structure = tesserae_io::poscar::Poscar::from_file(&poscar);
            let species_names = species
                .map(|s| parse_species(&s))
                .unwrap_or_else(|| default_species(comp.len()));
            if species_names.len() != comp.len() {
                eprintln!(
                    "error: --species has {} names but --composition has {} entries",
                    species_names.len(),
                    comp.len()
                );
                std::process::exit(1);
            }
            if output_dir.is_some() && multisite {
                eprintln!("error: --output-dir is not yet supported with --multisite");
                std::process::exit(1);
            }
            if let Some(ref dir) = output_dir {
                std::fs::create_dir_all(dir)
                    .unwrap_or_else(|e| panic!("failed to create {}: {e}", dir.display()));
            }
            let start = Instant::now();
            let mut count: u64 = 0;
            let cb = &mut |coloring: &[u32]| {
                count += 1;
                if print {
                    let s: Vec<String> = coloring.iter().map(|c| c.to_string()).collect();
                    println!("{}", s.join(" "));
                }
                if let Some(ref dir) = output_dir {
                    write_structure(
                        &output_format,
                        &structure.lattice,
                        &supercell_matrix,
                        &species_names,
                        coloring,
                        count,
                        dir,
                    );
                }
            };
            if multisite && filter_superperiodic {
                tesserae_io::enumerate_from_poscar_multisite_primitive(
                    &poscar,
                    &supercell_matrix,
                    &comp,
                    symprec,
                    cb,
                );
            } else if multisite {
                tesserae_io::enumerate_from_poscar_multisite(
                    &poscar,
                    &supercell_matrix,
                    &comp,
                    symprec,
                    cb,
                );
            } else if filter_superperiodic {
                tesserae_io::enumerate_from_poscar_primitive(
                    &poscar,
                    &supercell_matrix,
                    &comp,
                    symprec,
                    cb,
                );
            } else {
                tesserae_io::enumerate_from_poscar(&poscar, &supercell_matrix, &comp, symprec, cb);
            }
            eprintln!("{count} colorings in {:.3?}", start.elapsed());
        }

        Command::FromCifAtIndex {
            cif,
            index,
            composition,
            symprec,
            print,
            filter_superperiodic,
        } => {
            let comp = parse_composition(&composition);
            let start = Instant::now();
            let mut count: u64 = 0;
            let cb = &mut |coloring: &[u32]| {
                count += 1;
                if print {
                    let s: Vec<String> = coloring.iter().map(|c| c.to_string()).collect();
                    println!("{}", s.join(" "));
                }
            };
            if filter_superperiodic {
                tesserae_io::enumerate_from_cif_at_index_primitive(&cif, index, &comp, symprec, cb);
            } else {
                tesserae_io::enumerate_from_cif_at_index(&cif, index, &comp, symprec, cb);
            }
            eprintln!("{count} colorings in {:.3?}", start.elapsed());
        }

        Command::FromPoscarAtIndex {
            poscar,
            index,
            composition,
            symprec,
            print,
            filter_superperiodic,
        } => {
            let comp = parse_composition(&composition);
            let start = Instant::now();
            let mut count: u64 = 0;
            let cb = &mut |coloring: &[u32]| {
                count += 1;
                if print {
                    let s: Vec<String> = coloring.iter().map(|c| c.to_string()).collect();
                    println!("{}", s.join(" "));
                }
            };
            if filter_superperiodic {
                tesserae_io::enumerate_from_poscar_at_index_primitive(
                    &poscar, index, &comp, symprec, cb,
                );
            } else {
                tesserae_io::enumerate_from_poscar_at_index(&poscar, index, &comp, symprec, cb);
            }
            eprintln!("{count} colorings in {:.3?}", start.elapsed());
        }

        Command::AtIndex {
            index,
            nspecies,
            dim,
            composition,
            concentration_range,
            point_group,
            print,
            filter_superperiodic,
            sample,
            seed,
        } => {
            let pg = parse_point_group(&point_group, dim);
            let start = Instant::now();

            if let Some(k) = sample {
                let mut sampler = tesserae_core::enumerate::ReservoirSampler::new(k, seed);
                let cb = &mut |coloring: &[u32]| sampler.observe(coloring);
                run_at_index_enumeration(
                    &pg,
                    index,
                    dim,
                    nspecies,
                    composition.as_deref(),
                    concentration_range.as_deref(),
                    filter_superperiodic,
                    cb,
                );
                let samples = sampler.into_samples();
                for s in &samples {
                    let formatted: Vec<String> = s.iter().map(|c| c.to_string()).collect();
                    println!("{}", formatted.join(" "));
                }
                eprintln!(
                    "{} sampled colorings in {:.3?}",
                    samples.len(),
                    start.elapsed()
                );
            } else {
                let mut count: u64 = 0;
                let cb = &mut |coloring: &[u32]| {
                    count += 1;
                    if print {
                        let s: Vec<String> = coloring.iter().map(|c| c.to_string()).collect();
                        println!("{}", s.join(" "));
                    }
                };
                run_at_index_enumeration(
                    &pg,
                    index,
                    dim,
                    nspecies,
                    composition.as_deref(),
                    concentration_range.as_deref(),
                    filter_superperiodic,
                    cb,
                );
                eprintln!("{count} colorings in {:.3?}", start.elapsed());
            }
        }
    }
}
