# tesserae

> **Note:** This project is in active development.

A Rust library for symmetry-inequivalent crystal structure enumeration via the semidirect product decomposition G = T ⋊ P of crystallographic symmetry groups. Translation symmetry is canonicalized algebraically using necklace algorithms (Booth 1980), reducing McKay-style canonical augmentation to the small point group P (|P| ≤ 48).

## Algorithm

Crystallographic symmetry groups acting on a supercell always decompose as a semidirect product G = T ⋊ P, where T is the Abelian translation subgroup (isomorphic to ℤ_{d₁} × ℤ_{d₂} × ℤ_{d₃} via Smith Normal Form) and P is the point group (|P| ≤ 48 in 3D). tesserae exploits this structure in two phases:

1. **T-canonicalization** — Find the lex-minimum coloring under all translations using Booth's O(n) necklace algorithm. This absorbs the dominant factor |T| at sublinear cost.
2. **P-canonical augmentation** — Apply McKay's canonical augmentation restricted to the small point group P, building colorings site-by-site and pruning non-canonical branches.

tesserae treats the translation subgroup algebraically rather than as an opaque list of permutations.

## Features

- Binary and k-ary enumeration on arbitrary 3D supercells
- Multi-site Wyckoff enumeration (multiple crystallographic orbits with independent compositions)
- Partial occupancy / constrained enumeration (per-site species constraints from CIF `_atom_site_occupancy`)
- CIF and POSCAR input with automatic symmetry detection via [moyo](https://github.com/spglib/moyo)
- Structure output (write enumerated structures as POSCAR or CIF)
- Concentration range sweeps
- Superperiodic filtering (primitive-only structures, matching enumlib convention)
- Reservoir sampling for uniform random structure selection
- Orbit multiplicity reporting
- Pólya cross-check on every run (correctness guarantee)
- CLI, Python bindings (PyO3), and Rust library API

## Installation

### Prerequisites

- Rust toolchain (≥ 1.80): [rustup.rs](https://rustup.rs)
- For Python bindings: [maturin](https://github.com/PyO3/maturin) and a Python ≥ 3.9 environment

### CLI

```bash
cargo install --path crates/tesserae-cli
```

### Python bindings

```bash
# In a conda or venv environment
pip install maturin
maturin develop -m crates/tesserae-python/Cargo.toml
```

### Rust library

Add to your `Cargo.toml`:

```toml
[dependencies]
tesserae-core = { git = "https://github.com/archimono/tesserae" }
tesserae-io   = { git = "https://github.com/archimono/tesserae" }  # for CIF/POSCAR support
```

## Usage

### CLI

**Enumerate a fixed supercell from a CIF file:**

```bash
tesserae from-cif structure.cif \
  --supercell "2,0,0;0,2,0;0,0,2" \
  --composition "4,4"
```

**Sweep all inequivalent supercells at a volume index:**

```bash
tesserae from-cif-at-index structure.cif 8 \
  --composition "4,4" \
  --filter-superperiodic
```

**Multi-site enumeration (parent cells with multiple Wyckoff orbits):**

```bash
tesserae from-cif structure.cif \
  --supercell "2,0,0;0,2,0;0,0,2" \
  --composition "4,4" \
  --multisite
```

**Constrained enumeration from a CIF with partial occupancy:**

```bash
tesserae from-cif structure.cif \
  --supercell "2,0,0;0,2,0;0,0,2" \
  --composition "4,4" \
  --constrained
```

**Write enumerated structures to files:**

```bash
tesserae from-cif structure.cif \
  --supercell "2,0,0;0,2,0;0,0,2" \
  --composition "4,4" \
  --output-dir out/ \
  --output-format cif \
  --species "Na,Cl"
```

**Concentration range sweep:**

```bash
tesserae at-index 8 \
  --point-group "$OH" \
  --concentration-range "1-6,2-6" \
  --filter-superperiodic
```

**Reservoir sampling (uniform random selection):**

```bash
tesserae at-index 12 \
  --point-group "$OH" \
  --sample 100 \
  --seed 42
```

**Provide the point group manually (no CIF required):**

```bash
# SC Oh generators: C4x + C4z + inversion
OH="1,0,0/0,0,-1/0,1,0;0,-1,0/1,0,0/0,0,1;-1,0,0/0,-1,0/0,0,-1"

tesserae at-index 8 \
  --nspecies 3 \
  --point-group "$OH" \
  --filter-superperiodic
```

### Python

```python
import tesserae

# Fixed supercell from a CIF file
colorings = tesserae.enumerate_from_cif(
    "structure.cif",
    supercell_matrix=[[2,0,0],[0,2,0],[0,0,2]],
    composition=[4, 4],
)

# Sweep all supercell shapes at a volume index
colorings = tesserae.enumerate_from_cif_at_index(
    "structure.cif",
    index=8,
    composition=[4, 4],
)

# Primitive (non-superperiodic) structures only
colorings = tesserae.enumerate_from_cif_at_index_primitive(
    "structure.cif",
    index=8,
    composition=[4, 4],
)

# Provide point group directly (48 Oh matrices)
oh = [...]  # list of 3x3 integer matrices
colorings = tesserae.enumerate_at_index(oh, index=8, composition=[4, 4])
```

Each coloring is a list of integer species labels (0-indexed) of length equal to the number of sites in the supercell.

Multi-site and constrained variants are also available: `enumerate_from_cif_multisite`, `enumerate_from_cif_constrained`. Structure writers (`write_poscar`, `write_cif`) currently support single-site colorings only.

### Rust

```rust
use tesserae_core::supercell::{enumerate_at_index, enumerate_supercell};

// Enumerate a fixed supercell given the parent point group
let supercell = vec![vec![2,0,0], vec![0,2,0], vec![0,0,2]];
let composition = vec![4, 4];
let point_group: Vec<Vec<Vec<i64>>> = /* 48 Oh matrices */;

let mut count = 0;
enumerate_supercell(&supercell, &composition, &point_group, &mut |_coloring| {
    count += 1;
});

// Sweep all supercell shapes at a volume index
enumerate_at_index(&point_group, 8, &composition, &mut |coloring| {
    println!("{coloring:?}");
});
```

For CIF/POSCAR input, use `tesserae-io`:

```rust
use tesserae_io::{enumerate_from_cif, enumerate_from_cif_at_index};

enumerate_from_cif(
    "structure.cif".as_ref(),
    &supercell,
    &composition,
    1e-5,
    &mut |_| { count += 1; },
);
```

Multi-site and constrained variants: `enumerate_from_cif_multisite`, `enumerate_from_cif_constrained`, `enumerate_from_poscar_multisite`.

## Performance

tesserae is 8–191× faster than SHRY for fixed-supercell enumeration, and 1.3–3× faster than enumlib on at-index sweeps at larger problem sizes.

**vs SHRY (fixed supercell)**

| Benchmark | Sites | Structures | tesserae | SHRY | Speedup |
|-----------|------:|-----------:|---------:|-----:|--------:|
| SC 2×2×2 binary [4,4] | 8 | 6 | 10 ms | 1.3 s | 133× |
| SC 3×3×3 binary [13,14] | 27 | 16,384 | 261 ms | 2.0 s | 8× |
| FCC 2×2×2 binary [4,4] | 8 | 4 | 10 ms | 2.0 s | 191× |
| FCC 3×3×3 binary [13,14] | 27 | 16,384 | 326 ms | 5.9 s | 18× |

SHRY timings include Python startup (~1 s).

**vs enumlib (at-index sweep, SC Oh, both compiled with -O3)**

| Index | Structures | tesserae | enumlib | Speedup |
|------:|-----------:|---------:|--------:|--------:|
| 8 | 491 | 21 ms | 26 ms | 1.3× |
| 12 | 8,734 | 61 ms | 162 ms | 2.6× |
| 16 | 115,845 | 374 ms | 658 ms | 1.8× |

tesserae's advantage grows with problem size. Full results (including ternary and FCC sweeps) in [`bench_data/benchmark_results.md`](bench_data/benchmark_results.md).

