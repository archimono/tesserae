# Publication Benchmark Plan

## Context

The existing `bench_data/benchmark_results.md` has stale numbers from before the public release. The user wants publication-quality benchmarks for a journal paper, extended to larger problem sizes and covering new features (multi-site, constrained enumeration). All tools are available: hyperfine, SHRY 1.1.8, enumlib (`/tmp/enumlib/src/enum.x`), Criterion 0.5.

## Phase 0: Feasibility Verification (before any script changes)

Run untimed dry-runs to determine what's practical:

1. **Extended indices** — test `tesserae at-index 14` and `at-index 16` (SC binary, `--filter-superperiodic`), plus enumlib at idx 14/16. If enumlib exceeds 5 min, mark those rows "tesserae-only".
2. **Ternary idx 10/12** — test both tools. Ternary idx 12 has 91 compositions; enumlib may be too slow.
3. **SHRY multi-site** — test `shry nacl_prim.cif --from-species Na --to-species Na4Cu4 --scaling-matrix 2 2 2 --count-only --no-write`. Determines if SHRY is a valid multi-site competitor.
4. **SHRY constrained** — test with a partial-occupancy CIF. Determines if SHRY handles this workflow.
5. **SC 4×4×4** — test `tesserae from-cif sc.cif --supercell 4,0,0;0,4,0;0,0,4 --composition 32,32`. Check if feasible for SHRY comparison.

## Phase 1: New CIF Files

Create in `bench_data/`:

- **`nacl_prim.cif`** — NaCl rock salt (Fm-3m #225), FCC primitive, 2 atoms (Na + Cl). For multi-site benchmarks: binary substitution on Na sublattice, Cl fixed.
- **`perovskite.cif`** — SrTiO3 (Pm-3m #221), 5 atoms (Sr, Ti, O×3). For multi-site: binary on Sr site, Ti+O fixed.
- **`partial_occ.cif`** — Na/K partial occupancy on same site + Cl (Fm-3m #225). For constrained enumeration.

Verify each with `tesserae from-cif --multisite` / `--constrained` before benchmarking.

## Phase 2: Benchmark Groups

### Group A — Binary at-index sweep (tesserae vs enumlib)

Both tools enumerate all compositions per index with superperiodic filtering.

| Lattice | Indices |
|---------|---------|
| SC (Oh) | 2, 4, 6, 8, 10, 12, 14, 16 (14/16 if enumlib feasible) |
| FCC (Oh) | 2, 4, 6, 8, 10, 12, 14, 16 (same caveat) |

Use `at-index` with explicit point-group generators (pure algorithm, no CIF overhead). enumlib via `struct_enum.in`. Count agreement verified for every row.

### Group B — Ternary at-index sweep (tesserae vs enumlib)

| Lattice | Indices |
|---------|---------|
| SC | 4, 6, 8, 10 (12 if enumlib feasible) |
| FCC | 4, 6, 8, 10 (same caveat) |

### Group C — Fixed supercell (tesserae vs SHRY)

| System | Supercell | Composition | Sites |
|--------|-----------|-------------|-------|
| SC 2×2×2 | diagonal | [4,4] | 8 |
| SC 3×3×3 | diagonal | [13,14] | 27 |
| SC 4×4×4 | diagonal | [32,32] | 64 (if feasible) |
| FCC 2×2×2 | diagonal | [4,4] | 8 |
| FCC 3×3×3 | diagonal | [13,14] | 27 |

Both tools read CIF files. SHRY Python startup (~3.5s) noted in caption. Record SHRY import baseline.

### Group D — Multi-site enumeration (tesserae vs SHRY/enumlib where applicable)

| System | CIF | Supercell | Description |
|--------|-----|-----------|-------------|
| NaCl 2×2×2 | nacl_prim.cif | 2×2×2 | Binary on Na, Cl fixed |
| NaCl 3×3×3 | nacl_prim.cif | 3×3×3 | Binary on Na, Cl fixed |
| SrTiO3 2×2×2 | perovskite.cif | 2×2×2 | Binary on Sr, Ti+O fixed |

tesserae: `from-cif ... --multisite --composition ...`
Competitors: pending Phase 0 verification.

### Group E — Constrained/partial-occupancy enumeration

| System | CIF | Supercell | Description |
|--------|-----|-----------|-------------|
| Na/K+Cl 2×2×2 | partial_occ.cif | 2×2×2 | Na/K on cation, Cl fixed |
| Na/K+Cl 3×3×3 | partial_occ.cif | 3×3×3 | Na/K on cation, Cl fixed |

tesserae: `from-cif ... --constrained --composition ...`
Competitors: pending Phase 0 verification.

## Phase 3: Script Modifications (`bench_data/run_benchmarks.sh`)

Key changes:
1. **Count-agreement verification** — extract count from both tools, abort on mismatch
2. **Stddev reporting** — extract from hyperfine JSON alongside mean
3. **SHRY import baseline** — measure `python3 -c "import shry"` separately
4. **Adaptive hyperfine settings** — more warmup/runs for small indices, fewer for large
5. **Modular functions** — one function per benchmark group (A–E)
6. **Extended point groups** — keep existing OH_GENS/FCC_GENS for Groups A/B

## Phase 4: Results Format (`bench_data/benchmark_results.md`)

Structure:
1. **Environment** — CPU, RAM, OS, tool versions, date
2. **Methodology** — how timings measured, count verification, SHRY startup note, enumlib all-at-once note
3. **Table 1** — Binary at-index sweep (SC + FCC vs enumlib)
4. **Table 2** — Ternary at-index sweep (SC + FCC vs enumlib)
5. **Table 3** — Fixed supercell (vs SHRY)
6. **Table 4** — Multi-site enumeration
7. **Table 5** — Constrained enumeration
8. **Analysis** — speedup trends, T-canonicalization scaling explanation

Include ±stddev in timing columns.

## Phase 5: README Update

Add a brief "Performance" section after "Usage" with 3-4 highlight rows from Tables 1 and 3. Link to full `bench_data/benchmark_results.md`.

Note: `bench_data/` is gitignored in the private repo but should be included in the public repo. Verify this.

## Phase 6: Criterion Extensions (optional, lower priority)

- Extend `at_index` sweeps to idx 14/16 in `enumerate.rs`
- Add T-canonical n=128,256,512 for scaling plot data
- These produce data for paper figures but aren't publication tables

## Files to Modify

- `bench_data/run_benchmarks.sh` — rewrite with modular benchmark groups
- `bench_data/benchmark_results.md` — regenerated by script
- `bench_data/nacl_prim.cif` — new file
- `bench_data/perovskite.cif` — new file
- `bench_data/partial_occ.cif` — new file
- `README.md` — add Performance section
- `crates/tesserae-core/benches/enumerate.rs` — extend idx range (Phase 6)

## Verification

- Every comparative table row has count agreement between tools
- `cargo test` passes (no regressions)
- Criterion benchmarks compile and run: `cargo bench -p tesserae-core`
- Script runs end-to-end: `bash bench_data/run_benchmarks.sh`
- Results markdown renders correctly
