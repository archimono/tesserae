use std::time::Duration;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use tesserae_core::enumerate::{enumerate, enumerate_decomposed, enumerate_decomposed_cyclic};
use tesserae_core::hnf::hnf_3d;
use tesserae_core::perm::{Permutation, generate_group};
use tesserae_core::snf::smith_normal_form;
use tesserae_core::supercell::{
    enumerate_at_index, enumerate_at_index_primitive, enumerate_supercell,
    enumerate_supercell_primitive, generate_matrix_group, inequivalent_hnfs, transform_point_group,
};
use tesserae_core::t_canon::{
    cyclic_translations, product_translations, product_translations_2d, t_canonical,
    t_canonical_cyclic, t_canonical_product,
};

// ---------------------------------------------------------------------------
// Point group helpers
// ---------------------------------------------------------------------------

/// Oh point group in Cartesian (SC) basis: 48 signed permutation matrices.
fn oh_point_group() -> Vec<Vec<Vec<i64>>> {
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
                let sign: i64 = if signs & (1 << i) != 0 { -1 } else { 1 };
                m[i][perm[i]] = sign;
            }
            ops.push(m);
        }
    }
    ops
}

/// Oh point group in FCC primitive basis (48 elements).
/// Generators from benchmark_results.md / enumlib convention.
fn fcc_oh_point_group() -> Vec<Vec<Vec<i64>>> {
    let generators = vec![
        vec![vec![1, 1, 1], vec![0, 0, -1], vec![-1, 0, 0]],
        vec![vec![0, -1, 0], vec![1, 1, 1], vec![-1, 0, 0]],
        vec![vec![-1, 0, 0], vec![0, -1, 0], vec![0, 0, -1]],
    ];
    generate_matrix_group(&generators, 3)
}

fn all_binary_compositions(n: usize) -> Vec<Vec<usize>> {
    (0..=n).map(|k| vec![k, n - k]).collect()
}

fn all_ternary_compositions(n: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    for a in 0..=n {
        for b in 0..=(n - a) {
            result.push(vec![a, b, n - a - b]);
        }
    }
    result
}

// ---------------------------------------------------------------------------
// 1. Low-level: T-canonicalization
// ---------------------------------------------------------------------------

fn bench_t_canonical_cyclic(c: &mut Criterion) {
    let mut group = c.benchmark_group("t_canonical/cyclic_booth");
    for n in [8, 16, 32, 64] {
        let coloring: Vec<u32> = (0..n).map(|i| if i < n / 2 { 0 } else { 1 }).collect();
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| t_canonical_cyclic(&coloring));
        });
    }
    group.finish();
}

fn bench_t_canonical_brute_vs_booth(c: &mut Criterion) {
    let mut group = c.benchmark_group("t_canonical/cyclic_brute_vs_booth");
    for n in [8, 16, 32] {
        let translations = cyclic_translations(n);
        let coloring: Vec<u32> = (0..n).map(|i| if i < n / 2 { 0 } else { 1 }).collect();
        group.bench_with_input(BenchmarkId::new("brute", n), &n, |b, _| {
            b.iter(|| t_canonical(&coloring, &translations));
        });
        group.bench_with_input(BenchmarkId::new("booth", n), &n, |b, _| {
            b.iter(|| t_canonical_cyclic(&coloring));
        });
    }
    group.finish();
}

fn bench_t_canonical_product_2d(c: &mut Criterion) {
    let mut group = c.benchmark_group("t_canonical/product_2d");
    for (d0, d1) in [(2, 4), (3, 4), (4, 4), (4, 5)] {
        let n = d0 * d1;
        let dims = vec![d0, d1];
        let translations = product_translations_2d(d0, d1);
        let coloring: Vec<u32> = (0..n).map(|i| if i < n / 2 { 0 } else { 1 }).collect();
        let label = format!("{d0}x{d1}");
        group.bench_with_input(BenchmarkId::new("fast", &label), &label, |b, _| {
            b.iter(|| t_canonical_product(&coloring, &dims));
        });
        group.bench_with_input(BenchmarkId::new("brute", &label), &label, |b, _| {
            b.iter(|| t_canonical(&coloring, &translations));
        });
    }
    group.finish();
}

fn bench_t_canonical_product_3d(c: &mut Criterion) {
    let mut group = c.benchmark_group("t_canonical/product_3d");
    for dims in &[[2, 2, 2], [2, 2, 3], [2, 2, 4], [2, 3, 4], [3, 3, 3]] {
        let n: usize = dims.iter().product();
        let dim_vec = dims.to_vec();
        let translations = product_translations(&dim_vec);
        let coloring: Vec<u32> = (0..n).map(|i| if i < n / 2 { 0 } else { 1 }).collect();
        let label = format!("{}x{}x{}", dims[0], dims[1], dims[2]);
        group.bench_with_input(BenchmarkId::new("fast", &label), &label, |b, _| {
            b.iter(|| t_canonical_product(&coloring, &dim_vec));
        });
        group.bench_with_input(BenchmarkId::new("brute", &label), &label, |b, _| {
            b.iter(|| t_canonical(&coloring, &translations));
        });
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// 2. Mid-level: enumerate with explicit groups
// ---------------------------------------------------------------------------

fn bench_enumerate_cyclic(c: &mut Criterion) {
    let mut group = c.benchmark_group("enumerate/cyclic");
    for n in [8, 12, 16, 20] {
        let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
        let g = generate_group(n, &[rot]);
        let comp = vec![n / 2, n - n / 2];
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| enumerate(&comp, &g, &mut |_| {}));
        });
    }
    group.finish();
}

fn bench_enumerate_dihedral(c: &mut Criterion) {
    let mut group = c.benchmark_group("enumerate/dihedral");
    for n in [8, 12, 16, 20] {
        let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
        let refl = Permutation::new((0..n).map(|i| ((n - i) % n) as u32).collect());
        let g = generate_group(n, &[rot, refl]);
        let comp = vec![n / 2, n - n / 2];
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| enumerate(&comp, &g, &mut |_| {}));
        });
    }
    group.finish();
}

fn bench_enumerate_ternary_cyclic(c: &mut Criterion) {
    let mut group = c.benchmark_group("enumerate/ternary_cyclic");
    for n in [9, 12, 15] {
        let rot = Permutation::new((0..n).map(|i| ((i + 1) % n) as u32).collect());
        let g = generate_group(n, &[rot]);
        let comp = vec![n / 3, n / 3, n - 2 * (n / 3)];
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| enumerate(&comp, &g, &mut |_| {}));
        });
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// 3. Mid-level: decomposed enumeration
// ---------------------------------------------------------------------------

fn bench_decomposed_cyclic_dihedral(c: &mut Criterion) {
    let mut group = c.benchmark_group("decomposed/cyclic_dihedral");
    for n in [8, 12, 16, 20] {
        let refl = Permutation::new((0..n).map(|i| ((n - i) % n) as u32).collect());
        let point_group = vec![Permutation::identity(n), refl];
        let comp = vec![n / 2, n - n / 2];
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| enumerate_decomposed_cyclic(&comp, &point_group, &mut |_| {}));
        });
    }
    group.finish();
}

fn bench_decomposed_vs_full_c4_4x4(c: &mut Criterion) {
    let mut group = c.benchmark_group("decomposed/c4_4x4");
    let (d0, d1) = (4, 4);
    let n = d0 * d1;
    let translations = product_translations_2d(d0, d1);

    let rot90 = Permutation::new(
        (0..d0)
            .flat_map(|a| (0..d1).map(move |b| (b * d1 + (d0 - 1 - a)) as u32))
            .collect(),
    );
    let point_group = generate_group(n, &[rot90]);
    let full_group = {
        let mut g = Vec::with_capacity(translations.len() * point_group.len());
        for t in &translations {
            for p in &point_group {
                g.push(p.compose(t));
            }
        }
        g
    };
    let comp = vec![n / 2, n - n / 2];

    group.bench_function("full_group", |b| {
        b.iter(|| enumerate(&comp, &full_group, &mut |_| {}));
    });
    group.bench_function("decomposed", |b| {
        b.iter(|| enumerate_decomposed(&comp, &translations, &point_group, &mut |_| {}));
    });
    group.finish();
}

// ---------------------------------------------------------------------------
// 4. High-level: HNF reduction
// ---------------------------------------------------------------------------

fn bench_hnf_reduction(c: &mut Criterion) {
    let oh = oh_point_group();
    let fcc = fcc_oh_point_group();
    let mut group = c.benchmark_group("hnf_reduction");
    for idx in [4, 8, 12, 16] {
        let hnfs = hnf_3d(idx);
        group.bench_with_input(BenchmarkId::new("SC", idx), &idx, |b, _| {
            b.iter(|| inequivalent_hnfs(&hnfs, &oh));
        });
        group.bench_with_input(BenchmarkId::new("FCC", idx), &idx, |b, _| {
            b.iter(|| inequivalent_hnfs(&hnfs, &fcc));
        });
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// 5. High-level: setup cost (SNF + translations + transform_point_group)
// ---------------------------------------------------------------------------

fn bench_setup_cost(c: &mut Criterion) {
    let oh = oh_point_group();
    let mut group = c.benchmark_group("setup_cost");
    for idx in [4, 8, 12] {
        let hnfs = hnf_3d(idx);
        let inequiv = inequivalent_hnfs(&hnfs, &oh);
        group.bench_with_input(BenchmarkId::new("SC", idx), &idx, |b, _| {
            b.iter(|| {
                for s in &inequiv {
                    let snf = smith_normal_form(s);
                    let dims: Vec<usize> = snf.diagonal.iter().map(|&d| d as usize).collect();
                    let _translations = product_translations(&dims);
                    let _pg = transform_point_group(&oh, &snf);
                }
            });
        });
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// 6. End-to-end: at-index scaling sweep (binary, primitive/filtered)
// ---------------------------------------------------------------------------

fn bench_at_index_sc_binary_primitive(c: &mut Criterion) {
    let oh = oh_point_group();
    let mut group = c.benchmark_group("at_index/sc_binary_primitive");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    for idx in [2, 4, 6, 8, 10, 12] {
        let hnfs = hnf_3d(idx);
        let inequiv = inequivalent_hnfs(&hnfs, &oh);
        let compositions = all_binary_compositions(idx);
        group.bench_with_input(BenchmarkId::from_parameter(idx), &idx, |b, _| {
            b.iter(|| {
                for s in &inequiv {
                    for comp in &compositions {
                        enumerate_supercell_primitive(s, comp, &oh, &mut |_| {});
                    }
                }
            });
        });
    }
    group.finish();
}

fn bench_at_index_fcc_binary_primitive(c: &mut Criterion) {
    let fcc = fcc_oh_point_group();
    let mut group = c.benchmark_group("at_index/fcc_binary_primitive");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    for idx in [2, 4, 6, 8, 10, 12] {
        let hnfs = hnf_3d(idx);
        let inequiv = inequivalent_hnfs(&hnfs, &fcc);
        let compositions = all_binary_compositions(idx);
        group.bench_with_input(BenchmarkId::from_parameter(idx), &idx, |b, _| {
            b.iter(|| {
                for s in &inequiv {
                    for comp in &compositions {
                        enumerate_supercell_primitive(s, comp, &fcc, &mut |_| {});
                    }
                }
            });
        });
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// 7. End-to-end: at-index scaling sweep (ternary, primitive/filtered)
// ---------------------------------------------------------------------------

fn bench_at_index_sc_ternary_primitive(c: &mut Criterion) {
    let oh = oh_point_group();
    let mut group = c.benchmark_group("at_index/sc_ternary_primitive");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    for idx in [4, 6, 8] {
        let hnfs = hnf_3d(idx);
        let inequiv = inequivalent_hnfs(&hnfs, &oh);
        let compositions = all_ternary_compositions(idx);
        group.bench_with_input(BenchmarkId::from_parameter(idx), &idx, |b, _| {
            b.iter(|| {
                for s in &inequiv {
                    for comp in &compositions {
                        enumerate_supercell_primitive(s, comp, &oh, &mut |_| {});
                    }
                }
            });
        });
    }
    group.finish();
}

fn bench_at_index_fcc_ternary_primitive(c: &mut Criterion) {
    let fcc = fcc_oh_point_group();
    let mut group = c.benchmark_group("at_index/fcc_ternary_primitive");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    for idx in [4, 6, 8] {
        let hnfs = hnf_3d(idx);
        let inequiv = inequivalent_hnfs(&hnfs, &fcc);
        let compositions = all_ternary_compositions(idx);
        group.bench_with_input(BenchmarkId::from_parameter(idx), &idx, |b, _| {
            b.iter(|| {
                for s in &inequiv {
                    for comp in &compositions {
                        enumerate_supercell_primitive(s, comp, &fcc, &mut |_| {});
                    }
                }
            });
        });
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// 8. Decomposed vs full at 3D Oh (2×2×2 supercell)
// ---------------------------------------------------------------------------

fn bench_decomposed_vs_full_oh_2x2x2(c: &mut Criterion) {
    let oh = oh_point_group();
    let supercell = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 2]];
    let snf = smith_normal_form(&supercell);
    let dims: Vec<usize> = snf.diagonal.iter().map(|&d| d as usize).collect();
    let n: usize = dims.iter().product();

    let translations = product_translations(&dims);
    let point_group = transform_point_group(&oh, &snf);

    let full_group: Vec<Permutation> = {
        let mut g = Vec::with_capacity(translations.len() * point_group.len());
        for t in &translations {
            for p in &point_group {
                g.push(p.compose(t));
            }
        }
        g
    };

    let comp = vec![n / 2, n - n / 2];
    let mut group = c.benchmark_group("decomposed/oh_2x2x2");
    group.bench_function("full_group", |b| {
        b.iter(|| enumerate(&comp, &full_group, &mut |_| {}));
    });
    group.bench_function("decomposed", |b| {
        b.iter(|| enumerate_decomposed(&comp, &translations, &point_group, &mut |_| {}));
    });
    group.finish();
}

// ---------------------------------------------------------------------------
// 9. End-to-end: at-index sweep (non-primitive, all colorings)
// ---------------------------------------------------------------------------

fn bench_at_index_sc_binary(c: &mut Criterion) {
    let oh = oh_point_group();
    let mut group = c.benchmark_group("at_index/sc_binary");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    for idx in [2, 4, 6, 8, 10, 12] {
        let hnfs = hnf_3d(idx);
        let inequiv = inequivalent_hnfs(&hnfs, &oh);
        let compositions = all_binary_compositions(idx);
        group.bench_with_input(BenchmarkId::from_parameter(idx), &idx, |b, _| {
            b.iter(|| {
                for s in &inequiv {
                    for comp in &compositions {
                        enumerate_supercell(s, comp, &oh, &mut |_| {});
                    }
                }
            });
        });
    }
    group.finish();
}

fn bench_at_index_fcc_binary(c: &mut Criterion) {
    let fcc = fcc_oh_point_group();
    let mut group = c.benchmark_group("at_index/fcc_binary");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    for idx in [2, 4, 6, 8, 10, 12] {
        let hnfs = hnf_3d(idx);
        let inequiv = inequivalent_hnfs(&hnfs, &fcc);
        let compositions = all_binary_compositions(idx);
        group.bench_with_input(BenchmarkId::from_parameter(idx), &idx, |b, _| {
            b.iter(|| {
                for s in &inequiv {
                    for comp in &compositions {
                        enumerate_supercell(s, comp, &fcc, &mut |_| {});
                    }
                }
            });
        });
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// 10. Fixed supercell (non-primitive), matching SHRY comparison
// ---------------------------------------------------------------------------

fn bench_supercell_fixed(c: &mut Criterion) {
    let oh = oh_point_group();
    let fcc = fcc_oh_point_group();
    let mut group = c.benchmark_group("supercell/fixed");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));

    let sc_2x2x2 = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 2]];
    group.bench_function("SC_2x2x2_4_4", |b| {
        b.iter(|| enumerate_supercell(&sc_2x2x2, &[4, 4], &oh, &mut |_| {}));
    });

    let sc_3x3x3 = vec![vec![3, 0, 0], vec![0, 3, 0], vec![0, 0, 3]];
    group.bench_function("SC_3x3x3_13_14", |b| {
        b.iter(|| enumerate_supercell(&sc_3x3x3, &[13, 14], &oh, &mut |_| {}));
    });

    let fcc_2x2x2 = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 2]];
    group.bench_function("FCC_2x2x2_4_4", |b| {
        b.iter(|| enumerate_supercell(&fcc_2x2x2, &[4, 4], &fcc, &mut |_| {}));
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// 11. enumerate_at_index end-to-end (single call, black-box)
// ---------------------------------------------------------------------------

fn bench_at_index_e2e(c: &mut Criterion) {
    let oh = oh_point_group();
    let fcc = fcc_oh_point_group();
    let mut group = c.benchmark_group("at_index_e2e");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));

    group.bench_function("SC_idx8_primitive", |b| {
        b.iter(|| {
            let compositions = all_binary_compositions(8);
            for comp in &compositions {
                enumerate_at_index_primitive(&oh, 8, comp, &mut |_| {});
            }
        });
    });

    group.bench_function("FCC_idx8_primitive", |b| {
        b.iter(|| {
            let compositions = all_binary_compositions(8);
            for comp in &compositions {
                enumerate_at_index_primitive(&fcc, 8, comp, &mut |_| {});
            }
        });
    });

    group.bench_function("SC_idx8", |b| {
        b.iter(|| {
            let compositions = all_binary_compositions(8);
            for comp in &compositions {
                enumerate_at_index(&oh, 8, comp, &mut |_| {});
            }
        });
    });

    group.bench_function("FCC_idx8", |b| {
        b.iter(|| {
            let compositions = all_binary_compositions(8);
            for comp in &compositions {
                enumerate_at_index(&fcc, 8, comp, &mut |_| {});
            }
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    // T-canonicalization
    bench_t_canonical_cyclic,
    bench_t_canonical_brute_vs_booth,
    bench_t_canonical_product_2d,
    bench_t_canonical_product_3d,
    // Explicit-group enumeration
    bench_enumerate_cyclic,
    bench_enumerate_dihedral,
    bench_enumerate_ternary_cyclic,
    // Decomposed enumeration
    bench_decomposed_cyclic_dihedral,
    bench_decomposed_vs_full_c4_4x4,
    bench_decomposed_vs_full_oh_2x2x2,
    // HNF reduction
    bench_hnf_reduction,
    // Setup cost
    bench_setup_cost,
    // At-index primitive (filtered, matches enumlib comparison)
    bench_at_index_sc_binary_primitive,
    bench_at_index_fcc_binary_primitive,
    bench_at_index_sc_ternary_primitive,
    bench_at_index_fcc_ternary_primitive,
    // At-index non-primitive (all colorings)
    bench_at_index_sc_binary,
    bench_at_index_fcc_binary,
    // Fixed supercell (matches SHRY comparison)
    bench_supercell_fixed,
    // End-to-end black-box
    bench_at_index_e2e,
);
criterion_main!(benches);
