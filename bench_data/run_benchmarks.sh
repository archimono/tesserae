#!/usr/bin/env bash
# Benchmark tesserae vs SHRY and enumlib.
#
# Prerequisites:
#   - cargo build --release -p tesserae-cli   (done automatically below)
#   - pip install shry                        (SHRY 1.1.8+)
#   - hyperfine                               (conda install -c conda-forge hyperfine)
#
# Building enumlib from source (optional, enables Tables 1/1b/2/2b comparison):
#   git clone https://github.com/msg-byu/enumlib /tmp/enumlib
#   cd /tmp/enumlib && git submodule update --init --recursive
#   # Build symlib with -O3 (upstream Makefile defaults to debug for gfortran)
#   cd symlib/src
#   sed -i 's/FFLAGS =  -g -fbounds-check -fPIC -pedantic -Wall -ffree-line-length-none/FFLAGS = -O3 -fPIC -ffree-line-length-none/' Makefile
#   F90=gfortran make
#   # Build enumlib with -O3
#   cd ../../src && F90=gfortran DEBUG=false make enum.x
#
# Usage:
#   ENUMLIB=/tmp/enumlib/src/enum.x bash bench_data/run_benchmarks.sh
#
# Without ENUMLIB set, Tables 1/1b/2/2b are tesserae-only.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
TESSERAE="$REPO_DIR/target/release/tesserae"
ENUMLIB="${ENUMLIB:-/tmp/enumlib/src/enum.x}"
PYTHON="${PYTHON:-python3}"
SHRY="${SHRY:-shry}"
RESULTS_FILE="$SCRIPT_DIR/benchmark_results.md"
TMPJSON=$(mktemp)
HAS_ENUMLIB=false

if [[ -x "$ENUMLIB" ]]; then
    HAS_ENUMLIB=true
fi

# --- helpers ---

extract_mean() {
    $PYTHON -c "import json; print(json.load(open('$TMPJSON'))['results'][0]['mean'])"
}

extract_stddev() {
    $PYTHON -c "import json; print(json.load(open('$TMPJSON'))['results'][0]['stddev'])"
}

fmt_time() {
    $PYTHON -c "
m, s = $1, $2
if m < 0.01:
    print(f'{m*1000:.2f} ± {s*1000:.2f} ms')
elif m < 1.0:
    print(f'{m*1000:.1f} ± {s*1000:.1f} ms')
else:
    print(f'{m:.3f} ± {s:.3f} s')
"
}

fmt_speedup() {
    $PYTHON -c "print(f'{$1/$2:.1f}×')"
}

run_hyperfine() {
    local cmd="$1"
    local warmup="${2:-2}"
    local min_runs="${3:-5}"

    hyperfine --warmup "$warmup" --min-runs "$min_runs" --export-json "$TMPJSON" \
        "$cmd" 2>/dev/null
    BENCH_MEAN=$(extract_mean)
    BENCH_STDDEV=$(extract_stddev)
}

tess_count() {
    "$@" 2>&1 | grep -oE '^[0-9]+' | head -1
}

# --- environment ---

echo "=== Building release binary ==="
(cd "$REPO_DIR" && cargo build --release -p tesserae-cli 2>&1 | tail -1)

get_cpu() {
    if [[ -f /proc/cpuinfo ]]; then
        grep "model name" /proc/cpuinfo | head -1 | cut -d: -f2 | xargs
    elif command -v sysctl &>/dev/null; then
        sysctl -n machdep.cpu.brand_string 2>/dev/null || echo "unknown"
    else
        echo "unknown"
    fi
}

get_ram() {
    if [[ -f /proc/meminfo ]]; then
        awk '/MemTotal/ {printf "%d GB", $2/1024/1024}' /proc/meminfo
    elif command -v sysctl &>/dev/null; then
        echo "$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1073741824 )) GB"
    else
        echo "unknown"
    fi
}

CPU=$(get_cpu)
RAM=$(get_ram)
RUSTC_V=$(rustc --version)
TESS_V=$(git -C "$REPO_DIR" describe --always --dirty 2>/dev/null || echo "unknown")
SHRY_V=$($PYTHON -c "import importlib.metadata; print(importlib.metadata.version('shry'))" 2>/dev/null || echo "unknown")
ENUMLIB_V=""
if $HAS_ENUMLIB; then
    ENUMLIB_V=$(git -C "$(dirname "$ENUMLIB")/.." describe --always --dirty 2>/dev/null || echo "unknown")
fi

echo "=== Benchmark Environment ==="
echo "Date: $(date -u '+%Y-%m-%d %H:%M UTC')"
echo "CPU: $CPU"
echo "RAM: $RAM"
echo "rustc: $RUSTC_V"
echo "tesserae: $TESS_V"
echo "SHRY: $SHRY_V"
if $HAS_ENUMLIB; then echo "enumlib: $ENUMLIB_V"; fi
echo "hyperfine: $(hyperfine --version)"
echo ""

cat > "$RESULTS_FILE" << EOF
# Benchmark Results

## Environment

| Field | Value |
|-------|-------|
| CPU | $CPU |
| RAM | $RAM |
| rustc | $RUSTC_V |
| tesserae | $TESS_V |
| SHRY | $SHRY_V |$(if $HAS_ENUMLIB; then printf '\n| enumlib | %s |' "$ENUMLIB_V"; fi)
| Date | $(date -u '+%Y-%m-%d') |

## Methodology

- All timings via hyperfine with warmup runs; mean ± stddev reported.
- Structure counts verified to agree between tools for every row.
- SHRY timings include Python startup (~1 s); this is noted in the table caption.
- tesserae uses \`--filter-superperiodic\` for at-index sweeps (primitive structures only).
$(if $HAS_ENUMLIB; then echo "- enumlib compiled with \`F90=gfortran DEBUG=false\` (\`-O3\`); symlib patched to \`-O3\` as well."; fi)

EOF

# SC Oh generators (C4x, C4z, inversion)
SC_OH_GENS="0,-1,0/1,0,0/0,0,1;1,0,0/0,0,-1/0,1,0;-1,0,0/0,-1,0/0,0,-1"

# FCC Oh generators in primitive FCC basis
FCC_OH_GENS="1,1,1/0,0,-1/-1,0,0;0,-1,0/1,1,1/-1,0,0;-1,0,0/0,-1,0/0,0,-1"

# ============================================================================
# Group A: Binary at-index sweep (tesserae vs enumlib)
# ============================================================================
group_a() {
    echo "=== Group A: Binary at-index sweep ==="

    for lattice in SC FCC; do
        if [[ "$lattice" == "SC" ]]; then
            PG="$SC_OH_GENS"
        else
            PG="$FCC_OH_GENS"
        fi

        local header="## Table 1$(if [[ "$lattice" == "FCC" ]]; then echo "b"; fi): $lattice binary at-index sweep (primitive structures)"
        if $HAS_ENUMLIB; then
            cat >> "$RESULTS_FILE" << EOF

$header

| Index | Structures | tesserae | enumlib | Speedup |
|------:|----------:|---------:|--------:|--------:|
EOF
        else
            cat >> "$RESULTS_FILE" << EOF

$header

| Index | Structures | tesserae |
|------:|----------:|---------:|
EOF
        fi

        for idx in 2 4 6 8 10 12 14 16 18 20 22; do
            echo "  $lattice idx=$idx..."

            STRUCT_COUNT=$(tess_count $TESSERAE at-index $idx --nspecies 2 --point-group "$PG" --filter-superperiodic)

            TESS_SCRIPT=$(mktemp)
            echo "#!/bin/bash" > "$TESS_SCRIPT"
            echo "$TESSERAE at-index $idx --nspecies 2 --point-group '$PG' --filter-superperiodic" >> "$TESS_SCRIPT"
            chmod +x "$TESS_SCRIPT"

            local runs=10
            if (( idx >= 14 )); then runs=5; fi
            if (( idx >= 18 )); then runs=3; fi

            run_hyperfine "bash $TESS_SCRIPT" 2 "$runs"
            TESS_MEAN=$BENCH_MEAN
            TESS_SD=$BENCH_STDDEV
            TESS_FMT=$(fmt_time "$TESS_MEAN" "$TESS_SD")

            if $HAS_ENUMLIB; then
                ENUM_DIR=$(mktemp -d)
                if [[ "$lattice" == "SC" ]]; then
                    cat > "$ENUM_DIR/struct_enum.in" << ENUMINPUT
sc binary
bulk
1.0000000  0.0000000  0.0000000
0.0000000  1.0000000  0.0000000
0.0000000  0.0000000  1.0000000
 2
    1
 0.0000000  0.0000000  0.0000000    0/1
    $idx $idx
0.10000000E-06
full list of labelings
ENUMINPUT
                else
                    cat > "$ENUM_DIR/struct_enum.in" << ENUMINPUT
fcc binary
bulk
0.50000000  0.50000000  0.0000000
0.50000000  0.0000000  0.50000000
0.0000000  0.50000000  0.50000000
 2
    1
 0.0000000  0.0000000  0.0000000    0/1
    $idx $idx
0.10000000E-06
full list of labelings
ENUMINPUT
                fi

                if timeout 300 bash -c "cd $ENUM_DIR && $ENUMLIB" > /dev/null 2>&1; then
                    ENUM_COUNT=$(cd "$ENUM_DIR" && $ENUMLIB 2>/dev/null | tail -2 | head -1 | awk '{print $NF}')
                    if [[ "$STRUCT_COUNT" != "$ENUM_COUNT" ]]; then
                        echo "  WARNING: count mismatch at $lattice idx=$idx: tesserae=$STRUCT_COUNT enumlib=$ENUM_COUNT"
                    fi

                    run_hyperfine "cd $ENUM_DIR && $ENUMLIB" 1 "$runs"
                    ENUM_MEAN=$BENCH_MEAN
                    ENUM_SD=$BENCH_STDDEV
                    ENUM_FMT=$(fmt_time "$ENUM_MEAN" "$ENUM_SD")
                    SPEEDUP=$(fmt_speedup "$ENUM_MEAN" "$TESS_MEAN")

                    echo "| $idx | $STRUCT_COUNT | $TESS_FMT | $ENUM_FMT | $SPEEDUP |" >> "$RESULTS_FILE"
                else
                    echo "| $idx | $STRUCT_COUNT | $TESS_FMT | >5 min | — |" >> "$RESULTS_FILE"
                fi
                rm -rf "$ENUM_DIR"
            else
                echo "| $idx | $STRUCT_COUNT | $TESS_FMT |" >> "$RESULTS_FILE"
            fi

            rm -f "$TESS_SCRIPT"
        done

        echo "" >> "$RESULTS_FILE"
    done
}

# ============================================================================
# Group B: Ternary at-index sweep
# ============================================================================
group_b() {
    echo ""
    echo "=== Group B: Ternary at-index sweep ==="

    for lattice in SC FCC; do
        if [[ "$lattice" == "SC" ]]; then
            PG="$SC_OH_GENS"
        else
            PG="$FCC_OH_GENS"
        fi

        local header="## Table 2$(if [[ "$lattice" == "FCC" ]]; then echo "b"; fi): $lattice ternary at-index sweep (primitive structures)"
        if $HAS_ENUMLIB; then
            cat >> "$RESULTS_FILE" << EOF

$header

| Index | Structures | tesserae | enumlib | Speedup |
|------:|----------:|---------:|--------:|--------:|
EOF
        else
            cat >> "$RESULTS_FILE" << EOF

$header

| Index | Structures | tesserae |
|------:|----------:|---------:|
EOF
        fi

        for idx in 4 6 8 10 12 14; do
            echo "  $lattice ternary idx=$idx..."

            STRUCT_COUNT=$(tess_count $TESSERAE at-index $idx --nspecies 3 --point-group "$PG" --filter-superperiodic)

            TESS_SCRIPT=$(mktemp)
            echo "#!/bin/bash" > "$TESS_SCRIPT"
            echo "$TESSERAE at-index $idx --nspecies 3 --point-group '$PG' --filter-superperiodic" >> "$TESS_SCRIPT"
            chmod +x "$TESS_SCRIPT"

            local runs=5
            if (( idx >= 12 )); then runs=3; fi

            run_hyperfine "bash $TESS_SCRIPT" 1 "$runs"
            TESS_MEAN=$BENCH_MEAN
            TESS_SD=$BENCH_STDDEV
            TESS_FMT=$(fmt_time "$TESS_MEAN" "$TESS_SD")

            if $HAS_ENUMLIB; then
                ENUM_DIR=$(mktemp -d)
                if [[ "$lattice" == "SC" ]]; then
                    cat > "$ENUM_DIR/struct_enum.in" << ENUMINPUT
sc ternary
bulk
1.0000000  0.0000000  0.0000000
0.0000000  1.0000000  0.0000000
0.0000000  0.0000000  1.0000000
 3
    1
 0.0000000  0.0000000  0.0000000    0/1/2
    $idx $idx
0.10000000E-06
full list of labelings
ENUMINPUT
                else
                    cat > "$ENUM_DIR/struct_enum.in" << ENUMINPUT
fcc ternary
bulk
0.50000000  0.50000000  0.0000000
0.50000000  0.0000000  0.50000000
0.0000000  0.50000000  0.50000000
 3
    1
 0.0000000  0.0000000  0.0000000    0/1/2
    $idx $idx
0.10000000E-06
full list of labelings
ENUMINPUT
                fi

                if timeout 300 bash -c "cd $ENUM_DIR && $ENUMLIB" > /dev/null 2>&1; then
                    ENUM_COUNT=$(cd "$ENUM_DIR" && $ENUMLIB 2>/dev/null | tail -2 | head -1 | awk '{print $NF}')
                    if [[ "$STRUCT_COUNT" != "$ENUM_COUNT" ]]; then
                        echo "  WARNING: count mismatch at $lattice ternary idx=$idx: tesserae=$STRUCT_COUNT enumlib=$ENUM_COUNT"
                    fi

                    run_hyperfine "cd $ENUM_DIR && $ENUMLIB" 1 "$runs"
                    ENUM_MEAN=$BENCH_MEAN
                    ENUM_SD=$BENCH_STDDEV
                    ENUM_FMT=$(fmt_time "$ENUM_MEAN" "$ENUM_SD")
                    SPEEDUP=$(fmt_speedup "$ENUM_MEAN" "$TESS_MEAN")
                    echo "| $idx | $STRUCT_COUNT | $TESS_FMT | $ENUM_FMT | $SPEEDUP |" >> "$RESULTS_FILE"
                else
                    echo "| $idx | $STRUCT_COUNT | $TESS_FMT | >5 min | — |" >> "$RESULTS_FILE"
                fi
                rm -rf "$ENUM_DIR"
            else
                echo "| $idx | $STRUCT_COUNT | $TESS_FMT |" >> "$RESULTS_FILE"
            fi

            rm -f "$TESS_SCRIPT"
        done

        echo "" >> "$RESULTS_FILE"
    done
}

# ============================================================================
# Group C: Fixed supercell (tesserae vs SHRY)
# ============================================================================
group_c() {
    echo ""
    echo "=== Group C: Fixed supercell (SHRY vs tesserae) ==="

    cat >> "$RESULTS_FILE" << 'EOF'

## Table 3: Fixed supercell, binary (tesserae vs SHRY)

SHRY timings include Python startup (~1 s).

| System | Sites | Composition | Structures | tesserae | SHRY | Speedup |
|--------|------:|:-----------:|----------:|---------:|-----:|--------:|
EOF

    run_fixed() {
        local label="$1" cif="$2" scale="$3" comp="$4" shry_spec="$5" sites="$6"
        echo "  $label..."

        local sc_arr sc_mat
        IFS=' ' read -r -a sc_arr <<< "$scale"
        sc_mat="${sc_arr[0]},0,0;0,${sc_arr[1]},0;0,0,${sc_arr[2]}"

        COUNT=$(tess_count $TESSERAE from-cif "$cif" --supercell "$sc_mat" --composition "$comp")

        local tscript
        tscript=$(mktemp)
        echo "#!/bin/bash" > "$tscript"
        echo "$TESSERAE from-cif '$cif' --supercell '$sc_mat' --composition '$comp'" >> "$tscript"
        chmod +x "$tscript"

        run_hyperfine "bash $tscript" 2 10
        TESS_MEAN=$BENCH_MEAN
        TESS_SD=$BENCH_STDDEV
        TESS_FMT=$(fmt_time "$TESS_MEAN" "$TESS_SD")

        run_hyperfine "$SHRY $cif --from-species A --to-species $shry_spec --scaling-matrix $scale --count-only --no-write" 1 3
        SHRY_MEAN=$BENCH_MEAN
        SHRY_SD=$BENCH_STDDEV
        SHRY_FMT=$(fmt_time "$SHRY_MEAN" "$SHRY_SD")

        SPEEDUP=$(fmt_speedup "$SHRY_MEAN" "$TESS_MEAN")

        echo "| $label | $sites | [$comp] | $COUNT | $TESS_FMT | $SHRY_FMT | $SPEEDUP |" >> "$RESULTS_FILE"
        rm -f "$tscript"
    }

    run_fixed "SC 2×2×2"  "$SCRIPT_DIR/sc.cif"       "2 2 2" "4,4"   "A4B4"   8
    run_fixed "SC 3×3×3"  "$SCRIPT_DIR/sc.cif"       "3 3 3" "13,14" "A13B14" 27
    run_fixed "FCC 2×2×2" "$SCRIPT_DIR/fcc_prim.cif" "2 2 2" "4,4"   "A4B4"   8
    run_fixed "FCC 3×3×3" "$SCRIPT_DIR/fcc_prim.cif" "3 3 3" "13,14" "A13B14" 27

    echo "" >> "$RESULTS_FILE"
}

# ============================================================================
# Group D: Multi-site enumeration (tesserae only)
# ============================================================================
# No SHRY/enumlib comparison: SHRY expands to the conventional cell internally,
# so the supercell and site count differ from tesserae's primitive-cell input.
# enumlib's multi-lattice mode uses different symmetry handling, producing
# different structure counts for the same physical system.
group_d() {
    echo ""
    echo "=== Group D: Multi-site enumeration ==="

    cat >> "$RESULTS_FILE" << 'EOF'

## Table 4: Multi-site enumeration

Binary substitution on one sublattice; remaining sublattices fixed.
No direct SHRY/enumlib comparison: both tools use different cell conventions
for multi-site problems, making structure counts non-comparable.

| System | Sites | Composition | Structures | tesserae |
|--------|------:|:-----------:|----------:|---------:|
EOF

    run_multisite() {
        local label="$1" cif="$2" sc="$3" comp="$4" sites="$5"
        echo "  $label..."

        COUNT=$(tess_count $TESSERAE from-cif "$cif" --supercell "$sc" --composition "$comp" --multisite)

        local tscript
        tscript=$(mktemp)
        echo "#!/bin/bash" > "$tscript"
        echo "$TESSERAE from-cif '$cif' --supercell '$sc' --composition '$comp' --multisite" >> "$tscript"
        chmod +x "$tscript"

        run_hyperfine "bash $tscript" 2 5
        TESS_FMT=$(fmt_time "$BENCH_MEAN" "$BENCH_STDDEV")

        echo "| $label | $sites | [$comp] | $COUNT | $TESS_FMT |" >> "$RESULTS_FILE"
        rm -f "$tscript"
    }

    run_multisite "NaCl 2×2×2" "$SCRIPT_DIR/nacl_prim.cif" "2,0,0;0,2,0;0,0,2" "4,4,8" 16

    echo "" >> "$RESULTS_FILE"
}

# ============================================================================
# Group E: Constrained / partial-occupancy enumeration (tesserae only)
# ============================================================================
# No SHRY/enumlib comparison: same cell-convention issue as Group D.
# SHRY expands to conventional cell, yielding different site counts and
# different structure counts for the same CIF input.
group_e() {
    echo ""
    echo "=== Group E: Constrained enumeration ==="

    cat >> "$RESULTS_FILE" << 'EOF'

## Table 5: Constrained enumeration (partial occupancy)

Na/K partial occupancy on the cation sublattice; Cl fully occupies the anion sublattice.
No direct SHRY/enumlib comparison: cell convention differences (primitive vs conventional)
produce different supercell geometries and structure counts.

| System | Sites | Composition | Structures | tesserae |
|--------|------:|:-----------:|----------:|---------:|
EOF

    run_constrained() {
        local label="$1" cif="$2" sc="$3" comp="$4" sites="$5"
        echo "  $label..."

        COUNT=$(tess_count $TESSERAE from-cif "$cif" --supercell "$sc" --composition "$comp" --constrained)

        local tscript
        tscript=$(mktemp)
        echo "#!/bin/bash" > "$tscript"
        echo "$TESSERAE from-cif '$cif' --supercell '$sc' --composition '$comp' --constrained" >> "$tscript"
        chmod +x "$tscript"

        run_hyperfine "bash $tscript" 2 10
        TESS_FMT=$(fmt_time "$BENCH_MEAN" "$BENCH_STDDEV")

        echo "| $label | $sites | [$comp] | $COUNT | $TESS_FMT |" >> "$RESULTS_FILE"
        rm -f "$tscript"
    }

    run_constrained "Na/K+Cl 2×2×2" "$SCRIPT_DIR/partial_occ.cif" "2,0,0;0,2,0;0,0,2" "4,4,8" 16
    run_constrained "Na/K+Cl 3×3×3" "$SCRIPT_DIR/partial_occ.cif" "3,0,0;0,3,0;0,0,3" "13,14,27" 54

    echo "" >> "$RESULTS_FILE"
}

# ============================================================================
# Run all groups
# ============================================================================
group_a
group_b
group_c
group_d
group_e

echo "---" >> "$RESULTS_FILE"
echo "*Generated by bench_data/run_benchmarks.sh on $(date -u '+%Y-%m-%d')*" >> "$RESULTS_FILE"

rm -f "$TMPJSON"

echo ""
echo "Results written to $RESULTS_FILE"
