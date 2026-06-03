#!/usr/bin/env bash
# fix_all.sh
# Run from ~/transfuse/:  bash fix_all.sh
# Fixes all 4 failing tests and all 8 warnings in one pass.

set -e
cd "$(dirname "$0")" 2>/dev/null || true
[ -f src/main.rs ] || { echo "Run from the transfuse project root (where Cargo.toml is)"; exit 1; }

python3 << 'PYEOF'
import re, os

def fix(path, replacements):
    with open(path) as f:
        src = f.read()
    original = src
    for old, new in replacements:
        if old in src:
            src = src.replace(old, new)
            print(f"  fixed: {os.path.basename(path)}: {repr(old[:60])}")
        else:
            print(f"  skip : {os.path.basename(path)}: not found: {repr(old[:60])}")
    if src != original:
        with open(path, 'w') as f:
            f.write(src)

# ── src/main.rs ───────────────────────────────────────────────────────────────
fix("src/main.rs", [
    ("use log::{info, warn};",           "use log::info;"),
    ("use std::path::{Path, PathBuf};",  "use std::path::PathBuf;"),
    ("use deps::BinaryDep;\n",           ""),
])

# ── src/cluster.rs ────────────────────────────────────────────────────────────
fix("src/cluster.rs", [
    ("use std::process::{Command, Stdio};", "use std::process::Command;"),
])

# ── src/consensus.rs ─────────────────────────────────────────────────────────
fix("src/consensus.rs", [
    ("use anyhow::{bail, Context, Result};", "use anyhow::{Context, Result};"),
])

# ── src/fasta.rs ──────────────────────────────────────────────────────────────
fix("src/fasta.rs", [
    ("use anyhow::{bail, Context, Result};", "use anyhow::{Context, Result};"),
    ("use std::fs::{File, OpenOptions};",    "use std::fs::File;"),
])

# ── src/filter.rs ─────────────────────────────────────────────────────────────
fix("src/filter.rs", [
    ("use std::collections::HashMap;\n", ""),
])

# ── Failing test fixes ────────────────────────────────────────────────────────

# bam.rs: replace any 4-space sequence pretending to be \t in split calls
# The split character must be a real tab escape \t
import re

for fname in ["src/bam.rs", "src/cluster.rs", "src/score.rs"]:
    with open(fname) as f:
        src = f.read()
    original = src
    # Fix split('    ') -> split('\t')  (4 spaces inside single quotes)
    src = re.sub(r"split\('    '\)",    r"split('\\t')",     src)
    src = re.sub(r"split\('	'\)",     r"split('\\t')",     src)
    # Fix delimiter(b'    ') -> delimiter(b'\t')
    src = re.sub(r"delimiter\(b'    '\)", r"delimiter(b'\\t')", src)
    src = re.sub(r"delimiter\(b'	'\)",  r"delimiter(b'\\t')", src)
    if src != original:
        with open(fname, 'w') as f:
            f.write(src)
        print(f"  fixed: {os.path.basename(fname)}: tab escape in split/delimiter")
    else:
        print(f"  skip : {os.path.basename(fname)}: no tab escape issue found")

# score.rs: fix scores["t1"] key lookup — the test must use .get() not index
# Root cause: parse_salmon_quant returns keys from CSV Name column
# If tabs are spaces, CSV reader sees whole line as one field -> no "t1" key
# This is fixed by the delimiter fix above; verify by checking the test
with open("src/score.rs") as f:
    score_src = f.read()

# Ensure test uses .get() with expect rather than direct index to get better errors
old_assert = '''        assert!((scores["t1"] - 1.0).abs() < 1e-6);
        // t3 has 0 TPM so score should be ~0
        assert!(scores["t3"] < 0.01);
        // t2 should be between 0 and 1
        assert!(scores["t2"] > 0.0 && scores["t2"] < 1.0);'''

new_assert = '''        let s1 = *scores.get("t1").expect("t1 missing from scores");
        let s2 = *scores.get("t2").expect("t2 missing from scores");
        let s3 = *scores.get("t3").expect("t3 missing from scores");
        assert!((s1 - 1.0).abs() < 1e-6, "t1 score should be 1.0, got {}", s1);
        assert!(s3 < 0.01,               "t3 score should be ~0, got {}",   s3);
        assert!(s2 > 0.0 && s2 < 1.0,   "t2 score should be 0<x<1, got {}", s2);'''

if old_assert in score_src:
    score_src = score_src.replace(old_assert, new_assert)
    with open("src/score.rs", 'w') as f:
        f.write(score_src)
    print("  fixed: score.rs: improved test assertions with .get().expect()")

# bam.rs: fix test key lookup
with open("src/bam.rs") as f:
    bam_src = f.read()

old_bam = '''        let s = &stats["contig1"];'''
new_bam = '''        let s = stats.get("contig1").expect("contig1 missing — check tab escapes in test TSV");'''

if old_bam in bam_src:
    bam_src = bam_src.replace(old_bam, new_bam)
    with open("src/bam.rs", 'w') as f:
        f.write(bam_src)
    print("  fixed: bam.rs: improved test assertion with .get().expect()")

print("\nAll fixes applied.")
PYEOF

echo ""
echo "Running cargo test..."
cargo test 2>&1
