// src/score.rs
//
// Contig quality scoring — replaces the transrate Ruby gem.
//
// Two-phase scoring pipeline for each assembly FASTA:
//   Phase A — Expression (Salmon quasi-mapping)
//     Builds a Salmon index, quantifies reads, extracts TPM per contig.
//   Phase B — Alignment coverage (minimap2 + samtools)
//     Aligns reads with minimap2, computes per-contig coverage stats.
//   Final score = 0.4 * expression_score + 0.6 * coverage_score

use anyhow::{bail, Context, Result};
use csv::{ReaderBuilder, WriterBuilder};
use log::{debug, info, warn};
use rayon::prelude::*;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{Arc, Mutex};

use crate::aligner;
use crate::bam;

#[derive(Debug, Clone)]
pub struct ContigScore {
    pub score: f64,
    pub p_good: f64,
    pub p_bases_covered: f64,
    pub coverage: f64,
}

pub type ScoreMap = HashMap<String, ContigScore>;

#[derive(Debug, Deserialize)]
pub struct SalmonRecord {
    #[serde(rename = "Name")]      pub name: String,
    #[serde(rename = "Length")]    pub length: u64,
    #[serde(rename = "EffectiveLength")] pub effective_length: f64,
    #[serde(rename = "TPM")]       pub tpm: f64,
    #[serde(rename = "NumReads")]  pub num_reads: f64,
}

#[derive(Debug, Deserialize)]
struct TransrateRow {
    #[serde(rename = "contig_name")]
    name: String,

    #[serde(rename = "score", default)]
    score: Option<f64>,

    #[serde(rename = "p_good", default)]
    p_good: Option<f64>,

    #[serde(rename = "p_bases_covered", default)]
    p_bases_covered: Option<f64>,

    #[serde(rename = "coverage", default)]
    coverage: Option<f64>,
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Score all contigs using Salmon + minimap2/samtools in parallel via Rayon.
pub fn score_assemblies(
    assembly_files: &[PathBuf],
    left: &Path,
    right: &Path,
    threads: usize,
    verbose: bool,
) -> Result<ScoreMap> {
    let asm_threads = ((threads as f64 / assembly_files.len() as f64).ceil() as usize).max(1);
    info!("Scoring {} assemblies using {} thread(s) each", assembly_files.len(), asm_threads);

    let score_map: Arc<Mutex<ScoreMap>> = Arc::new(Mutex::new(HashMap::new()));

    let results: Vec<Result<()>> = assembly_files
        .par_iter()
        .map(|asm_path| {
            let prefix = asm_path
                .file_stem()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| "asm".into());

            info!("  Scoring assembly: {:?} (prefix={})", asm_path, prefix);

            let salmon_scores = run_salmon(asm_path, left, right, asm_threads, verbose)
                .unwrap_or_else(|e| {
                    warn!("  Salmon failed for {:?}: {}. Using zero expression scores.", asm_path, e);
                    HashMap::new()
                });

            let cov_scores = run_alignment_coverage(asm_path, left, right, asm_threads, verbose)
                .unwrap_or_else(|e| {
                    warn!("  Alignment failed for {:?}: {}. Using zero coverage scores.", asm_path, e);
                    HashMap::new()
                });

            // Combine: 40% expression, 60% coverage
            let mut batch: ScoreMap = HashMap::new();
            for (raw_id, &cov_score) in &cov_scores {
    let expr_score = salmon_scores
        .get(raw_id)
        .copied()
        .unwrap_or(0.0);

    let final_score =
        0.4 * expr_score + 0.6 * cov_score;

    batch.insert(
        format!("{}__{}", prefix, raw_id),
        ContigScore {
            score: final_score,

            // proxy for Transrate p_good
            p_good: expr_score,

            // proxy for p_bases_covered
            p_bases_covered: cov_score,

            // raw coverage estimate
            coverage: cov_score * 100.0,
        },
    );
}
            for (raw_id, &expr_score) in &salmon_scores {
                let prefixed = format!("{}__{}", prefix, raw_id);
                if !batch.contains_key(&prefixed) {
                    batch.insert(
    prefixed,
    ContigScore {
        score: 0.4 * expr_score,
        p_good: expr_score,
        p_bases_covered: 0.0,
        coverage: 0.0,
    },
);
                }
            }

            score_map.lock().unwrap().extend(batch);
            Ok(())
        })
        .collect();

    for r in results { r?; }

    Ok(Arc::try_unwrap(score_map).expect("Arc still has multiple owners")
        .into_inner().unwrap())
}

pub fn load_scores_from_csv(csv_files: &[PathBuf]) -> Result<ScoreMap> {
    let mut map = ScoreMap::new();

    for path in csv_files {
        info!("  Loading scores from {:?}", path);

        let file = File::open(path)
            .with_context(|| format!("Cannot open CSV {:?}", path))?;

        let mut rdr = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_reader(BufReader::new(file));

        for result in rdr.deserialize::<TransrateRow>() {
            match result {
                Ok(row) => {
                    map.insert(
                        row.name,
                        ContigScore {
                            score: row.score.unwrap_or(0.0),
                            p_good: row.p_good.unwrap_or(0.0),
                            p_bases_covered: row
                                .p_bases_covered
                                .unwrap_or(0.0),
                            coverage: row.coverage.unwrap_or(0.0),
                        },
                    );
                }
                Err(e) => warn!("Skipping malformed score row: {}", e),
            }
        }
    }

    info!("Loaded {} scores from file", map.len());

    Ok(map)
}

// ── CSV writer ────────────────────────────────────────────────────────────────

/// Serialise a ScoreMap to a sorted CSV file alongside output_path.
/// Filename: <output_stem>_scores.csv  (e.g. merged.fa -> merged_scores.csv)
/// Columns: contig_name, score  (sorted by score descending)
/// Compatible with load_scores_from_csv() for re-use in subsequent runs.
pub fn write_scores_csv(
    scores: &ScoreMap,
    output_path: &Path,
) -> Result<PathBuf> {
    let parent =
        output_path.parent().unwrap_or(Path::new("."));

    let stem = output_path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "scores".into());

    let csv_path =
        parent.join(format!("{stem}_scores.csv"));

    let mut rows: Vec<(&String, &ContigScore)> =
        scores.iter().collect();

    rows.sort_by(|a, b| {
        b.1.score
            .partial_cmp(&a.1.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut wtr = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(&csv_path)
        .with_context(|| {
            format!(
                "Cannot create scores file {:?}",
                csv_path
            )
        })?;

    wtr.write_record([
        "contig_name",
        "score",
        "p_good",
        "p_bases_covered",
        "coverage",
    ])?;

    for (name, metrics) in &rows {
        wtr.write_record([
            name.as_str(),
            &format!("{:.6}", metrics.score),
            &format!("{:.6}", metrics.p_good),
            &format!("{:.6}",
                metrics.p_bases_covered),
            &format!("{:.6}",
                metrics.coverage),
        ])?;
    }

    wtr.flush()?;

    info!(
        "Scores written to {:?} ({} contigs)",
        csv_path,
        rows.len()
    );

    Ok(csv_path)
}

// ── Phase A — Salmon ──────────────────────────────────────────────────────────

fn run_salmon(
    fasta: &Path, left: &Path, right: &Path, threads: usize, verbose: bool,
) -> Result<HashMap<String, f64>> {
    let stem   = fasta.file_stem().map(|s| s.to_string_lossy().into_owned())
                      .unwrap_or_else(|| "ref".into());
    let parent = fasta.parent().unwrap_or(Path::new("."));
    let index_dir = parent.join(format!("{stem}_salmon_index"));
    let quant_dir = parent.join(format!("{stem}_salmon_quant"));

    if !index_dir.exists() {
        debug!("  Building Salmon index for {:?}", fasta);
        let status = Command::new("salmon")
            .args(["index", "-t"]).arg(fasta)
            .args(["-i"]).arg(&index_dir)
            .args(["--threads", &threads.to_string()])
            .stderr(if verbose { Stdio::inherit() } else { Stdio::null() })
            .status().context("Failed to run salmon index")?;
        if !status.success() { bail!("salmon index failed for {:?}", fasta); }
    }

    if !quant_dir.exists() {
        debug!("  Running Salmon quant for {:?}", fasta);
        let status = Command::new("salmon")
            .args(["quant", "-i"]).arg(&index_dir)
            .args(["-l", "A", "-1"]).arg(left)
            .args(["-2"]).arg(right)
            .args(["-o"]).arg(&quant_dir)
            .args(["-p", &threads.to_string()])
            .args(["--validateMappings"])
            .stderr(if verbose { Stdio::inherit() } else { Stdio::null() })
            .status().context("Failed to run salmon quant")?;
        if !status.success() { bail!("salmon quant failed for {:?}", fasta); }
    }

    let quant_sf = quant_dir.join("quant.sf");
    if !quant_sf.exists() { bail!("quant.sf not found at {:?}", quant_sf); }
    parse_salmon_quant(&quant_sf)
}

fn parse_salmon_quant(quant_sf: &Path) -> Result<HashMap<String, f64>> {
    let file = File::open(quant_sf)
        .with_context(|| format!("Cannot open {:?}", quant_sf))?;
    let mut rdr = ReaderBuilder::new().delimiter(b'\t').has_headers(true)
        .from_reader(BufReader::new(file));

    let mut records: Vec<SalmonRecord> = Vec::new();
    for r in rdr.deserialize() {
        records.push(r.with_context(|| format!("Failed to parse {:?}", quant_sf))?);
    }

    let max_tpm = records.iter().map(|r| r.tpm).fold(0.0f64, f64::max);
    let log_max = (max_tpm + 1.0).log10().max(1e-9);

    let mut map = HashMap::new();
    for rec in &records {
        let score = (rec.tpm + 1.0).log10() / log_max;
        map.insert(rec.name.trim().to_string(), score.clamp(0.0, 1.0));
    }
    Ok(map)
}

// ── Phase B — minimap2 + samtools coverage ────────────────────────────────────

fn run_alignment_coverage(
    fasta: &Path, left: &Path, right: &Path, threads: usize, verbose: bool,
) -> Result<HashMap<String, f64>> {
    let bam_path  = aligner::align_reads(fasta, left, right, threads, verbose)?;
    let cov_stats = bam::compute_coverage(&bam_path)?;
    Ok(cov_stats.into_iter().map(|(id, stats)| (id, stats.score())).collect())
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_salmon_quant() {
        let tsv = "\
Name\tLength\tEffectiveLength\tTPM\tNumReads
t1\t500\t450.0\t1000.0\t200.0
t2\t300\t250.0\t100.0\t20.0
t3\t200\t150.0\t0.0\t0.0
";
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), tsv).unwrap();
        let scores = parse_salmon_quant(tmp.path()).unwrap();
        assert_eq!(scores.len(), 3);
        assert!((scores["t1"] - 1.0).abs() < 1e-6);
        assert!(scores["t3"] < 0.01);
        assert!(scores["t2"] > 0.0 && scores["t2"] < 1.0);
    }

    #[test]
    fn test_write_scores_csv_creates_file() {
        use tempfile::tempdir;
        let dir = tempdir().unwrap();
        let output = dir.path().join("merged.fa");
        let scores: ScoreMap = [
    (
        "k31__seq1".to_string(),
        ContigScore {
            score: 0.921,
            p_good: 0.921,
            p_bases_covered: 0.921,
            coverage: 2.0,
        },
    ),
    (
        "k41__seq1".to_string(),
        ContigScore {
            score: 0.654,
            p_good: 0.654,
            p_bases_covered: 0.654,
            coverage: 2.0,
        },
    ),
    (
        "k31__seq2".to_string(),
        ContigScore {
            score: 0.123,
            p_good: 0.123,
            p_bases_covered: 0.123,
            coverage: 2.0,
        },
    ),
]
.into_iter()
.collect();
        let csv_path = write_scores_csv(&scores, &output).unwrap();
        assert_eq!(csv_path.file_name().unwrap(), "merged_scores.csv");
        assert!(csv_path.exists());
        // Must round-trip through load_scores_from_csv
        let reloaded = load_scores_from_csv(&[csv_path]).unwrap();
        assert_eq!(reloaded.len(), 3);
        assert!((reloaded["k31__seq1"].score - 0.921).abs() < 1e-5);
    }

    #[test]
    fn test_write_scores_csv_sorted_descending() {
        use tempfile::tempdir;
        let dir = tempdir().unwrap();
        let output = dir.path().join("out.fa");
        let scores: ScoreMap = [
    (
        "a".to_string(),
        ContigScore {
            score: 0.3,
            p_good: 0.3,
            p_bases_covered: 0.3,
            coverage: 2.0,
        },
    ),
    (
        "b".to_string(),
        ContigScore {
            score: 0.9,
            p_good: 0.9,
            p_bases_covered: 0.9,
            coverage: 2.0,
        },
    ),
    (
        "c".to_string(),
        ContigScore {
            score: 0.6,
            p_good: 0.6,
            p_bases_covered: 0.6,
            coverage: 2.0,
        },
    ),
]
.into_iter()
.collect();
        let csv_path = write_scores_csv(&scores, &output).unwrap();
        let content = std::fs::read_to_string(&csv_path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert!(lines[1].contains("b"), "first data row should be highest score 'b'");
        assert!(lines[3].contains("a"), "last row should be lowest score 'a'");
    }

#[test]
fn test_load_scores_from_csv() {
    let csv = "contig_name\tscore\tp_good\tp_bases_covered\tcoverage
t1\t0.85\t0.80\t0.90\t2.0
t2\t0.40\t0.35\t0.50\t1.2
";

    let tmp = tempfile::NamedTempFile::new().unwrap();

    std::fs::write(tmp.path(), csv).unwrap();

    let scores =
        load_scores_from_csv(
            &[tmp.path().to_path_buf()]
        )
        .unwrap();

    assert_eq!(scores.len(), 2);

    assert!(
        (scores["t1"].score - 0.85)
            .abs()
            < 1e-6
    );

    assert!(
        (scores["t1"].coverage - 2.0)
            .abs()
            < 1e-6
    );
}
}
