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
use csv::ReaderBuilder;
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

pub type ScoreMap = HashMap<String, f64>;

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
    #[serde(rename = "contig_name")] name: String,
    #[serde(rename = "score", default)]     score: Option<f64>,
    #[serde(rename = "p_seq_true", default)] p_seq_true: Option<f64>,
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
                let expr_score = salmon_scores.get(raw_id).copied().unwrap_or(0.0);
                batch.insert(format!("{}__{}", prefix, raw_id), 0.4 * expr_score + 0.6 * cov_score);
            }
            for (raw_id, &expr_score) in &salmon_scores {
                let prefixed = format!("{}__{}", prefix, raw_id);
                if !batch.contains_key(&prefixed) {
                    batch.insert(prefixed, 0.4 * expr_score);
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
        let mut rdr = ReaderBuilder::new().has_headers(true)
            .from_reader(BufReader::new(file));
        for result in rdr.deserialize::<TransrateRow>() {
            match result {
                Ok(row) => {
                    let score = row.score.or(row.p_seq_true).unwrap_or(0.0);
                    map.insert(row.name, score);
                }
                Err(e) => warn!("  Skipping malformed CSV row: {}", e),
            }
        }
    }
    info!("  Loaded {} scores from CSV", map.len());
    Ok(map)
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
        map.insert(rec.name.clone(), score.clamp(0.0, 1.0));
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
        let tsv = "Name	Length	EffectiveLength	TPM	NumReads
                   t1	500	450.0	1000.0	200.0
                   t2	300	250.0	100.0	20.0
                   t3	200	150.0	0.0	0.0
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
    fn test_load_scores_from_csv() {
        let csv = "contig_name,score,length
t1,0.85,500
t2,0.40,300
";
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), csv).unwrap();
        let scores = load_scores_from_csv(&[tmp.path().to_path_buf()]).unwrap();
        assert_eq!(scores.len(), 2);
        assert!((scores["t1"] - 0.85).abs() < 1e-6);
    }
}
