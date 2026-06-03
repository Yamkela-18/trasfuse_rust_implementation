// src/filter.rs
//
// Contig filtering — replaces Transfuse#filter in transfuse.rb.
//
// Two public functions:
//   filter_assemblies()    -> keep assemblies with >= 1 contig above threshold
//   write_filtered_fasta() -> write subset of a FASTA above min_score

use anyhow::{Context, Result};
use log::{debug, info};
use std::path::{Path, PathBuf};

use crate::fasta::{self, FastaRecord};
use crate::score::ScoreMap;

// ── Public API ────────────────────────────────────────────────────────────────

/// Return the subset of assembly_files that contain at least one contig
/// whose score exceeds min_score.
/// Assemblies with no scored contigs are kept (unknown prefix scheme).
pub fn filter_assemblies(
    assembly_files: &[PathBuf],
    scores: &ScoreMap,
    min_score: f64,
) -> Result<Vec<PathBuf>> {
    if min_score <= 0.0 {
        return Ok(assembly_files.to_vec());
    }

    let mut kept = Vec::new();
    for path in assembly_files {
        let prefix = path
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "asm".into());

        let has_passing = scores.iter().any(|(id, &score)| {
            id.starts_with(&format!("{prefix}__")) && score >= min_score
        });

        if has_passing {
            kept.push(path.clone());
            debug!("  Keeping assembly {:?}", path);
        } else {
            info!("  Dropping assembly {:?} (no contigs above score threshold {})", path, min_score);
        }
    }
    Ok(kept)
}

/// Read input_fasta, keep only contigs with score >= min_score, write to output.
/// Unknown contigs (not in scores map) are kept.
pub fn write_filtered_fasta(
    input_fasta: &Path,
    scores: &ScoreMap,
    min_score: f64,
    output: &Path,
) -> Result<PathBuf> {
    let records = fasta::load_fasta_ordered(input_fasta)?;

    let kept: Vec<FastaRecord> = records
        .into_iter()
        .filter(|rec| {
            let score = scores.get(rec.id()).copied().unwrap_or(1.0);
            if score < min_score {
                debug!("  Dropping low-score contig {} (score={:.3})", rec.id(), score);
                false
            } else {
                true
            }
        })
        .collect();

    info!("  Final filter: kept {}/{} contigs", kept.len(), scores.len());

    fasta::write_fasta_file(output, &kept)
        .with_context(|| format!("Failed to write filtered FASTA to {:?}", output))?;

    Ok(output.to_path_buf())
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_filter_assemblies_keeps_all_when_threshold_zero() {
        let dir = tempdir().unwrap();
        let a = dir.path().join("a.fa");
        let b = dir.path().join("b.fa");
        std::fs::write(&a, ">s1
ACGT
").unwrap();
        std::fs::write(&b, ">s1
ACGT
").unwrap();
        let scores: ScoreMap = [("a__s1".into(), 0.1), ("b__s1".into(), 0.9)]
            .into_iter().collect();
        let result = filter_assemblies(&[a.clone(), b.clone()], &scores, 0.0).unwrap();
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_filter_assemblies_drops_below_threshold() {
        let dir = tempdir().unwrap();
        let a = dir.path().join("a.fa");
        let b = dir.path().join("b.fa");
        std::fs::write(&a, ">s1
ACGT
").unwrap();
        std::fs::write(&b, ">s1
ACGT
").unwrap();
        let scores: ScoreMap = [("a__s1".into(), 0.1), ("b__s1".into(), 0.9)]
            .into_iter().collect();
        let result = filter_assemblies(&[a, b.clone()], &scores, 0.5).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], b);
    }

    #[test]
    fn test_write_filtered_fasta() {
        let dir = tempdir().unwrap();
        let input  = dir.path().join("input.fa");
        let output = dir.path().join("output.fa");
        std::fs::write(&input, ">good
AAAA
>bad
CCCC
").unwrap();
        let scores: ScoreMap = [("good".into(), 0.9), ("bad".into(), 0.1)]
            .into_iter().collect();
        write_filtered_fasta(&input, &scores, 0.5, &output).unwrap();
        let records = fasta::load_fasta_ordered(&output).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id(), "good");
    }
}
