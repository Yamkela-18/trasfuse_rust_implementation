// src/bam.rs
//
// BAM coverage computation — replaces the bam-read C binary from transrate-tools.
//
// Uses samtools coverage (available since samtools 1.12) to extract per-contig
// alignment statistics from a sorted BAM file.
//
// samtools coverage output columns (TSV):
//   #rname  startpos  endpos  numreads  covbases  coverage  meandepth  meanbaseq  meanmapq
//
// ContigStats::score() approximates the Transrate contig score:
//   0.5 * coverage_fraction + 0.3 * normalised_depth + 0.2 * mapq_fraction

use anyhow::{bail, Context, Result};
use log::debug;
use std::collections::HashMap;
use std::path::Path;
use std::process::Command;

// ── Public types ──────────────────────────────────────────────────────────────

/// Per-contig alignment statistics extracted from the BAM file.
#[derive(Debug, Clone)]
pub struct ContigStats {
    pub num_reads: u64,
    pub covered_bases: u64,
    pub coverage: f64,         // fraction 0.0-1.0
    pub mean_depth: f64,
    pub mean_base_quality: f64,
    pub mean_map_quality: f64,
    pub length: u64,
}

impl ContigStats {
    /// Composite score in [0.0, 1.0] approximating the Transrate contig score.
    /// Weights: 0.5 x coverage + 0.3 x depth (capped at 10x) + 0.2 x mapq (capped at 30)
    pub fn score(&self) -> f64 {
        if self.length == 0 { return 0.0; }
        let cov_frac    = self.coverage.clamp(0.0, 1.0);
        let depth_score = (self.mean_depth / 10.0).min(1.0);
        let mapq_score  = (self.mean_map_quality / 30.0).min(1.0);
        0.5 * cov_frac + 0.3 * depth_score + 0.2 * mapq_score
    }
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Run samtools coverage on bam_path and return per-contig stats.
pub fn compute_coverage(bam_path: &Path) -> Result<HashMap<String, ContigStats>> {
    debug!("  Running samtools coverage on {:?}", bam_path);

    let output = Command::new("samtools")
        .args(["coverage", "-d", "0"])  // -d 0: no depth cap
        .arg(bam_path)
        .output()
        .context("Failed to run samtools coverage")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        bail!("samtools coverage failed: {}", stderr);
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    parse_coverage_output(&stdout)
}

// ── Parser ────────────────────────────────────────────────────────────────────

fn parse_coverage_output(text: &str) -> Result<HashMap<String, ContigStats>> {
    let mut stats = HashMap::new();

    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 { continue; }

        let name       = cols[0].to_string();
        let start: u64 = cols[1].parse().unwrap_or(0);
        let end: u64   = cols[2].parse().unwrap_or(0);
        let length     = end.saturating_sub(start).max(1);
        let num_reads  = cols[3].parse().unwrap_or(0u64);
        let cov_bases  = cols[4].parse().unwrap_or(0u64);
        let coverage   = cols[5].parse().unwrap_or(0.0f64) / 100.0; // samtools gives %
        let mean_depth = cols[6].parse().unwrap_or(0.0f64);
        let mean_baseq = cols[7].parse().unwrap_or(0.0f64);
        let mean_mapq  = cols[8].parse().unwrap_or(0.0f64);

        stats.insert(name, ContigStats {
            num_reads, covered_bases: cov_bases, coverage,
            mean_depth, mean_base_quality: mean_baseq,
            mean_map_quality: mean_mapq, length,
        });
    }
    Ok(stats)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_coverage_line() {
        let tsv = "name	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
                   contig1	1	1000	150	980	98.0	8.5	35.0	42.0
";
        let stats = parse_coverage_output(tsv).unwrap();
        assert_eq!(stats.len(), 1);
        let s = stats.get("contig1").expect("contig1 missing — check tab escapes in test TSV");
        assert_eq!(s.num_reads, 150);
        assert!((s.coverage - 0.98).abs() < 1e-6);
        assert!((s.mean_depth - 8.5).abs() < 1e-6);
    }

    #[test]
    fn test_score_range() {
        let s = ContigStats {
            num_reads: 100, covered_bases: 900, coverage: 0.90,
            mean_depth: 15.0, mean_base_quality: 35.0, mean_map_quality: 40.0,
            length: 1000,
        };
        let score = s.score();
        assert!(score > 0.0 && score <= 1.0, "score out of range: {score}");
    }

    #[test]
    fn test_zero_length_score() {
        let s = ContigStats {
            num_reads: 0, covered_bases: 0, coverage: 0.0,
            mean_depth: 0.0, mean_base_quality: 0.0, mean_map_quality: 0.0,
            length: 0,
        };
        assert_eq!(s.score(), 0.0);
    }
}
