// src/score.rs
// Rust-native Transfuse scoring replacement
// - Read-support based scoring

use anyhow::{bail, Result};
use csv::{ReaderBuilder, WriterBuilder};
use log::info;
use rayon::prelude::*;
use serde::Deserialize;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use crate::aligner;
use crate::bam;



#[derive(Debug, Clone)]
pub struct ContigScore {

    pub score:f64,
    pub p_good:f64,            // compatibility with TransFuse
    pub p_bases_covered:f64,   // fraction covered
    pub coverage:f64,          // mean depth
    pub depth_score:f64,       // diagnostics
    pub uniformity_score:f64,
    pub pair_score:f64,
}
impl Default for ContigScore {

    fn default()->Self {
        Self {
            score:0.0,
            p_good:0.0,
            p_bases_covered:0.0,
            coverage:0.0,
            depth_score:0.0,
            uniformity_score:0.0,
            pair_score:0.0,
        }
    }
}



pub type ScoreMap = HashMap<String, ContigScore>;
#[derive(Debug, Deserialize)]
struct ScoreRow {
    #[serde(rename="contig_name")]
    name:String,
    score:f64,
    p_good:f64,
    p_bases_covered:f64,
    coverage:f64,
}

// Main scoring

pub fn score_assemblies(
    assembly_files:&[PathBuf],
    left:&Path,
    right:&Path,
    threads:usize,
    verbose:bool,
    prefix_keys:bool,
)->Result<ScoreMap>{

    if assembly_files.is_empty(){
        bail!("No assemblies supplied");
    }

    let asm_threads =
        ((threads as f64 /
            assembly_files.len() as f64)
            .ceil() as usize)
            .max(1);

    info!(
        "Scoring {} assemblies using {} threads",
        assembly_files.len(),
        asm_threads
    );

    let results:Vec<Result<ScoreMap>> =
        assembly_files
            .par_iter()
            .enumerate()
            .map(|(idx, assembly)|{

                // Stable per-assembly index ("contig0", "contig1", ...), matching
                // the same scheme used by fasta::concatenate_assemblies. This
                // function is only ever called with prefix_keys=true against the
                // FULL original assembly list (see main.rs), so idx here lines up
                // exactly with the indices used later when the combined FASTA is
                // written, even if some assemblies get dropped by filtering in
                // between.
                let contig_prefix =
                    format!("contig{idx}");



                let metrics =
                    run_alignment_coverage(
                        assembly,
                        left,
                        right,
                        asm_threads,
                        verbose
                    )?;

                let mut scores =
                    ScoreMap::new();

                for (id,m) in metrics {
                    let final_score =
                        weighted_geometric_mean(&[

                            (m.p_bases_covered,0.45),

                            (m.depth_score,0.30),

                            (m.uniformity_score,0.15),

                            (m.pair_score,0.10),

                        ]);

                    let key =
                        if prefix_keys {
                            format!(
                                "{}_{}",
                                contig_prefix,
                                id
                            )
                        } else {
                            id.clone()
                        };
                    scores.insert(
                        key,
                        ContigScore {

                            score:final_score,

                            p_good:final_score,

                            p_bases_covered:
                            m.p_bases_covered,

                            coverage:
                            m.coverage,

                            depth_score:
                            m.depth_score,

                            uniformity_score:
                            m.uniformity_score,

                            pair_score:
                            m.pair_score,

                        }
                    );
                }


                Ok(scores)

            })
            .collect();



    let mut combined =
        ScoreMap::new();



    for r in results {

        combined.extend(r?);

    }



    Ok(combined)

}

// Alignment + coverage
fn run_alignment_coverage(

    fasta:&Path,

    left:&Path,

    right:&Path,

    threads:usize,

    verbose:bool,

)->Result<HashMap<String,ContigScore>>{
    let bam =
        aligner::align_reads(
            fasta,
            left,
            right,
            threads,
            verbose
        )?;

    let stats =
        bam::compute_coverage(&bam)?;


    let mut out =
        HashMap::new();

    for (id,s) in stats {


        let p_bases_covered =
            if s.length > 0 {
                s.covered_bases as f64 / s.length as f64
            } else {
                0.0

            };



        let coverage =
            s.mean_depth;



        let depth =
            depth_score(
                s.mean_depth
            );



        let pairs =
            pair_score(
                s.num_reads
            );

        let uniformity =
            1.0;

        let final_score =
            weighted_geometric_mean(&[

                (coverage,0.45),

                (depth,0.30),

                (uniformity,0.15),

                (pairs,0.10),

            ]);



        out.insert(id, ContigScore {

            score:final_score,

            p_good:final_score,

            p_bases_covered:
            p_bases_covered,

            coverage:
            s.mean_depth,


            depth_score:
            depth,

            uniformity_score:
            uniformity,

            pair_score:
            pairs,

        });

    }


    Ok(out)

}





// ============================================================
// Metrics
// ============================================================


fn depth_score(depth:f64)->f64{
    if depth <= 0.0 {
        return 0.0;
    }

    (
        1.0 -
            (-depth / 10.0).exp()
    )
        .clamp(0.0,1.0)
}

fn pair_score(reads:u64)->f64{
    if reads == 0 {
        return 0.0;
    }
    (
        1.0 -
            (-(reads as f64)/100.0).exp()
    )
        .clamp(0.0,1.0)

}

// Geometric mean
fn weighted_geometric_mean(
    values:&[(f64,f64)]
)->f64{
    let eps = 1e-6;
    let weight:f64 =
        values
            .iter()
            .map(|x|x.1)
            .sum();

    if weight == 0.0 {
        return 0.0;
    }

    let value =
        values
            .iter()
            .map(|(v,w)|{

                w *
                    v.clamp(0.0,1.0)
                        .max(eps)
                        .ln()

            })
            .sum::<f64>();

    (value / weight)
        .exp()
        .clamp(0.0,1.0)

}
// CSV output

pub fn write_scores_csv(

    scores:&ScoreMap,
    output:&Path
)->Result<PathBuf>{

    let stem =
        output
            .file_stem()
            .unwrap()
            .to_string_lossy();

    let path =
        output
            .with_file_name(
                format!("{}_scores.csv",stem)
            );

    let mut writer =
        WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(&path)?;

    writer.write_record([

        "contig_name",
        "score",
        "p_good",
        "p_bases_covered",
        "coverage",

    ])?;

    let mut rows:
        Vec<_> =
        scores.iter()
            .collect();

    rows.sort_by(|a,b|

        b.1.score
            .partial_cmp(&a.1.score)
            .unwrap()

    );
    for (id,s) in rows {
        writer.write_record([

            id,
            &format!("{:.6}",s.score),
            &format!("{:.6}",s.p_good),
            &format!("{:.6}",
                     s.p_bases_covered),
            &format!("{:.6}",
                     s.coverage),
        ])?;

    }


    writer.flush()?;
    Ok(path)
}
// Load scores

pub fn load_scores_from_csv(
    files:&[PathBuf]
)->Result<ScoreMap>{
    let mut map =
        ScoreMap::new();


    for file in files {
        let reader =
            ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(true)
                .from_reader(
                    BufReader::new(
                        File::open(file)?
                    )
                );

        for row in reader
            .into_deserialize::<ScoreRow>() {
            let row=row?;
            map.insert(
                row.name,
                ContigScore {
                    score:
                    row.score,
                    p_good:
                    row.p_good,
                    p_bases_covered:
                    row.p_bases_covered,
                    coverage:
                    row.coverage,
                    ..Default::default()

                }

            );

        }

    }


    Ok(map)

}