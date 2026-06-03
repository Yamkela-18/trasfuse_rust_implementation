// src/fasta.rs
//
// FASTA I/O — replaces BioRuby (bio gem) Bio::FlatFile / Bio::FastaFormat.
//
// Pure-Rust implementation; no external parser crate required.
// Handles:
//   - load_fasta()             -> HashMap<id, FastaRecord>
//   - concatenate_assemblies() -> merged FASTA with stem__id prefixing
//   - write_fasta_record()     -> 60-char line wrapping
//   - write_fasta_file()       -> write an ordered set of records

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

// ── Public types ──────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FastaRecord {
    pub header: String,
    pub sequence: String,
}

impl FastaRecord {
    pub fn id(&self) -> &str {
        self.header.split_whitespace().next().unwrap_or(&self.header)
    }
}

// ── load_fasta ────────────────────────────────────────────────────────────────

pub fn load_fasta(path: &Path) -> Result<HashMap<String, FastaRecord>> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open FASTA file {:?}", path))?;
    parse_fasta(BufReader::new(file))
}

pub fn parse_fasta<R: BufRead>(reader: R) -> Result<HashMap<String, FastaRecord>> {
    let mut map: HashMap<String, FastaRecord> = HashMap::new();
    let mut current_header: Option<String> = None;
    let mut current_seq = String::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("I/O error on line {}", line_no + 1))?;
        let line = line.trim_end();
        if line.is_empty() { continue; }

        if let Some(header) = line.strip_prefix('>') {
            if let Some(hdr) = current_header.take() {
                let rec = FastaRecord { header: hdr, sequence: current_seq.to_uppercase() };
                map.insert(rec.id().to_string(), rec);
            }
            current_header = Some(header.to_string());
            current_seq = String::new();
        } else {
            current_seq.push_str(line);
        }
    }
    if let Some(hdr) = current_header {
        let rec = FastaRecord { header: hdr, sequence: current_seq.to_uppercase() };
        map.insert(rec.id().to_string(), rec);
    }
    Ok(map)
}

pub fn load_fasta_ordered(path: &Path) -> Result<Vec<FastaRecord>> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open FASTA file {:?}", path))?;
    parse_fasta_ordered(BufReader::new(file))
}

pub fn parse_fasta_ordered<R: BufRead>(reader: R) -> Result<Vec<FastaRecord>> {
    let mut records = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim_end();
        if line.is_empty() { continue; }

        if let Some(header) = line.strip_prefix('>') {
            if let Some(hdr) = current_header.take() {
                records.push(FastaRecord { header: hdr, sequence: current_seq.to_uppercase() });
            }
            current_header = Some(header.to_string());
            current_seq = String::new();
        } else {
            current_seq.push_str(line);
        }
    }
    if let Some(hdr) = current_header {
        records.push(FastaRecord { header: hdr, sequence: current_seq.to_uppercase() });
    }
    Ok(records)
}

// ── concatenate_assemblies ────────────────────────────────────────────────────

pub fn concatenate_assemblies(
    assembly_files: &[PathBuf],
    output_path: &Path,
) -> Result<PathBuf> {
    let parent = output_path.parent().unwrap_or(Path::new("."));
    let stem = output_path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "combined".into());
    let cat_path = parent.join(format!("{stem}_combined.fa"));

    let mut writer = BufWriter::new(
        File::create(&cat_path)
            .with_context(|| format!("Cannot create {:?}", cat_path))?
    );

    for asm_path in assembly_files {
        let prefix = asm_path
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "asm".into());

        for rec in load_fasta_ordered(asm_path)? {
            let new_id = format!("{}__{}", prefix, rec.id());
            let description = rec.header[rec.id().len()..].trim();
            let new_header = if description.is_empty() {
                new_id.clone()
            } else {
                format!("{new_id} {description}")
            };
            write_fasta_record(&mut writer, &new_header, &rec.sequence)?;
        }
    }
    writer.flush()?;
    Ok(cat_path)
}

// ── Write helpers ─────────────────────────────────────────────────────────────

pub fn write_fasta_record<W: Write>(
    writer: &mut W, header: &str, sequence: &str,
) -> Result<()> {
    writeln!(writer, ">{header}")?;
    for chunk in sequence.as_bytes().chunks(60) {
        writer.write_all(chunk)?;
        writeln!(writer)?;
    }
    Ok(())
}

pub fn write_fasta_file(path: &Path, records: &[FastaRecord]) -> Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("Cannot create {:?}", path))?
    );
    for rec in records {
        write_fasta_record(&mut writer, &rec.header, &rec.sequence)?;
    }
    writer.flush()?;
    Ok(())
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_single_record() {
        let fasta = ">seq1 description
ACGT
ACGT
";
        let map = parse_fasta(BufReader::new(Cursor::new(fasta))).unwrap();
        assert_eq!(map.len(), 1);
        let rec = map.get("seq1").unwrap();
        assert_eq!(rec.sequence, "ACGTACGT");
        assert_eq!(rec.id(), "seq1");
    }

    #[test]
    fn test_parse_multiple_records() {
        let fasta = ">seq1
ACGT
>seq2
TTTT
>seq3
GGGG
";
        let map = parse_fasta(BufReader::new(Cursor::new(fasta))).unwrap();
        assert_eq!(map.len(), 3);
        assert_eq!(map["seq2"].sequence, "TTTT");
    }

    #[test]
    fn test_empty_lines_ignored() {
        let fasta = ">seq1

ACGT

ACGT

";
        let map = parse_fasta(BufReader::new(Cursor::new(fasta))).unwrap();
        assert_eq!(map["seq1"].sequence, "ACGTACGT");
    }

    #[test]
    fn test_lowercase_normalised() {
        let fasta = ">seq1
acgt
";
        let map = parse_fasta(BufReader::new(Cursor::new(fasta))).unwrap();
        assert_eq!(map["seq1"].sequence, "ACGT");
    }

    #[test]
    fn test_write_round_trip() {
        use tempfile::NamedTempFile;
        let tmp = NamedTempFile::new().unwrap();
        let records = vec![
            FastaRecord { header: "r1".into(), sequence: "A".repeat(120) },
            FastaRecord { header: "r2".into(), sequence: "GCGC".into() },
        ];
        write_fasta_file(tmp.path(), &records).unwrap();
        let loaded = load_fasta_ordered(tmp.path()).unwrap();
        assert_eq!(loaded.len(), 2);
        assert_eq!(loaded[0].sequence, "A".repeat(120));
        assert_eq!(loaded[1].sequence, "GCGC");
    }

    #[test]
    fn test_concatenate_prefixes_ids() {
        use tempfile::tempdir;
        let dir = tempdir().unwrap();
        let asm1 = dir.path().join("k31.fa");
        let asm2 = dir.path().join("k41.fa");
        std::fs::write(&asm1, ">t1
AAAA
>t2
CCCC
").unwrap();
        std::fs::write(&asm2, ">t1
GGGG
").unwrap();
        let output = dir.path().join("out.fa");
        let cat = concatenate_assemblies(&[asm1, asm2], &output).unwrap();
        let records = load_fasta_ordered(&cat).unwrap();
        assert_eq!(records.len(), 3);
        assert!(records.iter().any(|r| r.id() == "k31__t1"));
        assert!(records.iter().any(|r| r.id() == "k41__t1"));
    }
}
