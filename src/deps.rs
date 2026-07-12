// src/deps.rs
//
// Binary dependency management — replaces the Ruby bindeps gem.
//
// Tools managed:
//   vsearch   — sequence clustering
//   minimap2  — read alignment (replaces SNAP)
//   samtools  — BAM processing + coverage (replaces bam-read)
//   salmon    — RNA-seq quasi-mapping quantification (modern 1.x)

use anyhow::{bail, Context, Result};
use colored::Colorize;
use log::{debug, info, warn};
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Debug, Clone)]
pub struct BinaryDep {
    pub name: &'static str,
    pub version: &'static str,
    pub version_arg: &'static str,
    pub version_pattern: &'static str,
    pub download_url_linux: &'static str,
    pub download_url_macos: &'static str,
}

pub const REQUIRED_DEPS: &[BinaryDep] = &[
    BinaryDep {
        name: "vsearch", version: "2.0",
        version_arg: "--version", version_pattern: "vsearch v",
        download_url_linux: "https://github.com/torognes/vsearch/releases/download/v2.28.1/vsearch-2.28.1-linux-x86_64.tar.gz",
        download_url_macos: "https://github.com/torognes/vsearch/releases/download/v2.28.1/vsearch-2.28.1-macos-aarch64.tar.gz",
    },
    BinaryDep {
        name: "minimap2", version: "2.0",
        version_arg: "--version", version_pattern: "",
        download_url_linux: "https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2",
        download_url_macos: "",  // brew install minimap2
    },
    BinaryDep {
        name: "samtools", version: "1.0",
        version_arg: "--version", version_pattern: "samtools ",
        download_url_linux: "",  // apt / conda / brew
        download_url_macos: "",
    },
    BinaryDep {
        name: "salmon", version: "1.0",
        version_arg: "--version", version_pattern: "salmon ",
        download_url_linux: "https://github.com/COMBINE-lab/salmon/releases/download/v1.10.3/salmon-1.10.3_linux_x86_64.tar.gz",
        download_url_macos: "https://github.com/COMBINE-lab/salmon/releases/download/v1.10.3/salmon-1.10.3_mac_osx.tar.gz",
    },
];

// ── Public API ────────────────────────────────────────────────────────────────

pub fn check_dependencies() -> Result<Vec<BinaryDep>> {
    let mut missing = Vec::new();
    for dep in REQUIRED_DEPS {
        match probe_binary(dep) {
            Ok(ver)  => info!("  {} {} ✓  (found {})", dep.name.green(), dep.version, ver.dimmed()),
            Err(e)   => { warn!("  {} {} ✗  ({})", dep.name.red(), dep.version, e);
                          missing.push(dep.clone()); }
        }
    }
    Ok(missing)
}

pub fn install_dependencies(missing: &[BinaryDep]) -> Result<()> {
    let install_dir = default_install_dir();
    std::fs::create_dir_all(&install_dir)
        .with_context(|| format!("Cannot create install directory {:?}", install_dir))?;
    info!("Installing to {:?}", install_dir);

    for dep in missing {
        let url = if cfg!(target_os = "macos") { dep.download_url_macos }
                  else { dep.download_url_linux };
        if url.is_empty() {
            warn!("No automatic download for '{}'. Install via apt, brew, or conda.", dep.name);
            continue;
        }
        info!("Downloading {} from {}", dep.name, url);
        download_and_install(dep.name, url, &install_dir)?;
        info!("{} installed to {:?}", dep.name.green(), install_dir);
    }
    println!("
Add {:?} to your PATH if not already included.", install_dir);
    Ok(())
}

// ── Internals ─────────────────────────────────────────────────────────────────

fn probe_binary(dep: &BinaryDep) -> Result<String> {
    let output = Command::new(dep.name)
        .arg(dep.version_arg)
        .output()
        .with_context(|| format!("'{}' not found on PATH", dep.name))?;
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    let combined = format!("{}{}", stdout, stderr);
    let first_line = combined.lines().next().unwrap_or("").to_string();
    debug!("  {} version output: {}", dep.name, first_line.trim());
    Ok(first_line)
}

fn download_and_install(name: &str, url: &str, dir: &Path) -> Result<()> {
    let tmp = tempfile::tempdir().context("Failed to create temp dir")?;
    let archive_name = url.split('/').last().unwrap_or("archive");
    let archive_path = tmp.path().join(archive_name);

    let status = Command::new("curl").args(["-L", "-o"]).arg(&archive_path).arg(url)
        .status().context("Failed to run curl. Please install curl.")?;
    if !status.success() { bail!("curl failed for {}", url); }

    let ext = archive_name.to_lowercase();
    if ext.ends_with(".tar.gz") || ext.ends_with(".tgz") {
        Command::new("tar").args(["xzf"]).arg(&archive_path).arg("-C").arg(tmp.path())
            .status().context("Failed to run tar")?;
    } else if ext.ends_with(".tar.bz2") {
        Command::new("tar").args(["xjf"]).arg(&archive_path).arg("-C").arg(tmp.path())
            .status().context("Failed to run tar")?;
    } else if ext.ends_with(".zip") {
        Command::new("unzip").arg("-q").arg(&archive_path).arg("-d").arg(tmp.path())
            .status().context("Failed to run unzip")?;
    }

    let binary_path = find_binary_in_dir(tmp.path(), name)?;
    let dest = dir.join(name);
    std::fs::copy(&binary_path, &dest)
        .with_context(|| format!("Failed to copy {:?} to {:?}", binary_path, dest))?;

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut perms = std::fs::metadata(&dest)?.permissions();
        perms.set_mode(0o755);
        std::fs::set_permissions(&dest, perms)?;
    }
    Ok(())
}

fn find_binary_in_dir(dir: &Path, name: &str) -> Result<PathBuf> {
    for entry in walkdir(dir) {
        if entry.file_name().to_string_lossy() == name {
            return Ok(entry.path());
        }
    }
    bail!("Could not find '{}' in unpacked archive", name)
}

fn walkdir(dir: &Path) -> Vec<std::fs::DirEntry> {
    let mut results = Vec::new();
    if let Ok(entries) = std::fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() { results.extend(walkdir(&path)); }
            else { results.push(entry); }
        }
    }
    results
}

fn default_install_dir() -> PathBuf {
    std::env::var("HOME").map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/tmp"))
        .join(".local").join("bin")
}
