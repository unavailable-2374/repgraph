//! dtr-gate — bounded, auditable post-discovery gate for te-looker.
//!
//! This is deliberately a standalone binary while dtr's fragment-stitch stage is
//! being migrated.  It provides the two non-negotiable anti-inflation stages:
//! minimum recruited members/length and cross-seed CD-HIT-EST merging.  Its JSON
//! declares `mode=partial`, so Pan_TE will NOT mistake it for a complete gated dtr
//! result until fragment stitching is also internalized.

use std::env;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::{exit, Command};

const DEFAULT_MIN_MEMBERS: usize = 5;
const DEFAULT_MIN_LENGTH: usize = 100;
const DEFAULT_IDENTITY: &str = "0.80";

fn flag(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|arg| arg == name)
        .and_then(|i| args.get(i + 1))
        .cloned()
}

fn need(args: &[String], name: &str) -> String {
    flag(args, name).unwrap_or_else(|| {
        eprintln!("[dtr-gate] {name} required");
        exit(2);
    })
}

fn members(header: &str) -> Option<usize> {
    header
        .split_whitespace()
        .find_map(|token| token.strip_prefix("members=")?.parse::<usize>().ok())
}

fn fasta_stats(path: &str) -> (usize, usize) {
    let file = File::open(path).unwrap_or_else(|err| {
        eprintln!("[dtr-gate] cannot read {path}: {err}");
        exit(1);
    });
    let mut families = 0usize;
    let mut bases = 0usize;
    for line in BufReader::new(file).lines() {
        let line = line.unwrap_or_else(|err| {
            eprintln!("[dtr-gate] cannot read {path}: {err}");
            exit(1);
        });
        if line.starts_with('>') {
            families += 1;
        } else {
            bases += line.trim().len();
        }
    }
    (families, bases)
}

fn json_string(value: &str) -> String {
    value.replace('\\', "\\\\").replace('"', "\\\"")
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 || args[1] == "--help" || args[1] == "-h" {
        eprintln!("usage: dtr-gate --input <families.fa> --out <families.fa> --provenance <json> [--threads N] [--min-members N] [--min-length N] [--identity F]");
        exit(if args.len() < 2 { 2 } else { 0 });
    }
    let input = need(&args, "--input");
    let output = need(&args, "--out");
    let provenance = need(&args, "--provenance");
    let threads = flag(&args, "--threads").unwrap_or_else(|| "1".into());
    let min_members = flag(&args, "--min-members")
        .map(|v| v.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(DEFAULT_MIN_MEMBERS);
    let min_length = flag(&args, "--min-length")
        .map(|v| v.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(DEFAULT_MIN_LENGTH);
    let identity = flag(&args, "--identity").unwrap_or_else(|| DEFAULT_IDENTITY.into());
    let gated = format!("{output}.premerge.fa");

    let fin = File::open(&input).unwrap_or_else(|err| {
        eprintln!("[dtr-gate] cannot open {input}: {err}");
        exit(1);
    });
    let mut fout = BufWriter::new(File::create(&gated).unwrap_or_else(|err| {
        eprintln!("[dtr-gate] cannot create {gated}: {err}");
        exit(1);
    }));
    let mut header = String::new();
    let mut seq = String::new();
    let mut kept = 0usize;
    let mut dropped = 0usize;
    let mut missing_members = 0usize;
    let mut flush = |header: &str, seq: &str| {
        if header.is_empty() {
            return;
        }
        let n = match members(header) {
            Some(n) => n,
            None => {
                missing_members += 1;
                DEFAULT_MIN_MEMBERS
            }
        };
        if n >= min_members && seq.len() >= min_length {
            writeln!(fout, ">{header}").unwrap();
            for chunk in seq.as_bytes().chunks(80) {
                writeln!(fout, "{}", String::from_utf8_lossy(chunk)).unwrap();
            }
            kept += 1;
        } else {
            dropped += 1;
        }
    };
    for line in BufReader::new(fin).lines() {
        let line = line.unwrap_or_else(|err| {
            eprintln!("[dtr-gate] read failure: {err}");
            exit(1);
        });
        if let Some(next) = line.strip_prefix('>') {
            flush(&header, &seq);
            header = next.to_string();
            seq.clear();
        } else {
            seq.push_str(line.trim());
        }
    }
    flush(&header, &seq);
    drop(fout);

    if kept == 0 {
        File::create(&output).unwrap_or_else(|err| {
            eprintln!("[dtr-gate] cannot create {output}: {err}");
            exit(1);
        });
    } else {
        let status = Command::new("cd-hit-est")
            .args([
                "-i", &gated, "-o", &output, "-c", &identity, "-aS", "0.6", "-n", "5", "-r", "1",
                "-g", "1", "-d", "0", "-M", "0", "-T", &threads,
            ])
            .status()
            .unwrap_or_else(|err| {
                eprintln!("[dtr-gate] failed to start cd-hit-est: {err}");
                exit(3);
            });
        if !status.success() {
            eprintln!("[dtr-gate] cd-hit-est failed: {status}");
            exit(3);
        }
    }
    let (input_families, input_bases) = fasta_stats(&input);
    let (gated_families, gated_bases) = fasta_stats(&gated);
    let (output_families, output_bases) = fasta_stats(&output);
    let mut out = BufWriter::new(File::create(&provenance).unwrap_or_else(|err| {
        eprintln!("[dtr-gate] cannot create provenance: {err}");
        exit(1);
    }));
    writeln!(out, "{{").unwrap();
    writeln!(out, "  \"contract\": \"te-looker-discovery-gate-v1\",").unwrap();
    writeln!(out, "  \"mode\": \"partial\",").unwrap();
    writeln!(
        out,
        "  \"final_library\": \"{}\",",
        json_string(
            Path::new(&output)
                .file_name()
                .unwrap()
                .to_string_lossy()
                .as_ref()
        )
    )
    .unwrap();
    writeln!(
        out,
        "  \"input\": {{\"path\": \"{}\", \"families\": {}, \"bases\": {}}},",
        json_string(&input),
        input_families,
        input_bases
    )
    .unwrap();
    writeln!(out, "  \"copy_length_gate\": {{\"min_members\": {}, \"min_length_bp\": {}, \"kept\": {}, \"dropped\": {}, \"headers_missing_members_assumed_minimum\": {}}},", min_members, min_length, kept, dropped, missing_members).unwrap();
    writeln!(
        out,
        "  \"premerge\": {{\"families\": {}, \"bases\": {}}},",
        gated_families, gated_bases
    )
    .unwrap();
    writeln!(out, "  \"cross_seed_merge\": {{\"tool\": \"cd-hit-est\", \"identity\": {}, \"shorter_coverage\": 0.6, \"both_strands\": true}},", identity).unwrap();
    writeln!(out, "  \"fragment_stitch\": \"not_implemented\",").unwrap();
    writeln!(
        out,
        "  \"output\": {{\"path\": \"{}\", \"families\": {}, \"bases\": {}}}",
        json_string(&output),
        output_families,
        output_bases
    )
    .unwrap();
    writeln!(out, "}}").unwrap();
    out.flush().unwrap();
    let _ = fs::remove_file(&gated);
    println!("families\t{output_families}");
}

#[cfg(test)]
mod tests {
    use super::{json_string, members, DEFAULT_MIN_MEMBERS};

    #[test]
    fn member_parser_requires_a_well_formed_header_token() {
        assert_eq!(members("family_1 members=5 length=120"), Some(5));
        assert_eq!(members("family_1 members=0007"), Some(7));
        assert_eq!(members("family_1 notmembers=8"), None);
        assert_eq!(members("family_1 members=bad"), None);
        assert_eq!(DEFAULT_MIN_MEMBERS, 5);
    }

    #[test]
    fn json_escaping_keeps_provenance_valid_for_paths() {
        assert_eq!(json_string(r#"a\\b\"c"#), r#"a\\\\b\\\"c"#);
    }
}
