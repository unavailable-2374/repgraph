//! dtr-hmm — Track 1 profile-HMM recruitment for te-looker.
//!
//! Executes nhmmer against the genome, parses its stable --tblout schema, removes
//! redundant overlapping domains per (model, target, strand), applies reproducible
//! length/copy floors, and writes TE-refine-compatible member BED plus provenance.
//! The caller may then build a gapped POA consensus from these recruited copies.

use std::collections::BTreeMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::{exit, Command};

const DEFAULT_EVALUE: &str = "1e-5";
const DEFAULT_MIN_MEMBERS: usize = 5;
const DEFAULT_MIN_LENGTH: usize = 80;
const MAX_OVERLAP_FRACTION: f64 = 0.50;

#[derive(Clone, Debug, PartialEq)]
struct Hit {
    target: String,
    query: String,
    start: i64,
    end: i64,
    strand: char,
    evalue: f64,
    score: f64,
}

fn flag(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|arg| arg == name)
        .and_then(|i| args.get(i + 1))
        .cloned()
}

fn need(args: &[String], name: &str) -> String {
    flag(args, name).unwrap_or_else(|| {
        eprintln!("[dtr-hmm] {name} required");
        exit(2);
    })
}

fn sanitize_id(value: &str) -> String {
    let clean: String = value
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' {
                c
            } else {
                '_'
            }
        })
        .collect();
    if clean.is_empty() {
        "hmm".to_string()
    } else {
        format!("hmm_{clean}")
    }
}

fn parse_tbl_line(line: &str) -> Option<Hit> {
    if line.is_empty() || line.starts_with('#') {
        return None;
    }
    let fields: Vec<&str> = line.split_whitespace().collect();
    // HMMER 3 --tblout: target, accession, query, accession, hmmfrom, hmmto,
    // alifrom, alito, envfrom, envto, sq_len, strand, E-value, score, bias, ...
    if fields.len() < 15 {
        return None;
    }
    let a = fields[6].parse::<i64>().ok()?;
    let b = fields[7].parse::<i64>().ok()?;
    let start = a.min(b) - 1; // BED: 0-based, half-open
    let end = a.max(b);
    if start < 0 || end <= start {
        return None;
    }
    Some(Hit {
        target: fields[0].to_string(),
        query: fields[2].to_string(),
        start,
        end,
        strand: if fields[11] == "-" { '-' } else { '+' },
        evalue: fields[12].parse().ok()?,
        score: fields[13].parse().ok()?,
    })
}

fn overlap_fraction(a: &Hit, b: &Hit) -> f64 {
    let overlap = (a.end.min(b.end) - a.start.max(b.start)).max(0) as f64;
    overlap / ((a.end - a.start).min(b.end - b.start).max(1) as f64)
}

fn select_hits(mut hits: Vec<Hit>, min_length: usize, min_members: usize) -> Vec<Hit> {
    hits.retain(|h| (h.end - h.start) as usize >= min_length);
    let mut by_locus: BTreeMap<(String, String), Vec<Hit>> = BTreeMap::new();
    for hit in hits {
        // An HMMER hit on both strands at the same genomic interval is one physical
        // copy, not two members. Strand remains in the emitted BED after the strongest
        // overlapping domain has been selected.
        by_locus
            .entry((hit.query.clone(), hit.target.clone()))
            .or_default()
            .push(hit);
    }
    let mut retained = Vec::new();
    for (_, mut group) in by_locus {
        // Keep the strongest domain first; a second domain must be substantially distinct.
        group.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap()
                .then_with(|| a.start.cmp(&b.start))
        });
        let mut chosen: Vec<Hit> = Vec::new();
        for hit in group {
            if chosen
                .iter()
                .all(|previous| overlap_fraction(&hit, previous) <= MAX_OVERLAP_FRACTION)
            {
                chosen.push(hit);
            }
        }
        retained.extend(chosen);
    }
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for hit in &retained {
        *counts.entry(hit.query.clone()).or_default() += 1;
    }
    retained.retain(|hit| counts.get(&hit.query).copied().unwrap_or(0) >= min_members);
    retained.sort_by(|a, b| {
        (&a.query, &a.target, a.start, a.end, a.strand)
            .cmp(&(&b.query, &b.target, b.start, b.end, b.strand))
    });
    retained
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 || args[1] == "--help" || args[1] == "-h" {
        eprintln!("usage: dtr-hmm --genome <fa> --hmm <profiles.hmm> --out <members.bed> --provenance <json> [--threads N] [--evalue E] [--min-members N] [--min-length N]");
        exit(if args.len() < 2 { 2 } else { 0 });
    }
    let genome = need(&args, "--genome");
    let hmm = need(&args, "--hmm");
    let output = need(&args, "--out");
    let provenance = need(&args, "--provenance");
    let threads = flag(&args, "--threads").unwrap_or_else(|| "1".into());
    let evalue = flag(&args, "--evalue").unwrap_or_else(|| DEFAULT_EVALUE.into());
    let min_members = flag(&args, "--min-members")
        .map(|v| v.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(DEFAULT_MIN_MEMBERS);
    let min_length = flag(&args, "--min-length")
        .map(|v| v.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(DEFAULT_MIN_LENGTH);
    let raw_table = format!("{output}.nhmmer.tbl");

    let status = Command::new("nhmmer")
        .args([
            "--noali",
            "-o",
            "/dev/null",
            "--tblout",
            &raw_table,
            "-E",
            &evalue,
            "--cpu",
            &threads,
            &hmm,
            &genome,
        ])
        .status()
        .unwrap_or_else(|err| {
            eprintln!("[dtr-hmm] failed to start nhmmer: {err}");
            exit(3);
        });
    if !status.success() {
        eprintln!("[dtr-hmm] nhmmer failed: {status}");
        exit(3);
    }

    let raw: Vec<Hit> = BufReader::new(File::open(&raw_table).unwrap())
        .lines()
        .filter_map(|line| line.ok().and_then(|text| parse_tbl_line(&text)))
        .collect();
    let selected = select_hits(raw.clone(), min_length, min_members);
    let mut out = BufWriter::new(File::create(&output).unwrap_or_else(|err| {
        eprintln!("[dtr-hmm] cannot write {output}: {err}");
        exit(1);
    }));
    for hit in &selected {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{:.3}\t{}",
            hit.target,
            hit.start,
            hit.end,
            sanitize_id(&hit.query),
            hit.score,
            hit.strand
        )
        .unwrap();
    }
    out.flush().unwrap();
    let mut families: BTreeMap<String, usize> = BTreeMap::new();
    for hit in &selected {
        *families.entry(hit.query.clone()).or_default() += 1;
    }
    let mut prov = BufWriter::new(File::create(&provenance).unwrap_or_else(|err| {
        eprintln!("[dtr-hmm] cannot write provenance: {err}");
        exit(1);
    }));
    writeln!(prov, "{{\n  \"contract\": \"te-looker-hmm-track-v1\",\n  \"mode\": \"partial\",\n  \"tool\": \"nhmmer\",\n  \"hmm\": \"{}\",\n  \"genome\": \"{}\",\n  \"evalue\": \"{}\",\n  \"min_members\": {},\n  \"min_length_bp\": {},\n  \"raw_hits\": {},\n  \"selected_hits\": {},\n  \"families\": {}\n}}", hmm.replace('\\', "\\\\").replace('"', "\\\""), genome.replace('\\', "\\\\").replace('"', "\\\""), evalue, min_members, min_length, raw.len(), selected.len(), families.len()).unwrap();
    prov.flush().unwrap();
    println!("hmm_hits\t{}\nfamilies\t{}", selected.len(), families.len());
}

#[cfg(test)]
mod tests {
    use super::{parse_tbl_line, sanitize_id, select_hits, Hit};

    fn hit(query: &str, start: i64, end: i64, score: f64) -> Hit {
        Hit {
            target: "chr1".into(),
            query: query.into(),
            start,
            end,
            strand: '+',
            evalue: 1e-8,
            score,
        }
    }

    #[test]
    fn parses_nhmmer_reverse_coordinates_as_bed() {
        let line = "chr1 - seed - 2 38 67 27 67 27 92 - 1.6e-08 20.4 0.1 -";
        let parsed = parse_tbl_line(line).unwrap();
        assert_eq!((parsed.start, parsed.end, parsed.strand), (26, 67, '-'));
    }

    #[test]
    fn collapses_redundant_domains_before_copy_gate() {
        let hits = vec![
            hit("A", 0, 100, 30.0),
            hit("A", 10, 110, 29.0),
            hit("A", 200, 300, 28.0),
            hit("A", 400, 500, 27.0),
        ];
        let selected = select_hits(hits, 80, 3);
        assert_eq!(selected.len(), 3);
        assert_eq!(selected[0].start, 0);
    }

    #[test]
    fn collapses_same_locus_reverse_strand_duplicate() {
        let mut reverse = hit("A", 0, 100, 29.0);
        reverse.strand = '-';
        let hits = vec![
            hit("A", 0, 100, 30.0),
            reverse,
            hit("A", 200, 300, 28.0),
            hit("A", 400, 500, 27.0),
        ];
        let selected = select_hits(hits, 80, 3);
        assert_eq!(selected.len(), 3);
        assert_eq!(selected[0].strand, '+');
    }

    #[test]
    fn preserves_safe_identifier_for_te_refine() {
        assert_eq!(sanitize_id("Gypsy/1.2"), "hmm_Gypsy_1_2");
    }
}
