//! dtr-nonte-guard conservatively removes candidates explained by an explicit non-TE reference.
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::{exit, Command};

#[derive(Clone)]
struct Record {
    header: String,
    seq: String,
}
#[derive(Clone)]
struct Hit {
    query: String,
    subject: String,
    identity: f64,
    qlen: usize,
    start: usize,
    end: usize,
}

fn value(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|x| x == name)
        .and_then(|i| args.get(i + 1))
        .cloned()
}
fn required(args: &[String], name: &str) -> String {
    value(args, name).unwrap_or_else(|| {
        eprintln!("[dtr-nonte-guard] {name} is required");
        exit(2)
    })
}
fn records(path: &str) -> Vec<Record> {
    let mut out = Vec::new();
    let mut header = String::new();
    let mut seq = String::new();
    for line in BufReader::new(File::open(path).unwrap()).lines() {
        let line = line.unwrap();
        if let Some(h) = line.strip_prefix(">") {
            if !header.is_empty() {
                out.push(Record { header, seq });
                seq = String::new();
            }
            header = h.to_string();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !header.is_empty() {
        out.push(Record { header, seq });
    }
    out
}
fn name(header: &str) -> &str {
    header.split_whitespace().next().unwrap_or(header)
}
fn parse(line: &str) -> Option<Hit> {
    let f: Vec<&str> = line.split("\t").collect();
    if f.len() < 8 {
        return None;
    }
    let a: usize = f[5].parse().ok()?;
    let b: usize = f[6].parse().ok()?;
    Some(Hit {
        query: f[0].to_string(),
        subject: f[1].to_string(),
        identity: f[2].parse().ok()?,
        qlen: f[4].parse().ok()?,
        start: a.min(b),
        end: a.max(b),
    })
}
fn covered(mut iv: Vec<(usize, usize)>) -> usize {
    iv.sort_unstable();
    let mut total = 0;
    let mut cur: Option<(usize, usize)> = None;
    for (a, b) in iv {
        match cur {
            None => cur = Some((a, b)),
            Some((x, y)) if a <= y + 1 => cur = Some((x, y.max(b))),
            Some((x, y)) => {
                total += y - x + 1;
                cur = Some((a, b));
            }
        }
    }
    if let Some((x, y)) = cur {
        total += y - x + 1;
    }
    total
}
fn esc(s: &str) -> String {
    s.replace(char::from(92), "\\\\").replace("\"", "\\\"")
}

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 2 || a[1] == "--help" {
        eprintln!("usage: dtr-nonte-guard --query <candidates.fa> --reference <non_te.fa> --out-retained <fa> --out-rejected <fa> --report <tsv> --provenance <json> [--min-identity 80 --min-qcov .8 --threads 1]");
        exit(if a.len() < 2 { 2 } else { 0 });
    }
    let query = required(&a, "--query");
    let reference = required(&a, "--reference");
    let retained = required(&a, "--out-retained");
    let rejected = required(&a, "--out-rejected");
    let report = required(&a, "--report");
    let provenance = required(&a, "--provenance");
    let min_identity: f64 = value(&a, "--min-identity")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(80.0);
    let min_qcov: f64 = value(&a, "--min-qcov")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(0.8);
    let threads = value(&a, "--threads").unwrap_or_else(|| "1".to_string());
    let db = format!("{report}.non_te_db");
    let raw = format!("{report}.raw_hsps.tsv");
    let status = Command::new("makeblastdb")
        .args(["-in", &reference, "-dbtype", "nucl", "-out", &db])
        .status()
        .unwrap();
    if !status.success() {
        exit(3);
    }
    let status = Command::new("blastn")
        .args([
            "-task",
            "blastn",
            "-query",
            &query,
            "-db",
            &db,
            "-dust",
            "no",
            "-soft_masking",
            "false",
            "-evalue",
            "1e-10",
            "-num_threads",
            &threads,
            "-outfmt",
            "6 qseqid sseqid pident length qlen qstart qend bitscore",
            "-out",
            &raw,
        ])
        .status()
        .unwrap();
    if !status.success() {
        exit(3);
    }
    let mut hits: HashMap<String, Vec<Hit>> = HashMap::new();
    let mut raw_hits = 0usize;
    for line in BufReader::new(File::open(&raw).unwrap()).lines() {
        if let Some(hit) = parse(&line.unwrap()) {
            raw_hits += 1;
            if hit.identity >= min_identity {
                hits.entry(hit.query.clone()).or_default().push(hit);
            }
        }
    }
    let mut detail: HashMap<String, (f64, f64, String)> = HashMap::new();
    for (q, hs) in hits {
        let cov = covered(hs.iter().map(|h| (h.start, h.end)).collect()) as f64 / hs[0].qlen as f64;
        let best = hs
            .iter()
            .max_by(|x, y| x.identity.total_cmp(&y.identity))
            .unwrap();
        if cov >= min_qcov {
            detail.insert(q, (cov, best.identity, best.subject.clone()));
        }
    }
    let rejected_ids: HashSet<String> = detail.keys().cloned().collect();
    let all = records(&query);
    let mut keep = BufWriter::new(File::create(&retained).unwrap());
    let mut drop = BufWriter::new(File::create(&rejected).unwrap());
    let mut rep = BufWriter::new(File::create(&report).unwrap());
    writeln!(
        rep,
        "candidate\tdecision\tqcov\tbest_identity\tbest_non_te_reference"
    )
    .unwrap();
    let mut kept = 0usize;
    let mut dropped = 0usize;
    for r in all {
        let id = name(&r.header);
        if rejected_ids.contains(id) {
            let (cov, ident, subject) = detail.get(id).unwrap();
            writeln!(drop, ">{}\n{}", r.header, r.seq).unwrap();
            writeln!(rep, "{id}\treject_non_te\t{cov:.6}\t{ident:.3}\t{subject}").unwrap();
            dropped += 1;
        } else {
            writeln!(keep, ">{}\n{}", r.header, r.seq).unwrap();
            writeln!(rep, "{id}\tretain\t0\t0\t.").unwrap();
            kept += 1;
        }
    }
    keep.flush().unwrap();
    drop.flush().unwrap();
    rep.flush().unwrap();
    writeln!(File::create(&provenance).unwrap(), "{{\"contract\":\"te-looker-non-te-guard-v1\",\"mode\":\"partial\",\"family_call\":false,\"query\":\"{}\",\"reference\":\"{}\",\"raw_hsps\":{},\"retained\":{},\"rejected\":{},\"min_identity\":{},\"min_qcov\":{}}}", esc(&query), esc(&reference), raw_hits, kept, dropped, min_identity, min_qcov).unwrap();
    println!("retained\t{kept}\nrejected_non_te\t{dropped}");
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn coverage_merges_overlapping_intervals() {
        assert_eq!(covered(vec![(2, 10), (8, 15), (20, 22)]), 17);
    }
    #[test]
    fn parses_blast6_record() {
        let h = parse("q\tr\t99.0\t50\t100\t2\t51\t9").unwrap();
        assert_eq!((h.start, h.end, h.qlen), (2, 51, 100));
    }
}
