//! dtr-boundary — conservative, copy-supported adaptive boundary proposals.
//! Uses non-nested copies plus symmetric genomic flanks.  Selection is deferred
//! to dtr-boundary-gate, which independently requires a coverage gain.
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::{exit, Command};

#[derive(Clone)]
struct Rec {
    id: String,
    header: String,
    seq: String,
}
#[derive(Clone)]
struct Copy {
    chrom: String,
    start: usize,
    end: usize,
    strand: char,
}

fn flag(a: &[String], k: &str) -> Option<String> {
    a.iter()
        .position(|x| x == k)
        .and_then(|i| a.get(i + 1))
        .cloned()
}
fn need(a: &[String], k: &str) -> String {
    flag(a, k).unwrap_or_else(|| {
        eprintln!("[dtr-boundary] {k} is required");
        exit(2)
    })
}
fn fasta(path: &str) -> Vec<Rec> {
    let mut v = Vec::new();
    let (mut h, mut s) = (String::new(), String::new());
    for l in BufReader::new(File::open(path).unwrap()).lines() {
        let l = l.unwrap();
        if let Some(x) = l.strip_prefix('>') {
            if !h.is_empty() {
                let id = h.split_whitespace().next().unwrap().to_string();
                v.push(Rec {
                    id,
                    header: h,
                    seq: s,
                });
                s = String::new()
            }
            h = x.to_string()
        } else {
            s.push_str(l.trim())
        }
    }
    if !h.is_empty() {
        let id = h.split_whitespace().next().unwrap().to_string();
        v.push(Rec {
            id,
            header: h,
            seq: s,
        })
    }
    v
}
fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'N' => 'N',
            _ => 'N',
        })
        .collect()
}
fn writefa(w: &mut BufWriter<File>, h: &str, s: &str) {
    writeln!(w, ">{h}").unwrap();
    for x in s.as_bytes().chunks(80) {
        w.write_all(x).unwrap();
        w.write_all(b"\n").unwrap()
    }
}

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 2 || a[1] == "--help" {
        eprintln!("usage: dtr-boundary --genome <fa> --query <fa> --copies <clean.bed> --out <proposal.fa> --provenance <json> [--flank 100 --min-members 5 --max-members 20 --max-poa-bases 4000 --max-refinement-families N --max-length-ratio 1.5]");
        exit(if a.len() < 2 { 2 } else { 0 })
    }
    let genome = need(&a, "--genome");
    let query = need(&a, "--query");
    let copies = need(&a, "--copies");
    let out = need(&a, "--out");
    let prov = need(&a, "--provenance");
    let flank: usize = flag(&a, "--flank")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(100);
    let min: usize = flag(&a, "--min-members")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(5);
    let max: usize = flag(&a, "--max-members")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(20);
    let max_poa: usize = flag(&a, "--max-poa-bases")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(4000);
    let budget: usize = flag(&a, "--max-refinement-families")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(usize::MAX);
    let ratio: f64 = flag(&a, "--max-length-ratio")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(1.5);
    let gs: HashMap<String, String> = fasta(&genome).into_iter().map(|r| (r.id, r.seq)).collect();
    let refs = fasta(&query);
    let mut by: HashMap<String, Vec<Copy>> = HashMap::new();
    for l in BufReader::new(File::open(&copies).unwrap())
        .lines()
        .map_while(Result::ok)
    {
        let f: Vec<_> = l.split('\t').collect();
        if f.len() >= 4 {
            if let (Ok(start), Ok(end)) = (f[1].parse(), f[2].parse()) {
                by.entry(f[3].into()).or_default().push(Copy {
                    chrom: f[0].into(),
                    start,
                    end,
                    strand: f.get(5).and_then(|x| x.chars().next()).unwrap_or('+'),
                })
            }
        }
    }
    for v in by.values_mut() {
        v.sort_by(|x, y| {
            x.chrom
                .cmp(&y.chrom)
                .then(x.start.cmp(&y.start))
                .then(x.end.cmp(&y.end))
        })
    }
    let mut w = BufWriter::new(File::create(&out).unwrap());
    let tmp = format!("{}.spoa_in.fa", out);
    let (mut eligible, mut attempted, mut proposed, mut unchanged) =
        (0usize, 0usize, 0usize, 0usize);
    for r in refs {
        let mut candidate = None;
        if let Some(v) = by.get(&r.id) {
            if v.len() >= min && r.seq.len() <= max_poa {
                eligible += 1;
                if attempted < budget {
                    attempted += 1;
                    let mut t = BufWriter::new(File::create(&tmp).unwrap());
                    let mut used = 0usize;
                    for c in v.iter().take(max) {
                        if let Some(g) = gs.get(&c.chrom) {
                            let lo = c.start.saturating_sub(flank);
                            let hi = (c.end + flank).min(g.len());
                            if lo < hi {
                                let s = if c.strand == '-' {
                                    rc(&g[lo..hi])
                                } else {
                                    g[lo..hi].to_string()
                                };
                                writeln!(t, ">{}_{}\n{}", r.id, used, s).unwrap();
                                used += 1
                            }
                        }
                    }
                    t.flush().unwrap();
                    if used >= min {
                        let o =
                            Command::new(env::var("SPOA_BIN").unwrap_or_else(|_| "spoa".into()))
                                .args(["-r", "0", &tmp])
                                .output()
                                .unwrap_or_else(|e| {
                                    eprintln!("[dtr-boundary] spoa: {e}");
                                    exit(3)
                                });
                        if !o.status.success() {
                            exit(3)
                        }
                        let c: String = String::from_utf8_lossy(&o.stdout)
                            .lines()
                            .filter(|x| !x.starts_with('>'))
                            .collect();
                        let l = c.len() as f64;
                        if c.len() >= 80
                            && l >= r.seq.len() as f64 / ratio
                            && l <= r.seq.len() as f64 * ratio
                        {
                            candidate = Some(c)
                        }
                    }
                }
            }
        }
        if let Some(c) = candidate {
            proposed += 1;
            writefa(&mut w, &r.header, &c)
        } else {
            unchanged += 1;
            writefa(&mut w, &r.header, &r.seq)
        }
    }
    w.flush().unwrap();
    let _ = std::fs::remove_file(&tmp);
    writeln!(File::create(&prov).unwrap(), "{{\"contract\":\"te-looker-adaptive-boundary-v1\",\"mode\":\"partial\",\"family_call\":false,\"selection_deferred_to_coverage_gate\":true,\"flank\":{},\"min_members\":{},\"max_members\":{},\"max_poa_bases\":{},\"max_refinement_families\":{},\"attempted_families\":{},\"eligible_families\":{},\"proposed_families\":{},\"unchanged_families\":{}}}", flank, min, max, max_poa, budget, attempted, eligible, proposed, unchanged).unwrap();
    println!("eligible_families\t{eligible}\nproposed_families\t{proposed}\nunchanged_families\t{unchanged}")
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn revcomp() {
        assert_eq!(rc("AaCGn"), "NCGTT")
    }
}
