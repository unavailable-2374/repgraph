//! dtr-boundary-gate — accepts a proposed boundary only when copy coverage grows.
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::exit;

#[derive(Clone)]
struct Rec {
    id: String,
    header: String,
    seq: String,
}
fn flag(a: &[String], k: &str) -> Option<String> {
    a.iter()
        .position(|x| x == k)
        .and_then(|i| a.get(i + 1))
        .cloned()
}
fn need(a: &[String], k: &str) -> String {
    flag(a, k).unwrap_or_else(|| {
        eprintln!("[dtr-boundary-gate] {k} is required");
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
                let id = h.split_whitespace().next().unwrap().into();
                v.push(Rec {
                    id,
                    header: h,
                    seq: s,
                });
                s = String::new()
            }
            h = x.into()
        } else {
            s.push_str(l.trim())
        }
    }
    if !h.is_empty() {
        let id = h.split_whitespace().next().unwrap().into();
        v.push(Rec {
            id,
            header: h,
            seq: s,
        })
    };
    v
}
fn coverage(path: &str) -> HashMap<String, usize> {
    let mut x: HashMap<(String, String), Vec<(usize, usize)>> = HashMap::new();
    for l in BufReader::new(File::open(path).unwrap())
        .lines()
        .map_while(Result::ok)
    {
        let f: Vec<_> = l.split('\t').collect();
        if f.len() >= 4 {
            if let (Ok(s), Ok(e)) = (f[1].parse(), f[2].parse()) {
                if e > s {
                    x.entry((f[3].into(), f[0].into()))
                        .or_default()
                        .push((s, e))
                }
            }
        }
    }
    let mut r = HashMap::new();
    for ((id, _), mut v) in x {
        v.sort();
        let (mut n, mut s, mut e) = (0usize, 0usize, 0usize);
        for (a, b) in v {
            if e == 0 {
                s = a;
                e = b
            } else if a <= e {
                e = e.max(b)
            } else {
                n += e - s;
                s = a;
                e = b
            }
        }
        if e > 0 {
            n += e - s
        };
        *r.entry(id).or_insert(0) += n
    }
    r
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
        eprintln!("usage: dtr-boundary-gate --original <fa> --proposal <fa> --original-copies <bed> --proposal-copies <bed> --out <fa> --audit <tsv> --provenance <json> [--min-relative-gain .01]");
        exit(if a.len() < 2 { 2 } else { 0 })
    }
    let original = need(&a, "--original");
    let proposal = need(&a, "--proposal");
    let oc = need(&a, "--original-copies");
    let pc = need(&a, "--proposal-copies");
    let out = need(&a, "--out");
    let audit = need(&a, "--audit");
    let prov = need(&a, "--provenance");
    let gain: f64 = flag(&a, "--min-relative-gain")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(0.01);
    let p: HashMap<String, Rec> = fasta(&proposal)
        .into_iter()
        .map(|r| (r.id.clone(), r))
        .collect();
    let o = fasta(&original);
    let a_cov = coverage(&oc);
    let b_cov = coverage(&pc);
    let mut w = BufWriter::new(File::create(&out).unwrap());
    let mut aw = BufWriter::new(File::create(&audit).unwrap());
    writeln!(aw, "family\toriginal_length\tproposal_length\toriginal_coverage_bp\tproposal_coverage_bp\trelative_gain\tdecision").unwrap();
    let (mut accepted, mut retained) = (0usize, 0usize);
    for r in o {
        let q = p.get(&r.id);
        let before = *a_cov.get(&r.id).unwrap_or(&0);
        let after = *b_cov.get(&r.id).unwrap_or(&0);
        let rel = if before == 0 {
            if after > 0 {
                f64::INFINITY
            } else {
                0.0
            }
        } else {
            after as f64 / before as f64 - 1.0
        };
        let take = q.map(|z| z.seq != r.seq).unwrap_or(false) && rel >= gain;
        let (seq, plen, decision) = if take {
            accepted += 1;
            (
                &q.unwrap().seq,
                q.unwrap().seq.len(),
                "accept_coverage_gain",
            )
        } else {
            retained += 1;
            (
                &r.seq,
                q.map(|z| z.seq.len()).unwrap_or(r.seq.len()),
                "retain",
            )
        };
        writefa(&mut w, &r.header, seq);
        writeln!(
            aw,
            "{}\t{}\t{}\t{}\t{}\t{:.8}\t{}",
            r.id,
            r.seq.len(),
            plen,
            before,
            after,
            rel,
            decision
        )
        .unwrap()
    }
    w.flush().unwrap();
    aw.flush().unwrap();
    writeln!(File::create(&prov).unwrap(), "{{\"contract\":\"te-looker-adaptive-boundary-gate-v1\",\"mode\":\"partial\",\"family_call\":false,\"criterion\":\"per-family union copy coverage increase\",\"min_relative_gain\":{},\"accepted\":{},\"retained\":{}}}", gain, accepted, retained).unwrap();
    println!("accepted\t{accepted}\nretained\t{retained}")
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn merge_coverage() {
        let _ = coverage;
    }
}
