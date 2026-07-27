//! dtr-copy — full-genome copy catalog for provisional Track 2 consensus sequences.
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::{exit, Command};
#[derive(Clone)]
struct H {
    q: String,
    c: String,
    st: char,
    ql: usize,
    qs: usize,
    qe: usize,
    ss: usize,
    se: usize,
    id: f64,
    len: usize,
}
fn flag(a: &[String], n: &str) -> Option<String> {
    a.iter()
        .position(|x| x == n)
        .and_then(|i| a.get(i + 1))
        .cloned()
}
fn need(a: &[String], n: &str) -> String {
    flag(a, n).unwrap_or_else(|| {
        eprintln!("[dtr-copy] {n} required");
        exit(2)
    })
}
fn parse(s: &str) -> Option<H> {
    let f: Vec<_> = s.split('\t').collect();
    if f.len() < 10 {
        return None;
    }
    let (qs, qe, ss, se): (usize, usize, usize, usize) = (
        f[5].parse().ok()?,
        f[6].parse().ok()?,
        f[7].parse().ok()?,
        f[8].parse().ok()?,
    );
    Some(H {
        q: f[0].into(),
        c: f[1].into(),
        st: if se >= ss { '+' } else { '-' },
        ql: f[4].parse().ok()?,
        qs: qs.min(qe),
        qe: qs.max(qe),
        ss: ss.min(se),
        se: ss.max(se),
        id: f[2].parse().ok()?,
        len: f[3].parse().ok()?,
    })
}
fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 2 || a[1] == "--help" {
        eprintln!("usage: dtr-copy --genome <fa> --query <fa> --out <members.bed> --provenance <json> [--min-identity 75 --min-qcov .6 --threads 4]");
        exit(if a.len() < 2 { 2 } else { 0 })
    }
    let g = need(&a, "--genome");
    let q = need(&a, "--query");
    let out = need(&a, "--out");
    let prov = need(&a, "--provenance");
    let mi = flag(&a, "--min-identity")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(75f64);
    let mc = flag(&a, "--min-qcov")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(0.6f64);
    let th = flag(&a, "--threads").unwrap_or_else(|| "1".into());
    let db = format!("{out}.genome_db");
    let raw = format!("{out}.raw_hsps.tsv");
    let st = Command::new("makeblastdb")
        .args(["-in", &g, "-dbtype", "nucl", "-out", &db])
        .status()
        .unwrap();
    if !st.success() {
        exit(3)
    }
    let st = Command::new("blastn")
        .args([
            "-task",
            "blastn",
            "-query",
            &q,
            "-db",
            &db,
            "-dust",
            "no",
            "-soft_masking",
            "false",
            "-evalue",
            "1e-10",
            "-num_threads",
            &th,
            "-max_target_seqs",
            "5000",
            "-outfmt",
            "6 qseqid sseqid pident length qlen qstart qend sstart send bitscore",
            "-out",
            &raw,
        ])
        .status()
        .unwrap();
    if !st.success() {
        exit(3)
    }
    let mut m: HashMap<(String, String, char, usize), Vec<H>> = HashMap::new();
    let mut nraw = 0;
    for l in BufReader::new(File::open(&raw).unwrap()).lines() {
        if let Some(h) = parse(&l.unwrap()) {
            nraw += 1;
            m.entry((h.q.clone(), h.c.clone(), h.st, h.ql))
                .or_default()
                .push(h)
        }
    }
    let mut rows = Vec::new();
    for ((q, c, st, ql), mut v) in m {
        v.sort_by_key(|h| (h.qs, h.ss));
        let mut b: Vec<H> = Vec::new();
        for h in v {
            if let Some(x) = b.last_mut() {
                if h.qs <= x.qe + 200 && h.ss <= x.se + 500 {
                    x.qe = x.qe.max(h.qe);
                    x.se = x.se.max(h.se);
                    x.id = (x.id * x.len as f64 + h.id * h.len as f64) / (x.len + h.len) as f64;
                    x.len += h.len;
                    continue;
                }
            }
            b.push(h)
        }
        for h in b {
            let cov = (h.qe - h.qs + 1) as f64 / ql as f64;
            if cov >= mc && h.id >= mi {
                rows.push((q.clone(), c.clone(), st, h.ss, h.se, cov, h.id))
            }
        }
    }
    rows.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then(a.1.cmp(&b.1))
            .then(a.2.cmp(&b.2))
            .then(a.3.cmp(&b.3))
    });
    let mut w = BufWriter::new(File::create(&out).unwrap());
    let mut kept = 0;
    let mut last: (String, String, char, usize) = Default::default();
    for r in rows {
        if r.0 == last.0 && r.1 == last.1 && r.2 == last.2 && r.3 <= last.3 + 100 {
            continue;
        }
        last = (r.0.clone(), r.1.clone(), r.2, r.4);
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:.3}\t{}",
            r.1,
            r.3 - 1,
            r.4,
            r.0,
            r.6,
            r.2
        )
        .unwrap();
        kept += 1
    }
    w.flush().unwrap();
    writeln!(File::create(&prov).unwrap(),"{{\"contract\":\"te-looker-copy-catalog-v1\",\"mode\":\"partial\",\"raw_hsps\":{},\"instances\":{},\"min_identity\":{},\"min_qcov\":{}}}",nraw,kept,mi,mc).unwrap();
    println!("instances\t{kept}")
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn parses_reverse_hsp() {
        let h = parse("q\tc\t90\t100\t200\t2\t101\t300\t201\t9").unwrap();
        assert_eq!((h.st, h.ss, h.se), ('-', 201, 300));
    }
}
