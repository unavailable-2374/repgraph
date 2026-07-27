//! dtr-component-poa — provisional POA for bounded Track 2 graph components.
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::{exit, Command};
fn flag(a: &[String], n: &str) -> Option<String> {
    a.iter()
        .position(|x| x == n)
        .and_then(|i| a.get(i + 1))
        .cloned()
}
fn need(a: &[String], n: &str) -> String {
    flag(a, n).unwrap_or_else(|| {
        eprintln!("[dtr-component-poa] {n} required");
        exit(2)
    })
}
fn fasta(p: &str) -> HashMap<String, String> {
    let f = File::open(p).unwrap();
    let (mut id, mut seq) = (String::new(), String::new());
    let mut r = HashMap::new();
    for l in BufReader::new(f).lines() {
        let l = l.unwrap();
        if let Some(h) = l.strip_prefix('>') {
            if !id.is_empty() {
                r.insert(id, seq);
                seq = String::new()
            }
            id = h.split_whitespace().next().unwrap_or("").to_string()
        } else {
            seq.push_str(l.trim())
        }
    }
    if !id.is_empty() {
        r.insert(id, seq);
    }
    r
}
fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 2 || a[1] == "--help" {
        eprintln!("usage: dtr-component-poa --input <windows.fa> --membership <tsv> --out <fa> --provenance <json> [--min-members 5 --max-members 30]");
        exit(if a.len() < 2 { 2 } else { 0 })
    }
    let input = need(&a, "--input");
    let mem = need(&a, "--membership");
    let out = need(&a, "--out");
    let prov = need(&a, "--provenance");
    let min = flag(&a, "--min-members")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(5usize);
    let max = flag(&a, "--max-members")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(30usize);
    let seqs = fasta(&input);
    let mut groups: HashMap<String, Vec<String>> = HashMap::new();
    for (n, l) in BufReader::new(File::open(&mem).unwrap())
        .lines()
        .enumerate()
    {
        let l = l.unwrap();
        if n == 0 {
            continue;
        }
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() >= 2 {
            groups.entry(f[1].into()).or_default().push(f[0].into())
        }
    }
    let mut keys: Vec<_> = groups.keys().cloned().collect();
    keys.sort();
    let tmp = format!("{}.spoa_in.fa", out);
    let mut w = BufWriter::new(File::create(&out).unwrap());
    let (mut eligible, mut refined) = (0, 0);
    for k in keys {
        let mut ids = groups.remove(&k).unwrap();
        ids.sort();
        if ids.len() < min {
            continue;
        }
        eligible += 1;
        if ids.len() > max {
            ids.truncate(max)
        }
        let mut t = BufWriter::new(File::create(&tmp).unwrap());
        for id in &ids {
            if let Some(s) = seqs.get(id) {
                writeln!(t, ">{id}").unwrap();
                writeln!(t, "{s}").unwrap()
            }
        }
        t.flush().unwrap();
        let o = Command::new(std::env::var("SPOA_BIN").unwrap_or_else(|_| "spoa".into()))
            .args(["-r", "0", &tmp])
            .output()
            .unwrap_or_else(|e| {
                eprintln!("[dtr-component-poa] spoa: {e}");
                exit(3)
            });
        if !o.status.success() {
            exit(3)
        }
        let c: String = String::from_utf8_lossy(&o.stdout)
            .lines()
            .filter(|x| !x.starts_with('>'))
            .collect();
        if c.len() >= 80 {
            refined += 1;
            writeln!(
                w,
                ">track2_{}_partial members={} len={}",
                k,
                ids.len(),
                c.len()
            )
            .unwrap();
            for x in c.as_bytes().chunks(80) {
                w.write_all(x).unwrap();
                w.write_all(b"\n").unwrap()
            }
        }
    }
    w.flush().unwrap();
    let _ = std::fs::remove_file(&tmp);
    let mut p = File::create(&prov).unwrap();
    writeln!(p,"{{\"contract\":\"te-looker-component-poa-v1\",\"mode\":\"partial\",\"min_members\":{},\"max_members\":{},\"eligible_components\":{},\"refined_consensi\":{},\"family_call\":false}}",min,max,eligible,refined).unwrap();
    println!("refined_consensi\t{refined}")
}
