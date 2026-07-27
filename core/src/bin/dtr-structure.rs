//! dtr-structure audits instance-level TSD and conservative consensus terminal signals.
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::exit;

#[derive(Clone)]
struct Rec {
    id: String,
    seq: String,
}
#[derive(Clone)]
struct Copy {
    chrom: String,
    start: usize,
    end: usize,
    family: String,
}
fn value(a: &[String], k: &str) -> Option<String> {
    a.iter()
        .position(|x| x == k)
        .and_then(|i| a.get(i + 1))
        .cloned()
}
fn need(a: &[String], k: &str) -> String {
    value(a, k).unwrap_or_else(|| {
        eprintln!("[dtr-structure] {k} is required");
        exit(2)
    })
}
fn fasta(path: &str) -> Vec<Rec> {
    let mut out = Vec::new();
    let mut h = String::new();
    let mut s = String::new();
    for l in BufReader::new(File::open(path).unwrap()).lines() {
        let l = l.unwrap();
        if let Some(x) = l.strip_prefix(">") {
            if !h.is_empty() {
                out.push(Rec { id: h, seq: s });
                s = String::new();
            }
            h = x.split_whitespace().next().unwrap().to_string();
        } else {
            s.push_str(l.trim());
        }
    }
    if !h.is_empty() {
        out.push(Rec { id: h, seq: s })
    };
    out
}
fn genome(path: &str) -> HashMap<String, String> {
    fasta(path)
        .into_iter()
        .map(|r| (r.id, r.seq.to_ascii_uppercase()))
        .collect()
}
fn copies(path: &str) -> Vec<Copy> {
    BufReader::new(File::open(path).unwrap())
        .lines()
        .map_while(Result::ok)
        .filter_map(|l| {
            let f: Vec<&str> = l.split("\t").collect();
            if f.len() < 4 {
                return None;
            };
            Some(Copy {
                chrom: f[0].to_string(),
                start: f[1].parse().ok()?,
                end: f[2].parse().ok()?,
                family: f[3].to_string(),
            })
        })
        .collect()
}
fn rc(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|x| match *x {
            65 => 84,
            67 => 71,
            71 => 67,
            84 => 65,
            97 => 84,
            99 => 71,
            103 => 67,
            116 => 65,
            _ => 78,
        })
        .collect()
}
fn terminal_seed(a: &[u8], b: &[u8]) -> usize {
    let lim = a.len().min(b.len()).min(80);
    for k in (12..=lim).rev() {
        for i in 0..=a.len() - k {
            if b.windows(k).any(|x| x == &a[i..i + k]) {
                return k;
            }
        }
    }
    0
}
fn tail_run(s: &[u8]) -> usize {
    let mut n = 0;
    if s.is_empty() {
        return 0;
    };
    let z = s[s.len() - 1];
    if z != 65 && z != 84 {
        return 0;
    };
    for x in s.iter().rev().take(100) {
        if *x == z {
            n += 1
        } else {
            break;
        }
    }
    n
}
fn tsd(g: &HashMap<String, String>, c: &[Copy], min_k: usize) -> (usize, usize, usize) {
    let mut valid = 0;
    let mut best_k = 0;
    let mut best_n = 0;
    for k in min_k..=15 {
        let mut n = 0;
        let mut d = 0;
        for x in c {
            if let Some(s) = g.get(&x.chrom) {
                let b = s.as_bytes();
                if x.start >= k && x.end + k <= b.len() {
                    d += 1;
                    if b[x.start - k..x.start] == b[x.end..x.end + k] {
                        n += 1
                    }
                }
            }
        }
        if k == min_k {
            valid = d
        };
        if n > best_n || (n == best_n && k > best_k) {
            best_k = k;
            best_n = n
        }
    }
    (valid, best_k, best_n)
}
fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 2 || a[1] == "--help" {
        eprintln!("usage: dtr-structure --query <fa> --genome <fa> --copies <bed> --out <tsv> --provenance <json> [--min-tsd 3 --min-support 5]");
        exit(if a.len() < 2 { 2 } else { 0 })
    }
    let q = need(&a, "--query");
    let g = need(&a, "--genome");
    let b = need(&a, "--copies");
    let out = need(&a, "--out");
    let p = need(&a, "--provenance");
    let min_k: usize = value(&a, "--min-tsd")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(3);
    let min_support: usize = value(&a, "--min-support")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(5);
    if min_k < 3 {
        eprintln!("[dtr-structure] min-tsd must be >=3; TA is non-informative");
        exit(2)
    }
    let gs = genome(&g);
    let cp = copies(&b);
    let mut by: HashMap<String, Vec<Copy>> = HashMap::new();
    for x in cp {
        by.entry(x.family.clone()).or_default().push(x)
    }
    let mut w = BufWriter::new(File::create(&out).unwrap());
    writeln!(w,"family\tcopies\tvalid_flanks\ttsd_k\ttsd_support\ttsd_fraction\tterminal_direct_seed\tterminal_inverted_seed\tpoly_at_tail\tboundary_uncertain\tstructure_evidence").unwrap();
    let mut evidence = 0usize;
    let mut total = 0usize;
    for r in fasta(&q) {
        total += 1;
        let cs = by.remove(&r.id).unwrap_or_default();
        let (valid, k, n) = tsd(&gs, &cs, min_k);
        let frac = if valid == 0 {
            0.0
        } else {
            n as f64 / valid as f64
        };
        let z = r.seq.to_ascii_uppercase();
        let cap = z.len().min(500);
        let left = &z.as_bytes()[..cap];
        let right = &z.as_bytes()[z.len() - cap..];
        let direct = terminal_seed(left, right);
        let inv = terminal_seed(left, &rc(right));
        let tail = tail_run(z.as_bytes());
        let uncertain = cs.len() < 20;
        let has =
            (n >= min_support && frac >= 0.20) || (direct >= 12) || (inv >= 12) || (tail >= 12);
        if has {
            evidence += 1
        }
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}",
            r.id,
            cs.len(),
            valid,
            k,
            n,
            frac,
            direct,
            inv,
            tail,
            uncertain,
            has
        )
        .unwrap()
    }
    w.flush().unwrap();
    writeln!(File::create(&p).unwrap(),"{{\"contract\":\"te-looker-structure-audit-v1\",\"mode\":\"partial\",\"family_call\":false,\"families\":{},\"with_structure_evidence\":{},\"min_tsd\":{},\"min_support\":{}}}",total,evidence,min_k,min_support).unwrap();
    println!("audited\t{total}\nwith_structure_evidence\t{evidence}")
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rc_is_correct() {
        assert_eq!(String::from_utf8(rc(b"ACGT")).unwrap(), "ACGT")
    }
    #[test]
    fn finds_terminal_seed() {
        assert_eq!(terminal_seed(b"AAAACCCCGGGG", b"TTTAAAACCCCGGGG"), 12)
    }
}
