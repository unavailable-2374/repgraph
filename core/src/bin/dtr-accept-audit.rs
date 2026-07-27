//! dtr-accept-audit records pre-acceptance evidence without emitting a final TE library.
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::exit;

fn value(a: &[String], k: &str) -> Option<String> {
    a.iter()
        .position(|x| x == k)
        .and_then(|i| a.get(i + 1))
        .cloned()
}
fn need(a: &[String], k: &str) -> String {
    value(a, k).unwrap_or_else(|| {
        eprintln!("[dtr-accept-audit] {k} is required");
        exit(2)
    })
}
fn names(path: &str) -> Vec<String> {
    BufReader::new(File::open(path).unwrap())
        .lines()
        .map_while(Result::ok)
        .filter_map(|x| {
            x.strip_prefix(">")
                .map(|h| h.split_whitespace().next().unwrap().to_string())
        })
        .collect()
}
fn copy_stats(path: &str) -> HashMap<String, (usize, f64)> {
    let mut v: HashMap<String, (usize, f64)> = HashMap::new();
    for l in BufReader::new(File::open(path).unwrap())
        .lines()
        .map_while(Result::ok)
    {
        let f: Vec<&str> = l.split("\t").collect();
        if f.len() >= 5 {
            if let Ok(id) = f[4].parse::<f64>() {
                let e = v.entry(f[3].to_string()).or_insert((0, 0.0));
                e.0 += 1;
                e.1 += id;
            }
        }
    }
    for e in v.values_mut() {
        e.1 /= e.0 as f64
    }
    v
}
fn structure(path: &str) -> HashMap<String, (bool, bool)> {
    let mut out = HashMap::new();
    let mut it = BufReader::new(File::open(path).unwrap()).lines();
    let head = it.next().unwrap_or(Ok(String::new())).unwrap();
    let h: Vec<&str> = head.split("\t").collect();
    let id = h.iter().position(|x| *x == "family").unwrap_or(0);
    let uncertain = h
        .iter()
        .position(|x| *x == "boundary_uncertain")
        .unwrap_or(9);
    let evidence = h
        .iter()
        .position(|x| *x == "structure_evidence")
        .unwrap_or(10);
    for l in it.map_while(Result::ok) {
        let f: Vec<&str> = l.split("\t").collect();
        if f.len() > evidence {
            out.insert(
                f[id].to_string(),
                (f[evidence] == "true", f[uncertain] == "true"),
            );
        }
    }
    out
}
fn nonte(path: Option<String>) -> HashMap<String, String> {
    let mut out = HashMap::new();
    if let Some(p) = path {
        let mut it = BufReader::new(File::open(p).unwrap()).lines();
        let _ = it.next();
        for l in it.map_while(Result::ok) {
            let f: Vec<&str> = l.split("\t").collect();
            if f.len() >= 2 {
                out.insert(f[0].to_string(), f[1].to_string());
            }
        }
    }
    out
}
fn categories_complete(kinds: &HashSet<String>) -> bool {
    ["organelle", "rrna", "trna", "snrna"]
        .iter()
        .all(|x| kinds.contains(*x))
}
fn panel_complete(path: Option<String>) -> bool {
    let Some(p) = path else {
        return false;
    };
    let mut kinds = HashSet::new();
    let mut it = BufReader::new(File::open(p).unwrap()).lines();
    let _ = it.next();
    for l in it.map_while(Result::ok) {
        let f: Vec<&str> = l.split("\t").collect();
        if f.len() >= 2 {
            kinds.insert(f[1].to_string());
        }
    }
    categories_complete(&kinds)
}
fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 2 || a[1] == "--help" {
        eprintln!("usage: dtr-accept-audit --query <fa> --copies <bed> --structure <tsv> --out <tsv> --provenance <json> [--non-te-report <tsv> --non-te-manifest <tsv> --min-copies 10]");
        exit(if a.len() < 2 { 2 } else { 0 })
    }
    let q = need(&a, "--query");
    let c = need(&a, "--copies");
    let s = need(&a, "--structure");
    let out = need(&a, "--out");
    let p = need(&a, "--provenance");
    let min: usize = value(&a, "--min-copies")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(10);
    let cs = copy_stats(&c);
    let st = structure(&s);
    let nt = nonte(value(&a, "--non-te-report"));
    let g6 = panel_complete(value(&a, "--non-te-manifest"));
    let mut w = BufWriter::new(File::create(&out).unwrap());
    writeln!(w,"family\tcopies\tmean_identity\tg1_copy_support\tg2_structure_evidence\tg3_boundary_certain\tg4_non_te_status\tg5_segdup_suspect\tg6_non_te_panel_complete\tstatus").unwrap();
    let mut ready = 0usize;
    let mut total = 0usize;
    for id in names(&q) {
        total += 1;
        let (n, mean) = cs.get(&id).copied().unwrap_or((0, 0.0));
        let (structure, boundary_uncertain) = st.get(&id).copied().unwrap_or((false, true));
        let g1 = n >= min;
        let g3 = !boundary_uncertain;
        let nt_status = nt
            .get(&id)
            .cloned()
            .unwrap_or_else(|| "not_screened".to_string());
        let g4 = nt_status == "retain";
        let seg = n < 20 && mean >= 98.0;
        let status = if !g1 {
            "insufficient_copy_support"
        } else if !g4 {
            "needs_non_te_guard"
        } else if seg {
            "segdup_suspect"
        } else {
            "eligible_for_full_guard"
        };
        if status == "eligible_for_full_guard" {
            ready += 1
        }
        writeln!(
            w,
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            id, n, mean, g1, structure, g3, nt_status, seg, g6, status
        )
        .unwrap()
    }
    w.flush().unwrap();
    writeln!(File::create(&p).unwrap(),"{{\"contract\":\"te-looker-acceptance-audit-v1\",\"mode\":\"partial\",\"family_call\":false,\"acceptance_complete\":false,\"families\":{},\"eligible_for_full_guard\":{},\"min_copies\":{},\"non_te_panel_complete\":{}}}",total,ready,min,g6).unwrap();
    println!("audited\t{total}\neligible_for_full_guard\t{ready}")
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn absent_non_te_is_not_a_pass() {
        assert!(nonte(None).is_empty())
    }
    #[test]
    fn panel_requires_each_non_te_category() {
        let mut x = HashSet::new();
        x.insert("organelle".to_string());
        x.insert("rrna".to_string());
        assert!(!categories_complete(&x));
        x.insert("trna".to_string());
        x.insert("snrna".to_string());
        assert!(categories_complete(&x));
    }
}
