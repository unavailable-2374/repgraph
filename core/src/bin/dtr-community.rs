//! dtr-community — deterministic connected-component prepartition for Track 2.
//!
//! Components are graph containers, not TE calls. They preserve isolated/low-copy nodes
//! for later evidence, and provide a stable bounded input for a future Leiden split.

use std::collections::{BTreeMap, HashMap};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::exit;
fn flag(a: &[String], n: &str) -> Option<String> {
    a.iter()
        .position(|x| x == n)
        .and_then(|i| a.get(i + 1))
        .cloned()
}
fn need(a: &[String], n: &str) -> String {
    flag(a, n).unwrap_or_else(|| {
        eprintln!("[dtr-community] {n} required");
        exit(2)
    })
}
fn fasta_ids(p: &str) -> Vec<String> {
    let f = File::open(p).unwrap_or_else(|e| {
        eprintln!("[dtr-community] {e}");
        exit(1)
    });
    BufReader::new(f)
        .lines()
        .map_while(Result::ok)
        .filter_map(|l| {
            l.strip_prefix('>')
                .map(|h| h.split_whitespace().next().unwrap_or("").to_string())
        })
        .collect()
}
fn find(p: &mut [usize], x: usize) -> usize {
    if p[x] != x {
        let r = find(p, p[x]);
        p[x] = r
    }
    p[x]
}
fn join(p: &mut [usize], a: usize, b: usize) {
    let (a, b) = (find(p, a), find(p, b));
    if a != b {
        p[b] = a
    }
}
fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 2 || a[1] == "--help" {
        eprintln!("usage: dtr-community --nodes <candidates.fa> --edges <lsh.tsv> --out <membership.tsv> --provenance <json>");
        exit(if a.len() < 2 { 2 } else { 0 })
    }
    let nodes = need(&a, "--nodes");
    let edges = need(&a, "--edges");
    let out = need(&a, "--out");
    let prov = need(&a, "--provenance");
    let ids = fasta_ids(&nodes);
    let index: HashMap<String, usize> = ids
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, x)| (x, i))
        .collect();
    if index.len() != ids.len() {
        eprintln!("[dtr-community] duplicate FASTA identifiers");
        exit(2)
    }
    let mut parent: Vec<usize> = (0..ids.len()).collect();
    let mut degree = vec![0usize; ids.len()];
    let mut accepted = 0usize;
    for (n, line) in BufReader::new(File::open(&edges).unwrap())
        .lines()
        .enumerate()
    {
        let l = line.unwrap();
        if n == 0 && l.starts_with("node_a\t") {
            continue;
        }
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 2 {
            continue;
        }
        if let (Some(&x), Some(&y)) = (index.get(f[0]), index.get(f[1])) {
            join(&mut parent, x, y);
            degree[x] += 1;
            degree[y] += 1;
            accepted += 1
        }
    }
    let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for i in 0..ids.len() {
        let r = find(&mut parent, i);
        groups.entry(ids[r].clone()).or_default().push(i)
    }
    let mut ordered: Vec<Vec<usize>> = groups.into_values().collect();
    for g in &mut ordered {
        g.sort_by(|&x, &y| ids[x].cmp(&ids[y]))
    }
    ordered.sort_by(|x, y| ids[x[0]].cmp(&ids[y[0]]));
    let mut w = BufWriter::new(File::create(&out).unwrap());
    writeln!(w, "node\tcomponent\tcomponent_size\tdegree").unwrap();
    let mut nonsingle = 0usize;
    let mut largest = 0usize;
    for (i, g) in ordered.iter().enumerate() {
        largest = largest.max(g.len());
        if g.len() > 1 {
            nonsingle += 1
        }
        for &n in g {
            writeln!(w, "{}\tcc_{:06}\t{}\t{}", ids[n], i + 1, g.len(), degree[n]).unwrap()
        }
    }
    w.flush().unwrap();
    let mut p = File::create(&prov).unwrap();
    writeln!(p,"{{\"contract\":\"te-looker-community-prepartition-v1\",\"nodes\":{},\"edges\":{},\"components\":{},\"non_singleton_components\":{},\"largest_component\":{},\"algorithm\":\"deterministic_connected_components\",\"family_call\":false}}",ids.len(),accepted,ordered.len(),nonsingle,largest).unwrap();
    println!(
        "components\t{}\nlargest_component\t{}",
        ordered.len(),
        largest
    )
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn union_find_merges_transitive_edges() {
        let mut p = (0..4).collect::<Vec<_>>();
        join(&mut p, 0, 1);
        join(&mut p, 1, 2);
        assert_eq!(find(&mut p, 0), find(&mut p, 2));
        assert_ne!(find(&mut p, 0), find(&mut p, 3));
    }
}
