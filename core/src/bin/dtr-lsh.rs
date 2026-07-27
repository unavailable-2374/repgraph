//! dtr-lsh — bounded deterministic MinHash/LSH candidate graph for Track 2.
//!
//! Input records are candidate windows. This stage only proposes and verifies sparse
//! similarity edges; it never calls a component a TE family. Bucket and degree budgets
//! make its memory/time O(N * bands + retained_edges), not O(N^2).

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::exit;

const DEFAULT_K: usize = 15;
const DEFAULT_SKETCH: usize = 128;
const DEFAULT_BANDS: usize = 8;
const DEFAULT_BUCKET_CAP: usize = 128;
const DEFAULT_MAX_NEIGHBORS: usize = 16;
const DEFAULT_MIN_JACCARD: f64 = 0.08;

#[derive(Clone)]
struct Rec {
    id: String,
    seq: Vec<u8>,
}
fn flag(a: &[String], n: &str) -> Option<String> {
    a.iter()
        .position(|x| x == n)
        .and_then(|i| a.get(i + 1))
        .cloned()
}
fn need(a: &[String], n: &str) -> String {
    flag(a, n).unwrap_or_else(|| {
        eprintln!("[dtr-lsh] {n} required");
        exit(2)
    })
}
fn mix64(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ceb9fe1a85ec53);
    x ^ (x >> 33)
}
fn base(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}
fn read_fasta(p: &str) -> Vec<Rec> {
    let f = File::open(p).unwrap_or_else(|e| {
        eprintln!("[dtr-lsh] {e}");
        exit(1)
    });
    let mut out = Vec::new();
    let (mut id, mut seq) = (String::new(), Vec::new());
    for l in BufReader::new(f).lines() {
        let l = l.unwrap();
        if let Some(h) = l.strip_prefix('>') {
            if !id.is_empty() {
                out.push(Rec { id, seq });
                seq = Vec::new()
            }
            id = h.split_whitespace().next().unwrap_or("").to_string()
        } else {
            seq.extend(l.trim().bytes())
        }
    }
    if !id.is_empty() {
        out.push(Rec { id, seq })
    }
    out
}
fn sketch(seq: &[u8], k: usize, size: usize) -> Vec<u64> {
    let mut set = BTreeSet::new();
    if seq.len() < k {
        return Vec::new();
    }
    for i in 0..=seq.len() - k {
        let (mut f, mut r) = (0u64, 0u64);
        let mut ok = true;
        for j in 0..k {
            if let Some(x) = base(seq[i + j]) {
                f = (f << 2) | x
            } else {
                ok = false;
                break;
            }
        }
        if !ok {
            continue;
        }
        for j in (0..k).rev() {
            let x = base(seq[i + j]).unwrap();
            r = (r << 2) | (3 - x)
        }
        set.insert(mix64(f.min(r)));
    }
    set.into_iter().take(size).collect()
}
fn jac(a: &[u64], b: &[u64]) -> f64 {
    let (mut i, mut j, mut same) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            same += 1;
            i += 1;
            j += 1
        } else if a[i] < b[j] {
            i += 1
        } else {
            j += 1
        }
    }
    let u = a.len() + b.len() - same;
    if u == 0 {
        0.0
    } else {
        same as f64 / u as f64
    }
}
fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 2 || a[1] == "--help" {
        eprintln!("usage: dtr-lsh --input <candidates.fa> --edges <tsv> --provenance <json> [--k 15 --sketch 128 --bands 8 --bucket-cap 128 --max-neighbors 16 --min-jaccard .08]");
        exit(if a.len() < 2 { 2 } else { 0 })
    }
    let input = need(&a, "--input");
    let edges = need(&a, "--edges");
    let prov = need(&a, "--provenance");
    let k = flag(&a, "--k")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(DEFAULT_K);
    let ss = flag(&a, "--sketch")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(DEFAULT_SKETCH);
    let bands = flag(&a, "--bands")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(DEFAULT_BANDS);
    let cap = flag(&a, "--bucket-cap")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(DEFAULT_BUCKET_CAP);
    let degree = flag(&a, "--max-neighbors")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(DEFAULT_MAX_NEIGHBORS);
    let minj = flag(&a, "--min-jaccard")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(DEFAULT_MIN_JACCARD);
    if bands == 0 || ss < bands {
        eprintln!("[dtr-lsh] sketch must be >= bands > 0");
        exit(2)
    }
    let recs = read_fasta(&input);
    let sig: Vec<Vec<u64>> = recs.iter().map(|r| sketch(&r.seq, k, ss)).collect();
    let width = ss / bands;
    let mut buckets: BTreeMap<(usize, u64), Vec<usize>> = BTreeMap::new();
    for (i, s) in sig.iter().enumerate() {
        // Short TE candidates (notably MITEs/SINEs) may have fewer distinct k-mers
        // than the nominal sketch size. Use every complete band they do have instead
        // of silently excluding them from Track 2.
        for b in 0..(s.len() / width) {
            let mut h = 0u64;
            for &x in &s[b * width..(b + 1) * width] {
                h = mix64(h ^ x)
            }
            buckets.entry((b, h)).or_default().push(i)
        }
    }
    let mut pairs = HashSet::new();
    let mut oversized = 0usize;
    for (_, mut v) in buckets {
        if v.len() > cap {
            oversized += 1;
            v.sort_by_key(|&i| mix64(i as u64));
            v.truncate(cap)
        }
        for x in 0..v.len() {
            for y in x + 1..v.len() {
                pairs.insert((v[x].min(v[y]), v[x].max(v[y])));
            }
        }
    }
    let mut scored: Vec<(usize, usize, f64)> = pairs
        .into_iter()
        .filter_map(|(i, j)| {
            let s = jac(&sig[i], &sig[j]);
            if s >= minj {
                Some((i, j, s))
            } else {
                None
            }
        })
        .collect();
    scored.sort_by(|a, b| {
        b.2.partial_cmp(&a.2)
            .unwrap()
            .then_with(|| recs[a.0].id.cmp(&recs[b.0].id))
            .then_with(|| recs[a.1].id.cmp(&recs[b.1].id))
    });
    let mut used = vec![0usize; recs.len()];
    let mut kept = Vec::new();
    for e in scored {
        if used[e.0] < degree && used[e.1] < degree {
            used[e.0] += 1;
            used[e.1] += 1;
            kept.push(e)
        }
    }
    let mut out = BufWriter::new(File::create(&edges).unwrap());
    writeln!(out, "node_a\tnode_b\tjaccard").unwrap();
    for (i, j, s) in &kept {
        writeln!(out, "{}\t{}\t{:.6}", recs[*i].id, recs[*j].id, s).unwrap()
    }
    out.flush().unwrap();
    let mut p = File::create(&prov).unwrap();
    writeln!(p,"{{\"contract\":\"te-looker-lsh-graph-v1\",\"input\":\"{}\",\"nodes\":{},\"k\":{},\"sketch_size\":{},\"bands\":{},\"bucket_cap\":{},\"oversized_buckets\":{},\"max_neighbors\":{},\"min_jaccard\":{},\"edges\":{}}}",input.replace('"',"\\\""),recs.len(),k,ss,bands,cap,oversized,degree,minj,kept.len()).unwrap();
    println!("nodes\t{}\nedges\t{}", recs.len(), kept.len())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sketch_similarity_separates_related_and_unrelated() {
        let a = sketch(b"ACGTTGCAACGTTGCAACGTTGCAACGTTGCA", 5, 32);
        let b = sketch(b"ACGTTGCAACGTTGCAACGTTGCAACGTTGCG", 5, 32);
        let c = sketch(b"TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG", 5, 32);
        assert!(jac(&a, &b) > jac(&a, &c));
    }
    #[test]
    fn canonical_sketch_is_strand_invariant() {
        assert_eq!(sketch(b"ACGTTGCA", 4, 32), sketch(b"TGCAACGT", 4, 32));
    }
}
