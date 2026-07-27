//! dtr-community-refine performs deterministic first-level modularity splitting, not a final family call.
use std::collections::{BTreeMap, BTreeSet, HashMap};
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
        eprintln!("[dtr-community-refine] {k} is required");
        exit(2)
    })
}
fn ids(path: &str) -> Vec<String> {
    BufReader::new(File::open(path).unwrap())
        .lines()
        .filter_map(Result::ok)
        .filter_map(|x| {
            x.strip_prefix(">")
                .map(|h| h.split_whitespace().next().unwrap().to_string())
        })
        .collect()
}
fn refine(names: &[String], adj: &[Vec<usize>], gamma: f64) -> Vec<usize> {
    let n = names.len();
    let deg: Vec<f64> = adj.iter().map(|v| v.len() as f64).collect();
    let m2: f64 = deg.iter().sum();
    if m2 == 0.0 {
        return (0..n).collect();
    }
    let mut comm: Vec<usize> = (0..n).collect();
    let mut total = deg.clone();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|x, y| names[*x].cmp(&names[*y]));
    for _ in 0..20 {
        let mut moved = false;
        for &i in &order {
            let old = comm[i];
            let ki = deg[i];
            total[old] -= ki;
            let mut weights: BTreeMap<usize, f64> = BTreeMap::new();
            for &j in &adj[i] {
                *weights.entry(comm[j]).or_insert(0.0) += 1.0
            }
            let current = weights.get(&old).copied().unwrap_or(0.0) - gamma * ki * total[old] / m2;
            let mut best = old;
            let mut best_score = current;
            for (c, kin) in weights {
                let score = kin - gamma * ki * total[c] / m2;
                if score > best_score + 1e-12 || (score - best_score).abs() <= 1e-12 && c < best {
                    best = c;
                    best_score = score
                }
            }
            comm[i] = best;
            total[best] += ki;
            if best != old {
                moved = true
            }
        }
        if !moved {
            break;
        }
    }
    comm
}
fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 2 || a[1] == "--help" {
        eprintln!("usage: dtr-community-refine --nodes <fa> --edges <tsv> --out <membership.tsv> --provenance <json> [--resolution 1.0]");
        exit(if a.len() < 2 { 2 } else { 0 })
    }
    let nodes = need(&a, "--nodes");
    let edges = need(&a, "--edges");
    let out = need(&a, "--out");
    let prov = need(&a, "--provenance");
    let gamma: f64 = value(&a, "--resolution")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(1.0);
    if gamma <= 0.0 {
        eprintln!("[dtr-community-refine] resolution must be positive");
        exit(2)
    }
    let names = ids(&nodes);
    let index: HashMap<String, usize> = names
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, x)| (x, i))
        .collect();
    if index.len() != names.len() {
        eprintln!("[dtr-community-refine] duplicate FASTA identifiers");
        exit(2)
    }
    let mut sets: Vec<BTreeSet<usize>> = (0..names.len()).map(|_| BTreeSet::new()).collect();
    let mut nedges = 0usize;
    for (n, l) in BufReader::new(File::open(edges).unwrap())
        .lines()
        .enumerate()
    {
        let l = l.unwrap();
        if n == 0 && l.starts_with("node_a\t") {
            continue;
        }
        let f: Vec<&str> = l.split("\t").collect();
        if f.len() >= 2 {
            if let (Some(&x), Some(&y)) = (index.get(f[0]), index.get(f[1])) {
                if x != y && sets[x].insert(y) {
                    sets[y].insert(x);
                    nedges += 1
                }
            }
        }
    }
    let adj: Vec<Vec<usize>> = sets.into_iter().map(|x| x.into_iter().collect()).collect();
    let labels = refine(&names, &adj, gamma);
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (i, c) in labels.iter().enumerate() {
        groups.entry(*c).or_default().push(i)
    }
    let mut ordered: Vec<Vec<usize>> = groups.into_values().collect();
    for g in &mut ordered {
        g.sort_by(|x, y| names[*x].cmp(&names[*y]))
    }
    ordered.sort_by(|a, b| names[a[0]].cmp(&names[b[0]]));
    let mut w = BufWriter::new(File::create(&out).unwrap());
    writeln!(w, "node\tcommunity\tcommunity_size\tdegree").unwrap();
    let mut multi = 0usize;
    let mut largest = 0usize;
    for (i, g) in ordered.iter().enumerate() {
        largest = largest.max(g.len());
        if g.len() > 1 {
            multi += 1
        }
        for &x in g {
            writeln!(
                w,
                "{}\tcm_{:06}\t{}\t{}",
                names[x],
                i + 1,
                g.len(),
                adj[x].len()
            )
            .unwrap()
        }
    }
    w.flush().unwrap();
    writeln!(File::create(&prov).unwrap(),"{{\"contract\":\"te-looker-community-refinement-v1\",\"mode\":\"partial\",\"family_call\":false,\"algorithm\":\"deterministic_first_level_louvain_local_moving_not_leiden\",\"nodes\":{},\"edges\":{},\"communities\":{},\"non_singleton_communities\":{},\"largest_community\":{},\"resolution\":{}}}",names.len(),nedges,ordered.len(),multi,largest,gamma).unwrap();
    println!(
        "communities\t{}\nlargest_community\t{}",
        ordered.len(),
        largest
    )
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn separates_two_triangles_with_bridge() {
        let n = vec![
            "a".into(),
            "b".into(),
            "c".into(),
            "d".into(),
            "e".into(),
            "f".into(),
        ];
        let mut g = vec![vec![]; 6];
        for (a, b) in [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (2, 3)] {
            g[a].push(b);
            g[b].push(a)
        }
        let x = refine(&n, &g, 1.0);
        assert_ne!(x[0], x[5])
    }
}
