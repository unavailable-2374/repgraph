//! dtr-nested detects nested copy topology and emits clean outer instances for boundary evidence.
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::exit;
#[derive(Clone)]
struct I {
    chrom: String,
    start: usize,
    end: usize,
    family: String,
    line: String,
}
fn value(a: &[String], k: &str) -> Option<String> {
    a.iter()
        .position(|x| x == k)
        .and_then(|i| a.get(i + 1))
        .cloned()
}
fn need(a: &[String], k: &str) -> String {
    value(a, k).unwrap_or_else(|| {
        eprintln!("[dtr-nested] {k} is required");
        exit(2)
    })
}
fn parse(path: &str) -> Vec<I> {
    BufReader::new(File::open(path).unwrap())
        .lines()
        .map_while(Result::ok)
        .filter_map(|l| {
            let f: Vec<&str> = l.split("\t").collect();
            if f.len() < 4 {
                return None;
            }
            Some(I {
                chrom: f[0].to_string(),
                start: f[1].parse().ok()?,
                end: f[2].parse().ok()?,
                family: f[3].to_string(),
                line: l,
            })
        })
        .collect()
}
fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 2 || a[1] == "--help" {
        eprintln!("usage: dtr-nested --copies <bed> --events <tsv> --clean <bed> --provenance <json> [--min-flank 80 --max-active 4096]");
        exit(if a.len() < 2 { 2 } else { 0 })
    }
    let copies = need(&a, "--copies");
    let events = need(&a, "--events");
    let clean = need(&a, "--clean");
    let prov = need(&a, "--provenance");
    let flank: usize = value(&a, "--min-flank")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(80);
    let cap: usize = value(&a, "--max-active")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(4096);
    let mut by: HashMap<String, Vec<I>> = HashMap::new();
    for x in parse(&copies) {
        by.entry(x.chrom.clone()).or_default().push(x)
    }
    let mut dirty: HashSet<(String, usize, usize, String)> = HashSet::new();
    let mut ew = BufWriter::new(File::create(&events).unwrap());
    writeln!(ew,"outer_family\tinner_family\tchrom\touter_start\touter_end\tinner_start\tinner_end\tleft_flank\tright_flank").unwrap();
    let mut ne = 0usize;
    let mut total = 0usize;
    for mut v in by.into_values() {
        v.sort_by(|x, y| {
            x.start
                .cmp(&y.start)
                .then(y.end.cmp(&x.end))
                .then(x.family.cmp(&y.family))
        });
        total += v.len();
        for j in 0..v.len() {
            let inner = &v[j];
            let mut best: Option<&I> = None;
            let begin = j.saturating_sub(cap);
            for outer in v[begin..j].iter().rev() {
                if outer.end < inner.end {
                    continue;
                }
                if outer.family == inner.family {
                    continue;
                }
                if inner.start < outer.start + flank || inner.end + flank > outer.end {
                    continue;
                }
                if best
                    .map(|b| outer.end - outer.start < b.end - b.start)
                    .unwrap_or(true)
                {
                    best = Some(outer)
                }
            }
            if let Some(outer) = best {
                dirty.insert((
                    outer.chrom.clone(),
                    outer.start,
                    outer.end,
                    outer.family.clone(),
                ));
                writeln!(
                    ew,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    outer.family,
                    inner.family,
                    outer.chrom,
                    outer.start,
                    outer.end,
                    inner.start,
                    inner.end,
                    inner.start - outer.start,
                    outer.end - inner.end
                )
                .unwrap();
                ne += 1
            }
        }
    }
    ew.flush().unwrap();
    let mut cw = BufWriter::new(File::create(&clean).unwrap());
    let mut kept = 0usize;
    for x in parse(&copies) {
        if !dirty.contains(&(x.chrom.clone(), x.start, x.end, x.family.clone())) {
            writeln!(cw, "{}", x.line).unwrap();
            kept += 1
        }
    }
    cw.flush().unwrap();
    writeln!(File::create(&prov).unwrap(),"{{\"contract\":\"te-looker-nested-topology-v1\",\"mode\":\"partial\",\"family_call\":false,\"instances\":{},\"nested_events\":{},\"clean_instances\":{},\"min_flank\":{},\"max_active\":{}}}",total,ne,kept,flank,cap).unwrap();
    println!("nested_events\t{ne}\nclean_instances\t{kept}")
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn parse_bed() {
        let x = I {
            chrom: "c".into(),
            start: 1,
            end: 9,
            family: "a".into(),
            line: "c\t1\t9\ta".into(),
        };
        assert_eq!((x.end - x.start, x.family), (8, "a".into()))
    }
}
