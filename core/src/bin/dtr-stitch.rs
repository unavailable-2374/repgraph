//! dtr-stitch — reconstruct co-linear discovery fragments into full-length families.
//!
//! Given dtr's refined families and its member BED, this program detects fragment
//! families that recur in a consistent order, rebuilds their spanning occurrences,
//! and delegates gapped consensus construction to the existing `te-refine` binary.
//! It is intentionally conservative: a chain needs repeated, strand-consistent,
//! copy-balanced evidence; failure to reconstruct a chain leaves its source families
//! untouched.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::{exit, Command};

const GAP: i64 = 300;
const WINDOW_K: usize = 4;
const MIN_COOCCURRENCE: usize = 3;
const MIN_COOCCURRENCE_RATIO: f64 = 0.10;
const MIN_DIRECTION_CONSISTENCY: f64 = 0.70;
const MIN_COPY_RATIO: f64 = 0.10;
const MAX_OVERLAP_FRACTION: f64 = 0.50;
const MIN_SPAN: i64 = 80;
const MAX_SPAN: i64 = 25_000;

#[derive(Clone)]
struct FastaRecord {
    id: String,
    header: String,
    seq: String,
}
#[derive(Clone)]
struct Instance {
    chrom: String,
    start: i64,
    end: i64,
    family: String,
    strand: char,
}

fn flag(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|arg| arg == name)
        .and_then(|i| args.get(i + 1))
        .cloned()
}
fn need(args: &[String], name: &str) -> String {
    flag(args, name).unwrap_or_else(|| {
        eprintln!("[dtr-stitch] {name} required");
        exit(2);
    })
}
fn sibling(name: &str) -> String {
    if let Ok(exe) = env::current_exe() {
        if let Some(dir) = exe.parent() {
            let path = dir.join(name);
            if path.exists() {
                return path.to_string_lossy().into_owned();
            }
        }
    }
    name.to_string()
}
fn parse_id(header: &str) -> String {
    let token = header.split_whitespace().next().unwrap_or("");
    token.strip_suffix("_ref").unwrap_or(token).to_string()
}
fn read_fasta(path: &str) -> Vec<FastaRecord> {
    let file = File::open(path).unwrap_or_else(|err| {
        eprintln!("[dtr-stitch] cannot open {path}: {err}");
        exit(1);
    });
    let mut records = Vec::new();
    let mut header = String::new();
    let mut seq = String::new();
    for line in BufReader::new(file).lines() {
        let line = line.unwrap_or_else(|err| {
            eprintln!("[dtr-stitch] read failure: {err}");
            exit(1);
        });
        if let Some(next) = line.strip_prefix('>') {
            if !header.is_empty() {
                records.push(FastaRecord {
                    id: parse_id(&header),
                    header: header.clone(),
                    seq: seq.clone(),
                });
            }
            header = next.to_string();
            seq.clear();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !header.is_empty() {
        records.push(FastaRecord {
            id: parse_id(&header),
            header,
            seq,
        });
    }
    records
}
fn read_bed(path: &str, known: &HashSet<String>) -> Vec<Instance> {
    let file = File::open(path).unwrap_or_else(|err| {
        eprintln!("[dtr-stitch] cannot open {path}: {err}");
        exit(1);
    });
    let mut out = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line.unwrap_or_else(|err| {
            eprintln!("[dtr-stitch] read failure: {err}");
            exit(1);
        });
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 || !known.contains(fields[3]) {
            continue;
        }
        let (start, end) = match (fields[1].parse::<i64>(), fields[2].parse::<i64>()) {
            (Ok(s), Ok(e)) if e > s => (s, e),
            _ => continue,
        };
        let strand = if fields[5] == "-" { '-' } else { '+' };
        out.push(Instance {
            chrom: fields[0].to_string(),
            start,
            end,
            family: fields[3].to_string(),
            strand,
        });
    }
    out.sort_by(|a, b| {
        (&a.chrom, a.start, a.end, &a.family).cmp(&(&b.chrom, b.start, b.end, &b.family))
    });
    out
}
fn chain_families(instances: &[Instance], lengths: &HashMap<String, usize>) -> Vec<Vec<String>> {
    let mut directed: HashMap<(String, String), usize> = HashMap::new();
    for i in 0..instances.len() {
        let left = &instances[i];
        let mut seen = HashSet::new();
        for right in instances
            .iter()
            .take((i + 1 + WINDOW_K).min(instances.len()))
            .skip(i + 1)
        {
            if right.chrom != left.chrom {
                break;
            }
            let gap = right.start - left.end;
            if gap > GAP {
                break;
            }
            if right.strand != left.strand
                || right.family == left.family
                || !seen.insert(right.family.clone())
            {
                continue;
            }
            if gap < 0 {
                let short = lengths
                    .get(&left.family)
                    .copied()
                    .unwrap_or(0)
                    .min(lengths.get(&right.family).copied().unwrap_or(0));
                if short > 0 && (-gap as f64) > short as f64 * MAX_OVERLAP_FRACTION {
                    continue;
                }
            }
            let pair = if left.strand == '+' {
                (left.family.clone(), right.family.clone())
            } else {
                (right.family.clone(), left.family.clone())
            };
            *directed.entry(pair).or_insert(0) += 1;
        }
    }
    let mut undirected: BTreeMap<(String, String), (usize, usize)> = BTreeMap::new();
    for ((a, b), n) in directed {
        let (key, forward) = if a <= b {
            ((a, b), true)
        } else {
            ((b, a), false)
        };
        let entry = undirected.entry(key).or_insert((0, 0));
        if forward {
            entry.0 += n;
        } else {
            entry.1 += n;
        }
    }
    let mut edges = Vec::new();
    let copies: HashMap<String, usize> = instances.iter().fold(HashMap::new(), |mut acc, x| {
        *acc.entry(x.family.clone()).or_insert(0) += 1;
        acc
    });
    for ((a, b), (forward, reverse)) in undirected {
        let total = forward + reverse;
        if total < MIN_COOCCURRENCE {
            continue;
        }
        let ca = copies.get(&a).copied().unwrap_or(1);
        let cb = copies.get(&b).copied().unwrap_or(1);
        if total as f64 / (ca.min(cb).max(1) as f64) < MIN_COOCCURRENCE_RATIO {
            continue;
        }
        if forward.max(reverse) as f64 / (total as f64) < MIN_DIRECTION_CONSISTENCY {
            continue;
        }
        if ca.min(cb) as f64 / (ca.max(cb).max(1) as f64) < MIN_COPY_RATIO {
            continue;
        }
        let (left, right) = if forward >= reverse { (a, b) } else { (b, a) };
        edges.push((left, right, total as f64 / ca.min(cb).max(1) as f64));
    }
    edges.sort_by(|x, y| {
        y.2.partial_cmp(&x.2)
            .unwrap()
            .then_with(|| x.0.cmp(&y.0))
            .then_with(|| x.1.cmp(&y.1))
    });
    let mut right_of: BTreeMap<String, String> = BTreeMap::new();
    let mut left_of: BTreeMap<String, String> = BTreeMap::new();
    for (left, right, _) in edges {
        if !right_of.contains_key(&left) && !left_of.contains_key(&right) {
            right_of.insert(left.clone(), right.clone());
            left_of.insert(right, left);
        }
    }
    let mut chains = Vec::new();
    let mut visited = HashSet::new();
    for head in right_of.keys() {
        if left_of.contains_key(head) || visited.contains(head) {
            continue;
        }
        let mut chain = vec![head.clone()];
        visited.insert(head.clone());
        let mut current = head;
        while let Some(next) = right_of.get(current) {
            if !visited.insert(next.clone()) {
                break;
            }
            chain.push(next.clone());
            current = next;
        }
        if chain.len() >= 2 {
            chains.push(chain);
        }
    }
    chains
}
fn write_spans(
    chains: &[Vec<String>],
    instances: &[Instance],
    path: &str,
) -> BTreeMap<String, Vec<String>> {
    let mut out = BufWriter::new(File::create(path).unwrap_or_else(|err| {
        eprintln!("[dtr-stitch] cannot create {path}: {err}");
        exit(1);
    }));
    let mut mapping = BTreeMap::new();
    for (idx, chain) in chains.iter().enumerate() {
        let id = format!("stitched_{}", idx + 1);
        let set: HashSet<&String> = chain.iter().collect();
        let mut groups: BTreeMap<(String, char), Vec<&Instance>> = BTreeMap::new();
        for instance in instances {
            if set.contains(&instance.family) {
                groups
                    .entry((instance.chrom.clone(), instance.strand))
                    .or_default()
                    .push(instance);
            }
        }
        let mut wrote = false;
        for ((chrom, strand), mut rows) in groups {
            rows.sort_by_key(|x| (x.start, x.end));
            let mut start: Option<i64> = None;
            let mut end = 0i64;
            let mut members: BTreeSet<String> = BTreeSet::new();
            let flush = |start: Option<i64>,
                         end: i64,
                         members: &BTreeSet<String>,
                         out: &mut BufWriter<File>|
             -> bool {
                if let Some(s) = start {
                    if members.len() >= 2 && end - s >= MIN_SPAN && end - s <= MAX_SPAN {
                        writeln!(out, "{chrom}\t{s}\t{end}\t{id}\t.\t{strand}").unwrap();
                        return true;
                    }
                }
                false
            };
            for row in rows {
                match start {
                    None => {
                        start = Some(row.start);
                        end = row.end;
                        members.clear();
                        members.insert(row.family.clone());
                    }
                    Some(_) if row.start - end <= GAP => {
                        end = end.max(row.end);
                        members.insert(row.family.clone());
                    }
                    Some(_) => {
                        wrote |= flush(start, end, &members, &mut out);
                        start = Some(row.start);
                        end = row.end;
                        members.clear();
                        members.insert(row.family.clone());
                    }
                }
            }
            wrote |= flush(start, end, &members, &mut out);
        }
        if wrote {
            mapping.insert(id, chain.clone());
        }
    }
    out.flush().unwrap();
    mapping
}
fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 || args[1] == "--help" || args[1] == "-h" {
        eprintln!("usage: dtr-stitch --genome <fa> --input <families.fa> --members <bed> --out <families.fa> --provenance <json>");
        exit(if args.len() < 2 { 2 } else { 0 });
    }
    let genome = need(&args, "--genome");
    let input = need(&args, "--input");
    let members = need(&args, "--members");
    let output = need(&args, "--out");
    let provenance = need(&args, "--provenance");
    let original = read_fasta(&input);
    let by_id: HashMap<String, FastaRecord> = original
        .iter()
        .cloned()
        .map(|x| (x.id.clone(), x))
        .collect();
    let known: HashSet<String> = by_id.keys().cloned().collect();
    let lengths: HashMap<String, usize> = by_id
        .iter()
        .map(|(id, rec)| (id.clone(), rec.seq.len()))
        .collect();
    let instances = read_bed(&members, &known);
    let chains = chain_families(&instances, &lengths);
    let spans = format!("{output}.members.bed");
    let valid = write_spans(&chains, &instances, &spans);
    let refined = format!("{output}.refined.fa");
    if !valid.is_empty() {
        let status = Command::new(sibling("te-refine"))
            .args([genome.as_str(), spans.as_str(), refined.as_str()])
            .status()
            .unwrap_or_else(|err| {
                eprintln!("[dtr-stitch] failed to start te-refine: {err}");
                exit(3);
            });
        if !status.success() {
            eprintln!("[dtr-stitch] te-refine failed: {status}");
            exit(3);
        }
    } else {
        File::create(&refined).unwrap();
    }
    let stitched = read_fasta(&refined);
    let mut consumed = HashSet::new();
    let mut usable = Vec::new();
    for rec in stitched {
        if let Some(chain) = valid.get(&rec.id) {
            consumed.extend(chain.iter().cloned());
            usable.push(rec);
        }
    }
    let mut out = BufWriter::new(File::create(&output).unwrap_or_else(|err| {
        eprintln!("[dtr-stitch] cannot create {output}: {err}");
        exit(1);
    }));
    for rec in &usable {
        writeln!(out, ">{}", rec.header).unwrap();
        for chunk in rec.seq.as_bytes().chunks(80) {
            writeln!(out, "{}", String::from_utf8_lossy(chunk)).unwrap();
        }
    }
    for rec in &original {
        if !consumed.contains(&rec.id) {
            writeln!(out, ">{}", rec.header).unwrap();
            for chunk in rec.seq.as_bytes().chunks(80) {
                writeln!(out, "{}", String::from_utf8_lossy(chunk)).unwrap();
            }
        }
    }
    out.flush().unwrap();
    let mut prov = BufWriter::new(File::create(&provenance).unwrap());
    writeln!(prov, "{{\n  \"contract\": \"te-looker-fragment-stitch-v1\",\n  \"mode\": \"partial\",\n  \"input_families\": {},\n  \"member_instances\": {},\n  \"chains_detected\": {},\n  \"chains_with_spans\": {},\n  \"stitched_families\": {},\n  \"source_families_consumed\": {},\n  \"output_families\": {}\n}}", original.len(), instances.len(), chains.len(), valid.len(), usable.len(), consumed.len(), usable.len() + original.len() - consumed.len()).unwrap();
    prov.flush().unwrap();
    println!("stitched\t{}", usable.len());
}

#[cfg(test)]
mod tests {
    use super::{chain_families, Instance};
    use std::collections::HashMap;

    fn item(start: i64, end: i64, family: &str, strand: char) -> Instance {
        Instance {
            chrom: "chr1".to_string(),
            start,
            end,
            family: family.to_string(),
            strand,
        }
    }

    fn lengths() -> HashMap<String, usize> {
        [("A".to_string(), 100usize), ("B".to_string(), 100usize)]
            .into_iter()
            .collect()
    }

    #[test]
    fn detects_repeated_plus_strand_collinear_chain() {
        let instances = vec![
            item(0, 100, "A", '+'),
            item(120, 220, "B", '+'),
            item(1000, 1100, "A", '+'),
            item(1120, 1220, "B", '+'),
            item(2000, 2100, "A", '+'),
            item(2120, 2220, "B", '+'),
        ];
        assert_eq!(
            chain_families(&instances, &lengths()),
            vec![vec![String::from("A"), String::from("B")]]
        );
    }

    #[test]
    fn normalizes_minus_strand_genomic_order_to_te_order() {
        let instances = vec![
            item(0, 100, "B", '-'),
            item(120, 220, "A", '-'),
            item(1000, 1100, "B", '-'),
            item(1120, 1220, "A", '-'),
            item(2000, 2100, "B", '-'),
            item(2120, 2220, "A", '-'),
        ];
        assert_eq!(
            chain_families(&instances, &lengths()),
            vec![vec![String::from("A"), String::from("B")]]
        );
    }

    #[test]
    fn rejects_pairs_with_excessive_overlap() {
        let instances = vec![
            item(0, 100, "A", '+'),
            item(0, 20, "B", '+'),
            item(1000, 1100, "A", '+'),
            item(1000, 1020, "B", '+'),
            item(2000, 2100, "A", '+'),
            item(2000, 2020, "B", '+'),
        ];
        assert!(chain_families(&instances, &lengths()).is_empty());
    }
}
