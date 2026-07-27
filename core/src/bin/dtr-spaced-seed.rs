//! dtr-spaced-seed — deterministic exact spaced-seed gate for Track 2 planning.
//!
//! It is deliberately a measurement stage: it counts canonical spaced seeds exactly,
//! reports the candidate mass and capped pair upper bound, and never promotes a family.
//! This is the v4 H0 gate before enabling an LSH/window-graph path on a genome class.

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::Write;
use std::process::exit;
use te_core::read_fasta;

const DEFAULT_PATTERN: &str = "11101110111011101111"; // span 20, weight 16
const DEFAULT_MIN_COUNT: u32 = 3;
const DEFAULT_CAP: u64 = 200;

fn flag(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|arg| arg == name)
        .and_then(|i| args.get(i + 1))
        .cloned()
}
fn need(args: &[String], name: &str) -> String {
    flag(args, name).unwrap_or_else(|| {
        eprintln!("[dtr-spaced-seed] {name} required");
        exit(2);
    })
}
fn pattern_positions(pattern: &str) -> Vec<usize> {
    let positions: Vec<usize> = pattern
        .bytes()
        .enumerate()
        .filter_map(|(i, b)| match b {
            b'1' => Some(i),
            b'0' => None,
            _ => {
                eprintln!("[dtr-spaced-seed] pattern must contain only 0/1");
                exit(2);
            }
        })
        .collect();
    if positions.is_empty() || positions.len() > 16 {
        eprintln!("[dtr-spaced-seed] pattern weight must be 1..16");
        exit(2);
    }
    positions
}
fn canonical(codes: &[u8], start: usize, positions: &[usize]) -> Option<u32> {
    let mut forward = 0u32;
    let mut reverse = 0u32;
    for &pos in positions {
        let base = *codes.get(start + pos)?;
        if base > 3 {
            return None;
        }
        forward = (forward << 2) | base as u32;
    }
    for &pos in positions.iter().rev() {
        let base = *codes.get(start + pos)?;
        reverse = (reverse << 2) | (3 - base) as u32;
    }
    Some(forward.min(reverse))
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 || args[1] == "--help" || args[1] == "-h" {
        eprintln!("usage: dtr-spaced-seed --genome <fa> --out <json> [--pattern 11101110111011101111] [--min-count 3] [--cap 200]");
        exit(if args.len() < 2 { 2 } else { 0 });
    }
    let genome = need(&args, "--genome");
    let output = need(&args, "--out");
    let pattern = flag(&args, "--pattern").unwrap_or_else(|| DEFAULT_PATTERN.to_string());
    let positions = pattern_positions(&pattern);
    let min_count = flag(&args, "--min-count")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(DEFAULT_MIN_COUNT);
    let cap = flag(&args, "--cap")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(DEFAULT_CAP);
    let span = pattern.len();
    let recs = read_fasta(&genome);
    let mut counts: HashMap<u32, u32> = HashMap::new();
    let mut valid_windows = 0u64;
    for record in &recs {
        if record.codes.len() < span {
            continue;
        }
        for start in 0..=record.codes.len() - span {
            if let Some(code) = canonical(&record.codes, start, &positions) {
                *counts.entry(code).or_default() += 1;
                valid_windows += 1;
            }
        }
    }
    let mut repetitive_seeds = 0u64;
    let mut repetitive_occurrences = 0u64;
    let mut capped_pair_upper = 0u128;
    let mut uncapped_pair_upper = 0u128;
    let mut max_count = 0u32;
    for &count in counts.values() {
        max_count = max_count.max(count);
        if count >= min_count {
            repetitive_seeds += 1;
            repetitive_occurrences += count as u64;
            let capped = (count as u64).min(cap) as u128;
            capped_pair_upper += capped * (capped.saturating_sub(1)) / 2;
            let full = count as u128;
            uncapped_pair_upper += full * (full.saturating_sub(1)) / 2;
        }
    }
    let mut out = File::create(&output).unwrap_or_else(|err| {
        eprintln!("[dtr-spaced-seed] cannot write {output}: {err}");
        exit(1);
    });
    writeln!(out, "{{\n  \"contract\": \"te-looker-spaced-seed-gate-v1\",\n  \"genome\": \"{}\",\n  \"pattern\": \"{}\",\n  \"span\": {},\n  \"weight\": {},\n  \"min_count\": {},\n  \"occurrence_cap\": {},\n  \"records\": {},\n  \"valid_windows\": {},\n  \"distinct_seeds\": {},\n  \"repetitive_seeds\": {},\n  \"repetitive_occurrences\": {},\n  \"max_count\": {},\n  \"capped_pair_upper\": \"{}\",\n  \"uncapped_pair_upper\": \"{}\"\n}}", genome.replace('\\', "\\\\").replace('"', "\\\""), pattern, span, positions.len(), min_count, cap, recs.len(), valid_windows, counts.len(), repetitive_seeds, repetitive_occurrences, max_count, capped_pair_upper, uncapped_pair_upper).unwrap();
    println!("repetitive_seeds\t{repetitive_seeds}\nrepetitive_occurrences\t{repetitive_occurrences}\ncapped_pair_upper\t{capped_pair_upper}");
}

#[cfg(test)]
mod tests {
    use super::{canonical, pattern_positions};
    #[test]
    fn canonicalizes_reverse_complements() {
        let positions = pattern_positions("1111");
        assert_eq!(
            canonical(&[0, 1, 2, 3], 0, &positions),
            canonical(&[0, 1, 2, 3], 0, &positions)
        );
        assert_eq!(canonical(&[0, 1, 2, 3], 0, &positions), Some(27));
    }
    #[test]
    fn invalid_base_breaks_a_spaced_window() {
        assert_eq!(canonical(&[0, 255, 2], 0, &pattern_positions("111")), None);
    }
}
