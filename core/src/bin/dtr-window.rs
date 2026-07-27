//! dtr-window — bounded deterministic spaced-seed candidate windows for Track 2.
use std::collections::{BinaryHeap, HashMap};
use std::env;
use std::fs::File;
use std::io::Write;
use std::process::exit;
use te_core::read_fasta;
const PAT: &str = "11101110111011101111";
fn flag(a: &[String], n: &str) -> Option<String> {
    a.iter()
        .position(|x| x == n)
        .and_then(|i| a.get(i + 1))
        .cloned()
}
fn need(a: &[String], n: &str) -> String {
    flag(a, n).unwrap_or_else(|| {
        eprintln!("[dtr-window] {n} required");
        exit(2)
    })
}
fn mix(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ceb9fe1a85ec53);
    x ^ (x >> 33)
}
fn pos(p: &str) -> Vec<usize> {
    let v: Vec<_> = p
        .bytes()
        .enumerate()
        .filter_map(|(i, b)| match b {
            b'1' => Some(i),
            b'0' => None,
            _ => {
                eprintln!("bad pattern");
                exit(2)
            }
        })
        .collect();
    if v.is_empty() || v.len() > 16 {
        exit(2)
    }
    v
}
fn code(c: &[u8], s: usize, p: &[usize]) -> Option<(u32, bool)> {
    let (mut f, mut r) = (0u32, 0u32);
    for &i in p {
        let b = *c.get(s + i)?;
        if b > 3 {
            return None;
        }
        f = (f << 2) | b as u32
    }
    for &i in p.iter().rev() {
        let b = *c.get(s + i)?;
        r = (r << 2) | (3 - b) as u32
    }
    Some((f.min(r), f <= r))
}
fn rev(mut v: Vec<u8>) -> Vec<u8> {
    for b in &mut v {
        *b = match *b {
            0 => 3,
            1 => 2,
            2 => 1,
            3 => 0,
            _ => 255,
        }
    }
    v.reverse();
    v
}
fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 2 || a[1] == "--help" {
        eprintln!("usage: dtr-window --genome <fa> --out <fa> --provenance <json> [--pattern P --min-count 3 --max-windows 50000 --flank 1000]");
        exit(if a.len() < 2 { 2 } else { 0 })
    }
    let g = need(&a, "--genome");
    let out = need(&a, "--out");
    let prov = need(&a, "--provenance");
    let pat = flag(&a, "--pattern").unwrap_or_else(|| PAT.into());
    let pp = pos(&pat);
    let min = flag(&a, "--min-count")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(3u32);
    let max = flag(&a, "--max-windows")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(50000usize);
    let flank = flag(&a, "--flank")
        .map(|x| x.parse().unwrap_or_else(|_| exit(2)))
        .unwrap_or(1000usize);
    let rec = read_fasta(&g);
    let mut counts = HashMap::new();
    for r in &rec {
        if r.codes.len() >= pat.len() {
            for s in 0..=r.codes.len() - pat.len() {
                if let Some((c, _)) = code(&r.codes, s, &pp) {
                    *counts.entry(c).or_insert(0u32) += 1
                }
            }
        }
    }
    let mut heap: BinaryHeap<(u64, usize, usize, u32, bool)> = BinaryHeap::new();
    let mut eligible = 0u64;
    for (ri, r) in rec.iter().enumerate() {
        if r.codes.len() < pat.len() {
            continue;
        }
        for s in 0..=r.codes.len() - pat.len() {
            if let Some((c, st)) = code(&r.codes, s, &pp) {
                if counts[&c] >= min {
                    eligible += 1;
                    let h = mix((c as u64) ^ ((ri as u64) << 32) ^ s as u64);
                    let x = (h, ri, s, c, st);
                    if heap.len() < max {
                        heap.push(x)
                    } else if h < heap.peek().unwrap().0 {
                        heap.pop();
                        heap.push(x)
                    }
                }
            }
        }
    }
    let mut v = heap.into_vec();
    v.sort_by_key(|x| x.0);
    let mut f = File::create(&out).unwrap();
    for (_, ri, s, c, st) in &v {
        let r = &rec[*ri];
        let lo = s.saturating_sub(flank);
        let hi = (s + pat.len() + flank).min(r.codes.len());
        let seq = if *st {
            r.codes[lo..hi].to_vec()
        } else {
            rev(r.codes[lo..hi].to_vec())
        };
        let dna: String = seq
            .iter()
            .map(|&b| {
                if b < 4 {
                    b"ACGT"[b as usize] as char
                } else {
                    'N'
                }
            })
            .collect();
        writeln!(
            f,
            ">win_{}_{}_{} seed={:08x} source={}:{}-{} strand={}",
            ri,
            s,
            c,
            c,
            r.name,
            lo,
            hi,
            if *st { '+' } else { '-' }
        )
        .unwrap();
        for x in dna.as_bytes().chunks(80) {
            f.write_all(x).unwrap();
            f.write_all(b"\n").unwrap()
        }
    }
    let mut p = File::create(&prov).unwrap();
    writeln!(p,"{{\"contract\":\"te-looker-spaced-window-v1\",\"genome\":\"{}\",\"pattern\":\"{}\",\"min_count\":{},\"flank\":{},\"eligible_occurrences\":{},\"max_windows\":{},\"selected_windows\":{},\"sampling\":\"deterministic_hash_bottom\"}}",g.replace('"',"\\\""),pat,min,flank,eligible,max,v.len()).unwrap();
    println!(
        "eligible_occurrences\t{eligible}\nselected_windows\t{}",
        v.len()
    )
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn canonical_orientation_is_reported() {
        let p = pos("1111");
        assert_eq!(code(&[0, 1, 2, 3], 0, &p), Some((27, true)));
    }
    #[test]
    fn reverse_complement_restores_forward() {
        assert_eq!(rev(vec![3, 2, 1, 0]), vec![3, 2, 1, 0]);
    }
}
