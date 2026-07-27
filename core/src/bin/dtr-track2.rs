//! dtr-track2 orchestrates the bounded Track 2 workflow and never emits a final family library.
use std::env;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{exit, Command};

fn value(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .position(|x| x == key)
        .and_then(|i| args.get(i + 1))
        .cloned()
}
fn required(args: &[String], key: &str) -> String {
    value(args, key).unwrap_or_else(|| {
        eprintln!("[dtr-track2] {key} is required");
        exit(2)
    })
}
fn sibling(name: &str) -> PathBuf {
    env::current_exe().unwrap().parent().unwrap().join(name)
}
fn run(name: &str, args: &[String]) {
    let status = Command::new(sibling(name))
        .args(args)
        .status()
        .unwrap_or_else(|e| {
            eprintln!("[dtr-track2] cannot execute {name}: {e}");
            exit(3)
        });
    if !status.success() {
        eprintln!("[dtr-track2] {name} failed");
        exit(3);
    }
}
fn fasta_count(path: &Path) -> usize {
    std::io::BufRead::lines(std::io::BufReader::new(File::open(path).unwrap()))
        .map_while(Result::ok)
        .filter(|line| line.starts_with(">"))
        .count()
}
fn json(s: &str) -> String {
    s.replace(char::from(92), "\\\\").replace("\"", "\\\"")
}

fn main() {
    let a: Vec<String> = env::args().collect();
    if a.len() < 2 || a[1] == "--help" {
        eprintln!("usage: dtr-track2 --genome <fa> --out <dir> [--non-te-reference <fa> --non-te-manifest <tsv> --max-windows 50000 --threads 1]");
        exit(if a.len() < 2 { 2 } else { 0 });
    }
    let genome = required(&a, "--genome");
    let out = PathBuf::from(required(&a, "--out"));
    let threads = value(&a, "--threads").unwrap_or_else(|| "1".to_string());
    let max_windows = value(&a, "--max-windows").unwrap_or_else(|| "50000".to_string());
    let non_te = value(&a, "--non-te-reference");
    let non_te_manifest = value(&a, "--non-te-manifest");
    fs::create_dir_all(&out).unwrap();
    let windows = out.join("windows.fa");
    let window_json = out.join("window_provenance.json");
    run(
        "dtr-window",
        &[
            "--genome".into(),
            genome.clone(),
            "--out".into(),
            windows.display().to_string(),
            "--provenance".into(),
            window_json.display().to_string(),
            "--max-windows".into(),
            max_windows,
        ],
    );
    let edges = out.join("lsh_edges.tsv");
    let lsh_json = out.join("lsh_provenance.json");
    run(
        "dtr-lsh",
        &[
            "--input".into(),
            windows.display().to_string(),
            "--edges".into(),
            edges.display().to_string(),
            "--provenance".into(),
            lsh_json.display().to_string(),
        ],
    );
    let members = out.join("components.tsv");
    let community_json = out.join("community_provenance.json");
    run(
        "dtr-community",
        &[
            "--nodes".into(),
            windows.display().to_string(),
            "--edges".into(),
            edges.display().to_string(),
            "--out".into(),
            members.display().to_string(),
            "--provenance".into(),
            community_json.display().to_string(),
        ],
    );
    let communities = out.join("communities.tsv");
    let community_refine_json = out.join("community_refinement_provenance.json");
    run(
        "dtr-community-refine",
        &[
            "--nodes".into(),
            windows.display().to_string(),
            "--edges".into(),
            edges.display().to_string(),
            "--out".into(),
            communities.display().to_string(),
            "--provenance".into(),
            community_refine_json.display().to_string(),
        ],
    );
    let provisional = out.join("provisional.fa");
    let poa_json = out.join("poa_provenance.json");
    run(
        "dtr-component-poa",
        &[
            "--input".into(),
            windows.display().to_string(),
            "--membership".into(),
            communities.display().to_string(),
            "--out".into(),
            provisional.display().to_string(),
            "--provenance".into(),
            poa_json.display().to_string(),
        ],
    );
    let query = if let Some(reference) = non_te.as_ref() {
        let retained = out.join("non_te_retained.fa");
        run(
            "dtr-nonte-guard",
            &[
                "--query".into(),
                provisional.display().to_string(),
                "--reference".into(),
                reference.clone(),
                "--out-retained".into(),
                retained.display().to_string(),
                "--out-rejected".into(),
                out.join("non_te_rejected.fa").display().to_string(),
                "--report".into(),
                out.join("non_te_report.tsv").display().to_string(),
                "--provenance".into(),
                out.join("non_te_provenance.json").display().to_string(),
                "--threads".into(),
                threads.clone(),
            ],
        );
        retained
    } else {
        provisional.clone()
    };
    let copies = out.join("copy_instances.bed");
    let copy_json = out.join("copy_provenance.json");
    if fasta_count(&query) > 0 {
        run(
            "dtr-copy",
            &[
                "--genome".into(),
                genome.clone(),
                "--query".into(),
                query.display().to_string(),
                "--out".into(),
                copies.display().to_string(),
                "--provenance".into(),
                copy_json.display().to_string(),
                "--threads".into(),
                threads,
            ],
        );
    } else {
        File::create(&copies).unwrap();
        writeln!(File::create(&copy_json).unwrap(), "{{\"contract\":\"te-looker-copy-catalog-v1\",\"mode\":\"partial\",\"instances\":0,\"skipped_empty_query\":true}}").unwrap();
    }
    let structure = out.join("structure_audit.tsv");
    let structure_json = out.join("structure_provenance.json");
    if fasta_count(&query) > 0 {
        run(
            "dtr-structure",
            &[
                "--query".into(),
                query.display().to_string(),
                "--genome".into(),
                genome.clone(),
                "--copies".into(),
                copies.display().to_string(),
                "--out".into(),
                structure.display().to_string(),
                "--provenance".into(),
                structure_json.display().to_string(),
            ],
        );
    } else {
        writeln!(File::create(&structure).unwrap(), "family\tcopies\tvalid_flanks\ttsd_k\ttsd_support\ttsd_fraction\tterminal_direct_seed\tterminal_inverted_seed\tpoly_at_tail\tboundary_uncertain\tstructure_evidence").unwrap();
        writeln!(File::create(&structure_json).unwrap(), "{{\"contract\":\"te-looker-structure-audit-v1\",\"mode\":\"partial\",\"family_call\":false,\"families\":0,\"skipped_empty_query\":true}}").unwrap();
    }
    let acceptance = out.join("acceptance_audit.tsv");
    let acceptance_json = out.join("acceptance_provenance.json");
    let mut accept_args = vec![
        "--query".into(),
        query.display().to_string(),
        "--copies".into(),
        copies.display().to_string(),
        "--structure".into(),
        structure.display().to_string(),
        "--out".into(),
        acceptance.display().to_string(),
        "--provenance".into(),
        acceptance_json.display().to_string(),
    ];
    if non_te.is_some() {
        accept_args.push("--non-te-report".into());
        accept_args.push(out.join("non_te_report.tsv").display().to_string());
    }
    if let Some(manifest) = non_te_manifest.as_ref() {
        accept_args.push("--non-te-manifest".into());
        accept_args.push(manifest.clone());
    }
    run("dtr-accept-audit", &accept_args);
    let summary = out.join("track2_provenance.json");
    writeln!(File::create(&summary).unwrap(), "{{\"contract\":\"te-looker-track2-v1\",\"mode\":\"partial\",\"family_call\":false,\"genome\":\"{}\",\"windows\":{},\"provisional_consensi\":{},\"eligible_consensi\":{},\"non_te_reference\":{},\"final_library\":null,\"stages\":[\"window\",\"lsh\",\"community_prepartition\",\"community_refinement\",\"component_poa\",\"non_te_guard_optional\",\"copy_catalog\",\"structure_audit\",\"acceptance_audit\"]}}", json(&genome), fasta_count(&windows), fasta_count(&provisional), fasta_count(&query), non_te.as_ref().map(|x| format!("\"{}\"", json(x))).unwrap_or_else(|| "null".to_string())).unwrap();
    println!(
        "mode\tpartial\nprovisional_consensi\t{}\neligible_consensi\t{}",
        fasta_count(&provisional),
        fasta_count(&query)
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn json_escapes_paths() {
        assert_eq!(json("a\\b\"c"), "a\\\\b\\\"c");
    }
}
