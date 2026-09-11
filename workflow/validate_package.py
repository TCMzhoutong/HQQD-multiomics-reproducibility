#!/usr/bin/env python3
"""Validate the compact package against the current manuscript results."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FAILURES: list[str] = []
PASSES: list[str] = []


def read_rows(relative: str, delimiter: str = ",") -> list[dict[str, str]]:
    path = ROOT / relative
    if not path.exists():
        FAILURES.append(f"Missing file: {relative}")
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def expect(label: str, observed, expected) -> None:
    if observed == expected:
        PASSES.append(f"{label}: {observed}")
    else:
        FAILURES.append(f"{label}: observed {observed!r}; expected {expected!r}")


# Metabolomics and multi-omics
metabolites_32 = read_rows(
    "03_serum_metabolomics/results_reference/HQQD_reverse_candidates_32.csv"
)
expect("Retained HQQD-reversed metabolites", len(metabolites_32), 32)

candidate_5 = read_rows(
    "04_multiomics_correlation/results_reference/candidate_metabolites_5.csv"
)
expect("Prioritised metabolites", len(candidate_5), 5)
expected_metabolite_ids = {
    "pos_7892", "pos_8061", "pos_4360", "neg_12260", "neg_12234"
}
expect("Prioritised metabolite IDs", {row.get("ID") for row in candidate_5}, expected_metabolite_ids)

# Network pharmacology
common_targets = read_rows(
    "06_network_pharmacology/results_reference/three_way_common_targets_30.csv"
)
expect("Three-way common targets", len(common_targets), 30)

nodes = read_rows(
    "06_network_pharmacology/06.hub_genes/string_nodes_for_cytoscape.tsv", "\t"
)
expect("STRING input nodes", len(nodes), 30)
connected = [
    row for row in nodes
    if row.get("connected_in_network", "").strip().lower() in {"true", "1", "yes"}
]
expect("Connected STRING nodes", len(connected), 28)

edges = read_rows(
    "06_network_pharmacology/06.hub_genes/string_edges_for_cytoscape.tsv", "\t"
)
expect("STRING edges", len(edges), 105)

hubs = read_rows(
    "06_network_pharmacology/06.hub_genes/03_hub_genes_intersection_all.csv"
)
expected_hubs = {"AKT1", "BCL2", "CASP1", "MAPK14", "MMP9", "PPARG", "PTGS2", "PTPRC", "SIRT1"}
expect("Hub genes", {row.get("name") for row in hubs}, expected_hubs)

# Docking and molecular dynamics
docking_pairs = read_rows(
    "09_molecular_docking_dynamics/01_inputs/docking_pairs.csv"
)
expect("Docking pairs", len(docking_pairs), 44)
expect("Docking ligands", len({row.get("ligand_id") for row in docking_pairs}), 11)
expect("Docking targets", len({row.get("target") for row in docking_pairs}), 4)

docking_scores = read_rows(
    "09_molecular_docking_dynamics/09_results_tables/final_mean_sd_long.csv"
)
expect("Docking score summaries", len(docking_scores), 44)
expect("Docking repetitions", {row.get("n_runs") for row in docking_scores}, {"3"})

md_pairs = read_rows(
    "09_molecular_docking_dynamics/12_md_candidates/01_selection/md_retained_4_pairs.csv"
)
expect("Retained molecular-dynamics systems", len(md_pairs), 4)

active_names = "\n".join(
    str(value) for row in docking_pairs + docking_scores + md_pairs for value in row.values()
).lower()
excluded_annotations = read_rows(
    "03_serum_metabolomics/excluded_metabolite_annotations.csv"
)
excluded_tokens = {
    str(row.get(field, "")).strip().lower()
    for row in excluded_annotations
    for field in ("ID", "Reported_Annotation")
    if str(row.get(field, "")).strip()
}
active_excluded_tokens = sorted(token for token in excluded_tokens if token in active_names)
expect(
    "Prespecified excluded annotations present in active docking and MD manifests",
    active_excluded_tokens,
    [],
)

# Portability and repository size
personal_path_pattern = re.compile(
    r"(?:[A-Za-z]:[\\/]+Users[\\/]|/mnt/[a-z]/Users/|/home/[^/]+/(?:mini|ana)conda)",
    re.I,
)
text_suffixes = {
    ".r", ".py", ".sh", ".md", ".txt", ".csv", ".tsv", ".yml", ".yaml",
    ".json", ".toml", ".ini", ".cfg", ".ps1", ".bat", ".mdp", ".top", ".itp"
}
personal_hits: list[str] = []
oversized: list[str] = []
for path in ROOT.rglob("*"):
    if not path.is_file():
        continue
    if path.resolve() == Path(__file__).resolve():
        continue
    relative = path.relative_to(ROOT).as_posix()
    if path.stat().st_size > 100 * 1024 * 1024:
        oversized.append(relative)
    if path.suffix.lower() in text_suffixes:
        try:
            content = path.read_text(encoding="utf-8-sig", errors="ignore")
        except OSError:
            continue
        if personal_path_pattern.search(content):
            personal_hits.append(relative)

expect("Files larger than 100 MB", oversized, [])
expect("Personal absolute paths in distributable text files", personal_hits, [])

print("HQQD reproducibility-package validation")
for item in PASSES:
    print(f"PASS  {item}")
for item in FAILURES:
    print(f"FAIL  {item}")
print(f"Summary: {len(PASSES)} passed; {len(FAILURES)} failed")
sys.exit(1 if FAILURES else 0)
