"""Normalize and validate the three manual Cytoscape hub-gene exports.

Expected manual inputs in ``08.cytoscape``:

* NCC.csv: cytoHubba Top-10 selected-node table (normalized to MCC.csv)
* cytonca.txt: all 28 connected nodes ranked by Subgragh centrality
* mcode.csv: nodes from the highest-scoring MCODE cluster

The script applies the original analysis rule to the only unfiltered export:
CytoNCA is truncated to the first 20 ranked nodes.  It then writes the three
standardized inputs consumed by ``venn_hub_genes.r``.
"""

from __future__ import annotations

import hashlib
import re
import shutil
from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent
INPUT_DIR = BASE / "08.cytoscape"
OUTPUT_DIR = BASE / "06.hub_genes"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_cytonca_report(path: Path) -> tuple[pd.DataFrame, str]:
    raw = path.read_bytes()
    text = None
    selected_encoding = None
    for encoding in ("utf-8", "gb18030", "cp936", "latin1"):
        try:
            candidate = raw.decode(encoding)
        except UnicodeDecodeError:
            continue
        if "Results ranked by Subgragh" in candidate:
            text = candidate
            selected_encoding = encoding
            break
    if text is None or selected_encoding is None:
        raise ValueError(f"Unable to decode a CytoNCA Subgragh report: {path}")

    records: list[dict[str, object]] = []
    for line in text.splitlines():
        fields = [field.strip() for field in line.split("\t") if field.strip()]
        if len(fields) < 3 or not fields[0].isdigit():
            continue
        record: dict[str, object] = {"Rank": int(fields[0]), "Gene": fields[1]}
        for field in fields[2:]:
            if ":" not in field:
                continue
            key, value = field.split(":", 1)
            record[key.strip()] = float(value.strip())
        records.append(record)

    expected_columns = [
        "Rank",
        "Gene",
        "Subgragh",
        "Degree",
        "Eigenvector",
        "Information",
        "LAC",
        "Betweenness",
        "Closeness",
        "Network",
    ]
    result = pd.DataFrame(records)
    missing = set(expected_columns) - set(result.columns)
    if missing:
        raise ValueError(f"CytoNCA report is missing fields: {sorted(missing)}")
    result = result[expected_columns].sort_values("Rank").reset_index(drop=True)
    return result, selected_encoding


def unique_gene_set(frame: pd.DataFrame, column: str, label: str) -> set[str]:
    if column not in frame.columns:
        raise ValueError(f"{label} is missing the {column!r} column")
    genes = frame[column].dropna().astype(str).str.strip()
    if (genes == "").any() or genes.duplicated().any():
        raise ValueError(f"{label} contains blank or duplicated gene symbols")
    return set(genes)


def main() -> None:
    source_mcc = INPUT_DIR / "NCC.csv"
    source_mcode = INPUT_DIR / "mcode.csv"
    source_cytonca = INPUT_DIR / "cytonca.txt"
    for path in (source_mcc, source_mcode, source_cytonca):
        if not path.is_file():
            raise FileNotFoundError(path)

    edges = pd.read_csv(OUTPUT_DIR / "string_edges_for_cytoscape.tsv", sep="\t")
    connected_genes = set(edges["source"]).union(edges["target"])
    if len(connected_genes) != 28:
        raise AssertionError(f"Expected 28 connected STRING nodes, found {len(connected_genes)}")

    mcc = pd.read_csv(source_mcc)
    mcc_genes = unique_gene_set(mcc, "name", "cytoHubba Top-10 export")
    if len(mcc_genes) != 10:
        raise AssertionError(f"Expected 10 cytoHubba genes, found {len(mcc_genes)}")

    mcode = pd.read_csv(source_mcode)
    mcode_genes = unique_gene_set(mcode, "name", "MCODE export")
    if "MCODE::Clusters (1)" in mcode.columns and mcode["MCODE::Clusters (1)"].nunique() != 1:
        raise AssertionError("The MCODE export contains more than one cluster")

    cytonca_all, cytonca_encoding = read_cytonca_report(source_cytonca)
    cytonca_all_genes = unique_gene_set(cytonca_all, "Gene", "CytoNCA full ranking")
    if len(cytonca_all) != 28 or list(cytonca_all["Rank"]) != list(range(1, 29)):
        raise AssertionError("CytoNCA report must contain ranks 1-28 exactly once")
    if cytonca_all_genes != connected_genes:
        missing = sorted(connected_genes - cytonca_all_genes)
        unexpected = sorted(cytonca_all_genes - connected_genes)
        raise AssertionError(
            f"CytoNCA genes do not match the 28-node STRING network; missing={missing}, unexpected={unexpected}"
        )

    for label, genes in (("cytoHubba", mcc_genes), ("MCODE", mcode_genes)):
        unexpected = sorted(genes - connected_genes)
        if unexpected:
            raise AssertionError(f"{label} contains genes outside the current STRING network: {unexpected}")

    cytonca_top20 = cytonca_all.head(20).copy()
    cytonca_top20_genes = set(cytonca_top20["Gene"])

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    mcc.to_csv(OUTPUT_DIR / "MCC.csv", index=False, encoding="utf-8")
    mcode.to_csv(OUTPUT_DIR / "MCODE.csv", index=False, encoding="utf-8")
    cytonca_all.to_csv(OUTPUT_DIR / "cytoNCA_all_28.csv", index=False, encoding="utf-8")
    cytonca_top20.to_csv(OUTPUT_DIR / "cytoNCA.csv", index=False, encoding="utf-8")
    shutil.copy2(source_cytonca, OUTPUT_DIR / "cytoNCA_source_report.txt")

    intersection = mcc_genes & cytonca_top20_genes & mcode_genes
    manifest = pd.DataFrame(
        [
            {
                "method": "cytoHubba",
                "source_file": source_mcc.name,
                "source_count": len(mcc),
                "selection_rule": "provided Top 10",
                "selected_count": len(mcc_genes),
                "standardized_output": "06.hub_genes/MCC.csv",
                "source_sha256": file_sha256(source_mcc),
            },
            {
                "method": "CytoNCA",
                "source_file": source_cytonca.name,
                "source_count": len(cytonca_all),
                "selection_rule": "top 20 by reported Subgragh rank",
                "selected_count": len(cytonca_top20),
                "standardized_output": "06.hub_genes/cytoNCA.csv",
                "source_sha256": file_sha256(source_cytonca),
            },
            {
                "method": "MCODE",
                "source_file": source_mcode.name,
                "source_count": len(mcode),
                "selection_rule": "all nodes in exported highest-scoring cluster",
                "selected_count": len(mcode_genes),
                "standardized_output": "06.hub_genes/MCODE.csv",
                "source_sha256": file_sha256(source_mcode),
            },
        ]
    )
    manifest.to_csv(OUTPUT_DIR / "cytoscape_hub_export_manifest.tsv", sep="\t", index=False)

    audit = f"""# Cytoscape hub-export preprocessing audit

- Current STRING network: 28 connected nodes and {len(edges)} edges.
- cytoHubba input: 10 provided nodes; no further numerical truncation applied.
- CytoNCA input: all 28 nodes read from `{source_cytonca.name}` using `{cytonca_encoding}`; ranks 1-20 retained by the original rule.
- MCODE input: all {len(mcode_genes)} nodes in the exported highest-scoring cluster retained.
- Three-way intersection expected before Venn rendering: {len(intersection)} genes.

Three-way genes: {', '.join(sorted(intersection))}.
"""
    (OUTPUT_DIR / "cytoscape_hub_preprocessing_audit.md").write_text(audit, encoding="utf-8")

    print(f"cytoHubba: {len(mcc_genes)} genes")
    print(f"CytoNCA: {len(cytonca_all)} ranked genes -> top {len(cytonca_top20)}")
    print(f"MCODE: {len(mcode_genes)} genes")
    print(f"Expected three-way intersection: {len(intersection)} genes")


if __name__ == "__main__":
    main()
