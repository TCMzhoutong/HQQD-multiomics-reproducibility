from io import StringIO
from pathlib import Path

import pandas as pd
import requests


BASE = Path(__file__).resolve().parent
TARGET_FILE = BASE / "05.disease_targets" / "04_disease_ingredient_metabolite_target.csv"
OUT_DIR = BASE / "06.hub_genes"
FIG_DIR = BASE / "figures"
API_ROOT = "https://string-db.org/api"
SPECIES = 9606
REQUIRED_SCORE = 400
CALLER_IDENTITY = "HQQD_KP_pneumonia_revision"


def post_string(output_format, method, genes):
    response = requests.post(
        f"{API_ROOT}/{output_format}/{method}",
        data={
            "identifiers": "\r".join(genes),
            "species": SPECIES,
            "required_score": REQUIRED_SCORE,
            "network_type": "functional",
            "add_nodes": 0,
            "show_query_node_labels": 1,
            "caller_identity": CALLER_IDENTITY,
        },
        timeout=300,
    )
    response.raise_for_status()
    return response


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    genes = (
        pd.read_csv(TARGET_FILE, encoding="utf-8-sig")["targets"]
        .dropna()
        .astype(str)
        .str.strip()
    )
    genes = sorted(g for g in genes.unique() if g)
    if not genes:
        raise RuntimeError("The three-way target intersection is empty")

    mapping_response = requests.post(
        f"{API_ROOT}/tsv/get_string_ids",
        data={
            "identifiers": "\r".join(genes),
            "species": SPECIES,
            "echo_query": 1,
            "caller_identity": CALLER_IDENTITY,
        },
        timeout=300,
    )
    mapping_response.raise_for_status()
    mapping = pd.read_csv(StringIO(mapping_response.text), sep="\t")
    mapping.to_csv(OUT_DIR / "string_id_mapping.tsv", sep="\t", index=False)

    missing = sorted(set(genes) - set(mapping["queryItem"].astype(str)))
    if missing:
        raise RuntimeError(f"STRING identifier mapping failed for: {', '.join(missing)}")

    network_response = post_string("tsv", "network", genes)
    raw = pd.read_csv(StringIO(network_response.text), sep="\t")
    raw.to_csv(OUT_DIR / "string_network_raw.tsv", sep="\t", index=False)

    edges = raw.rename(
        columns={
            "preferredName_A": "source",
            "preferredName_B": "target",
            "score": "combined_score",
            "nscore": "neighborhood_on_chromosome",
            "fscore": "gene_fusion",
            "pscore": "phylogenetic_cooccurrence",
            "ascore": "coexpression",
            "escore": "experimentally_determined_interaction",
            "dscore": "database_annotated",
            "tscore": "automated_textmining",
        }
    )
    edge_columns = [
        "source",
        "target",
        "combined_score",
        "neighborhood_on_chromosome",
        "gene_fusion",
        "phylogenetic_cooccurrence",
        "coexpression",
        "experimentally_determined_interaction",
        "database_annotated",
        "automated_textmining",
    ]
    edges = edges[edge_columns].sort_values(
        ["combined_score", "source", "target"], ascending=[False, True, True]
    )
    edges.to_csv(OUT_DIR / "string_edges_for_cytoscape.tsv", sep="\t", index=False)

    connected = set(edges["source"]).union(edges["target"])
    map_lookup = mapping.drop_duplicates("queryItem").set_index("queryItem")
    nodes = pd.DataFrame({"name": genes})
    nodes["string_id"] = nodes["name"].map(map_lookup["stringId"])
    nodes["connected_in_network"] = nodes["name"].isin(connected)
    nodes["type"] = "Shared target"
    nodes.to_csv(OUT_DIR / "string_nodes_for_cytoscape.tsv", sep="\t", index=False)

    compatible = pd.DataFrame(
        {
            "#node1": raw["preferredName_A"],
            "node2": raw["preferredName_B"],
            "node1_string_id": raw["stringId_A"],
            "node2_string_id": raw["stringId_B"],
            "neighborhood_on_chromosome": raw["nscore"],
            "gene_fusion": raw["fscore"],
            "phylogenetic_cooccurrence": raw["pscore"],
            "homology": pd.NA,
            "coexpression": raw["ascore"],
            "experimentally_determined_interaction": raw["escore"],
            "database_annotated": raw["dscore"],
            "automated_textmining": raw["tscore"],
            "combined_score": raw["score"],
        }
    )
    compatible.to_csv(OUT_DIR / "string_interactions_short.tsv", sep="\t", index=False)

    svg_response = post_string("svg", "network", genes)
    (FIG_DIR / "string_ppi_network.svg").write_bytes(svg_response.content)

    image_response = post_string("image", "network", genes)
    (FIG_DIR / "string_ppi_network.png").write_bytes(image_response.content)

    link_response = post_string("tsv-no-header", "get_link", genes)
    (OUT_DIR / "string_network_link.txt").write_text(link_response.text.strip() + "\n", encoding="utf-8")

    print(f"STRING input genes: {len(genes)}")
    print(f"Mapped genes: {mapping['queryItem'].nunique()}")
    print(f"Connected nodes: {len(connected)}")
    print(f"Edges at required score >= {REQUIRED_SCORE}: {len(edges)}")
    print("Required Cytoscape network file:")
    print(OUT_DIR / "string_edges_for_cytoscape.tsv")
    print("Audit-only input-node status table (do not import as a second network):")
    print(OUT_DIR / "string_nodes_for_cytoscape.tsv")
    print("STRING SVG:")
    print(FIG_DIR / "string_ppi_network.svg")


if __name__ == "__main__":
    main()
