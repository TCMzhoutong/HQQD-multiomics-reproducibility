"""Build Cytoscape-ready HQQD-metabolite-disease-target network tables.

The active network is derived from the five retained metabolites and the
30-way target intersection. Cytoscape styling/layout remains a manual step.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent
OUTPUT = BASE / "08.cytoscape"
DISEASE = "Bacterial pneumonia"


def edge_frame(source: str, targets: pd.DataFrame, source_col: str, target_col: str) -> pd.DataFrame:
    source_nodes = pd.DataFrame(
        {"source": source, "interaction": "connects", "target": targets[source_col].drop_duplicates()}
    )
    target_edges = targets[[source_col, target_col]].drop_duplicates().rename(
        columns={source_col: "source", target_col: "target"}
    )
    target_edges.insert(1, "interaction", "connects")
    return pd.concat([source_nodes, target_edges], ignore_index=True)


def node_frame(edges: pd.DataFrame, ingredients: set[str], metabolites: set[str], genes: set[str], disease: bool) -> pd.DataFrame:
    node_type = {"HQQD": (1, "Formula"), "Metabolite": (2, "Metabolite_Source")}
    if disease:
        node_type[DISEASE] = (0, "Disease")
    node_type.update({name: (3, "Ingredient") for name in ingredients})
    node_type.update({name: (4, "Metabolite") for name in metabolites})
    node_type.update({name: (5, "Target") for name in genes})

    degree = pd.concat([edges["source"], edges["target"]]).value_counts()
    rows = [
        {"name": name, "type": node_type[name][0], "type_name": node_type[name][1], "degree": int(degree.get(name, 0))}
        for name in sorted(node_type, key=lambda x: (node_type[x][0], x))
    ]
    return pd.DataFrame(rows)


def write_pair(edges: pd.DataFrame, nodes: pd.DataFrame, edge_stem: str, node_stem: str) -> None:
    edges.to_csv(OUTPUT / f"{edge_stem}.tsv", sep="\t", index=False, encoding="utf-8-sig")
    nodes.to_csv(OUTPUT / f"{node_stem}.tsv", sep="\t", index=False, encoding="utf-8-sig")
    edges.to_excel(OUTPUT / f"{edge_stem}.xlsx", index=False)
    nodes.to_excel(OUTPUT / f"{node_stem}.xlsx", index=False)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)

    ing = pd.read_csv(BASE / "03.ingredients_targets" / "02_hplc-ms_target.csv", encoding="utf-8-sig")
    met = pd.read_csv(BASE / "04.GutMicrobe_metabolite_targets" / "03_metabolite_target.csv", encoding="utf-8-sig")
    common = pd.read_csv(
        BASE / "05.disease_targets" / "04_disease_ingredient_metabolite_target.csv", encoding="utf-8-sig"
    )
    genes = set(common["targets"].dropna().astype(str).str.strip())

    ing = ing.dropna(subset=["Ingredient_name", "Common name"])
    met = met.dropna(subset=["CM_Name", "Common name"])
    ing = ing[ing["Common name"].isin(genes)].copy()
    met = met[met["Common name"].isin(genes)].copy()

    ingredients = set(ing["Ingredient_name"].astype(str))
    metabolites = set(met["CM_Name"].astype(str))
    edges = pd.concat(
        [
            edge_frame("HQQD", ing, "Ingredient_name", "Common name"),
            edge_frame("Metabolite", met, "CM_Name", "Common name"),
        ],
        ignore_index=True,
    ).drop_duplicates()
    nodes = node_frame(edges, ingredients, metabolites, genes, disease=False)

    disease_edges = pd.DataFrame(
        {"source": DISEASE, "interaction": "connects", "target": sorted(genes)}
    )
    edges_with_disease = pd.concat([edges, disease_edges], ignore_index=True).drop_duplicates()
    nodes_with_disease = node_frame(edges_with_disease, ingredients, metabolites, genes, disease=True)

    if len(metabolites) != 5 or len(genes) != 30:
        raise ValueError(f"Unexpected active network: {len(metabolites)} metabolites, {len(genes)} targets")
    if len(ingredients) != 37 or len(edges) != 254 or len(edges_with_disease) != 284:
        raise ValueError(
            f"Unexpected network dimensions: {len(ingredients)} ingredients, "
            f"{len(edges)} edges, {len(edges_with_disease)} disease-inclusive edges"
        )

    write_pair(edges, nodes, "network_all_edges", "type_all_nodes")
    write_pair(edges_with_disease, nodes_with_disease, "network_all_with_disease_edges", "type_all_with_disease_nodes")

    ingredient_audit = nodes_with_disease[nodes_with_disease["type_name"].eq("Ingredient")].sort_values(
        ["degree", "name"], ascending=[False, True]
    )
    ingredient_audit.to_csv(OUTPUT / "ingredient_degree_audit.csv", index=False, encoding="utf-8-sig")

    manifest = pd.DataFrame(
        [
            {"item": "retained_metabolites", "value": len(metabolites)},
            {"item": "shared_targets", "value": len(genes)},
            {"item": "ingredient_nodes", "value": len(ingredients)},
            {"item": "nodes_without_disease", "value": len(nodes)},
            {"item": "edges_without_disease", "value": len(edges)},
            {"item": "nodes_with_disease", "value": len(nodes_with_disease)},
            {"item": "edges_with_disease", "value": len(edges_with_disease)},
        ]
    )
    manifest.to_csv(OUTPUT / "cytoscape_network_manifest.tsv", sep="\t", index=False, encoding="utf-8-sig")

    print(f"Cytoscape network: {len(nodes_with_disease)} nodes, {len(edges_with_disease)} edges")
    print(f"Top ingredient degree: {ingredient_audit['degree'].max()}")


if __name__ == "__main__":
    main()
