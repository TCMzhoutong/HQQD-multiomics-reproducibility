from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent
EXPECTED_METABOLITES = {
    "Thyroxine",
    "N-Oleoyl Alanine",
    "8-Amino-7-oxononanoate",
    "PE(P-18:0/18:1(12Z)-2OH(9,10))",
    "12,13-DiHODE",
}
EXPECTED_COMMON_TARGETS = {
    "ADAM17", "AKR1B1", "AKT1", "ANPEP", "APEX1", "BCHE", "BCL2", "CA1",
    "CA9", "CASP1", "CDC25B", "CES1", "CTSD", "EPHX2", "F2", "FUCA1", "FYN",
    "GPBAR1", "GSTP1", "LCK", "MAPK14", "MME", "MMP8", "MMP9", "PGR", "PPARG",
    "PTGER4", "PTGS2", "PTPRC", "SIRT1",
}


def main():
    inputs = pd.read_csv(
        BASE / "04.GutMicrobe_metabolite_targets" / "02_metabolite_SMILES.csv",
        encoding="utf-8-sig",
    )
    common = pd.read_csv(
        BASE / "05.disease_targets" / "04_disease_ingredient_metabolite_target.csv",
        encoding="utf-8-sig",
    )
    edges = pd.read_csv(
        BASE / "06.hub_genes" / "string_edges_for_cytoscape.tsv",
        sep="\t",
        encoding="utf-8-sig",
    )
    nodes = pd.read_csv(
        BASE / "06.hub_genes" / "string_nodes_for_cytoscape.tsv",
        sep="\t",
        encoding="utf-8-sig",
    )

    assert set(inputs["CM_Name"]) == EXPECTED_METABOLITES
    assert len(common) == 30
    assert set(common["targets"]) == EXPECTED_COMMON_TARGETS
    assert "MGAM" not in set(common["targets"])
    assert "PTPRC" in set(common["targets"])
    assert len(nodes) == 30
    assert int(nodes["connected_in_network"].sum()) == 28
    assert set(nodes.loc[~nodes["connected_in_network"], "name"]) == {"CES1", "FUCA1"}
    assert len(edges) == 105

    print("PASS: 5 candidate metabolites and 30 three-way common targets")
    print("PASS: STRING network contains 28 connected nodes and 105 edges")


if __name__ == "__main__":
    main()
