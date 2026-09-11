from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent
HUB_DIR = BASE / "06.hub_genes"
EXPECTED_HUBS = {
    "AKT1",
    "BCL2",
    "CASP1",
    "MAPK14",
    "MMP9",
    "PPARG",
    "PTGS2",
    "PTPRC",
    "SIRT1",
}


def main() -> None:
    mcc = pd.read_csv(HUB_DIR / "MCC.csv")
    cytonca_all = pd.read_csv(HUB_DIR / "cytoNCA_all_28.csv")
    cytonca = pd.read_csv(HUB_DIR / "cytoNCA.csv")
    mcode = pd.read_csv(HUB_DIR / "MCODE.csv")
    intersection = pd.read_csv(HUB_DIR / "03_hub_genes_intersection_all.csv")

    assert len(mcc) == 10 and mcc["name"].nunique() == 10
    assert len(cytonca_all) == 28 and cytonca_all["Gene"].nunique() == 28
    assert list(cytonca_all["Rank"]) == list(range(1, 29))
    assert len(cytonca) == 20 and list(cytonca["Rank"]) == list(range(1, 21))
    assert len(mcode) == 9 and mcode["name"].nunique() == 9
    assert set(intersection["name"]) == EXPECTED_HUBS

    print("PASS: cytoHubba 10, CytoNCA 28->20, MCODE 9")
    print("PASS: exact nine-gene intersection")


if __name__ == "__main__":
    main()
