from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent
EXPECTED_CANDIDATES = {
    "Thyroxine",
    "N-Oleoyl Alanine",
    "8-Amino-7-oxononanoate",
    "PE(P-18:0/18:1(12Z)-2OH(9,10))",
    "12,13-DiHODE",
}
EXPECTED_RESTORED_NEGATIVE_IDS = {
    "neg_9297",
    "neg_9101",
    "neg_8898",
    "neg_8900",
    "neg_9896",
    "neg_9514",
    "neg_9416",
    "neg_7945",
    "neg_8249",
    "neg_4991",
}


def main():
    reference = BASE / "results_reference"
    metabolites = pd.read_csv(
        reference / "differential_metabolites_32.csv",
        encoding="utf-8-sig",
    )
    rho = pd.read_csv(
        reference / "spearman_rho_matrix.csv",
        index_col=0,
        encoding="utf-8-sig",
    )
    p_value = pd.read_csv(
        reference / "spearman_pvalue_matrix.csv",
        index_col=0,
        encoding="utf-8-sig",
    )
    taxa_mantel = pd.read_csv(
        reference / "taxa_mantel_results.csv",
        encoding="utf-8-sig",
    )
    metabolite_mantel = pd.read_csv(
        reference / "metabolite_mantel_results.csv",
        encoding="utf-8-sig",
    )
    candidates = pd.read_csv(
        reference / "candidate_metabolites_5.csv",
        encoding="utf-8-sig",
    )

    assert len(metabolites) == 32
    assert EXPECTED_RESTORED_NEGATIVE_IDS.issubset(set(metabolites["ID"]))
    assert rho.shape == (18, 32)
    assert p_value.shape == (18, 32)
    assert int((p_value < 0.05).sum().sum()) == 237
    assert len(taxa_mantel) == 90
    assert int((taxa_mantel["p"] < 0.05).sum()) == 15
    assert len(metabolite_mantel) == 160
    assert int((metabolite_mantel["p"] < 0.05).sum()) == 27
    assert set(candidates["Name"]) == EXPECTED_CANDIDATES

    print("PASS: 18 taxa x 32 metabolites; 237 significant Spearman associations")
    print("PASS: all 10 previously omitted negative-ion features entered the statistics")
    print("PASS: 15 taxa-phenotype and 27 metabolite-phenotype Mantel associations")
    print("PASS: exact five-metabolite candidate set matches the current manuscript")


if __name__ == "__main__":
    main()
