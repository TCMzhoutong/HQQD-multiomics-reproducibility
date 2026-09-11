from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent
def main():
    reference = BASE / "results_reference"
    summary = pd.read_csv(reference / "HQQD_reverse_summary.csv").iloc[0]
    candidates = pd.read_csv(
        reference / "HQQD_reverse_candidates_32.csv",
        encoding="utf-8-sig",
    )
    excluded = pd.read_csv(
        reference / "excluded_annotations_6.csv",
        encoding="utf-8-sig",
    )

    assert int(summary["ReverseCandidatesBeforeExclusion"]) == 38
    assert int(summary["ExcludedAnnotations"]) == 6
    assert int(summary["ReverseCandidates"]) == 32
    assert len(candidates) == 32
    assert excluded["ID"].nunique() == 6
    assert not set(candidates["ID"]).intersection(set(excluded["ID"]))

    # Guard against silent locale/encoding regressions in both ion modes.
    pos_features = pd.read_csv(BASE / "metabolites_exp_pos.csv")
    neg_features = pd.read_csv(BASE / "metabolites_exp_neg.csv")
    assert pos_features.loc[pos_features["ID"] == "pos_8351", "name"].item() == (
        "3-Oxo-5β-chola-8(14),11-dien-24-oic Acid"
    )
    assert neg_features.loc[neg_features["ID"] == "neg_5611", "name"].item() == (
        "17α,20β-Hydroxyprogesterone sulfate"
    )

    print("PASS: 38 raw reverse candidates - 6 excluded annotations = 32 retained candidates")
    print("PASS: UTF-8 metabolite labels preserved in both ion modes")


if __name__ == "__main__":
    main()
