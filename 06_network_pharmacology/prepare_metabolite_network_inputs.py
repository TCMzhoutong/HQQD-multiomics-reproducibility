from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent
CANDIDATES = (
    BASE.parent
    / "04_correlation_analysis"
    / "network_pharmacology_candidates"
    / "candidate_metabolites_intersection.csv"
)
TARGET_DIR = BASE / "04.GutMicrobe_metabolite_targets"
EXISTING_STRUCTURES = TARGET_DIR / "02_metabolite_SMILES.csv"


CURATED_STRUCTURES = {
    "pos_4360": {
        "HMDB_ID": "",
        "KEGG_annotation": "",
        "PubChem_CID": "5460197",
        "SMILES": "CC(C(=O)CCCCCC(=O)[O-])N",
        "Search_Source": "CM_Name",
        "Matched_Query": "8-Amino-7-oxononanoate",
    },
}


def clean_string(value):
    if pd.isna(value):
        return ""
    return str(value).strip()


def main():
    candidates = pd.read_csv(CANDIDATES, encoding="utf-8-sig")
    candidates["CM_Name"] = candidates["Name"]

    old = pd.read_csv(EXISTING_STRUCTURES, encoding="utf-8-sig")
    structure_columns = [
        "HMDB_ID",
        "KEGG_annotation",
        "PubChem_CID",
        "SMILES",
        "Search_Source",
        "Matched_Query",
    ]
    old_lookup = old.set_index("ID")[structure_columns].to_dict(orient="index")

    structure_rows = []
    for feature_id in candidates["ID"]:
        if feature_id in CURATED_STRUCTURES:
            structure_rows.append({"ID": feature_id, **CURATED_STRUCTURES[feature_id]})
        elif feature_id in old_lookup:
            structure_rows.append({"ID": feature_id, **old_lookup[feature_id]})
        else:
            raise RuntimeError(f"No curated structure record is available for {feature_id}")

    structures = pd.DataFrame(structure_rows)
    merged = candidates.merge(structures, on="ID", how="left", validate="one_to_one")

    required = ["ID", "Name", "PubChem_CID", "SMILES"]
    missing = {
        col: merged.loc[merged[col].map(clean_string).eq(""), "ID"].tolist()
        for col in required
        if merged[col].map(clean_string).eq("").any()
    }
    if missing:
        raise RuntimeError(f"Incomplete metabolite structure metadata: {missing}")

    target_columns = list(old.columns)
    for col in target_columns:
        if col not in merged.columns:
            merged[col] = pd.NA
    merged = merged[target_columns]

    candidates.to_csv(TARGET_DIR / "HQQD_reverse_diff.csv", index=False, encoding="utf-8-sig")
    candidates.merge(
        structures[["ID", "HMDB_ID", "KEGG_annotation"]],
        on="ID",
        how="left",
        validate="one_to_one",
    ).to_csv(
        TARGET_DIR / "HQQD_reverse_diff_with_HMDB_ID.csv",
        index=False,
        encoding="utf-8-sig",
    )
    merged.to_csv(EXISTING_STRUCTURES, index=False, encoding="utf-8-sig")

    manifest = merged[
        ["ID", "CM_Name", "HMDB_ID", "PubChem_CID", "SMILES", "Search_Source", "Matched_Query"]
    ].copy()
    manifest["Included_In_Network_Pharmacology"] = True
    manifest.to_csv(
        TARGET_DIR / "network_metabolite_input_manifest.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print(f"Prepared {len(merged)} network-pharmacology metabolite inputs")
    print(merged[["ID", "CM_Name", "HMDB_ID", "PubChem_CID"]].to_string(index=False))


if __name__ == "__main__":
    main()
