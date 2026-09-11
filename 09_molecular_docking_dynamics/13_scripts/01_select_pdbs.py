from pathlib import Path
import json
import urllib.request

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

UNIPROT = {
    "PTPRC": "P08575",
    "CASP1": "P29466",
    "BCL2": "P10415",
    "MAPK14": "Q16539",
}

PREFERRED = {
    "BCL2": "8HTS",
    "CASP1": "6PZP",
    "PTPRC": "1YGU",
    "MAPK14": "3HLL",
}

GRAPHQL = """query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    struct{title}
    exptl{method}
    rcsb_entry_info{resolution_combined}
    polymer_entities{
      rcsb_polymer_entity_container_identifiers{auth_asym_ids uniprot_ids}
      rcsb_entity_source_organism{scientific_name}
      rcsb_polymer_entity{pdbx_description}
    }
    nonpolymer_entities{
      pdbx_entity_nonpoly{comp_id name}
      rcsb_nonpolymer_entity_container_identifiers{auth_asym_ids}
    }
  }
}"""


def post_json(url: str, payload: dict) -> dict:
    request = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(request, timeout=90) as response:
        return json.load(response)


def search_uniprot(uniprot: str) -> list[str]:
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot,
            },
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 100},
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }
    data = post_json("https://search.rcsb.org/rcsbsearch/v2/query", query)
    return [item["identifier"] for item in data.get("result_set", [])]


def fetch_entries(ids: list[str]) -> list[dict]:
    if not ids:
        return []
    payload = {"query": GRAPHQL, "variables": {"ids": ids}}
    data = post_json("https://data.rcsb.org/graphql", payload)
    return data.get("data", {}).get("entries", [])


def summarize_entry(target: str, entry: dict) -> dict:
    method = "; ".join(x.get("method", "") for x in entry.get("exptl") or [])
    resolutions = entry.get("rcsb_entry_info", {}).get("resolution_combined") or []
    resolution = min(resolutions) if resolutions else None
    human = False
    uniprot_match = False
    protein_chains = []
    polymer_descriptions = []
    for entity in entry.get("polymer_entities") or []:
        ids = entity.get("rcsb_polymer_entity_container_identifiers") or {}
        orgs = entity.get("rcsb_entity_source_organism") or []
        desc = (entity.get("rcsb_polymer_entity") or {}).get("pdbx_description")
        polymer_descriptions.append(desc)
        if UNIPROT[target] in (ids.get("uniprot_ids") or []):
            uniprot_match = True
            protein_chains.extend(ids.get("auth_asym_ids") or [])
        if any((org.get("scientific_name") or "").lower() == "homo sapiens" for org in orgs):
            human = True
    ligands = []
    for entity in entry.get("nonpolymer_entities") or []:
        non = entity.get("pdbx_entity_nonpoly") or {}
        comp_id = non.get("comp_id")
        if comp_id and comp_id not in {"HOH", "WAT"}:
            ligands.append(f"{comp_id}:{non.get('name')}")
    return {
        "target": target,
        "pdb_id": entry["rcsb_id"],
        "title": (entry.get("struct") or {}).get("title"),
        "method": method,
        "resolution_A": resolution,
        "human": human,
        "uniprot_match": uniprot_match,
        "protein_chains_for_uniprot": ",".join(protein_chains),
        "has_nonwater_ligand": bool(ligands),
        "nonwater_ligands": "; ".join(ligands),
        "polymer_descriptions": "; ".join(x for x in polymer_descriptions if x),
    }


def main() -> None:
    rows = []
    for target, uniprot in UNIPROT.items():
        ids = search_uniprot(uniprot)
        if PREFERRED[target] not in ids:
            ids.insert(0, PREFERRED[target])
        entries = fetch_entries(ids[:100])
        for entry in entries:
            if entry:
                rows.append(summarize_entry(target, entry))
    candidates = pd.DataFrame(rows)
    candidates["meets_resolution_le_2_5"] = candidates["resolution_A"].le(2.5)
    candidates["meets_core_rule"] = (
        candidates["human"]
        & candidates["uniprot_match"]
        & candidates["method"].str.contains("X-RAY", case=False, na=False)
        & candidates["meets_resolution_le_2_5"]
        & candidates["has_nonwater_ligand"]
    )
    candidates = candidates.sort_values(
        ["target", "meets_core_rule", "resolution_A", "has_nonwater_ligand"],
        ascending=[True, False, True, False],
    )
    out = ROOT / "02_pdb_selection" / "pdb_candidates_from_rcsb.csv"
    candidates.to_csv(out, index=False, encoding="utf-8-sig")

    selected = pd.DataFrame(
        [
            {
                "target": "PTPRC",
                "pdb_id": "1YGU",
                "decision": "Use user-prior structure with pTyr peptide. No <=2.5 A human PTPRC ligand-bound PDB was found in current RCSB candidate scan; document limitation.",
            },
            {"target": "CASP1", "pdb_id": "6PZP", "decision": "Human, X-ray, 1.94 A, ligand-bound VX-765/P7S; matches prior choice."},
            {"target": "BCL2", "pdb_id": "8HTS", "decision": "Human, X-ray, 1.25 A, ligand-bound; better resolution than previous 6O0K."},
            {"target": "MAPK14", "pdb_id": "3HLL", "decision": "Human, X-ray, 1.95 A, ligand-bound PH-797804/I45."},
        ]
    )
    selected = selected.merge(
        candidates.drop_duplicates(["target", "pdb_id"]),
        on=["target", "pdb_id"],
        how="left",
    )
    selected.to_csv(ROOT / "02_pdb_selection" / "selected_pdbs.csv", index=False, encoding="utf-8-sig")
    print(selected[["target", "pdb_id", "resolution_A", "has_nonwater_ligand", "decision"]].to_string(index=False))


if __name__ == "__main__":
    main()
