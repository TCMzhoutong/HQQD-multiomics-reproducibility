#!/usr/bin/env python3
"""QC atom compatibility between topology SDF and docked ligand PDB."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdFMCS


def read_pdb_elements(path: Path) -> list[str]:
    elems = []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            elem = line[76:78].strip()
            if not elem:
                name = line[12:16].strip()
                elem = "".join(c for c in name if c.isalpha())[:1]
            elems.append(elem.upper())
    return elems


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--topology-sdf", required=True)
    ap.add_argument("--docked-pdb", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    mol = Chem.SDMolSupplier(args.topology_sdf, removeHs=False)[0]
    if mol is None:
        raise SystemExit(f"Cannot read SDF: {args.topology_sdf}")
    sdf_elems = [a.GetSymbol().upper() for a in mol.GetAtoms()]
    sdf_heavy = [e for e in sdf_elems if e != "H"]
    pdb_elems = read_pdb_elements(Path(args.docked_pdb))
    pdb_heavy = [e for e in pdb_elems if e != "H"]
    pdb_mol = Chem.MolFromPDBFile(args.docked_pdb, removeHs=False, sanitize=False)
    if pdb_mol is None:
        raise SystemExit(f"Cannot read docked PDB: {args.docked_pdb}")
    sdf_noh = Chem.RemoveHs(mol)
    pdb_noh = Chem.RemoveHs(pdb_mol)

    rows = []
    ok = len(sdf_heavy) == len(pdb_heavy)
    if ok:
        mcs = rdFMCS.FindMCS(
            [sdf_noh, pdb_noh],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            timeout=20,
        )
        query = Chem.MolFromSmarts(mcs.smartsString)
        sdf_match = sdf_noh.GetSubstructMatch(query) if query else ()
        pdb_match = pdb_noh.GetSubstructMatch(query) if query else ()
        ok = mcs.numAtoms == len(sdf_heavy) == len(pdb_heavy)
    else:
        sdf_match = ()
        pdb_match = ()

    match_by_sdf = {}
    for s_idx, p_idx in zip(sdf_match, pdb_match):
        match_by_sdf[s_idx] = p_idx

    for i, s in enumerate(sdf_heavy):
        p_idx = match_by_sdf.get(i, "")
        p = pdb_heavy[p_idx] if isinstance(p_idx, int) and p_idx < len(pdb_heavy) else ""
        rows.append(
            {
                "sdf_heavy_index_0based": i,
                "sdf_element": s,
                "docked_heavy_index_0based": p_idx,
                "docked_pdb_element": p,
                "mapped": isinstance(p_idx, int) and s == p,
            }
        )

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "sdf_heavy_index_0based",
                "sdf_element",
                "docked_heavy_index_0based",
                "docked_pdb_element",
                "mapped",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    status = "PASS" if ok else "FAIL"
    print(
        f"ligand heavy-atom mapping QC: {status}; "
        f"sdf_atoms={len(sdf_elems)} sdf_heavy={len(sdf_heavy)} "
        f"docked_atoms={len(pdb_elems)} docked_heavy={len(pdb_heavy)}"
    )
    if not ok:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
