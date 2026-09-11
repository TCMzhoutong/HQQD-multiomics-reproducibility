#!/usr/bin/env python3
"""Create an SDF with original bond orders and docked heavy-atom coordinates.

The resulting file is used as the ACPYPE input, so ligand topology and ligand
coordinates come from one consistent molecule. Hydrogens are regenerated from
the docked heavy-atom geometry instead of being left at the original PubChem
conformer positions.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdFMCS


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--topology-sdf", required=True)
    ap.add_argument("--docked-pdb", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    topo = Chem.SDMolSupplier(args.topology_sdf, removeHs=False)[0]
    if topo is None:
        raise SystemExit(f"Cannot read topology SDF: {args.topology_sdf}")
    docked = Chem.MolFromPDBFile(args.docked_pdb, removeHs=False, sanitize=False)
    if docked is None:
        raise SystemExit(f"Cannot read docked PDB: {args.docked_pdb}")

    topo_hvy = Chem.RemoveHs(topo)
    docked_hvy = Chem.RemoveHs(docked)
    mcs = rdFMCS.FindMCS(
        [topo_hvy, docked_hvy],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        timeout=20,
    )
    query = Chem.MolFromSmarts(mcs.smartsString)
    topo_match = topo_hvy.GetSubstructMatch(query) if query else ()
    docked_match = docked_hvy.GetSubstructMatch(query) if query else ()
    if mcs.numAtoms != topo_hvy.GetNumAtoms() or mcs.numAtoms != docked_hvy.GetNumAtoms():
        raise SystemExit(
            f"MCS does not cover all heavy atoms: mcs={mcs.numAtoms} "
            f"topo_heavy={topo_hvy.GetNumAtoms()} docked_heavy={docked_hvy.GetNumAtoms()}"
        )

    dock_conf = docked_hvy.GetConformer()
    conf = topo_hvy.GetConformer()
    for t_idx, d_idx in zip(topo_match, docked_match):
        conf.SetAtomPosition(t_idx, dock_conf.GetAtomPosition(d_idx))

    # Add hydrogens with coordinates around the docked heavy-atom conformer.
    docked_topology = Chem.AddHs(topo_hvy, addCoords=True)
    Chem.SanitizeMol(docked_topology)
    writer = Chem.SDWriter(args.out)
    writer.write(docked_topology)
    writer.close()
    print(
        f"Wrote {args.out}; heavy_atoms={topo_hvy.GetNumAtoms()} "
        f"total_atoms={docked_topology.GetNumAtoms()}"
    )


if __name__ == "__main__":
    main()
