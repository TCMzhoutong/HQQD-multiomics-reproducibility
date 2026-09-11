#!/usr/bin/env python3
"""Make a docked thyroxine SDF chemically consistent with net charge -1.

The PubChem-like neutral thyroxine SDF can contain a neutral carboxylic acid.
If ACPYPE is asked to parameterize that neutral structure with total charge -1,
AmberTools sees an odd electron count and AM1-BCC fails. For the MD workflow we
use a closed-shell monoanion by deprotonating the carboxylic acid oxygen.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def find_carboxyl_oh_hydrogen(mol: Chem.Mol) -> tuple[int, int]:
    for carbon in mol.GetAtoms():
        if carbon.GetSymbol() != "C":
            continue
        oxygens = []
        for oxygen in carbon.GetNeighbors():
            if oxygen.GetSymbol() != "O":
                continue
            bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), oxygen.GetIdx())
            oxygens.append((oxygen, bond.GetBondType()))
        has_carbonyl = any(bond_type == Chem.BondType.DOUBLE for _, bond_type in oxygens)
        if not has_carbonyl:
            continue
        for oxygen, bond_type in oxygens:
            if bond_type != Chem.BondType.SINGLE:
                continue
            for neighbor in oxygen.GetNeighbors():
                if neighbor.GetSymbol() == "H":
                    return oxygen.GetIdx(), neighbor.GetIdx()
    raise SystemExit("Could not find a carboxylic-acid OH hydrogen in thyroxine SDF.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", required=True)
    parser.add_argument("--outfile", required=True)
    parser.add_argument("--site", choices=["carboxylate"], default="carboxylate")
    args = parser.parse_args()

    infile = Path(args.infile)
    mol = Chem.SDMolSupplier(str(infile), removeHs=False)[0]
    if mol is None:
        raise SystemExit(f"Could not read SDF: {infile}")

    oxygen_idx, hydrogen_idx = find_carboxyl_oh_hydrogen(mol)
    rw = Chem.RWMol(mol)
    rw.RemoveAtom(hydrogen_idx)
    if hydrogen_idx < oxygen_idx:
        oxygen_idx -= 1
    rw.GetAtomWithIdx(oxygen_idx).SetFormalCharge(-1)
    fixed = rw.GetMol()
    Chem.SanitizeMol(fixed)

    charge = sum(atom.GetFormalCharge() for atom in fixed.GetAtoms())
    formula = rdMolDescriptors.CalcMolFormula(fixed)
    if charge != -1:
        raise SystemExit(f"Expected thyroxine monoanion charge -1, got {charge}.")

    writer = Chem.SDWriter(str(args.outfile))
    writer.write(fixed)
    writer.close()
    print(f"Wrote thyroxine monoanion SDF: formula={formula} formal_charge={charge}")


if __name__ == "__main__":
    main()
