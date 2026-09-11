#!/usr/bin/env python3
"""Build a protein-ligand GROMACS complex after ligand atom-order QC passes.

The script replaces the ACPYPE ligand coordinates with docked coordinates,
merges protein and ligand GRO files, and patches topology includes in a
reproducible way.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from pathlib import Path


def parse_gro(path: Path):
    lines = path.read_text().splitlines()
    title = lines[0]
    natoms = int(lines[1].strip())
    atoms = lines[2:2 + natoms]
    box = lines[2 + natoms]
    return title, atoms, box


def parse_pdb_coords(path: Path):
    coords = []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            x = float(line[30:38]) / 10.0
            y = float(line[38:46]) / 10.0
            z = float(line[46:54]) / 10.0
            coords.append((x, y, z))
    return coords


def gro_atom_element(line: str) -> str:
    name = line[10:15].strip()
    letters = "".join(c for c in name if c.isalpha())
    if not letters:
        return ""
    if letters[0].upper() == "H":
        return "H"
    return letters[0].upper()


def read_mapping(path: Path) -> dict[int, int]:
    mapping = {}
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if str(row.get("mapped", "")).lower() != "true":
                continue
            mapping[int(row["sdf_heavy_index_0based"])] = int(row["docked_heavy_index_0based"])
    return mapping


def replace_gro_heavy_coords(gro_atom_lines, coords, mapping):
    heavy_idx = [i for i, line in enumerate(gro_atom_lines) if gro_atom_element(line) != "H"]
    if len(heavy_idx) != len(mapping):
        raise SystemExit(f"Heavy atom mapping mismatch: gro_heavy={len(heavy_idx)} mapped={len(mapping)}")
    out = []
    coord_by_gro_idx = {}
    for sdf_heavy_pos, docked_heavy_pos in mapping.items():
        coord_by_gro_idx[heavy_idx[sdf_heavy_pos]] = coords[docked_heavy_pos]
    for i, line in enumerate(gro_atom_lines):
        if i in coord_by_gro_idx:
            x, y, z = coord_by_gro_idx[i]
            out.append(f"{line[:20]}{x:8.3f}{y:8.3f}{z:8.3f}{line[44:] if len(line) > 44 else ''}")
        else:
            out.append(line)
    return out


def patch_topology(topol: Path, ligand_itp_rel: str) -> None:
    txt = topol.read_text()
    include_block = (
        f'\n; Include ligand topology\n#include "{ligand_itp_rel}"\n\n'
        '; Include ligand position restraints\n#ifdef POSRES\n#include "posre_lig.itp"\n#endif\n'
    )
    if "Include ligand topology" not in txt:
        marker = '#include "amber99sb-ildn.ff/forcefield.itp"'
        if marker in txt:
            txt = txt.replace(marker, marker + include_block, 1)
        else:
            txt = include_block + "\n" + txt
    if re.search(r"^UNL\\s+\\d+", txt, flags=re.M) is None:
        txt = txt.rstrip() + "\nUNL                 1\n"
    topol.write_text(txt)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair-dir", required=True)
    args = ap.parse_args()
    pair_dir = Path(args.pair_dir)
    prep = pair_dir / "01_prepare"
    run = pair_dir / "02_run"
    run.mkdir(exist_ok=True)

    qc = prep / "ligand_mapping_qc.csv"
    if not qc.exists() or "False" in qc.read_text():
        raise SystemExit("Ligand atom-order QC has not passed; refusing to merge system.")

    protein_title, protein_atoms, protein_box = parse_gro(prep / "protein.gro")
    ligand_title, ligand_atoms, _ = parse_gro(prep / "UNL.acpype" / "UNL_GMX.gro")
    # UNL_GMX.gro was generated from ligand_docked_for_acpype.sdf, which already
    # contains docked heavy-atom coordinates and regenerated hydrogen positions.

    complex_atoms = protein_atoms + ligand_atoms
    complex_gro = ["Protein-ligand complex", str(len(complex_atoms)), *complex_atoms, protein_box]
    (run / "complex.gro").write_text("\n".join(complex_gro) + "\n")

    topol = run / "topol.top"
    topol.write_text((prep / "topol.top").read_text())
    for itp in prep.glob("*.itp"):
        shutil.copy2(itp, run / itp.name)
    patch_topology(topol, "../01_prepare/UNL.acpype/UNL_GMX.itp")
    print(f"Wrote {run / 'complex.gro'} and {topol}")


if __name__ == "__main__":
    main()
