#!/usr/bin/env python3
"""Align SWISS-MODEL repaired receptors back to the docking receptor frame."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ORIGINAL_DIR = ROOT / "03_receptor_preparation" / "02_clean_single_receptor_pdb"
SWISSMODEL_DIR = ROOT / "12_md_candidates" / "07_swissmodel_repaired_outputs"
OUT_DIR = ROOT / "12_md_candidates" / "09_aligned_repaired_receptors"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


@dataclass
class Atom:
    line: str
    resname: str
    atomname: str
    coord: np.ndarray


@dataclass
class Residue:
    resname: str
    ca: np.ndarray


PAIRS = {
    "MAPK14_3HLL": (
        ORIGINAL_DIR / "MAPK14_3HLL_clean.pdb",
        SWISSMODEL_DIR / "MAPK14.pdb",
        OUT_DIR / "MAPK14_3HLL_aligned_repaired_chainA.pdb",
    ),
    "BCL2_8HTS": (
        ORIGINAL_DIR / "BCL2_8HTS_clean.pdb",
        SWISSMODEL_DIR / "BCL2.pdb",
        OUT_DIR / "BCL2_8HTS_aligned_repaired_chainA.pdb",
    ),
    "PTPRC_1YGU": (
        ORIGINAL_DIR / "PTPRC_1YGU_clean.pdb",
        SWISSMODEL_DIR / "PTPRC.pdb",
        OUT_DIR / "PTPRC_1YGU_aligned_repaired_chainA.pdb",
    ),
}


def read_atoms(path: Path, chain: str | None = None) -> list[Atom]:
    atoms = []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if chain is not None and line[21].strip() != chain:
            continue
        resname = line[17:20].strip()
        atomname = line[12:16].strip()
        coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        atoms.append(Atom(line=line, resname=resname, atomname=atomname, coord=coord))
    return atoms


def ca_residues(atoms: list[Atom]) -> list[Residue]:
    return [Residue(a.resname, a.coord) for a in atoms if a.atomname == "CA"]


def needleman_wunsch(a: str, b: str) -> tuple[str, str]:
    match = 2
    mismatch = -1
    gap = -2
    n, m = len(a), len(b)
    score = np.zeros((n + 1, m + 1), dtype=int)
    trace = np.zeros((n + 1, m + 1), dtype=np.int8)
    for i in range(1, n + 1):
        score[i, 0] = score[i - 1, 0] + gap
        trace[i, 0] = 1
    for j in range(1, m + 1):
        score[0, j] = score[0, j - 1] + gap
        trace[0, j] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1, j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            up = score[i - 1, j] + gap
            left = score[i, j - 1] + gap
            best = max(diag, up, left)
            score[i, j] = best
            trace[i, j] = 0 if best == diag else 1 if best == up else 2
    aa, bb = [], []
    i, j = n, m
    while i > 0 or j > 0:
        t = trace[i, j]
        if i > 0 and j > 0 and t == 0:
            aa.append(a[i - 1]); bb.append(b[j - 1])
            i -= 1; j -= 1
        elif i > 0 and (j == 0 or t == 1):
            aa.append(a[i - 1]); bb.append("-")
            i -= 1
        else:
            aa.append("-"); bb.append(b[j - 1])
            j -= 1
    return "".join(reversed(aa)), "".join(reversed(bb))


def matched_ca(original: list[Residue], repaired: list[Residue]) -> tuple[np.ndarray, np.ndarray, int]:
    seq_o = "".join(THREE_TO_ONE.get(r.resname, "X") for r in original)
    seq_r = "".join(THREE_TO_ONE.get(r.resname, "X") for r in repaired)
    aln_o, aln_r = needleman_wunsch(seq_o, seq_r)
    oi = ri = 0
    o_pts, r_pts = [], []
    identities = 0
    for co, cr in zip(aln_o, aln_r):
        current_o = original[oi] if co != "-" else None
        current_r = repaired[ri] if cr != "-" else None
        if current_o is not None and current_r is not None and co == cr:
            o_pts.append(current_o.ca)
            r_pts.append(current_r.ca)
            identities += 1
        if co != "-":
            oi += 1
        if cr != "-":
            ri += 1
    return np.array(o_pts), np.array(r_pts), identities


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    mobile_centroid = mobile.mean(axis=0)
    target_centroid = target.mean(axis=0)
    p = mobile - mobile_centroid
    q = target - target_centroid
    c = p.T @ q
    v, _, wt = np.linalg.svd(c)
    d = np.sign(np.linalg.det(v @ wt))
    u = v @ np.diag([1.0, 1.0, d]) @ wt
    aligned = p @ u + target_centroid
    rmsd = float(np.sqrt(((aligned - target) ** 2).sum(axis=1).mean()))
    return u, target_centroid - mobile_centroid @ u, rmsd


def transform_line(line: str, coord: np.ndarray) -> str:
    return f"{line[:30]}{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}{line[54:]}"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, (original_path, repaired_path, out_path) in PAIRS.items():
        original_atoms = read_atoms(original_path, chain="A")
        repaired_atoms = read_atoms(repaired_path, chain="A")
        target_ca, mobile_ca, identities = matched_ca(ca_residues(original_atoms), ca_residues(repaired_atoms))
        if len(target_ca) < 20:
            raise RuntimeError(f"{name}: too few matched CA atoms ({len(target_ca)})")
        rot, shift, rmsd = kabsch(mobile_ca, target_ca)
        out_lines = []
        for atom in repaired_atoms:
            aligned = atom.coord @ rot + shift
            out_lines.append(transform_line(atom.line, aligned))
        out_lines.extend(["TER", "END"])
        out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8", newline="\n")
        rows.append((name, len(target_ca), identities, rmsd, out_path))
    summary = OUT_DIR / "aligned_repaired_receptors_summary.csv"
    summary.write_text(
        "target_pdb,matched_ca,sequence_identities_for_fit,fit_ca_rmsd_A,aligned_receptor_pdb\n"
        + "\n".join(f"{name},{matched},{identities},{rmsd:.4f},{out_path}" for name, matched, identities, rmsd, out_path in rows)
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(summary)
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
