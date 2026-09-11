#!/usr/bin/env python3
"""Run standard MD analyses and draw publication-style figures.

This script expects completed `02_run/md.tpr` and `02_run/md.xtc`.
It writes raw XVG/XPM outputs plus PNG/SVG figures to `03_analysis`.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def run(cmd: str, cwd: Path) -> None:
    print(cmd)
    subprocess.run(cmd, cwd=cwd, shell=True, check=True, executable="/bin/bash")


def read_xvg(path: Path) -> pd.DataFrame:
    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        if not line.strip() or line.startswith(("#", "@")):
            continue
        rows.append([float(x) for x in line.split()])
    return pd.DataFrame(rows)


def plot_xvg(path: Path, out_stem: Path, xlabel: str, ylabel: str, title: str) -> None:
    df = read_xvg(path)
    fig, ax = plt.subplots(figsize=(5.2, 3.2), dpi=300)
    ax.plot(df.iloc[:, 0], df.iloc[:, 1], lw=1.2, color="#2F6F9F")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, color="#e7e7e7", lw=0.5)
    fig.tight_layout()
    fig.savefig(out_stem.with_suffix(".png"))
    fig.savefig(out_stem.with_suffix(".svg"))
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair-dir", required=True)
    args = ap.parse_args()
    pair_dir = Path(args.pair_dir)
    run_dir = pair_dir / "02_run"
    out = pair_dir / "03_analysis"
    out.mkdir(exist_ok=True)

    if not (run_dir / "md.tpr").exists() or not (run_dir / "md.xtc").exists():
        raise SystemExit("md.tpr/md.xtc not found. Run production MD first.")

    run("printf 'System\\n' | gmx trjconv -s md.tpr -f md.xtc -o ../03_analysis/step1_nojump.xtc -pbc nojump -n index.ndx", run_dir)
    run("printf 'Protein_UNL\\nSystem\\n' | gmx trjconv -s md.tpr -f ../03_analysis/step1_nojump.xtc -o ../03_analysis/md_fixed.xtc -center -fit rot+trans -n index.ndx", run_dir)
    run("printf 'Protein\\nProtein\\n' | gmx rms -s md.tpr -f ../03_analysis/md_fixed.xtc -o ../03_analysis/rmsd_protein.xvg -tu ns -n index.ndx", run_dir)
    run("printf 'Protein\\nUNL\\n' | gmx rms -s md.tpr -f ../03_analysis/md_fixed.xtc -o ../03_analysis/rmsd_ligand.xvg -tu ns -n index.ndx", run_dir)
    run("printf 'Protein\\n' | gmx rmsf -s md.tpr -f ../03_analysis/md_fixed.xtc -o ../03_analysis/rmsf_protein.xvg -res -n index.ndx", run_dir)
    run("printf 'Protein\\n' | gmx gyrate -s md.tpr -f ../03_analysis/md_fixed.xtc -o ../03_analysis/rg_protein.xvg -n index.ndx", run_dir)
    run("printf 'Protein\\nProtein\\n' | gmx sasa -s md.tpr -f ../03_analysis/md_fixed.xtc -o ../03_analysis/sasa_total.xvg -or ../03_analysis/sasa_residue.xvg -n index.ndx", run_dir)
    run("printf 'Protein\\nUNL\\n' | gmx hbond -s md.tpr -f ../03_analysis/md_fixed.xtc -num ../03_analysis/hbond_num.xvg -n index.ndx", run_dir)
    run("printf 'Protein\\nUNL\\n' | gmx mindist -s md.tpr -f ../03_analysis/md_fixed.xtc -n index.ndx -d 0.35 -on ../03_analysis/contacts_0p35nm.xvg", run_dir)

    plot_xvg(out / "rmsd_protein.xvg", out / "RMSD_Protein", "Time (ns)", "RMSD (nm)", "Protein RMSD")
    plot_xvg(out / "rmsd_ligand.xvg", out / "RMSD_Ligand", "Time (ns)", "RMSD (nm)", "Ligand RMSD")
    plot_xvg(out / "rmsf_protein.xvg", out / "RMSF_Protein", "Residue", "RMSF (nm)", "Protein RMSF")
    plot_xvg(out / "rg_protein.xvg", out / "Rg_Protein", "Time (ns)", "Rg (nm)", "Radius of gyration")
    plot_xvg(out / "sasa_total.xvg", out / "SASA_Protein", "Time (ns)", "SASA (nm²)", "Protein SASA")
    plot_xvg(out / "hbond_num.xvg", out / "Hbond_Complex", "Time (ns)", "H-bonds", "Protein-ligand H-bonds")
    plot_xvg(out / "contacts_0p35nm.xvg", out / "Contacts_0p35nm", "Time (ns)", "Contacts", "Close contacts within 0.35 nm")


if __name__ == "__main__":
    main()
