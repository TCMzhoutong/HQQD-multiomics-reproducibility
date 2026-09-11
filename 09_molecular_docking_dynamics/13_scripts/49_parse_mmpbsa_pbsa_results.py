#!/usr/bin/env python3
"""Parse unified MM-PBSA inp=1 outputs and write paper-style tables."""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "12_md_candidates" / "06_md_analysis"
PBSA_ROOT = OUT / "09_mmpbsa_pbsa_final20ns"

SYSTEMS = {
    "MAPK14_3HLL__L01_Pelargonidin": ("MAPK14-Pelargonidin", "#0072B2"),
    "MAPK14_3HLL__L07_Thyroxine": ("MAPK14-Thyroxine", "#D55E00"),
    "PTPRC_1YGU__L03_5_7_dihydroxy_2_phenyl_4H_chromen_4_one": ("PTPRC-Flavone", "#009E73"),
    "BCL2_8HTS__L03_5_7_dihydroxy_2_phenyl_4H_chromen_4_one": ("BCL2-Flavone", "#CC79A7"),
}

PAPER_ROWS = [
    ("Delta VDWAALS", r"$\Delta E_{vdw}$", "VDWAALS"),
    ("Delta EEL", r"$\Delta E_{elec}$", "EEL"),
    ("Delta EPB", r"$\Delta G_{PB}$", "EPB"),
    ("Delta ENPOLAR", r"$\Delta G_{nonpolar}$", "ENPOLAR"),
    ("Delta TOTAL", r"$\Delta G_{bind}$ (MM/PBSA)", "TOTAL"),
]


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.labelsize": 9,
            "axes.titlesize": 9.5,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 7.5,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.dpi": 300,
            "savefig.dpi": 600,
        }
    )


def save_figure(fig: plt.Figure, out_stem: Path) -> None:
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in [".png", ".svg", ".pdf"]:
        fig.savefig(out_stem.with_suffix(ext), bbox_inches="tight")


def parse_results(path: Path, pair_id: str) -> tuple[pd.DataFrame, int | None]:
    text = path.read_text(errors="ignore")
    frames = None
    m = re.search(r"Calculations performed using\s+(\d+)\s+complex frames", text)
    if m:
        frames = int(m.group(1))

    in_delta = False
    rows = []
    for line in text.splitlines():
        if "Delta (Complex - Receptor - Ligand)" in line:
            in_delta = True
            continue
        if not in_delta:
            continue
        stripped = line.strip()
        if not stripped or stripped.startswith("-") or "Energy Component" in stripped:
            continue
        clean = stripped.lstrip("Δ")
        parts = clean.split()
        if len(parts) < 6:
            continue
        comp = parts[0]
        if not re.search(r"[A-Z]", comp):
            continue
        try:
            vals = [float(x) for x in parts[1:6]]
        except ValueError:
            continue
        rows.append(
            {
                "pair_id": pair_id,
                "system": SYSTEMS[pair_id][0],
                "component": comp,
                "average_kcal_mol": vals[0],
                "sd_prop": vals[1],
                "sd": vals[2],
                "sem_prop": vals[3],
                "sem": vals[4],
            }
        )
    return pd.DataFrame(rows), frames


def make_paper_table(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for raw_name, paper_name, component in PAPER_ROWS:
        row = {
            "Original output (gmx_MMPBSA)": raw_name,
            "Paper-style term": paper_name,
        }
        for pair_id, (label, _) in SYSTEMS.items():
            sub = summary[(summary["pair_id"] == pair_id) & (summary["component"] == component)]
            if sub.empty:
                row[label] = "NA"
            else:
                avg = float(sub.iloc[0]["average_kcal_mol"])
                sem = float(sub.iloc[0]["sem"])
                row[label] = f"{avg:.2f} ± {sem:.2f}"
        rows.append(row)
    return pd.DataFrame(rows)


def plot_total(summary: pd.DataFrame) -> None:
    total = summary[summary["component"] == "TOTAL"].copy()
    total["order"] = total["pair_id"].map({p: i for i, p in enumerate(SYSTEMS)})
    total = total.sort_values("order")
    x = np.arange(len(total))
    colors = [SYSTEMS[p][1] for p in total["pair_id"]]
    labels = [SYSTEMS[p][0].replace("-", "\n", 1) for p in total["pair_id"]]

    fig, ax = plt.subplots(figsize=(6.1, 3.4))
    ax.bar(x, total["average_kcal_mol"], yerr=total["sem"], capsize=3, color=colors, edgecolor="#333333", linewidth=0.5)
    ax.axhline(0, color="#222222", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("MM-PBSA ΔGbind (kcal/mol)")
    ax.set_title("MM-PBSA binding free energy, 80-100 ns")
    ax.grid(True, axis="y", color="#E6E8EB", linewidth=0.55)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save_figure(fig, OUT / "04_figures" / "09_mmpbsa_pbsa_total_binding_energy")
    plt.close(fig)


def plot_components(summary: pd.DataFrame) -> None:
    components = [row[2] for row in PAPER_ROWS]
    df = summary[summary["component"].isin(components)].copy()
    systems = list(SYSTEMS)
    x = np.arange(len(components))
    width = 0.15
    fig, ax = plt.subplots(figsize=(8.2, 3.9))
    for i, pair_id in enumerate(systems):
        sub = df[df["pair_id"] == pair_id].set_index("component").reindex(components)
        ax.bar(
            x + (i - (len(systems) - 1) / 2) * width,
            sub["average_kcal_mol"],
            width=width,
            label=SYSTEMS[pair_id][0],
            color=SYSTEMS[pair_id][1],
        )
    ax.axhline(0, color="#222222", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(["VDWAALS", "EEL", "EPB", "ENPOLAR", "TOTAL"], rotation=25, ha="right")
    ax.set_ylabel("Energy component (kcal/mol)")
    ax.set_title("MM-PBSA energy components")
    ax.legend(frameon=False, ncol=2)
    ax.grid(True, axis="y", color="#E6E8EB", linewidth=0.55)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save_figure(fig, OUT / "04_figures" / "10_mmpbsa_pbsa_energy_components")
    plt.close(fig)


def main() -> None:
    configure_style()
    tables_dir = OUT / "03_tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    status_rows = []
    parts = []
    for pair_id, (label, _) in SYSTEMS.items():
        pair_dir = PBSA_ROOT / pair_id
        result_file = pair_dir / "FINAL_RESULTS_MMPBSA.dat"
        failed = (pair_dir / "FAILED.flag").exists()
        frames = None
        if result_file.exists() and not failed:
            df, frames = parse_results(result_file, pair_id)
            if not df.empty:
                parts.append(df)
        status_rows.append(
            {
                "pair_id": pair_id,
                "system": label,
                "status": "failed" if failed else ("complete" if result_file.exists() else "missing"),
                "frames": frames,
            }
        )

    summary = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    paper = make_paper_table(summary)
    status = pd.DataFrame(status_rows)

    summary.to_csv(tables_dir / "mmpbsa_pbsa_energy_components.csv", index=False)
    status.to_csv(tables_dir / "mmpbsa_pbsa_job_status.csv", index=False)
    paper.to_csv(tables_dir / "mmpbsa_pbsa_paper_style_energy_table.csv", index=False, encoding="utf-8-sig")
    with pd.ExcelWriter(tables_dir / "mmpbsa_pbsa_analysis_tables.xlsx", engine="openpyxl") as writer:
        status.to_excel(writer, sheet_name="job_status", index=False)
        summary.to_excel(writer, sheet_name="energy_components", index=False)
        paper.to_excel(writer, sheet_name="paper_style_table", index=False)

    if not summary.empty:
        plot_total(summary)
        plot_components(summary)

    report = OUT / "05_reports" / "mmpbsa_pbsa_analysis_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "\n".join(
            [
                "# MM-PBSA inp=1 analysis report",
                "",
                "- Method: gmx_MMPBSA v1.5.0.3, Poisson-Boltzmann solver in sander.",
                "- PB nonpolar model: `inp=1`.",
                "- Snapshots: 21 frames from 80-100 ns at 1 ns intervals.",
                "- Ionic strength: `istrng=0.150`.",
                "- This unified PBSA rerun avoids mixing GBSA and PBSA while avoiding the default PBSA dispersion term dominating total Delta G.",
                "- Units: kcal/mol. More negative Delta Gbind indicates stronger predicted binding.",
                "",
                "## Job status",
                "",
                "```",
                status.to_string(index=False),
                "```",
                "",
                "## Paper-style table",
                "",
                "```",
                paper.to_string(index=False),
                "```",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(f"Finished MM-PBSA inp=1 parsing and tables under {OUT}")


if __name__ == "__main__":
    main()
