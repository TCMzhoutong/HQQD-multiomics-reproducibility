from pathlib import Path
import argparse
import math
import os
import re
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
PYMOL = os.environ.get("PYMOL_EXECUTABLE") or shutil.which("pymol") or "pymol"
OBABEL = os.environ.get("OBABEL_EXECUTABLE") or shutil.which("obabel") or "obabel"
BABEL_DATADIR = os.environ.get("BABEL_DATADIR")
VINA = os.environ.get("VINA_EXECUTABLE") or shutil.which("vina") or "vina"


def current_project_path(value):
    """Rebase historical absolute paths onto the current project root."""
    if pd.isna(value):
        return value
    text = str(value).replace("\\", "/")
    markers = [
        "03_receptor_preparation/",
        "06_vina_run1/",
        "07_vina_run2/",
        "08_vina_run3/",
    ]
    lower = text.lower()
    for marker in markers:
        index = lower.find(marker.lower())
        if index >= 0:
            return str(ROOT / Path(text[index:]))
    return value


def run_cmd(cmd: list[str], log_path: Path | None = None) -> int:
    env = os.environ.copy()
    if BABEL_DATADIR:
        env["BABEL_DATADIR"] = BABEL_DATADIR
    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("w", encoding="utf-8") as handle:
            return subprocess.run(cmd, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, env=env).returncode
    return subprocess.run(cmd, cwd=ROOT, env=env).returncode


def parse_score(log_path: Path) -> float | None:
    pattern = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+")
    if not log_path.exists():
        return None
    for line in log_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = pattern.match(line)
        if m:
            return float(m.group(1))
    return None


def atom_xyz_from_pdbqt(path: Path) -> list[tuple[str, float, float, float]]:
    atoms = []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")):
            elem = line[77:79].strip() or line[12:16].strip()[0]
            atoms.append((elem.upper(), float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return atoms


def atom_xyz_from_pdb(path: Path) -> list[tuple[str, float, float, float]]:
    atoms = []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            elem = line[76:78].strip() or line[12:16].strip()[0]
            if elem.upper() != "H":
                atoms.append((elem.upper(), float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return atoms


def kabsch_rmsd(native_atoms, docked_atoms) -> float | None:
    native = [a for a in native_atoms if a[0] != "H"]
    docked = [a for a in docked_atoms if a[0] != "H"]
    n = min(len(native), len(docked))
    if n < 3:
        return None
    p = np.array([[a[1], a[2], a[3]] for a in native[:n]], dtype=float)
    q = np.array([[a[1], a[2], a[3]] for a in docked[:n]], dtype=float)
    # Redocking poses are already in the receptor coordinate frame, so a direct
    # heavy-atom RMSD is the most transparent automated QC here.
    return float(np.sqrt(((p - q) ** 2).sum() / n))


def receptor_atoms(path: Path) -> list[dict]:
    atoms = []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("ATOM"):
            elem = line[76:78].strip() or line[12:16].strip()[0]
            atoms.append(
                {
                    "elem": elem.upper(),
                    "resn": line[17:20].strip(),
                    "chain": line[21].strip(),
                    "resi": line[22:26].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                }
            )
    return atoms


def dist(a, b) -> float:
    return math.sqrt((a[1] - b["x"]) ** 2 + (a[2] - b["y"]) ** 2 + (a[3] - b["z"]) ** 2)


def contact_residues(receptor_pdb: Path, ligand_pdbqt: Path, max_res=12) -> list[dict]:
    lig = atom_xyz_from_pdbqt(ligand_pdbqt)
    rec = receptor_atoms(receptor_pdb)
    by_res = {}
    for la in lig:
        for ra in rec:
            d = dist(la, ra)
            if d <= 4.2:
                key = (ra["resn"], ra["chain"], ra["resi"])
                item = by_res.setdefault(key, {"resn": ra["resn"], "chain": ra["chain"], "resi": ra["resi"], "min_dist": d, "types": set()})
                item["min_dist"] = min(item["min_dist"], d)
                if d <= 3.5 and la[0] in {"N", "O", "S"} and ra["elem"] in {"N", "O", "S"}:
                    item["types"].add("H-bond")
                elif la[0] == "C" and ra["elem"] == "C":
                    item["types"].add("Hydrophobic")
                else:
                    item["types"].add("van der Waals")
    rows = sorted(by_res.values(), key=lambda x: x["min_dist"])[:max_res]
    for r in rows:
        if not r["types"]:
            r["types"].add("van der Waals")
        r["types"] = sorted(r["types"])
    return rows


def make_2d_panel(pair_id: str, ligand_name: str, contacts: list[dict], out: Path) -> None:
    img = Image.new("RGB", (900, 620), "white")
    draw = ImageDraw.Draw(img)
    try:
        font_big = ImageFont.truetype("arial.ttf", 34)
        font = ImageFont.truetype("arial.ttf", 22)
        font_small = ImageFont.truetype("arial.ttf", 18)
    except Exception:
        font_big = font = font_small = None
    draw.text((30, 20), "2D interaction map", fill=(30, 30, 30), font=font_big)
    draw.rounded_rectangle((310, 170, 590, 420), radius=20, outline=(170, 170, 170), width=3)
    draw.text((345, 275), ligand_name[:22], fill=(120, 120, 120), font=font)
    center = (450, 295)
    colors = {"H-bond": (255, 190, 230), "Hydrophobic": (180, 235, 180), "van der Waals": (205, 245, 205)}
    n = max(len(contacts), 1)
    for i, c in enumerate(contacts):
        angle = 2 * math.pi * i / n - math.pi / 2
        x = int(center[0] + 320 * math.cos(angle))
        y = int(center[1] + 210 * math.sin(angle))
        ctype = c["types"][0]
        color = colors.get(ctype, (210, 235, 210))
        draw.line((center[0], center[1], x, y), fill=(230, 150, 220) if ctype == "H-bond" else (175, 220, 175), width=2)
        draw.ellipse((x - 42, y - 42, x + 42, y + 42), fill=color, outline=(160, 160, 160), width=2)
        label = f"{c['resn']}\\n{c['chain']}:{c['resi']}"
        draw.multiline_text((x - 35, y - 24), label, fill=(40, 80, 40), font=font_small, align="center")
    y0 = 505
    for j, (name, color) in enumerate(colors.items()):
        draw.rectangle((40 + j * 250, y0, 70 + j * 250, y0 + 25), fill=color, outline=(120, 120, 120))
        draw.text((80 + j * 250, y0 - 2), name, fill=(30, 30, 30), font=font_small)
    out.parent.mkdir(parents=True, exist_ok=True)
    img.save(out)


def render_pymol_pair(pair: dict, contacts: list[dict], out_prefix: Path) -> tuple[Path, Path]:
    overview = out_prefix.with_name(out_prefix.name + "_overview.png")
    zoom = out_prefix.with_name(out_prefix.name + "_zoom.png")
    res_select = " or ".join([f"(chain {c['chain']} and resi {c['resi']})" for c in contacts[:8]]) or "none"
    pml = out_prefix.with_suffix(".pml")
    pml.parent.mkdir(parents=True, exist_ok=True)
    pml.write_text(
        f"""
load {pair['receptor_pdb']}, receptor
load {pair['out_pdbqt']}, ligand
hide everything
bg_color white
show cartoon, receptor
color aquamarine, receptor
show sticks, ligand
color cyan, ligand
set cartoon_transparency, 0.0
orient receptor
zoom receptor, 5
ray 850, 620
png {overview}, dpi=220
hide sticks, receptor
select pocket_res, {res_select}
show sticks, pocket_res
color red, pocket_res
show cartoon, receptor
set cartoon_transparency, 0.45
show sticks, ligand
color cyan, ligand
label pocket_res and name CA, "%s-%s" % (resn, resi)
set label_size, 18
set label_color, black
orient ligand
zoom ligand or pocket_res, 8
ray 900, 620
png {zoom}, dpi=220
quit
""",
        encoding="utf-8",
    )
    run_cmd([PYMOL, "-cq", str(pml)])
    return overview, zoom


def compose_triptych(pair: dict, overview: Path, zoom: Path, panel2d: Path, out: Path) -> None:
    imgs = [Image.open(p).convert("RGB") for p in [overview, zoom, panel2d]]
    h = 620
    resized = [im.resize((900, h)) for im in imgs]
    canvas = Image.new("RGB", (2700, h), "white")
    x = 0
    for im in resized:
        canvas.paste(im, (x, 0))
        x += 900
    draw = ImageDraw.Draw(canvas)
    try:
        font_big = ImageFont.truetype("arial.ttf", 36)
    except Exception:
        font_big = None
    draw.text((90, 555), pair["target"], fill=(0, 0, 0), font=font_big)
    draw.text((1120, 25), pair["ligand_name"][:32], fill=(0, 0, 0), font=font_big)
    out.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out)


def redocking_validation() -> pd.DataFrame:
    grids = pd.read_csv(ROOT / "03_receptor_preparation" / "grid_boxes.csv")
    rows = []
    native_dir = ROOT / "05_redocking_validation" / "01_native_ligands"
    for native_sdf in native_dir.glob("*_native.sdf"):
        target = native_sdf.stem.split("_")[0]
        grid = grids[grids["target"].eq(target)].iloc[0].to_dict()
        pdbqt = ROOT / "05_redocking_validation" / "02_native_pdbqt" / f"{native_sdf.stem}.pdbqt"
        config = ROOT / "05_redocking_validation" / "03_configs" / f"{native_sdf.stem}.txt"
        out_pose = ROOT / "05_redocking_validation" / "04_poses" / f"{native_sdf.stem}_redock.pdbqt"
        log = ROOT / "05_redocking_validation" / "05_logs" / f"{native_sdf.stem}.log"
        if not pdbqt.exists():
            pdbqt.parent.mkdir(parents=True, exist_ok=True)
            run_cmd([OBABEL, str(native_sdf), "-O", str(pdbqt), "-h", "--partialcharge", "gasteiger"])
        if not pdbqt.exists() or pdbqt.stat().st_size == 0:
            rows.append({"target": target, "native_ligand": native_sdf.name, "redocking_score": None, "rmsd_A": None, "note": "Native ligand PDBQT conversion failed; redocking skipped."})
            continue
        config.parent.mkdir(parents=True, exist_ok=True)
        config.write_text(
            f"""receptor = {grid['receptor_pdbqt']}
ligand = {pdbqt}
center_x = {grid['center_x']}
center_y = {grid['center_y']}
center_z = {grid['center_z']}
size_x = {grid['size_x']}
size_y = {grid['size_y']}
size_z = {grid['size_z']}
exhaustiveness = 10
num_modes = 9
seed = 2026062199
""",
            encoding="utf-8",
        )
        if not out_pose.exists() or parse_score(log) is None:
            out_pose.parent.mkdir(parents=True, exist_ok=True)
            run_cmd([VINA, "--config", str(config), "--out", str(out_pose)], log)
        if not out_pose.exists():
            rows.append({"target": target, "native_ligand": native_sdf.name, "redocking_score": parse_score(log), "rmsd_A": None, "note": "Redocking pose missing."})
            continue
        site_ref = Path(grid["site_reference"])
        rmsd = kabsch_rmsd(atom_xyz_from_pdb(site_ref), atom_xyz_from_pdbqt(out_pose))
        rows.append({"target": target, "native_ligand": native_sdf.name, "redocking_score": parse_score(log), "rmsd_A": rmsd, "note": "Direct heavy-atom RMSD in receptor coordinates; automated redocking QC."})
    rows.append({"target": "PTPRC", "native_ligand": "", "redocking_score": None, "rmsd_A": None, "note": "No small-molecule co-crystal ligand in 1YGU; pTyr peptide used only for grid definition."})
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "05_redocking_validation" / "redocking_validation.csv", index=False, encoding="utf-8-sig")
    return df


def make_tables_and_heatmap() -> tuple[pd.DataFrame, pd.DataFrame]:
    long = pd.read_csv(ROOT / "09_results_tables" / "all_3_runs_long.csv")
    stats = pd.read_csv(ROOT / "09_results_tables" / "final_mean_sd_long.csv")
    active_pair_ids = set(pd.read_csv(ROOT / "01_inputs" / "docking_pairs.csv")["pair_id"])
    long = long[long["pair_id"].isin(active_pair_ids)].copy()
    stats = stats[stats["pair_id"].isin(active_pair_ids)].copy()
    for column in ["out_pdbqt", "receptor_pdb", "ligand_pdbqt", "config", "log"]:
        if column in long.columns:
            long[column] = long[column].map(current_project_path)
    stats.to_csv(ROOT / "09_results_tables" / "final_mean_sd_long.csv", index=False, encoding="utf-8-sig")
    long.to_csv(ROOT / "09_results_tables" / "all_3_runs_long.csv", index=False, encoding="utf-8-sig")
    stats["mean_sd"] = stats.apply(lambda r: f"{r['mean_affinity_kcal_mol']:.2f}±{r['sd_affinity_kcal_mol']:.2f}", axis=1)
    ligand_order = pd.read_csv(ROOT / "01_inputs" / "ligands.csv")[["name", "source"]].rename(columns={"name": "ligand_name"})
    final = ligand_order.merge(stats[["target", "ligand_name", "mean_sd", "mean_affinity_kcal_mol", "sd_affinity_kcal_mol"]], on="ligand_name", how="left")
    wide = final.pivot_table(index=["source", "ligand_name"], columns="target", values="mean_sd", aggfunc="first").reset_index()
    wide = wide.rename(columns={"ligand_name": "compound"})
    for col in ["PTPRC", "CASP1", "BCL2", "MAPK14"]:
        if col not in wide.columns:
            wide[col] = ""
    wide = wide[["source", "compound", "PTPRC", "CASP1", "BCL2", "MAPK14"]]
    redock_path = ROOT / "05_redocking_validation" / "redocking_validation.csv"
    redock = pd.read_csv(redock_path) if redock_path.exists() else redocking_validation()
    with pd.ExcelWriter(ROOT / "09_results_tables" / "final_docking_results_mean_sd.xlsx", engine="openpyxl") as writer:
        wide.to_excel(writer, sheet_name="Final_mean_SD", index=False)
        for run_id in [1, 2, 3]:
            run_df = pd.read_csv(ROOT / "09_results_tables" / f"run{run_id}_raw_scores.csv")
            run_df = run_df[run_df["pair_id"].isin(active_pair_ids)].copy()
            run_df.to_csv(ROOT / "09_results_tables" / f"run{run_id}_raw_scores.csv", index=False, encoding="utf-8-sig")
            run_df.to_excel(writer, sheet_name=f"Run_{run_id}_raw", index=False)
        stats.to_excel(writer, sheet_name="Final_long", index=False)
        redock.to_excel(writer, sheet_name="Redocking_RMSD", index=False)
        pd.read_csv(ROOT / "02_pdb_selection" / "selected_pdbs.csv").to_excel(writer, sheet_name="PDB_selection", index=False)
        pd.read_csv(ROOT / "03_receptor_preparation" / "receptor_preparation_audit.csv").to_excel(writer, sheet_name="Receptor_prep", index=False)
        ligand_prep = pd.read_csv(ROOT / "04_ligand_preparation" / "ligand_preparation_audit.csv")
        active_ligand_ids = set(pd.read_csv(ROOT / "01_inputs" / "ligands.csv")["ligand_id"])
        ligand_prep = ligand_prep[ligand_prep["ligand_id"].isin(active_ligand_ids)].copy()
        ligand_prep.to_excel(writer, sheet_name="Ligand_prep", index=False)

    return stats, long


def make_figures(stats: pd.DataFrame, long: pd.DataFrame) -> None:
    best = long.sort_values("best_affinity_kcal_mol").groupby("pair_id", as_index=False).head(1)
    rows = stats.merge(best[["pair_id", "out_pdbqt", "receptor_pdb", "best_affinity_kcal_mol"]], on="pair_id", how="left")
    for _, row in rows.iterrows():
        out = ROOT / "10_figures_48_pairs" / f"{row['pair_id']}.png"
        if out.exists():
            continue
        contacts = contact_residues(Path(row["receptor_pdb"]), Path(row["out_pdbqt"]))
        prefix = ROOT / "10_figures_48_pairs" / "work" / row["pair_id"]
        overview, zoom = render_pymol_pair(row.to_dict(), contacts, prefix)
        panel2d = prefix.with_name(prefix.name + "_2d.png")
        make_2d_panel(row["pair_id"], row["ligand_name"], contacts, panel2d)
        compose_triptych(row.to_dict(), overview, zoom, panel2d, out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables-only", action="store_true", help="Regenerate active result tables and heatmap only.")
    args = parser.parse_args()
    stats, long = make_tables_and_heatmap()
    subprocess.run([sys.executable, str(ROOT / "13_scripts" / "61_restore_heatmap_matplotlib.py")], check=True)
    if not args.tables_only:
        make_figures(stats, long)
    print(f"Final tables and publication-style heatmap completed for {len(stats)} active docking pairs.")


if __name__ == "__main__":
    main()
