from pathlib import Path
import argparse
import math
import re
import shutil
import subprocess
import urllib.error
import urllib.request
import os

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OBABEL = os.environ.get("OBABEL_EXECUTABLE") or shutil.which("obabel") or "obabel"
BABEL_DATADIR = os.environ.get("BABEL_DATADIR")

CHAIN_RULES = {
    "PTPRC": {
        "pdb_id": "1YGU",
        "receptor_chains": ["A"],
        "site_polymer_chains": ["C"],
        "native_ligand_resnames": [],
        "note": "Keep one human CD45/PTPRC chain A; use pTyr peptide chain C for pocket definition. Chain B/D are symmetry-related copy.",
    },
    "CASP1": {
        "pdb_id": "6PZP",
        "receptor_chains": ["A", "B"],
        "site_polymer_chains": [],
        "native_ligand_resnames": ["P7S"],
        "note": "Keep A+B because caspase-1 active site is represented by the p20/p10 enzyme assembly, not a simple duplicate chain.",
    },
    "BCL2": {
        "pdb_id": "8HTS",
        "receptor_chains": ["A"],
        "site_polymer_chains": [],
        "native_ligand_resnames": ["N2L"],
        "note": "Single human BCL2 chain with co-crystal ligand N2L.",
    },
    "MAPK14": {
        "pdb_id": "3HLL",
        "receptor_chains": ["A"],
        "site_polymer_chains": [],
        "native_ligand_resnames": ["I45"],
        "note": "Single human MAPK14 chain; use PH-797804/I45 as the primary pocket ligand.",
    },
}

EXCLUDED_HET = {"HOH", "WAT", "DOD", "SO4", "PO4", "PO2", "CL", "NA", "K", "CA", "MG", "ZN", "MN", "ACT", "EDO", "GOL"}


def slugify(value: str) -> str:
    value = str(value).strip()
    value = re.sub(r"[^A-Za-z0-9]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value[:80] or "unknown"


def run(cmd: list[str], log: Path | None = None) -> None:
    print(" ".join(str(x) for x in cmd))
    env = os.environ.copy()
    if BABEL_DATADIR:
        env["BABEL_DATADIR"] = BABEL_DATADIR
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("w", encoding="utf-8") as handle:
            subprocess.run(cmd, cwd=ROOT, check=True, stdout=handle, stderr=subprocess.STDOUT, env=env)
    else:
        subprocess.run(cmd, cwd=ROOT, check=True, env=env)


def download(url: str, dest: Path) -> bool:
    if dest.exists() and dest.stat().st_size > 0:
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url, timeout=90) as response:
            data = response.read()
        if data:
            dest.write_bytes(data)
            return True
    except urllib.error.URLError as exc:
        print(f"Download failed: {url} ({exc})")
    return False


def pdb_xyz(line: str) -> tuple[float, float, float]:
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


def normalize_atom_line(line: str) -> str:
    resn = line[17:20]
    if line.startswith("HETATM") and resn.strip() == "MSE":
        line = "ATOM  " + line[6:]
        line = line[:17] + "MET" + line[20:]
        if line[12:16].strip() == "SE":
            line = line[:12] + " SD " + line[16:76] + " S  " + line[78:]
    return line


def prepare_tables() -> None:
    ligands = pd.read_csv(ROOT / "01_inputs" / "ligands.csv")
    selected = pd.read_csv(ROOT / "01_inputs" / "targets.csv")
    required_ligand_columns = {"ligand_id", "name", "source", "PubChem_CID", "SMILES"}
    required_target_columns = {"target", "pdb_id"}
    if not required_ligand_columns.issubset(ligands.columns):
        raise ValueError("The canonical ligand table is missing required columns")
    if not required_target_columns.issubset(selected.columns):
        raise ValueError("The canonical target table is missing required columns")
    pairs = []
    for _, target in selected.iterrows():
        for _, ligand in ligands.iterrows():
            pairs.append(
                {
                    "pair_id": f"{target['target']}_{target['pdb_id']}__{ligand['ligand_id']}",
                    "target": target["target"],
                    "pdb_id": target["pdb_id"],
                    "ligand_id": ligand["ligand_id"],
                    "ligand_name": ligand["name"],
                    "source": ligand["source"],
                    "pubchem_cid": ligand.get("PubChem_CID", ""),
                    "smiles": ligand.get("SMILES", ""),
                }
            )
    pd.DataFrame(pairs).to_csv(ROOT / "01_inputs" / "docking_pairs.csv", index=False, encoding="utf-8-sig")


def download_ligands() -> None:
    ligands = pd.read_csv(ROOT / "01_inputs" / "ligands.csv")
    audit_rows = []
    for _, row in ligands.iterrows():
        ligand_id = row["ligand_id"]
        raw_sdf = ROOT / "04_ligand_preparation" / "01_pubchem_sdf" / f"{ligand_id}.sdf"
        cid = row.get("PubChem_CID", "")
        ok = False
        if pd.notna(cid) and str(cid).strip():
            cid_int = str(int(float(cid)))
            for record_type in ["3d", "2d"]:
                ok = download(
                    f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_int}/record/SDF/?record_type={record_type}",
                    raw_sdf,
                )
                if ok:
                    break
        if not ok:
            smi = ROOT / "04_ligand_preparation" / "01_pubchem_sdf" / f"{ligand_id}.smi"
            smi.write_text(str(row["SMILES"]).strip() + f" {ligand_id}\n", encoding="utf-8")
            run([OBABEL, str(smi), "-O", str(raw_sdf), "--gen3d"])

        min_sdf = ROOT / "04_ligand_preparation" / "02_minimized_sdf" / f"{ligand_id}.sdf"
        pdbqt = ROOT / "04_ligand_preparation" / "03_pdbqt" / f"{ligand_id}.pdbqt"
        is_large_flexible = ligand_id.startswith("L11_")
        min_steps = "200" if is_large_flexible else "2500"
        if not min_sdf.exists():
            run([OBABEL, str(raw_sdf), "-O", str(min_sdf), "--gen3d", "--minimize", "--ff", "MMFF94", "--steps", min_steps, "-h"])
        if not pdbqt.exists():
            run([OBABEL, str(min_sdf), "-O", str(pdbqt), "-h", "--partialcharge", "gasteiger"])
        audit_rows.append(
            {
                "ligand_id": ligand_id,
                "ligand_name": row["name"],
                "source": row["source"],
                "pubchem_cid": row.get("PubChem_CID", ""),
                "input_sdf": str(raw_sdf),
                "minimized_sdf": str(min_sdf),
                "pdbqt": str(pdbqt),
                "forcefield": "MMFF94",
                "minimization_steps": int(min_steps),
                "charge": "Gasteiger via OpenBabel",
                "note": "Reduced minimization steps for large flexible ligand." if is_large_flexible else "",
            }
        )
    pd.DataFrame(audit_rows).to_csv(ROOT / "04_ligand_preparation" / "ligand_preparation_audit.csv", index=False, encoding="utf-8-sig")


def download_targets() -> None:
    selected = pd.read_csv(ROOT / "02_pdb_selection" / "selected_pdbs.csv")
    for _, row in selected.iterrows():
        pdb_id = str(row["pdb_id"]).upper()
        target = row["target"]
        dest = ROOT / "03_receptor_preparation" / "01_raw_pdb" / f"{target}_{pdb_id}.pdb"
        download(f"https://files.rcsb.org/download/{pdb_id}.pdb", dest)


def prepare_receptors() -> None:
    selected = pd.read_csv(ROOT / "02_pdb_selection" / "selected_pdbs.csv")
    audit_rows = []
    grid_rows = []
    for _, row in selected.iterrows():
        target = row["target"]
        pdb_id = str(row["pdb_id"]).upper()
        rule = CHAIN_RULES[target]
        raw = ROOT / "03_receptor_preparation" / "01_raw_pdb" / f"{target}_{pdb_id}.pdb"
        clean = ROOT / "03_receptor_preparation" / "02_clean_single_receptor_pdb" / f"{target}_{pdb_id}_clean.pdb"
        site_ref = ROOT / "03_receptor_preparation" / "03_site_reference_ligands" / f"{target}_{pdb_id}_site_reference.pdb"
        native_sdf = ROOT / "05_redocking_validation" / "01_native_ligands" / f"{target}_{pdb_id}_native.sdf"
        receptor_pdbqt = ROOT / "03_receptor_preparation" / "04_receptor_pdbqt" / f"{target}_{pdb_id}.pdbqt"
        lines = raw.read_text(encoding="utf-8", errors="ignore").splitlines(True)
        atom_lines = []
        site_lines = []
        missing = [l.rstrip() for l in lines if l.startswith("REMARK 465")]
        for line in lines:
            norm = normalize_atom_line(line)
            if norm.startswith("ATOM") and norm[21] in rule["receptor_chains"]:
                atom_lines.append(norm)
            elif norm.startswith("ATOM") and norm[21] in rule["site_polymer_chains"]:
                site_lines.append(norm)
            elif norm.startswith("HETATM"):
                resn = norm[17:20].strip()
                chain = norm[21]
                if resn in rule["native_ligand_resnames"]:
                    site_lines.append(norm)
        clean.parent.mkdir(parents=True, exist_ok=True)
        site_ref.parent.mkdir(parents=True, exist_ok=True)
        clean.write_text("".join(atom_lines) + "END\n", encoding="utf-8")
        site_ref.write_text("".join(site_lines) + ("END\n" if site_lines else ""), encoding="utf-8")
        if site_lines and rule["native_ligand_resnames"]:
            run([OBABEL, str(site_ref), "-O", str(native_sdf), "-h"])
        if not receptor_pdbqt.exists():
            run([OBABEL, str(clean), "-O", str(receptor_pdbqt), "-xr", "-h", "--partialcharge", "gasteiger"])

        points = [pdb_xyz(l) for l in site_lines if l.startswith(("ATOM", "HETATM"))]
        if not points:
            points = [pdb_xyz(l) for l in atom_lines if l.startswith("ATOM")]
            grid_source = "protein_center_fallback"
        elif rule["native_ligand_resnames"]:
            grid_source = "native_cocrystal_ligand"
        else:
            grid_source = "native_phosphopeptide_site_reference"
        xs, ys, zs = zip(*points)
        cx, cy, cz = sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)
        sx = min(max(max(xs) - min(xs) + 18.0, 22.0), 38.0 if grid_source != "protein_center_fallback" else 70.0)
        sy = min(max(max(ys) - min(ys) + 18.0, 22.0), 38.0 if grid_source != "protein_center_fallback" else 70.0)
        sz = min(max(max(zs) - min(zs) + 18.0, 22.0), 38.0 if grid_source != "protein_center_fallback" else 70.0)
        grid_rows.append(
            {
                "target": target,
                "pdb_id": pdb_id,
                "receptor_chains": ",".join(rule["receptor_chains"]),
                "site_reference": str(site_ref),
                "receptor_pdb": str(clean),
                "receptor_pdbqt": str(receptor_pdbqt),
                "grid_source": grid_source,
                "center_x": round(cx, 3),
                "center_y": round(cy, 3),
                "center_z": round(cz, 3),
                "size_x": round(sx, 3),
                "size_y": round(sy, 3),
                "size_z": round(sz, 3),
            }
        )
        audit_rows.append(
            {
                "target": target,
                "pdb_id": pdb_id,
                "kept_receptor_chains": ",".join(rule["receptor_chains"]),
                "site_reference_chains_or_ligands": ",".join(rule["site_polymer_chains"] + rule["native_ligand_resnames"]),
                "chain_decision_note": rule["note"],
                "raw_missing_residue_remark465_lines": len(missing),
                "clean_receptor_atom_lines": len(atom_lines),
                "clean_receptor_has_HETATM": False,
                "hydrogenation_charge_tool": "OpenBabel -h --partialcharge gasteiger; note Kollman/ADT not available in current environment.",
                "repair_decision": "Missing/disordered residues audited from REMARK 465. Full homology loop modelling not applied in this automated pass; must be reviewed near pocket before manuscript use.",
            }
        )
    pd.DataFrame(grid_rows).to_csv(ROOT / "03_receptor_preparation" / "grid_boxes.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(audit_rows).to_csv(ROOT / "03_receptor_preparation" / "receptor_preparation_audit.csv", index=False, encoding="utf-8-sig")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tables-only",
        action="store_true",
        help="Regenerate active ligand/target/pair tables without downloading or preparing structures.",
    )
    args = parser.parse_args()
    for folder in [
        "04_ligand_preparation/01_pubchem_sdf",
        "04_ligand_preparation/02_minimized_sdf",
        "04_ligand_preparation/03_pdbqt",
        "03_receptor_preparation/01_raw_pdb",
        "03_receptor_preparation/02_clean_single_receptor_pdb",
        "03_receptor_preparation/03_site_reference_ligands",
        "03_receptor_preparation/04_receptor_pdbqt",
        "05_redocking_validation/01_native_ligands",
    ]:
        (ROOT / folder).mkdir(parents=True, exist_ok=True)
    prepare_tables()
    if args.tables_only:
        print("Active input tables completed.")
        return
    download_ligands()
    download_targets()
    prepare_receptors()
    print("Strict input preparation completed.")


if __name__ == "__main__":
    main()
