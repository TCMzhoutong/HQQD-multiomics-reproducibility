#!/usr/bin/env python3
"""Prepare audited GROMACS MD work folders for selected docking poses.

This script intentionally stops before long production MD. It creates one
folder per selected pair, writes reproducible WSL run scripts, and records the
charge/protonation decisions used for ligand parameterization.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SELECTED_CSV = ROOT / "12_md_candidates" / "01_selection" / "md_retained_4_pairs.csv"
CHARGE_CSV = ROOT / "12_md_candidates" / "01_selection" / "ligand_charge_protonation_decisions.csv"
SYSTEMS_DIR = ROOT / "12_md_candidates" / "05_gromacs_systems"
MDP_DIR = ROOT / "mdp"
ALIGNED_REPAIRED_DIR = ROOT / "12_md_candidates" / "09_aligned_repaired_receptors"

REPAIRED_RECEPTORS = {
    "MAPK14_3HLL": ALIGNED_REPAIRED_DIR / "MAPK14_3HLL_aligned_repaired_chainA.pdb",
    "BCL2_8HTS": ALIGNED_REPAIRED_DIR / "BCL2_8HTS_aligned_repaired_chainA.pdb",
    "PTPRC_1YGU": ALIGNED_REPAIRED_DIR / "PTPRC_1YGU_aligned_repaired_chainA.pdb",
}

MD_GMX_ENV = """if [[ -n "${MD_GMX_ENV:-}" ]]; then
  export PATH="$MD_GMX_ENV/bin:$MD_GMX_ENV/bin.AVX2_256:$PATH"
  export LD_LIBRARY_PATH="$MD_GMX_ENV/lib:${LD_LIBRARY_PATH:-}"
fi
"""

MD_CUDA_ENV = """if [[ -n "${MD_CUDA_ENV:-}" ]]; then
  if ! command -v conda >/dev/null 2>&1; then
    echo "MD_CUDA_ENV is set but conda is not available on PATH." >&2
    exit 1
  fi
  eval "$(conda shell.bash hook)"
  conda activate "$MD_CUDA_ENV"
fi
"""


def copy_required(src: Path, dst: Path) -> None:
    if not src.exists():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def assert_repaired_receptor_is_aligned(target_pdb: str, receptor_src: Path) -> None:
    if target_pdb not in REPAIRED_RECEPTORS:
        return
    receptor_resolved = receptor_src.resolve()
    aligned_root = ALIGNED_REPAIRED_DIR.resolve()
    if aligned_root not in receptor_resolved.parents:
        raise RuntimeError(
            f"{target_pdb} is a modeled/repaired receptor, but the selected input is not aligned "
            f"to the original docking frame: {receptor_resolved}"
        )
    summary = ALIGNED_REPAIRED_DIR / "aligned_repaired_receptors_summary.csv"
    if not summary.exists():
        raise FileNotFoundError(
            f"Missing alignment summary: {summary}. Run 13_scripts/32_align_repaired_receptors_to_docking_frame.py first."
        )


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def main() -> None:
    selected = pd.read_csv(SELECTED_CSV)
    charges = pd.read_csv(CHARGE_CSV).set_index("ligand_id")
    systems_rows = []

    for row in selected.to_dict("records"):
        pair_id = row["pair_id"]
        ligand_id = row["ligand_id"]
        ligand_name = row["ligand_name"]
        target_pdb = f"{row['target']}_{row['pdb_id']}"
        charge = charges.loc[ligand_id].to_dict()

        pair_dir = SYSTEMS_DIR / pair_id
        inputs = pair_dir / "00_inputs"
        prep = pair_dir / "01_prepare"
        run = pair_dir / "02_run"
        analysis = pair_dir / "03_analysis"
        for d in [inputs, prep, run, analysis]:
            d.mkdir(parents=True, exist_ok=True)

        default_receptor_src = ROOT / "12_md_candidates" / "04_clean_receptor_pdb" / f"{target_pdb}_clean.pdb"
        receptor_src = REPAIRED_RECEPTORS.get(target_pdb, default_receptor_src)
        receptor_source = "swissmodel_repaired_clean_chainA_aligned_to_docking_frame" if target_pdb in REPAIRED_RECEPTORS else "original_clean_receptor"
        assert_repaired_receptor_is_aligned(target_pdb, receptor_src)
        pose_src = ROOT / "12_md_candidates" / "03_representative_docked_pdbqt" / f"{pair_id}_run{int(row['representative_run_id'])}_representative_out.pdbqt"
        ligand_sdf_src = ROOT / str(charge["topology_source"])

        receptor_dst = inputs / "protein_clean.pdb"
        pose_dst = inputs / "ligand_docked.pdbqt"
        ligand_sdf_dst = inputs / "ligand_topology_source.sdf"
        copy_required(receptor_src, receptor_dst)
        copy_required(pose_src, pose_dst)
        copy_required(ligand_sdf_src, ligand_sdf_dst)
        for mdp in ["em.mdp", "nvt.mdp", "npt.mdp", "md.mdp"]:
            copy_required(MDP_DIR / mdp, inputs / mdp)

        charge_row = {
            "pair_id": pair_id,
            "ligand_id": ligand_id,
            "ligand_name": row["ligand_name"],
            "net_charge": int(charge["net_charge"]),
            "protonation_state": charge["protonation_state"],
            "charge_method": charge["charge_method"],
            "force_field": charge["force_field"],
            "status": charge["status"],
            "note": charge["note"],
        }
        with (inputs / "ligand_charge_decision.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(charge_row))
            writer.writeheader()
            writer.writerow(charge_row)

        prepare_sh = f"""#!/usr/bin/env bash
set -eo pipefail
{MD_GMX_ENV}

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
PAIR_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT="$(cd "$PAIR_DIR/../../.." && pwd)"
cd "$PAIR_DIR"

mkdir -p 01_prepare 02_run 03_analysis
cd 01_prepare

echo "== Docked ligand coordinates =="
python3 "$ROOT/13_scripts/26_extract_first_pdbqt_model.py" \\
  --infile ../00_inputs/ligand_docked.pdbqt \\
  --outfile ligand_docked_model1.pdbqt
obabel ligand_docked_model1.pdbqt -O ligand_docked.pdb -h

echo "== Ligand topology from audited SDF plus docked coordinates =="
python3 "$ROOT/13_scripts/31_make_docked_ligand_sdf.py" \\
  --topology-sdf ../00_inputs/ligand_topology_source.sdf \\
  --docked-pdb ligand_docked.pdb \\
  --out ligand_docked_for_acpype.sdf
if [[ "{ligand_name}" == "Thyroxine" && "{int(charge['net_charge'])}" == "-1" ]]; then
  python3 "$ROOT/13_scripts/33_fix_thyroxine_monoanion_sdf.py" \\
    --infile ligand_docked_for_acpype.sdf \\
    --outfile ligand_docked_for_acpype.sdf \\
    --site carboxylate
fi
rm -rf UNL.acpype UNL_*
acpype -i ligand_docked_for_acpype.sdf -c bcc -n {int(charge['net_charge'])} -a gaff -b UNL

echo "== Protein topology =="
pdb4amber -i ../00_inputs/protein_clean.pdb -o protein_repaired.pdb -y --add-missing-atoms -l protein_repair_pdb4amber.log
gmx pdb2gmx -f protein_repaired.pdb -o protein.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh

echo "== Ligand atom-order QC =="
python3 "$ROOT/13_scripts/22_qc_ligand_mapping.py" \\
  --topology-sdf ../00_inputs/ligand_topology_source.sdf \\
  --docked-pdb ligand_docked.pdb \\
  --out ligand_mapping_qc.csv

echo "Preparation stage finished. Review 01_prepare/ligand_mapping_qc.csv before system merge."
"""
        write_text(prep / "prepare_01_topology_and_qc.sh", prepare_sh)

        build_sh = f"""#!/usr/bin/env bash
set -eo pipefail
{MD_GMX_ENV}

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
PAIR_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT="$(cd "$PAIR_DIR/../../.." && pwd)"
cd "$PAIR_DIR"

python3 "$ROOT/13_scripts/23_build_gromacs_complex.py" \\
  --pair-dir "$PAIR_DIR"

cd 02_run
cp ../00_inputs/*.mdp .

echo "== Box, solvation, ions =="
gmx editconf -f complex.gro -o newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp -f em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 -rmin 0.5

echo "== EM/NVT/NPT/MD commands prepared =="
"""
        write_text(run / "prepare_02_box_solvate_ions.sh", build_sh)

        run_short_sh = f"""#!/usr/bin/env bash
set -eo pipefail
{MD_GMX_ENV}
SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
cd "$SCRIPT_DIR"

gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -v

gmx make_ndx -f em.gro -o index.ndx <<'EOF'
13
name 21 UNL
1 | 21
name 22 Protein_UNL
q
EOF

gmx genrestr -f ../01_prepare/UNL.acpype/UNL_GMX.gro -o posre_lig.itp -fc 1000 1000 1000 <<'EOF'
0
EOF

run_gpu_equilibration() {{
  local deffnm="$1"
  echo "== Running $deffnm with GPU acceleration =="
  if gmx mdrun -deffnm "$deffnm" -nb gpu -pme gpu -bonded gpu -update gpu -v > "${{deffnm}}_gpu_update_attempt.log" 2>&1; then
    cat "${{deffnm}}_gpu_update_attempt.log"
    return 0
  fi

  if grep -Eiq 'Update groups can not|update groups cannot|GPU update|update gpu|not supported with.*update' "${{deffnm}}_gpu_update_attempt.log"; then
    echo "GPU update was not accepted for $deffnm; retrying with GPU nonbonded/PME/bonded and CPU update."
    rm -f "${{deffnm}}.log" "${{deffnm}}.edr" "${{deffnm}}.trr" "${{deffnm}}.xtc" "${{deffnm}}.gro" "${{deffnm}}.cpt" "${{deffnm}}_prev.cpt"
    gmx mdrun -deffnm "$deffnm" -nb gpu -pme gpu -bonded gpu -v
  else
    cat "${{deffnm}}_gpu_update_attempt.log"
    echo "GPU equilibration failed for a reason other than unsupported GPU update."
    return 1
  fi
}}

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
{MD_CUDA_ENV}
run_gpu_equilibration nvt

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr
run_gpu_equilibration npt

echo "Equilibration finished. Production MD is not started by this script."
"""
        write_text(run / "run_03_em_nvt_npt.sh", run_short_sh)

        prod_sh = f"""#!/usr/bin/env bash
set -eo pipefail
{MD_CUDA_ENV}
SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
cd "$SCRIPT_DIR"

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr
gmx mdrun -deffnm md -nb gpu -pme gpu -bonded gpu -update gpu -ntmpi 1 -ntomp 12 -pin on -v
"""
        write_text(run / "run_04_production_100ns.sh", prod_sh)

        analysis_sh = f"""#!/usr/bin/env bash
set -eo pipefail
{MD_GMX_ENV}
SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
PAIR_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT="$(cd "$PAIR_DIR/../../.." && pwd)"
cd "$PAIR_DIR"

python3 "$ROOT/13_scripts/24_analyze_md_and_plot.py" --pair-dir "$PAIR_DIR"
"""
        write_text(analysis / "run_05_analysis_and_figures.sh", analysis_sh)

        systems_rows.append({
            "pair_id": pair_id,
            "pair_dir": pair_dir.relative_to(ROOT).as_posix(),
            "receptor_source": receptor_source,
            "receptor_input": receptor_src.relative_to(ROOT).as_posix(),
            "net_charge": int(charge["net_charge"]),
            "prepare_script": (prep / "prepare_01_topology_and_qc.sh").relative_to(ROOT).as_posix(),
            "build_script": (run / "prepare_02_box_solvate_ions.sh").relative_to(ROOT).as_posix(),
            "equilibration_script": (run / "run_03_em_nvt_npt.sh").relative_to(ROOT).as_posix(),
            "production_script": (run / "run_04_production_100ns.sh").relative_to(ROOT).as_posix(),
            "analysis_script": (analysis / "run_05_analysis_and_figures.sh").relative_to(ROOT).as_posix(),
        })

    summary = pd.DataFrame(systems_rows)
    summary.to_csv(SYSTEMS_DIR / "md_systems_manifest.csv", index=False, encoding="utf-8-sig")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
