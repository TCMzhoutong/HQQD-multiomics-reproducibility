#!/usr/bin/env bash
set -eo pipefail
if [[ -n "${MD_GMX_ENV:-}" ]]; then
  export PATH="$MD_GMX_ENV/bin:$MD_GMX_ENV/bin.AVX2_256:$PATH"
  export LD_LIBRARY_PATH="$MD_GMX_ENV/lib:${LD_LIBRARY_PATH:-}"
fi


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAIR_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT="$(cd "$PAIR_DIR/../../.." && pwd)"
cd "$PAIR_DIR"

mkdir -p 01_prepare 02_run 03_analysis
cd 01_prepare

echo "== Docked ligand coordinates =="
python3 "$ROOT/13_scripts/26_extract_first_pdbqt_model.py" \
  --infile ../00_inputs/ligand_docked.pdbqt \
  --outfile ligand_docked_model1.pdbqt
obabel ligand_docked_model1.pdbqt -O ligand_docked.pdb -h

echo "== Ligand topology from audited SDF plus docked coordinates =="
python3 "$ROOT/13_scripts/31_make_docked_ligand_sdf.py" \
  --topology-sdf ../00_inputs/ligand_topology_source.sdf \
  --docked-pdb ligand_docked.pdb \
  --out ligand_docked_for_acpype.sdf
if [[ "Thyroxine" == "Thyroxine" && "-1" == "-1" ]]; then
  python3 "$ROOT/13_scripts/33_fix_thyroxine_monoanion_sdf.py" \
    --infile ligand_docked_for_acpype.sdf \
    --outfile ligand_docked_for_acpype.sdf \
    --site carboxylate
fi
rm -rf UNL.acpype UNL_*
acpype -i ligand_docked_for_acpype.sdf -c bcc -n -1 -a gaff -b UNL

echo "== Protein topology =="
pdb4amber -i ../00_inputs/protein_clean.pdb -o protein_repaired.pdb -y --add-missing-atoms -l protein_repair_pdb4amber.log
gmx pdb2gmx -f protein_repaired.pdb -o protein.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh

echo "== Ligand atom-order QC =="
python3 "$ROOT/13_scripts/22_qc_ligand_mapping.py" \
  --topology-sdf ../00_inputs/ligand_topology_source.sdf \
  --docked-pdb ligand_docked.pdb \
  --out ligand_mapping_qc.csv

echo "Preparation stage finished. Review 01_prepare/ligand_mapping_qc.csv before system merge."
