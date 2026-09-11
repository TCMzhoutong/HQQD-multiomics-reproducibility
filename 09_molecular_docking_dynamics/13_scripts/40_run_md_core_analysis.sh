#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -n "${MD_GMX_ENV:-}" ]]; then
  export PATH="$MD_GMX_ENV/bin:$MD_GMX_ENV/bin.AVX2_256:$PATH"
  export LD_LIBRARY_PATH="$MD_GMX_ENV/lib:${LD_LIBRARY_PATH:-}"
fi

SYSTEMS_DIR="$ROOT/12_md_candidates/05_gromacs_systems"
ANALYSIS_ROOT="$ROOT/12_md_candidates/06_md_analysis"
LOG_DIR="$ANALYSIS_ROOT/00_logs"
mkdir -p "$LOG_DIR"

pairs=(
  "MAPK14_3HLL__L01_Pelargonidin"
  "MAPK14_3HLL__L07_Thyroxine"
  "PTPRC_1YGU__L03_5_7_dihydroxy_2_phenyl_4H_chromen_4_one"
  "BCL2_8HTS__L03_5_7_dihydroxy_2_phenyl_4H_chromen_4_one"
)

group_number() {
  local index_file="$1"
  local group_name="$2"
  awk -v target="$group_name" '
    /^\[.*\]$/ {
      name=$0
      gsub(/^\[[[:space:]]*/, "", name)
      gsub(/[[:space:]]*\]$/, "", name)
      if (name == target) {
        print n
        exit
      }
      n++
    }
  ' "$index_file"
}

run_gmx() {
  local log_file="$1"
  shift
  echo "[$(date '+%F %T')] $*" | tee -a "$log_file"
  "$@" >>"$log_file" 2>&1
}

run_pipe_if_missing() {
  local output_file="$1"
  local log_file="$2"
  local input_text="$3"
  shift 3
  if [[ -s "$output_file" ]]; then
    echo "[$(date '+%F %T')] skip existing $output_file" | tee -a "$log_file"
    return
  fi
  printf "%b" "$input_text" | run_gmx "$log_file" "$@"
}

for pair in "${pairs[@]}"; do
  pair_dir="$SYSTEMS_DIR/$pair"
  run_dir="$pair_dir/02_run"
  out_dir="$pair_dir/03_analysis"
  traj_dir="$out_dir/01_trajectory_preprocessed"
  xvg_dir="$out_dir/02_xvg_metrics"
  hbond_dir="$out_dir/03_hbond_strict"
  qc_dir="$out_dir/04_energy_qc"
  dssp_dir="$out_dir/05_secondary_structure"
  mkdir -p "$traj_dir" "$xvg_dir" "$hbond_dir" "$qc_dir" "$dssp_dir"
  log_file="$LOG_DIR/${pair}.analysis.log"
  : > "$log_file"

  echo "== $pair =="
  if [[ ! -s "$run_dir/md.tpr" || ! -s "$run_dir/md.xtc" || ! -s "$run_dir/index.ndx" ]]; then
    echo "Missing md.tpr/md.xtc/index.ndx for $pair" | tee -a "$log_file"
    exit 1
  fi

  protein_unl_group="$(group_number "$run_dir/index.ndx" "Protein_UNL")"
  unl_group="$(group_number "$run_dir/index.ndx" "UNL")"
  protein_group="$(group_number "$run_dir/index.ndx" "Protein")"
  if [[ -z "$protein_unl_group" || -z "$unl_group" || -z "$protein_group" ]]; then
    echo "Could not find Protein/UNL/Protein_UNL groups in index.ndx for $pair" | tee -a "$log_file"
    exit 1
  fi

  cd "$run_dir"

  run_pipe_if_missing "$traj_dir/md_nojump.xtc" "$log_file" "System\n" \
      gmx trjconv -s md.tpr -f md.xtc \
        -o "$traj_dir/md_nojump.xtc" -pbc nojump -n index.ndx

  run_pipe_if_missing "$traj_dir/md_centered.xtc" "$log_file" "${protein_unl_group}\n0\n" \
      gmx trjconv -s md.tpr -f "$traj_dir/md_nojump.xtc" \
        -o "$traj_dir/md_centered.xtc" -pbc mol -ur compact -center -n index.ndx

  run_pipe_if_missing "$traj_dir/md_fit_backbone.xtc" "$log_file" "Backbone\nSystem\n" \
      gmx trjconv -s md.tpr -f "$traj_dir/md_centered.xtc" \
        -o "$traj_dir/md_fit_backbone.xtc" -fit rot+trans -n index.ndx

  run_pipe_if_missing "$xvg_dir/rmsd_protein_backbone.xvg" "$log_file" "Backbone\nBackbone\n" \
      gmx rms -s md.tpr -f "$traj_dir/md_fit_backbone.xtc" \
        -o "$xvg_dir/rmsd_protein_backbone.xvg" -tu ns -n index.ndx

  run_pipe_if_missing "$xvg_dir/rmsd_complex_protein_ligand.xvg" "$log_file" "${protein_unl_group}\n${protein_unl_group}\n" \
      gmx rms -s md.tpr -f "$traj_dir/md_fit_backbone.xtc" \
        -o "$xvg_dir/rmsd_complex_protein_ligand.xvg" -tu ns -n index.ndx

  run_pipe_if_missing "$xvg_dir/rmsd_ligand_backbone_aligned.xvg" "$log_file" "Backbone\n${unl_group}\n" \
      gmx rms -s md.tpr -f "$traj_dir/md_fit_backbone.xtc" \
        -o "$xvg_dir/rmsd_ligand_backbone_aligned.xvg" -tu ns -n index.ndx

  run_pipe_if_missing "$xvg_dir/rmsf_calpha_by_residue.xvg" "$log_file" "C-alpha\n" \
      gmx rmsf -s md.tpr -f "$traj_dir/md_fit_backbone.xtc" \
        -o "$xvg_dir/rmsf_calpha_by_residue.xvg" -res -n index.ndx

  run_pipe_if_missing "$xvg_dir/rg_protein.xvg" "$log_file" "${protein_group}\n" \
      gmx gyrate -s md.tpr -f "$traj_dir/md_fit_backbone.xtc" \
        -o "$xvg_dir/rg_protein.xvg" -n index.ndx

  run_pipe_if_missing "$xvg_dir/sasa_protein_total.xvg" "$log_file" "${protein_group}\n${protein_group}\n" \
      gmx sasa -s md.tpr -f "$traj_dir/md_fit_backbone.xtc" \
        -o "$xvg_dir/sasa_protein_total.xvg" -or "$xvg_dir/sasa_protein_by_residue.xvg" -n index.ndx

  run_pipe_if_missing "$hbond_dir/hbond_protein_ligand_num.xvg" "$log_file" "${protein_group}\n${unl_group}\n" \
      gmx hbond-legacy -s md.tpr -f "$traj_dir/md_fit_backbone.xtc" -n index.ndx \
        -num "$hbond_dir/hbond_protein_ligand_num.xvg" \
        -g "$hbond_dir/hbond_protein_ligand.log" \
        -hbn "$hbond_dir/hbond_protein_ligand.ndx" \
        -tu ns -dt 0.1 -r 0.35 -a 30

  if [[ ! -s "$dssp_dir/secondary_structure_counts.xvg" ]]; then
    run_gmx "$log_file" \
      gmx dssp -s md.tpr -f "$traj_dir/md_fit_backbone.xtc" -n index.ndx \
        -tu ns -dt 0.2 \
        -num "$dssp_dir/secondary_structure_counts.xvg" \
        -o "$dssp_dir/secondary_structure_assignments.dat" \
        -hmode dssp -clear
  else
    echo "[$(date '+%F %T')] skip existing $dssp_dir/secondary_structure_counts.xvg" | tee -a "$log_file"
  fi

  run_pipe_if_missing "$qc_dir/temperature_pressure_density.xvg" "$log_file" "Temperature\nPressure\nDensity\n0\n" \
      gmx energy -f md.edr -o "$qc_dir/temperature_pressure_density.xvg"
done

echo "Core GROMACS analysis finished: $ANALYSIS_ROOT"
