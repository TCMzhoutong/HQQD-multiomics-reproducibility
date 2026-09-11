#!/usr/bin/env bash
set -eo pipefail
if [[ -n "${MD_GMX_ENV:-}" ]]; then
  export PATH="$MD_GMX_ENV/bin:$MD_GMX_ENV/bin.AVX2_256:$PATH"
  export LD_LIBRARY_PATH="$MD_GMX_ENV/lib:${LD_LIBRARY_PATH:-}"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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

run_gpu_equilibration() {
  local deffnm="$1"
  echo "== Running $deffnm with GPU acceleration =="
  if gmx mdrun -deffnm "$deffnm" -nb gpu -pme gpu -bonded gpu -update gpu -v > "${deffnm}_gpu_update_attempt.log" 2>&1; then
    cat "${deffnm}_gpu_update_attempt.log"
    return 0
  fi

  if grep -Eiq 'Update groups can not|update groups cannot|GPU update|update gpu|not supported with.*update' "${deffnm}_gpu_update_attempt.log"; then
    echo "GPU update was not accepted for $deffnm; retrying with GPU nonbonded/PME/bonded and CPU update."
    rm -f "${deffnm}.log" "${deffnm}.edr" "${deffnm}.trr" "${deffnm}.xtc" "${deffnm}.gro" "${deffnm}.cpt" "${deffnm}_prev.cpt"
    gmx mdrun -deffnm "$deffnm" -nb gpu -pme gpu -bonded gpu -v
  else
    cat "${deffnm}_gpu_update_attempt.log"
    echo "GPU equilibration failed for a reason other than unsupported GPU update."
    return 1
  fi
}

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
if [[ -n "${MD_CUDA_ENV:-}" ]]; then
  if ! command -v conda >/dev/null 2>&1; then
    echo "MD_CUDA_ENV is set but conda is not available on PATH." >&2
    exit 1
  fi
  eval "$(conda shell.bash hook)"
  conda activate "$MD_CUDA_ENV"
fi

run_gpu_equilibration nvt

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr
run_gpu_equilibration npt

echo "Equilibration finished. Production MD is not started by this script."
