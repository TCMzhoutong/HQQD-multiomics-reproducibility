#!/usr/bin/env bash
set -eo pipefail
if [[ -n "${MD_CUDA_ENV:-}" ]]; then
  if ! command -v conda >/dev/null 2>&1; then
    echo "MD_CUDA_ENV is set but conda is not available on PATH." >&2
    exit 1
  fi
  eval "$(conda shell.bash hook)"
  conda activate "$MD_CUDA_ENV"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr
gmx mdrun -deffnm md -nb gpu -pme gpu -bonded gpu -update gpu -ntmpi 1 -ntomp 12 -pin on -v
