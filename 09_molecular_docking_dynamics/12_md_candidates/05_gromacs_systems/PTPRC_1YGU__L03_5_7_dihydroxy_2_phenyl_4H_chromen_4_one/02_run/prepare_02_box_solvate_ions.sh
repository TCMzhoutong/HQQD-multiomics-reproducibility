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

python3 "$ROOT/13_scripts/23_build_gromacs_complex.py" \
  --pair-dir "$PAIR_DIR"

cd 02_run
cp ../00_inputs/*.mdp .

echo "== Box, solvation, ions =="
gmx editconf -f complex.gro -o newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp -f em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 -rmin 0.5

echo "== EM/NVT/NPT/MD commands prepared =="
