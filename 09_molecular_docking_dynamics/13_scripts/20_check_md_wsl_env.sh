#!/usr/bin/env bash
set -eo pipefail

if [[ -n "${MD_GMX_ENV:-}" ]]; then
  export PATH="$MD_GMX_ENV/bin:$MD_GMX_ENV/bin.AVX2_256:$PATH"
  export LD_LIBRARY_PATH="$MD_GMX_ENV/lib:${LD_LIBRARY_PATH:-}"
fi

echo "== Conda environment =="
if command -v conda >/dev/null 2>&1; then
  conda info --envs | sed -n '1,20p'
else
  echo "Conda is not available on PATH."
fi

echo
echo "== GPU =="
nvidia-smi | sed -n '1,20p' || true

echo
echo "== GROMACS =="
gmx --version | sed -n '1,40p'

echo
echo "== Small molecule / MM-PBSA tools =="
obabel -V
(antechamber -h 2>&1 | sed -n '1,5p') || true
(parmchk2 -h 2>&1 | sed -n '1,5p') || true
(acpype -h 2>&1 | sed -n '1,12p') || true
(gmx_MMPBSA -h 2>&1 | sed -n '1,20p') || true

echo
echo "== Python analysis packages =="
python3 - <<'PY'
import MDAnalysis, pandas, numpy, scipy, matplotlib, seaborn, rdkit
print("MDAnalysis", MDAnalysis.__version__)
print("pandas", pandas.__version__)
print("numpy", numpy.__version__)
print("scipy", scipy.__version__)
print("matplotlib", matplotlib.__version__)
print("seaborn", seaborn.__version__)
print("rdkit", rdkit.__version__)
PY

echo
echo "== DuIvyTools =="
command -v dit
(dit -h 2>&1 | sed -n '1,12p') || true
