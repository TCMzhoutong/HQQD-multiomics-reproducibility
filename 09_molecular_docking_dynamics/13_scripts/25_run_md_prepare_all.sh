#!/usr/bin/env bash
set -eo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -n "${MD_GMX_ENV:-}" ]]; then
  export PATH="$MD_GMX_ENV/bin:$MD_GMX_ENV/bin.AVX2_256:$PATH"
  export LD_LIBRARY_PATH="$MD_GMX_ENV/lib:${LD_LIBRARY_PATH:-}"
fi
cd "$ROOT"

python3 13_scripts/21_prepare_md_systems.py

while IFS= read -r script; do
  echo "== Running $script =="
  bash "$script"
done < <(find "$ROOT/12_md_candidates/05_gromacs_systems" -path "*/01_prepare/prepare_01_topology_and_qc.sh" | sort)

echo "All prepare/QC stages finished."
