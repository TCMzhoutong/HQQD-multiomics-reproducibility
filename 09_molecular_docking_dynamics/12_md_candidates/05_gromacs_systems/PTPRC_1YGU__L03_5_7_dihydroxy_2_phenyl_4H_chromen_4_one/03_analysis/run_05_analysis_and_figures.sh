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

python3 "$ROOT/13_scripts/24_analyze_md_and_plot.py" --pair-dir "$PAIR_DIR"
