#!/usr/bin/env bash
set -eo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "Starting production MD for the four retained systems."

while IFS= read -r script; do
  echo "== Running $script =="
  bash "$script"
done < <(find "$ROOT/12_md_candidates/05_gromacs_systems" -path "*/02_run/run_04_production_100ns.sh" | sort)

echo "All production MD stages finished."
