#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SYSTEMS_DIR="$ROOT/12_md_candidates/05_gromacs_systems"
OUT_ROOT="$ROOT/12_md_candidates/06_md_analysis/09_mmpbsa_pbsa_final20ns"

if [[ -n "${MD_GMX_ENV:-}" ]]; then
  export PATH="$PATH:$MD_GMX_ENV/bin:$MD_GMX_ENV/bin.AVX2_256"
  export LD_LIBRARY_PATH="$MD_GMX_ENV/lib:${LD_LIBRARY_PATH:-}"
fi

pairs=(
  "MAPK14_3HLL__L01_Pelargonidin"
  "MAPK14_3HLL__L07_Thyroxine"
  "PTPRC_1YGU__L03_5_7_dihydroxy_2_phenyl_4H_chromen_4_one"
  "BCL2_8HTS__L03_5_7_dihydroxy_2_phenyl_4H_chromen_4_one"
)

mkdir -p "$OUT_ROOT"

for pair in "${pairs[@]}"; do
  echo "== MM-PBSA inp=1 $pair =="
  pair_dir="$SYSTEMS_DIR/$pair"
  run_dir="$pair_dir/02_run"
  analysis_dir="$pair_dir/03_analysis"
  out_dir="$OUT_ROOT/$pair"
  mkdir -p "$out_dir"
  cd "$out_dir"

  if [[ -s FINAL_RESULTS_MMPBSA.dat && ! -s FAILED.flag ]]; then
    echo "Skip existing MM-PBSA inp=1 result for $pair"
    continue
  fi

  rm -rf _GMXMMPBSA_* FINAL_RESULTS_MMPBSA.dat gmx_MMPBSA.log mmpbsa_pbsa.in FAILED.flag run_gmx_MMPBSA.stdout.log mmpbsa_80_100ns_dt1ns.xtc

  printf "System\n" | gmx trjconv \
    -s "$run_dir/md.tpr" \
    -f "$analysis_dir/01_trajectory_preprocessed/md_fit_backbone.xtc" \
    -o "$out_dir/mmpbsa_80_100ns_dt1ns.xtc" \
    -b 80000 -e 100000 -dt 1000 \
    -n "$run_dir/index.ndx"

  cat > mmpbsa_pbsa.in <<EOF
&general
  sys_name="$pair",
  startframe=1,
  endframe=999999,
  interval=1,
  verbose=1,
  keep_files=0,
  temperature=300,
/
&pb
  istrng=0.150,
  inp=1,
/
EOF

  if ! gmx_MMPBSA -O -i mmpbsa_pbsa.in \
    -cs "$run_dir/md.tpr" \
    -ci "$run_dir/index.ndx" \
    -cg 1 21 \
    -ct "$out_dir/mmpbsa_80_100ns_dt1ns.xtc" \
    -cp "$run_dir/topol.top" \
    -o "$out_dir/FINAL_RESULTS_MMPBSA.dat" \
    -nogui | tee "$out_dir/run_gmx_MMPBSA.stdout.log"; then
    echo "MM-PBSA inp=1 failed for $pair; see $out_dir/run_gmx_MMPBSA.stdout.log" | tee "$out_dir/FAILED.flag"
    continue
  fi
done

echo "All MM-PBSA inp=1 jobs finished: $OUT_ROOT"
