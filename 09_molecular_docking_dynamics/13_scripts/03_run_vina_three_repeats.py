from pathlib import Path
import os
import re
import shutil
import subprocess

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
VINA = os.environ.get("VINA_EXECUTABLE") or shutil.which("vina") or "vina"
BABEL_DATADIR = os.environ.get("BABEL_DATADIR")

RUNS = {
    1: {"folder": "06_vina_run1", "seed_base": 2026062101},
    2: {"folder": "07_vina_run2", "seed_base": 2026062202},
    3: {"folder": "08_vina_run3", "seed_base": 2026062303},
}


def run_cmd(cmd: list[str], log_path: Path | None = None) -> int:
    env = os.environ.copy()
    if BABEL_DATADIR:
        env["BABEL_DATADIR"] = BABEL_DATADIR
    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("w", encoding="utf-8") as handle:
            return subprocess.run(cmd, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, env=env).returncode
    return subprocess.run(cmd, cwd=ROOT, env=env).returncode


def parse_score(log_path: Path) -> float | None:
    if not log_path.exists():
        return None
    pattern = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+")
    for line in log_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = pattern.match(line)
        if match:
            return float(match.group(1))
    return None


def write_config(path: Path, receptor: str, ligand: str, grid: dict, seed: int, exhaustiveness: int = 10) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"""receptor = {receptor}
ligand = {ligand}

center_x = {grid['center_x']}
center_y = {grid['center_y']}
center_z = {grid['center_z']}

size_x = {grid['size_x']}
size_y = {grid['size_y']}
size_z = {grid['size_z']}

exhaustiveness = {exhaustiveness}
num_modes = 9
energy_range = 4
seed = {seed}
""",
        encoding="utf-8",
    )


def build_job_manifest() -> pd.DataFrame:
    pairs = pd.read_csv(ROOT / "01_inputs" / "docking_pairs.csv")
    grids = pd.read_csv(ROOT / "03_receptor_preparation" / "grid_boxes.csv")
    jobs = []
    for i, pair in pairs.iterrows():
        grid = grids[grids["target"].eq(pair["target"])].iloc[0].to_dict()
        ligand = ROOT / "04_ligand_preparation" / "03_pdbqt" / f"{pair['ligand_id']}.pdbqt"
        ligand_relative = ligand.relative_to(ROOT).as_posix()
        large_flexible = pair["ligand_id"].startswith("L11_")
        for run_id, run_info in RUNS.items():
            folder = ROOT / run_info["folder"]
            config = folder / "01_configs" / f"{pair['pair_id']}.txt"
            out_pose = folder / "02_poses" / f"{pair['pair_id']}_out.pdbqt"
            log = folder / "03_logs" / f"{pair['pair_id']}.log"
            seed = run_info["seed_base"] + i
            write_config(config, grid["receptor_pdbqt"], ligand_relative, grid, seed, exhaustiveness=10)
            jobs.append(
                {
                    **pair.to_dict(),
                    "run_id": run_id,
                    "seed": seed,
                    "exhaustiveness": 10,
                    "large_flexible_ligand": large_flexible,
                    "receptor_pdb": grid["receptor_pdb"],
                    "receptor_pdbqt": grid["receptor_pdbqt"],
                    "site_reference": grid["site_reference"],
                    "grid_source": grid["grid_source"],
                    "center_x": grid["center_x"],
                    "center_y": grid["center_y"],
                    "center_z": grid["center_z"],
                    "size_x": grid["size_x"],
                    "size_y": grid["size_y"],
                    "size_z": grid["size_z"],
                    "config_path": config.relative_to(ROOT).as_posix(),
                    "out_pdbqt": out_pose.relative_to(ROOT).as_posix(),
                    "log_path": log.relative_to(ROOT).as_posix(),
                }
            )
    manifest = pd.DataFrame(jobs)
    manifest.to_csv(ROOT / "09_results_tables" / "vina_three_run_job_manifest.csv", index=False, encoding="utf-8-sig")
    return manifest


def run_docking_jobs(manifest: pd.DataFrame) -> None:
    for _, job in manifest.iterrows():
        out_pose = ROOT / Path(job["out_pdbqt"])
        log = ROOT / Path(job["log_path"])
        if out_pose.exists() and out_pose.stat().st_size > 0 and parse_score(log) is not None:
            continue
        print(f"Run {job['run_id']} {job['pair_id']}")
        out_pose.parent.mkdir(parents=True, exist_ok=True)
        run_cmd([VINA, "--config", str(ROOT / Path(job["config_path"])), "--out", str(out_pose)], log)


def summarize_runs(manifest: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, job in manifest.iterrows():
        rows.append({**job.to_dict(), "best_affinity_kcal_mol": parse_score(ROOT / Path(job["log_path"]))})
    scored = pd.DataFrame(rows)
    scored.to_csv(ROOT / "09_results_tables" / "all_3_runs_long.csv", index=False, encoding="utf-8-sig")
    for run_id in RUNS:
        scored[scored["run_id"].eq(run_id)].to_csv(
            ROOT / "09_results_tables" / f"run{run_id}_raw_scores.csv", index=False, encoding="utf-8-sig"
        )
    stats = (
        scored.groupby(["pair_id", "target", "pdb_id", "ligand_id", "ligand_name", "source"], as_index=False)
        .agg(
            mean_affinity_kcal_mol=("best_affinity_kcal_mol", "mean"),
            sd_affinity_kcal_mol=("best_affinity_kcal_mol", "std"),
            n_runs=("best_affinity_kcal_mol", "count"),
        )
        .sort_values(["target", "mean_affinity_kcal_mol"])
    )
    stats.to_csv(ROOT / "09_results_tables" / "final_mean_sd_long.csv", index=False, encoding="utf-8-sig")
    return scored


def main() -> None:
    manifest = build_job_manifest()
    run_docking_jobs(manifest)
    scored = summarize_runs(manifest)
    done = scored["best_affinity_kcal_mol"].notna().sum()
    print(f"Scored {done}/{len(scored)} docking runs.")


if __name__ == "__main__":
    main()
