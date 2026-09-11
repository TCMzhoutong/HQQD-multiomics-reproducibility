"""Process manual Cytoscape exports, calculate the intersection, and publish it."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path(__file__).resolve().parent
PYTHON = Path(sys.executable)
RSCRIPT = shutil.which("Rscript") or "Rscript"


def run(command: list[str]) -> None:
    print("RUN:", " ".join(command), flush=True)
    subprocess.run(command, cwd=BASE, check=True)


def main() -> None:
    run([str(PYTHON), "prepare_cytoscape_hub_exports.py"])
    run([str(RSCRIPT), "venn_hub_genes.r"])
    run([str(PYTHON), "validate_hub_selection_outputs.py"])

    print("Processed Cytoscape exports and validated the hub-gene intersection.")


if __name__ == "__main__":
    main()
