from pathlib import Path
import shutil
import subprocess
import sys


BASE = Path(__file__).resolve().parent
PYTHON = Path(sys.executable)
RSCRIPT = shutil.which("Rscript") or "Rscript"


def run(command):
    print("RUN:", " ".join(map(str, command)))
    subprocess.run(command, cwd=BASE, check=True)


def main():
    run([PYTHON, BASE / "prepare_metabolite_network_inputs.py"])
    run([PYTHON, BASE / "rerun_metabolite_swisstarget.py"])
    run([RSCRIPT, BASE / "venn_disease_ingredient_metabolite.r"])
    run([PYTHON, BASE / "fetch_string_ppi.py"])
    run([PYTHON, BASE / "validate_pre_cytoscape_outputs.py"])


if __name__ == "__main__":
    main()
