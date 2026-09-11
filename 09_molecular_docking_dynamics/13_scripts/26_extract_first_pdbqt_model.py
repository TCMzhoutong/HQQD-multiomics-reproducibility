#!/usr/bin/env python3
"""Extract the first MODEL from a Vina PDBQT output."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--outfile", required=True)
    args = ap.parse_args()
    lines = Path(args.infile).read_text(errors="ignore").splitlines()
    out = []
    in_first = False
    seen_model = False
    for line in lines:
        if line.startswith("MODEL"):
            if seen_model:
                break
            seen_model = True
            in_first = True
            out.append(line)
            continue
        if line.startswith("ENDMDL") and in_first:
            out.append(line)
            break
        if in_first or not seen_model:
            if line.startswith(("ATOM", "HETATM", "ROOT", "BRANCH", "ENDBRANCH", "TORSDOF", "REMARK")):
                out.append(line)
    Path(args.outfile).write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
