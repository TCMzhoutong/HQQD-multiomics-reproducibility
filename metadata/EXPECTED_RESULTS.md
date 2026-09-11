# Expected results for the manuscript-aligned release

These checks define the active analysis state. A rerun should reproduce these identities and counts, allowing only numerical tolerance appropriate to software versions and stochastic algorithms.

## Serum metabolomics

- Retained HQQD-reversed serum metabolites: **32**.
- Prioritised metabolites: **5**:
  - Thyroxine (`pos_7892`)
  - N-Oleoyl Alanine (`pos_8061`)
  - 8-Amino-7-oxononanoate (`pos_4360`)
  - PE(P-18:0/18:1(12Z)-2OH(9,10)) (`neg_12260`)
  - 12,13-DiHODE (`neg_12234`)

## Multi-omics correlation analysis

- HQQD-responsive taxa included in the principal correlation matrix: **18**.
- Retained metabolite annotations included in the matrix: **32**.
- Significant taxon–metabolite Spearman associations: **237** (**164 positive**, **73 negative**).
- Metabolite–phenotype Mantel tests: **160**, of which **27** are significant.
- Significant taxon–phenotype Mantel associations: **15**.

## Network pharmacology

- Metabolite targets: **317**.
- Three-way common targets shared by HQQD-constituent targets, prioritised-metabolite targets, and bacterial-pneumonia targets: **30**.
- STRING network: **28 connected nodes** and **105 edges**.
- Disconnected targets: `CES1` and `FUCA1`.
- Hub genes: `AKT1`, `BCL2`, `CASP1`, `MAPK14`, `MMP9`, `PPARG`, `PTGS2`, `PTPRC`, and `SIRT1`.
- Targets prioritised by downstream immune and single-cell analyses: `PTPRC`, `CASP1`, `BCL2`, and `MAPK14`.

## Molecular docking and dynamics

- Docking ligands: **11** (six HQQD constituents and five prioritised serum metabolites).
- Docking targets: **4** (`PTPRC`, `CASP1`, `BCL2`, and `MAPK14`).
- Docking combinations: **44**, with **three Vina repeats per combination**.
- Retained molecular-dynamics systems: **4**:
  - MAPK14–Pelargonidin
  - MAPK14–Thyroxine
  - PTPRC–5,7-dihydroxy-2-phenyl-4H-chromen-4-one
  - BCL2–5,7-dihydroxy-2-phenyl-4H-chromen-4-one

## Version note

The former 31-target list differs from the current 30-target list by removal of `MGAM` only. The nine hub genes are unchanged.
