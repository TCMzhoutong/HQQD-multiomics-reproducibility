# HQQD reproducibility package

This repository contains the analysis code, compact analysis inputs, and reference result tables supporting the manuscript **“Haoqin Qingdan Decoction ameliorates *Klebsiella pneumoniae*-induced pneumonia: An integrated study based on multi-omics and network pharmacology.”**

## Scope

The package is organised by analytical module:

1. animal experiment;
2. shotgun metagenomics;
3. untargeted serum metabolomics;
4. multi-omics correlation analysis;
5. Mendelian randomisation;
6. network pharmacology;
7. immune-infiltration analysis;
8. single-cell transcriptomic analysis; and
9. molecular docking and molecular dynamics simulation.

Each module contains only the files needed to understand or reproduce the reported analysis. Manuscript-layout files, publication figures, Chinese narrative reports, software binaries, temporary files, backups, and superseded analyses are intentionally excluded.

## Data organisation

Each numbered module contains a `README.md` describing its entry scripts, bundled inputs, reference outputs, and external-data requirements. Some modules preserve the original relative directory layout because their scripts depend on it; the package therefore does not force every module into an artificial common folder pattern.

- Compact distributable inputs and reference result tables are bundled with the relevant module.
- `metadata/EXTERNAL_DATA_MANIFEST.tsv` records primary or public datasets that are not bundled.
- `metadata/EXPECTED_RESULTS.md` defines the numerical and identity checks for the current manuscript.
- `environment/` records the required R, Python, and external software dependencies.
- `workflow/validate_package.py` checks the principal manuscript-aligned outputs and portability constraints.
- `audit/` documents inclusion, exclusion, unresolved author inputs, and file integrity.

This repository is intentionally GitHub-friendly. Large primary data and public reference datasets must be deposited in, or downloaded from, suitable domain repositories rather than committed to GitHub.

## Reproduction order

Run modules in numerical order. Modules 04 and 06 consume outputs from earlier modules. Module-specific README files describe their inputs, outputs, and any external software requirements.

## Important current-analysis constraints

- Thirty-two serum metabolite annotations were retained after the prespecified annotation exclusions.
- Five serum metabolites were prioritised by the multi-omics correlation workflow.
- Network pharmacology identified 30 three-way common targets; 28 formed the connected STRING network, and nine hub genes were retained.
- Molecular docking comprises 44 combinations of 11 ligands and four targets, each run three times.
- Molecular dynamics results are reported for four retained systems.

## Repositories and accessions

The processed shotgun metagenomic data and associated analysis files have been deposited in OMIX under accession **OMIX020530**. Other external datasets and repository details requiring confirmation are listed in `metadata/EXTERNAL_DATA_MANIFEST.tsv` and `audit/AUTHOR_INPUT_NEEDED.md`.

## Licence

No licence has been assigned in this draft package. Add an author-approved open-source licence before public release.
