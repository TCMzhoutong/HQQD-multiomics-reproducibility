# Validation report

Validation date: 2026-09-10

## Manuscript-aligned checks

All 17 automated checks in `workflow/validate_package.py` passed:

- 32 retained HQQD-reversed metabolite annotations;
- five prioritised metabolites with the expected feature identifiers;
- 30 three-way common targets;
- 30 STRING input nodes, 28 connected nodes, and 105 edges;
- the expected nine hub genes;
- 44 docking combinations comprising 11 ligands and four targets;
- three docking repetitions per combination;
- four retained molecular-dynamics systems;
- no prespecified excluded annotation in the active docking or molecular-dynamics manifests;
- no file larger than 100 MB; and
- no personal absolute path in distributable text files.

## Static syntax checks

- Python: 27 files parsed successfully with Python 3; no syntax failures.
- R: 87 files parsed successfully with R 4.5.2; no syntax failures.
- Shell: 28 files passed `bash -n` using Git Bash; no syntax failures.

## Module-level checks

- Serum metabolomics: 38 initial reverse-candidate annotations, six prespecified exclusions, and 32 retained annotations; UTF-8 labels verified in both ion modes.
- Multi-omics correlation: 18 taxa by 32 metabolites, 237 significant Spearman associations, 15 taxon-phenotype and 27 metabolite-phenotype Mantel associations, and the expected five prioritised metabolites.
- Network pharmacology: five metabolite inputs, 30 three-way common targets, a 28-node/105-edge connected STRING network, CytoNCA 28-to-20 selection, and the expected nine-gene hub intersection.

## Repository-hygiene checks

- All 487 code and text-data files are readable as UTF-8.
- No CJK text remains in distributable code, documentation, or text-data files.
- No relative filename contains non-ASCII characters.
- No empty, backup, temporary, `publication_package`, or `pocket` directory remains.
- No software executable, Cytoscape session, Adobe Illustrator file, molecular-dynamics trajectory, or checkpoint is bundled.
- The largest bundled file is below GitHub's 100 MB per-file limit.

The source negative-ion metabolite table contained six malformed byte sequences confined to two long, non-prioritised annotation labels. In the clean package these were converted to question-mark placeholders during UTF-8 normalisation; identifiers and quantitative measurements were not changed.

## Scope of validation

This report confirms package structure, static syntax, portability checks, and agreement of compact reference outputs with the current manuscript. It is not a claim that every analysis was rerun end-to-end during package curation. Full reruns require the external primary/public datasets and software listed in `metadata/EXTERNAL_DATA_MANIFEST.tsv` and `environment/SOFTWARE.md`.
