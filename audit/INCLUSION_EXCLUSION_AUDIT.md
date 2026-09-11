# Inclusion and exclusion audit

## Source protection

The original analysis project was treated as read-only. No source file was deleted, renamed, or edited while constructing this package.

## Inclusion rule

Files were included only when they serve one of four functions:

1. executable analysis code used by a manuscript-reported analysis;
2. a compact, distributable input needed by that code;
3. a compact reference table used to verify the rerun; or
4. metadata needed to locate externally hosted primary/public data.

## Global exclusions

- root and subproject `publication_package` directories;
- manuscript-layout assets, assembled publication figures, Adobe Illustrator files, and figure-production helpers;
- Chinese narrative result reports and administrative spreadsheets;
- pathology-report documents and extracted report images;
- backups, archives, temporary files, IDE settings, caches, and compiled bytecode;
- vendored R/Python libraries and unrelated colocalisation analyses;
- software binaries such as Vina executables and vendor DLLs;
- duplicated raw/processed files already represented by the canonical input;
- large public GWAS/GEO datasets and genomic reference panels that should be retrieved from their original repositories;
- molecular-dynamics trajectories, checkpoints, and bulky visualisation exports;
- annotations excluded during metabolite curation and their downstream analysis artefacts.

## Module decisions

| Module | Included | Excluded |
|---|---|---|
| Animal experiment | Canonical six-sample-per-group phenotype table, analysis script, compact statistical reference tables | Duplicate summary workbooks, pathology report, extracted images, FCS files, publication figures |
| Metagenomics | Analysis scripts, compact abundance matrices, key statistical/contribution tables | Primary reads, 200 MB gene-abundance table, large gene lists, publication assets |
| Serum metabolomics | Revised pipeline, processed intensity/annotation tables, exclusion list, 32-metabolite reference table | Primary instrument files, plots, publication package, superseded outputs |
| Multi-omics correlation | Revised scripts and final 32-metabolite/5-candidate tables | Duplicated upstream inputs, vendored package copy, figures and narrative reports |
| Mendelian randomisation | Reusable functions, compact final intermediates/results, public-data manifest | Multi-gigabyte GWAS files, LD-reference binaries, mediation analyses not reported in the manuscript, backups, colocalisation workflow |
| Network pharmacology | Current pipeline scripts, processed constituent/target inputs, 30-target/28-node/9-hub reference tables | Vendor raw files, notebooks with temporary paths, Cytoscape sessions, colocalisation and publication assets |
| Immune infiltration | Analysis scripts, batch-corrected expression matrix and compact result tables | R workspace objects, publication figures, narrative reports, unlicensed LM22 redistribution |
| Single-cell transcriptomics | Current scripts, compact result tables and GEO manifest | Raw GEO matrices, multi-gigabyte RDS/MTX objects, archived reruns and publication figures |
| Docking and dynamics | Active input definitions, prepared structures, active score tables, reproducibility scripts, four-system compact summaries | Vina binary, publication visualisations, excluded-analysis artefacts, raw MD trajectories/checkpoints |

`FILE_MANIFEST_SHA256.tsv` records every included file after final validation, except the manifest itself.
