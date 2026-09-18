# Mendelian randomisation

Run `Rscript main.R` from this directory after setting `ACTIVE_EXPOSURES` and `ACTIVE_OUTCOME` in `config.R` for each analysis listed in `analysis_jobs.csv`, and recording the public data locations in `data_sources_config.R`. PLINK is resolved from `PLINK_EXECUTABLE` or `PATH`; set `PLINK_REFERENCE_PREFIX` and `LD_REFERENCE_DIRECTORY` for local linkage-disequilibrium resources.

Large GWAS, pQTL, FinnGen, and linkage-disequilibrium reference files are not bundled. Compact manuscript-relevant result tables are retained in `results_reference/`.

The active manuscript-aligned jobs comprise gut-microbiota exposures from MiBioGen, Rühlemann, and Qin; thyroxine from the METSIM, Chen, and Rhee plasma-metabolite GWAS resources; FinnGen R12 bacterial pneumonia as the primary outcome; FinnGen R7 bacterial pneumonia not elsewhere classified as the supplementary outcome; and UKB-PPP plasma proteins for pQTL-MR. Dataset citations and accessions are recorded in `../metadata/EXTERNAL_DATA_MANIFEST.tsv`.
