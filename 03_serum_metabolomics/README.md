# Untargeted serum metabolomics

Run `Rscript run_revised_pipeline_utf8.R` from this directory. The revised workflow applies the current annotation-exclusion list, retains 32 HQQD-reversed metabolite features, and writes the compact outputs used by downstream multi-omics analysis.

`excluded_metabolite_annotations.csv` records the six excluded annotations. Raw mass-spectrometry files and vendor project files are not bundled and must be deposited separately. Run `python validate_revised_outputs.py` to check the revised outputs.

For cross-platform portability, `metabolites_exp_neg.csv` was normalised to UTF-8. Six malformed byte sequences confined to two long, non-prioritised annotation labels were replaced with question-mark placeholders; feature identifiers and all quantitative values were preserved.
