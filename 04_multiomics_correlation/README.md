# Multi-omics correlation analysis

Run `Rscript run_correlation_pipeline_revised.R` from this directory. The workflow reads the canonical animal table and the revised metabolomics outputs directly from modules 01 and 03, avoiding duplicate input copies.

The analysis prioritises five metabolites associated with HQQD-responsive microbial taxa and phenotype modules. Run `python validate_revised_outputs.py` to verify the expected candidate set.
