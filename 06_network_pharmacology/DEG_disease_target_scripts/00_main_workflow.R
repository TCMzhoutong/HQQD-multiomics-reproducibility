# ==============================================================================
# ==============================================================================
# 
# 
# 
#      - GSE40012_series_matrix.txt.gz
#      - GSE20346_series_matrix.txt.gz
#   
#      source('00_main_workflow.R')
#   
#      source('scripts/01_load_data.R')
#      source('scripts/02_sample_selection.R')
#      ...
# 
# 
# ==============================================================================

cat("
==============================================================================
        GEO differential-expression analysis - workflow started
==============================================================================

Analysis scope:
  - Datasets: GSE40012 & GSE20346
  - Contrast: Bacterial pneumonia vs healthy control
  - Platform: GPL6947 (Illumina HumanHT-12 V3.0)

Analysis steps:
  1. Load data
  2. Select samples
  3. Map probes to genes
  4. Correct batch effects
  5. Differential-expression analysis
  6. PCA visualisation
  7. Volcano-plot visualisation
  8. Heatmap visualisation

==============================================================================
\n")

start_time <- Sys.time()
cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

# ==============================================================================
# ==============================================================================
cat(">>> Step1: Load data...\n")
cat("-------------------------------------\n")
source('scripts/01_load_data.R')
cat("\n")

# ==============================================================================
# ==============================================================================
cat(">>> Step2: Select samples...\n")
cat("-------------------------------------\n")
source('scripts/02_sample_selection.R')
cat("\n")

# ==============================================================================
# ==============================================================================
cat(">>> Step3: Map probes to genes...\n")
cat("-------------------------------------\n")
source('scripts/03_probe_to_gene.R')
cat("\n")

# ==============================================================================
# ==============================================================================
cat(">>> Step4: Correct batch effects...\n")
cat("-------------------------------------\n")
source('scripts/04_batch_correction.R')
cat("\n")

# ==============================================================================
# ==============================================================================
cat(">>> Step5: Differential-expression analysis...\n")
cat("-------------------------------------\n")
source('scripts/05_DEG_analysis.R')
cat("\n")

# ==============================================================================
# ==============================================================================
cat(">>> Step6: PCA visualisation...\n")
cat("-------------------------------------\n")
source('scripts/06_PCA_plots.R')
cat("\n")

# ==============================================================================
# ==============================================================================
cat(">>> Step7: Volcano-plot visualisation...\n")
cat("-------------------------------------\n")
source('scripts/07_volcano_plot.R')
cat("\n")

# ==============================================================================
# ==============================================================================
cat(">>> Step8: Heatmap visualisation...\n")
cat("-------------------------------------\n")
source('scripts/08_heatmap.R')
cat("\n")

# ==============================================================================
# ==============================================================================
end_time <- Sys.time()
elapsed_time <- end_time - start_time

cat("
==============================================================================
                    All analyses complete.
==============================================================================

Elapsed time: ", format(elapsed_time), "

Generated outputs:

[Data files]
  01_raw_data/
    - GSE40012_data.RData
    - GSE20346_data.RData
  
  02_processed_data/
    - GSE40012_selected.RData
    - GSE20346_selected.RData
    - data_ready_for_merge.RData
  
  03_batch_corrected/
    - data_before_combat.RData
    - data_after_combat.RData
    - expression_matrix_combat.csv
    - sample_metadata.csv
  
  04_DEG_results/
    - all_genes_results.csv
    - significant_DEGs.csv
    - DEG_results.RData

[Figure files]
  05_figures/
    - Figure_AB_PCA_batch_correction.pdf/png
    - Figure_A_PCA_before_correction.pdf/png
    - Figure_B_PCA_after_correction.pdf/png
    - Supplementary_PCA_by_group.pdf/png
    - Figure_C_Volcano_plot.pdf/png
    - Figure_C_Volcano_plot_simple.pdf/png
    - Figure_D_Heatmap.pdf/png
    - Supplementary_Heatmap_clustered.pdf/png
    - Top50_DEGs_expression.csv
    - Top50_DEGs_statistics.csv
    - Top50_DEGs_for_volcano_labels.csv

==============================================================================
", sep = "")

cat("\nEnd time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nAnalysis complete. All results were saved.\n\n")
