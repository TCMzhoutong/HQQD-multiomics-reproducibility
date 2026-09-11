# ============================================================================
# ============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("sva", quietly = TRUE)) {
  BiocManager::install("sva")
}

library(sva)

script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0 || script_dir == "") {
  project_root <- getwd()
  if (basename(project_root) == "scripts") {
    project_root <- dirname(project_root)
  }
} else {
  project_root <- dirname(script_dir)
}

processed_data_dir <- file.path(project_root, "02_processed_data")
batch_corrected_dir <- file.path(project_root, "03_batch_corrected")

if (!dir.exists(batch_corrected_dir)) {
  dir.create(batch_corrected_dir, recursive = TRUE)
}

cat("=== Load data ===\n")

load(file.path(processed_data_dir, "data_ready_for_merge.RData"))

cat("\n=== Merge expression matrices ===\n")

expr_40012_common <- expr_40012_common[common_genes, ]
expr_20346_common <- expr_20346_common[common_genes, ]

expr_combined <- cbind(expr_40012_common, expr_20346_common)

cat("Merged expression-matrix dimensions:", dim(expr_combined), "\n")
cat("- Genes:", nrow(expr_combined), "\n")
cat("- Samples:", ncol(expr_combined), "\n")

cat("\n=== Create metadata ===\n")

batch <- c(
  rep("GSE40012", ncol(expr_40012_common)),
  rep("GSE20346", ncol(expr_20346_common))
)

group <- c(group_40012_selected, group_20346_selected)

sample_ids <- colnames(expr_combined)

metadata <- data.frame(
  SampleID = sample_ids,
  Batch = batch,
  Group = group,
  stringsAsFactors = FALSE
)

cat("Metadata:\n")
print(table(metadata$Batch, metadata$Group))

cat("\n=== Save data before batch correction ===\n")

expr_before_combat <- expr_combined

save(expr_before_combat, metadata, 
     file = file.path(batch_corrected_dir, "data_before_combat.RData"))

cat("Pre-correction data saved.\n")

cat("\n=== Run ComBat batch correction ===\n")

mod <- model.matrix(~ Group, data = metadata)

cat("Design matrix:\n")
print(head(mod))

cat("\nRunning ComBat...\n")

expr_combat <- ComBat(
  dat = expr_combined,
  batch = batch,
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

cat("\nComBat complete.\n")
cat("Corrected expression-matrix dimensions:", dim(expr_combat), "\n")

cat("\n=== Check corrected data ===\n")

cat("\nExpression-value range before correction:\n")
cat("Min:", round(min(expr_before_combat), 3), "\n")
cat("Max:", round(max(expr_before_combat), 3), "\n")
cat("Mean:", round(mean(expr_before_combat), 3), "\n")

cat("\nExpression-value range after correction:\n")
cat("Min:", round(min(expr_combat), 3), "\n")
cat("Max:", round(max(expr_combat), 3), "\n")
cat("Mean:", round(mean(expr_combat), 3), "\n")

cat("\n=== Save batch-corrected data ===\n")

expr_after_combat <- expr_combat

save(expr_after_combat, metadata,
     file = file.path(batch_corrected_dir, "data_after_combat.RData"))

write.csv(expr_after_combat, 
          file.path(batch_corrected_dir, "expression_matrix_combat.csv"))
write.csv(metadata, 
          file.path(batch_corrected_dir, "sample_metadata.csv"), 
          row.names = FALSE)

cat("Batch-corrected data saved to 03_batch_corrected.\n")

cat("\n
==========================================================
              Batch-correction summary
==========================================================

Input data:
  - GSE40012Samples: ", ncol(expr_40012_common), "
  - GSE20346Samples: ", ncol(expr_20346_common), "
  - Shared genes: ", length(common_genes), "

Sample groups:
")
print(table(metadata$Group))
cat("
Batch information:
")
print(table(metadata$Batch))
cat("
Corrected-data range:
  - Min: ", round(min(expr_combat), 3), "
  - Max: ", round(max(expr_combat), 3), "
  - Mean: ", round(mean(expr_combat), 3), "

Next step: source('scripts/05_DEG_analysis.R')
==========================================================
", sep = "")
