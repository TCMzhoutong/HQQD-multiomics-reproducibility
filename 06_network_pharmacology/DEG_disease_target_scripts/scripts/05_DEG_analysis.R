# ============================================================================
# ============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma")
}

library(limma)

script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0 || script_dir == "") {
  project_root <- getwd()
  if (basename(project_root) == "scripts") {
    project_root <- dirname(project_root)
  }
} else {
  project_root <- dirname(script_dir)
}

batch_corrected_dir <- file.path(project_root, "03_batch_corrected")
deg_results_dir <- file.path(project_root, "04_DEG_results")

if (!dir.exists(deg_results_dir)) {
  dir.create(deg_results_dir, recursive = TRUE)
}

cat("=== Load data ===\n")

load(file.path(batch_corrected_dir, "data_after_combat.RData"))

cat("Expression-matrix dimensions:", dim(expr_after_combat), "\n")
cat("- Genes:", nrow(expr_after_combat), "\n")
cat("- Samples:", ncol(expr_after_combat), "\n")

cat("\nSample groups:\n")
print(table(metadata$Group))

cat("\n=== Create design matrix ===\n")

group <- factor(metadata$Group, levels = c("Control", "Bacterial_Pneumonia"))

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

cat("Design matrix:\n")
print(head(design))
cat("\nSamples per group:\n")
print(colSums(design))

cat("\n=== Fit linear model ===\n")

fit <- lmFit(expr_after_combat, design)

cat("\n=== Create contrast matrix ===\n")

contrast_matrix <- makeContrasts(
  Bacterial_Pneumonia_vs_Control = Bacterial_Pneumonia - Control,
  levels = design
)

cat("Contrast matrix:\n")
print(contrast_matrix)

cat("\n=== Run contrast analysis ===\n")

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

cat("\n=== Obtain differential-expression results ===\n")

deg_results <- topTable(
  fit2, 
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

cat("Result-column definitions:\n")
cat("- logFC: log2 fold change (Bacterial_Pneumonia vs Control)\n")
cat("- AveExpr: average expression\n")
cat("- t: t statistic\n")
cat("- P.Value: raw p value\n")
cat("- adj.P.Val: adjusted p value (BH method)\n")
cat("- B: B statistic (log-odds)\n")

cat("\nDifferential-expression result preview:\n")
print(head(deg_results))

cat("\n=== Filter DEGs ===\n")

logFC_threshold <- 0.8
pval_threshold <- 0.05

deg_results$Significant <- ifelse(
  abs(deg_results$logFC) > logFC_threshold & deg_results$adj.P.Val < pval_threshold,
  ifelse(deg_results$logFC > 0, "Up", "Down"),
  "Not Sig"
)

cat("\nDEG selection criteria:\n")
cat("- |log2FC| >", logFC_threshold, "\n")
cat("- adj.P.Val <", pval_threshold, "\n")

cat("\nDEG summary:\n")
print(table(deg_results$Significant))

deg_up <- deg_results[deg_results$Significant == "Up", ]
deg_down <- deg_results[deg_results$Significant == "Down", ]

cat("\nUpregulated genes:", nrow(deg_up), "\n")
cat("Downregulated genes:", nrow(deg_down), "\n")
cat("Total DEGs:", nrow(deg_up) + nrow(deg_down), "\n")

cat("\n=== Top 10 upregulated genes ===\n")
print(head(deg_up[order(deg_up$logFC, decreasing = TRUE), 
                   c("logFC", "adj.P.Val")], 10))

cat("\n=== Top 10 downregulated genes ===\n")
print(head(deg_down[order(deg_down$logFC), 
                     c("logFC", "adj.P.Val")], 10))

cat("\n=== Save results ===\n")

write.csv(deg_results, 
          file.path(deg_results_dir, "all_genes_results.csv"))

deg_significant <- deg_results[deg_results$Significant != "Not Sig", ]
write.csv(deg_significant, 
          file.path(deg_results_dir, "significant_DEGs.csv"))

save(deg_results, deg_significant, deg_up, deg_down,
     logFC_threshold, pval_threshold,
     file = file.path(deg_results_dir, "DEG_results.RData"))

save(expr_after_combat, metadata, deg_results,
     file = file.path(deg_results_dir, "data_for_visualization.RData"))

cat("Results saved to 04_DEG_results.\n")

cat("\n
==========================================================
            Differential-expression summary
==========================================================

Analysis method: limma
Contrast: Bacterial_Pneumonia vs Control

Selection criteria:
  - |log2FC| > ", logFC_threshold, "
  - adj.P.Val < ", pval_threshold, "

Results:
  - Total genes: ", nrow(deg_results), "
  - upregulated genes: ", nrow(deg_up), "
  - downregulated genes: ", nrow(deg_down), "
  - Total DEGs: ", nrow(deg_up) + nrow(deg_down), "

Top 5 upregulated genes:
", sep = "")
print(head(deg_up[order(deg_up$logFC, decreasing = TRUE), "logFC", drop = FALSE], 5))
cat("
Top 5 downregulated genes:
")
print(head(deg_down[order(deg_down$logFC), "logFC", drop = FALSE], 5))
cat("
Next step: source('scripts/06_PCA_plots.R')
==========================================================
")
