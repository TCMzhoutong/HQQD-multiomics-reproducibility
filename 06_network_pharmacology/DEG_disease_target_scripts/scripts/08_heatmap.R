# ============================================================================
# ============================================================================

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("dplyr")) install.packages("dplyr")

library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0 || script_dir == "") {
  project_root <- getwd()
  if (basename(project_root) == "scripts") {
    project_root <- dirname(project_root)
  }
} else {
  project_root <- dirname(script_dir)
}

deg_results_dir <- file.path(project_root, "04_DEG_results")
figures_dir <- file.path(project_root, "05_figures")

cat("=== Load data ===\n")

load(file.path(deg_results_dir, "data_for_visualization.RData"))

cat("\n=== Select top 50 DEGs ===\n")

logFC_threshold <- 0.8
pval_threshold <- 0.05

deg_significant <- deg_results %>%
  filter(abs(logFC) > logFC_threshold & adj.P.Val < pval_threshold)

cat("Significant DEGs:", nrow(deg_significant), "\n")

n_top <- min(50, nrow(deg_significant))

top_degs <- deg_significant %>%
  arrange(desc(abs(logFC))) %>%
  head(n_top)

cat("Selected top DEGs:", nrow(top_degs), "\n")

top_genes <- rownames(top_degs)

cat("\n=== Extract expression matrices ===\n")

expr_top <- expr_after_combat[top_genes, ]

cat("Heatmap expression-matrix dimensions:", dim(expr_top), "\n")

expr_scaled <- t(scale(t(expr_top)))

cat("Range after Z-score standardisation:\n")
cat("Min:", round(min(expr_scaled), 2), "\n")
cat("Max:", round(max(expr_scaled), 2), "\n")

cat("\n=== Prepare annotations ===\n")

annotation_col <- data.frame(
  Group = metadata$Group,
  Batch = metadata$Batch,
  row.names = metadata$SampleID
)

annotation_row <- data.frame(
  Regulation = ifelse(top_degs$logFC > 0, "Up", "Down"),
  row.names = rownames(top_degs)
)

annotation_colors <- list(
  Group = c(Control = "#00A087", Bacterial_Pneumonia = "#E64B35"),
  Batch = c(GSE40012 = "#3C5488", GSE20346 = "#F39B7F"),
  Regulation = c(Up = "#E64B35", Down = "#4DBBD5")
)

cat("\n=== Plot heatmap ===\n")

heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

sample_order <- order(metadata$Group)
expr_scaled_ordered <- expr_scaled[, sample_order]
annotation_col_ordered <- annotation_col[sample_order, ]

p_heatmap <- pheatmap(
  expr_scaled_ordered,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col_ordered,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  fontsize_row = 6,
  fontsize = 10,
  main = "Heatmap of Top 50 DEGs",
  border_color = NA,
  scale = "none",
  cellwidth = NA,
  cellheight = 8
)

pdf(file.path(figures_dir, "Figure_D_Heatmap.pdf"), width = 12, height = 14)
print(p_heatmap)
dev.off()

png(file.path(figures_dir, "Figure_D_Heatmap.png"), 
    width = 12, height = 14, units = "in", res = 300)
print(p_heatmap)
dev.off()

svg(file.path(figures_dir, "Figure_D_Heatmap.svg"), 
    width = 12, height = 14, onefile = FALSE)
print(p_heatmap)
dev.off()

cat("Heatmap saved to 05_figures.\n")

p_heatmap_clustered <- pheatmap(
  expr_scaled,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  fontsize_row = 6,
  fontsize = 10,
  main = "Heatmap of Top 50 DEGs (Clustered)",
  border_color = NA,
  scale = "none",
  cellheight = 8
)

pdf(file.path(figures_dir, "Supplementary_Heatmap_clustered.pdf"), 
    width = 12, height = 14)
print(p_heatmap_clustered)
dev.off()

png(file.path(figures_dir, "Supplementary_Heatmap_clustered.png"), 
    width = 12, height = 14, units = "in", res = 300)
print(p_heatmap_clustered)
dev.off()

svg(file.path(figures_dir, "Supplementary_Heatmap_clustered.svg"), 
    width = 12, height = 14, onefile = FALSE)
print(p_heatmap_clustered)
dev.off()

cat("Clustered heatmap saved.\n")

write.csv(expr_top, 
          file.path(figures_dir, "Top50_DEGs_expression.csv"))
write.csv(top_degs, 
          file.path(figures_dir, "Top50_DEGs_statistics.csv"))

cat("\n=== Heatmap visualisation complete ===\n")
cat("
Generated files (saved to 05_figures):
  - Figure_D_Heatmap.pdf/png: main heatmap ordered by group
  - Supplementary_Heatmap_clustered.pdf/png: clustered heatmap
  - Top50_DEGs_expression.csv: top-50-gene expression matrix
  - Top50_DEGs_statistics.csv: top-50-gene statistics

All visualisations complete.
==========================================================
")
