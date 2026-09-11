# ============================================================================
# ============================================================================

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("cowplot")) install.packages("cowplot")

library(ggplot2)
library(ggrepel)
library(cowplot)

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
figures_dir <- file.path(project_root, "05_figures")

if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

cat("=== Load data ===\n")

load(file.path(batch_corrected_dir, "data_before_combat.RData"))
load(file.path(batch_corrected_dir, "data_after_combat.RData"))

perform_pca <- function(expr_matrix) {
  expr_t <- t(expr_matrix)
  
  pca_result <- prcomp(expr_t, center = TRUE, scale. = TRUE)
  
  pca_data <- as.data.frame(pca_result$x[, 1:2])
  colnames(pca_data) <- c("PC1", "PC2")
  
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  
  return(list(
    data = pca_data,
    var_explained = var_explained
  ))
}

cat("\n=== Run PCA ===\n")

pca_before <- perform_pca(expr_before_combat)
pca_before_data <- cbind(pca_before$data, metadata)

cat("\nPCA variance explained before batch correction:\n")
cat("PC1:", round(pca_before$var_explained[1], 1), "%\n")
cat("PC2:", round(pca_before$var_explained[2], 1), "%\n")

pca_after <- perform_pca(expr_after_combat)
pca_after_data <- cbind(pca_after$data, metadata)

cat("\nPCA variance explained after batch correction:\n")
cat("PC1:", round(pca_after$var_explained[1], 1), "%\n")
cat("PC2:", round(pca_after$var_explained[2], 1), "%\n")

cat("\n=== Plot PCA ===\n")

theme_publication <- theme_bw(base_family = "sans") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

batch_colors <- c("GSE40012" = "#E64B35", "GSE20346" = "#4DBBD5")
group_colors <- c("Control" = "#00A087", "Bacterial_Pneumonia" = "#E64B35")

p_before_batch <- ggplot(pca_before_data, aes(x = PC1, y = PC2, color = Batch, shape = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = batch_colors) +
  labs(
    title = "Before Batch Correction",
    x = paste0("PC1 (", round(pca_before$var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_before$var_explained[2], 1), "%)")
  ) +
  theme_publication +
  stat_ellipse(aes(group = Batch), type = "norm", linetype = 2, alpha = 0.5)

p_after_batch <- ggplot(pca_after_data, aes(x = PC1, y = PC2, color = Batch, shape = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = batch_colors) +
  labs(
    title = "After Batch Correction",
    x = paste0("PC1 (", round(pca_after$var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_after$var_explained[2], 1), "%)")
  ) +
  theme_publication +
  stat_ellipse(aes(group = Batch), type = "norm", linetype = 2, alpha = 0.5)

pca_combined <- plot_grid(
  p_before_batch, 
  p_after_batch,
  ncol = 2,
  align = "h"
)

ggsave(file.path(figures_dir, "Figure_AB_PCA_batch_correction.pdf"), 
       pca_combined, 
       width = 14, height = 6)

ggsave(file.path(figures_dir, "Figure_AB_PCA_batch_correction.png"), 
       pca_combined, 
       width = 14, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "Figure_AB_PCA_batch_correction.svg"), 
       pca_combined, 
       width = 14, height = 6)

cat("PCA panels A-B saved to 05_figures.\n")

p_before_group <- ggplot(pca_before_data, aes(x = PC1, y = PC2, color = Group, shape = Batch)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Before Correction (by Disease Status)",
    x = paste0("PC1 (", round(pca_before$var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_before$var_explained[2], 1), "%)")
  ) +
  theme_publication +
  stat_ellipse(aes(group = Group), type = "norm", linetype = 2, alpha = 0.5)

p_after_group <- ggplot(pca_after_data, aes(x = PC1, y = PC2, color = Group, shape = Batch)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "After Correction (by Disease Status)",
    x = paste0("PC1 (", round(pca_after$var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_after$var_explained[2], 1), "%)")
  ) +
  theme_publication +
  stat_ellipse(aes(group = Group), type = "norm", linetype = 2, alpha = 0.5)

pca_group_combined <- plot_grid(
  p_before_group, 
  p_after_group,
  ncol = 2,
  align = "h"
)

ggsave(file.path(figures_dir, "Supplementary_PCA_by_group.pdf"), 
       pca_group_combined, 
       width = 14, height = 6)

ggsave(file.path(figures_dir, "Supplementary_PCA_by_group.png"), 
       pca_group_combined, 
       width = 14, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "Supplementary_PCA_by_group.svg"), 
       pca_group_combined, 
       width = 14, height = 6)

cat("Supplementary PCA by disease state saved.\n")

ggsave(file.path(figures_dir, "Figure_A_PCA_before_correction.pdf"), 
       p_before_batch, width = 7, height = 6)
ggsave(file.path(figures_dir, "Figure_B_PCA_after_correction.pdf"), 
       p_after_batch, width = 7, height = 6)

ggsave(file.path(figures_dir, "Figure_A_PCA_before_correction.png"), 
       p_before_batch, width = 7, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "Figure_B_PCA_after_correction.png"), 
       p_after_batch, width = 7, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "Figure_A_PCA_before_correction.svg"), 
       p_before_batch, width = 7, height = 6)
ggsave(file.path(figures_dir, "Figure_B_PCA_after_correction.svg"), 
       p_after_batch, width = 7, height = 6)

cat("\n=== PCA visualisation complete ===\n")
cat("
Generated files (saved to 05_figures):
  - Figure_AB_PCA_batch_correction.pdf/png: combined PCA panels (A-B)
  - Figure_A_PCA_before_correction.pdf/png: PCA before correction
  - Figure_B_PCA_after_correction.pdf/png: PCA after correction
  - Supplementary_PCA_by_group.pdf/png: PCA coloured by disease state

Next step: source('scripts/07_volcano_plot.R')
")
