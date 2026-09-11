# ============================================================================
# ============================================================================

# Clear environment
rm(list = ls())

old_lifecycle_verbosity <- getOption("lifecycle_verbosity")
options(lifecycle_verbosity = "quiet")
on.exit(options(lifecycle_verbosity = old_lifecycle_verbosity), add = TRUE)

library(dplyr)
library(tibble)
library(RColorBrewer)
library(ggplot2)

load_linket <- function() {
  if (requireNamespace("linkET", quietly = TRUE)) {
    suppressPackageStartupMessages(library(linkET))
    return("installed linkET package")
  }
  candidates <- c(
    file.path("..", "04_correlation_analysis", "linkET-master", "R"),
    file.path("04_correlation_analysis", "linkET-master", "R")
  )
  linket_dir <- candidates[dir.exists(candidates)][1]
  if (is.na(linket_dir)) {
    stop("Hy4m/linkET source directory not found.")
  }
  files <- list.files(linket_dir, pattern = "\\.R$", full.names = TRUE)
  for (f in files) source(f, encoding = "UTF-8")
  "local Hy4m/linkET source"
}
linket_source <- load_linket()

cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(cmd_file) && nzchar(cmd_file)) dirname(normalizePath(cmd_file)) else ""
if (length(script_dir) == 0 || script_dir == "") {
  project_root <- getwd()
  if (basename(project_root) == "scripts") {
    project_root <- dirname(project_root)
  }
} else {
  project_root <- dirname(script_dir)
}

raw_dir <- file.path(project_root, "01_raw_data")
batch_corrected_dir <- file.path(project_root, "03_batch_corrected")
immune_dir <- file.path(project_root, "06_immune_infiltration")
output_dir <- file.path(project_root, "06_immune_infiltration", "hubgene_immune_correlation")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(paste0("Created output directory: ", output_dir, "\n"))
}

cat("\n=== Data Import ===\n")

hubgene_candidates <- c(
  file.path(dirname(project_root), "06_network_pharmacology", "06.hub_genes", "03_hub_genes_intersection_all.csv"),
  file.path(raw_dir, "03_hubgenes_disease_ingredient_metabolite.csv")
)
hubgene_file <- hubgene_candidates[file.exists(hubgene_candidates)][1]
expr_file <- file.path(batch_corrected_dir, "expression_matrix_combat.csv")
meta_file <- file.path(batch_corrected_dir, "sample_metadata.csv")
cibersort_wide_file <- file.path(immune_dir, "cibersort_LM22_results_with_metadata.csv")

if (is.na(hubgene_file) || !file.exists(hubgene_file)) {
  stop("Hub-gene file not found: ../06_network_pharmacology/06.hub_genes/03_hub_genes_intersection_all.csv or 01_raw_data/03_hubgenes_disease_ingredient_metabolite.csv")
}
if (!file.exists(expr_file)) {
  stop("Expression matrix not found: 03_batch_corrected/expression_matrix_combat.csv")
}
if (!file.exists(meta_file)) {
  stop("Sample metadata not found: 03_batch_corrected/sample_metadata.csv")
}
if (!file.exists(cibersort_wide_file)) {
  stop("CIBERSORT wide-format result not found: 06_immune_infiltration/cibersort_LM22_results_with_metadata.csv. Run scripts/09_immune_infiltration_CIBERSORT.R first.")
}

hubgene_tbl <- read.csv(hubgene_file, stringsAsFactors = FALSE)
if (!("name" %in% colnames(hubgene_tbl))) {
  stop("The hub-gene file is missing the 'name' column")
}

hubgenes_raw <- hubgene_tbl$name
hubgenes_raw <- hubgenes_raw[!is.na(hubgenes_raw) & nzchar(hubgenes_raw)]
hubgenes <- unique(trimws(hubgenes_raw))

expr_mat <- read.csv(expr_file, row.names = 1, check.names = FALSE)
expr_mat <- as.matrix(expr_mat)
metadata <- read.csv(meta_file, stringsAsFactors = FALSE)

required_cols <- c("SampleID", "Batch", "Group")
if (!all(required_cols %in% colnames(metadata))) {
  stop("sample_metadata.csv is missing required columns: SampleID, Batch, Group")
}

common_samples <- intersect(colnames(expr_mat), metadata$SampleID)
if (length(common_samples) < 2) {
  stop("Too few shared samples between the expression matrix and metadata")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]
metadata <- metadata[match(common_samples, metadata$SampleID), ]

cat(paste0("Hub genes requested: ", length(hubgenes), "\n"))
cat(paste0("Hub gene source: ", normalizePath(hubgene_file, winslash = "/", mustWork = TRUE), "\n"))
cat(paste0("Expression matrix: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples\n"))

expr_genes <- rownames(expr_mat)
hubgenes_upper <- toupper(hubgenes)
expr_genes_upper <- toupper(expr_genes)

present_idx <- match(hubgenes_upper, expr_genes_upper)
present_flag <- !is.na(present_idx)

present_hubgenes <- hubgenes[present_flag]
missing_hubgenes <- hubgenes[!present_flag]
present_expr_genes <- expr_genes[present_idx[present_flag]]

cat("\n=== Hub Gene Matching ===\n")
cat(paste0("Genes found in GSE40012 + GSE20346 merged matrix: ", length(present_hubgenes), "\n"))
cat(paste0("Genes NOT found: ", length(missing_hubgenes), "\n"))

if (length(missing_hubgenes) > 0) {
  cat("Missing genes:\n")
  print(missing_hubgenes)
}

write.csv(
  data.frame(name = hubgenes, in_dataset = present_flag, stringsAsFactors = FALSE),
  file.path(output_dir, "hubgenes_in_dataset_status.csv"),
  row.names = FALSE
)

if (length(present_expr_genes) < 1) {
  stop("None of the hub genes are present in the expression matrix")
}

cat("\n=== Loading Existing CIBERSORT LM22 Wide Table ===\n")

cibersort_df <- read.csv(cibersort_wide_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!("SampleID" %in% colnames(cibersort_df))) {
  stop("cibersort_LM22_results_with_metadata.csv is missing the SampleID column")
}

lm22_cells <- c(
  "B cells naive",
  "B cells memory",
  "Plasma cells",
  "T cells CD8",
  "T cells CD4 naive",
  "T cells CD4 memory resting",
  "T cells CD4 memory activated",
  "T cells follicular helper",
  "T cells regulatory (Tregs)",
  "T cells gamma delta",
  "NK cells resting",
  "NK cells activated",
  "Monocytes",
  "Macrophages M0",
  "Macrophages M1",
  "Macrophages M2",
  "Dendritic cells resting",
  "Dendritic cells activated",
  "Mast cells resting",
  "Mast cells activated",
  "Eosinophils",
  "Neutrophils"
)

cell_cols <- intersect(lm22_cells, colnames(cibersort_df))
missing_cells <- setdiff(lm22_cells, cell_cols)

if (length(cell_cols) != 22) {
  warning(paste0("The number of LM22 cell types is not 22; observed: ", length(cell_cols)))
  if (length(missing_cells) > 0) {
    cat("Missing LM22 cell columns:\n")
    print(missing_cells)
  }
}

if (length(cell_cols) < 1) {
  stop("No LM22 cell columns were found in cibersort_LM22_results_with_metadata.csv")
}

ciber_samples <- cibersort_df$SampleID
aligned_samples <- intersect(common_samples, ciber_samples)

if (length(aligned_samples) < 3) {
  stop("Too few samples remain after aligning the expression matrix with the CIBERSORT table")
}

immune_matrix <- cibersort_df[match(aligned_samples, cibersort_df$SampleID), cell_cols, drop = FALSE]
immune_matrix <- as.data.frame(immune_matrix, stringsAsFactors = FALSE)
immune_matrix[] <- lapply(immune_matrix, function(x) suppressWarnings(as.numeric(x)))
immune_matrix <- as.matrix(immune_matrix)

if (any(is.na(immune_matrix))) {
  warning("The immune-cell matrix contains NA values; corresponding missing values will be ignored in correlation analysis")
}

gene_matrix <- t(expr_mat[present_expr_genes, aligned_samples, drop = FALSE])
gene_matrix <- apply(gene_matrix, 2, as.numeric)
rownames(gene_matrix) <- aligned_samples

cat(paste0("Aligned samples: ", length(aligned_samples), "\n"))
cat(paste0("Immune matrix: ", nrow(immune_matrix), " samples x ", ncol(immune_matrix), " cells\n"))
cat(paste0("Gene matrix: ", nrow(gene_matrix), " samples x ", ncol(gene_matrix), " hub genes\n"))

cat("\n=== Spearman Correlation via linkET ===\n")

immune_df <- as.data.frame(immune_matrix, check.names = FALSE)
gene_df_for_cor <- as.data.frame(gene_matrix, check.names = FALSE)

cor_obj <- correlate(immune_df, gene_df_for_cor, method = "spearman", use = "pairwise.complete.obs")
cor_matrix <- cor_obj$r
pval_matrix <- cor_obj$p

write.csv(cor_matrix, file.path(output_dir, "spearman_correlation_matrix_immune_vs_hubgenes.csv"))
write.csv(pval_matrix, file.path(output_dir, "spearman_pvalue_matrix_immune_vs_hubgenes.csv"))

cat("Correlation analysis complete.\n")

n_sig_05 <- sum(pval_matrix < 0.05, na.rm = TRUE)
n_sig_01 <- sum(pval_matrix < 0.01, na.rm = TRUE)
n_sig_001 <- sum(pval_matrix < 0.001, na.rm = TRUE)

cat(paste0("Significant correlations (p < 0.05): ", n_sig_05, "\n"))
cat(paste0("Significant correlations (p < 0.01): ", n_sig_01, "\n"))
cat(paste0("Significant correlations (p < 0.001): ", n_sig_001, "\n"))

color_palette <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                                    "#FFFFFF",
                                    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100)

build_and_save_linket_heatmap <- function(cor_mat, p_mat, pdf_file, png_file, svg_file = NULL,
                                          width = 11, height = 6.5) {
  cor_plot <- cor_mat
  p_plot <- p_mat
  cor_for_cluster <- cor_plot
  cor_for_cluster[!is.finite(cor_for_cluster)] <- 0

  if (nrow(cor_plot) > 1) {
    row_ord <- hclust(dist(cor_for_cluster), method = "complete")$order
    cor_plot <- cor_plot[row_ord, , drop = FALSE]
    p_plot <- p_plot[row_ord, , drop = FALSE]
    cor_for_cluster <- cor_for_cluster[row_ord, , drop = FALSE]
  }
  if (ncol(cor_plot) > 1) {
    col_ord <- hclust(dist(t(cor_for_cluster)), method = "complete")$order
    cor_plot <- cor_plot[, col_ord, drop = FALSE]
    p_plot <- p_plot[, col_ord, drop = FALSE]
  }

  sig_df <- expand.grid(
    .rownames = rownames(p_plot),
    .colnames = colnames(p_plot),
    stringsAsFactors = FALSE
  )
  sig_df$p <- as.vector(p_plot)
  sig_df$star <- ifelse(
    is.na(sig_df$p),
    "",
    ifelse(sig_df$p < 0.001, "***", ifelse(sig_df$p < 0.01, "**", ifelse(sig_df$p < 0.05, "*", "")))
  )
  sig_df <- sig_df[sig_df$star != "", , drop = FALSE]

  p <- qcorrplot(as_correlate(cor_plot, p = p_plot, is_corr = TRUE), is_corr = TRUE) +
    geom_tile(colour = "grey90", linewidth = 0.2, width = 1, height = 1) +
    geom_text(
      data = sig_df,
      aes(x = .colnames, y = .rownames, label = star),
      inherit.aes = FALSE,
      size = 3.2,
      color = "black",
      family = "sans"
    ) +
    scale_fill_gradientn(
      colours = color_palette,
      limits = c(-1, 1),
      name = "Spearman rho",
      guide = guide_colorbar(display = "rectangles")
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_family = "sans", base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10)
    )

  ggsave(pdf_file, plot = p, width = width, height = height,
         units = "in")
  ggsave(png_file, plot = p, width = width, height = height,
         units = "in", dpi = 300)
  if (!is.null(svg_file)) {
    ggsave(svg_file, plot = p, width = width, height = height,
           units = "in", device = "svg")
  }

  return(list(width = width, height = height))
}

cat("\n=== Generating Correlation Heatmap ===\n")

plot_size <- build_and_save_linket_heatmap(
  cor_matrix,
  pval_matrix,
  file.path(output_dir, "immune_cell_vs_hubgene_correlation_heatmap.pdf"),
  file.path(output_dir, "immune_cell_vs_hubgene_correlation_heatmap.png"),
  file.path(output_dir, "immune_cell_vs_hubgene_correlation_heatmap.svg"),
  width = 11,
  height = 6.5
)

cat(paste0("Heatmap dimensions: ", round(plot_size$width, 2), " x ",
           round(plot_size$height, 2), " inches\n"))

cat("\n=== Generating Filtered Heatmap (p < 0.05) ===\n")

sig_rows <- apply(pval_matrix, 1, function(x) any(x < 0.05, na.rm = TRUE))
sig_cols <- apply(pval_matrix, 2, function(x) any(x < 0.05, na.rm = TRUE))

if (sum(sig_rows) > 0 && sum(sig_cols) > 0) {
  cor_matrix_filtered <- cor_matrix[sig_rows, sig_cols, drop = FALSE]
  pval_matrix_filtered <- pval_matrix[sig_rows, sig_cols, drop = FALSE]

  cat(paste0("Filtered matrix: ", nrow(cor_matrix_filtered), " immune cells x ",
             ncol(cor_matrix_filtered), " hub genes\n"))

  build_and_save_linket_heatmap(
    cor_matrix_filtered,
    pval_matrix_filtered,
    file.path(output_dir, "immune_cell_vs_hubgene_correlation_heatmap_significant.pdf"),
    file.path(output_dir, "immune_cell_vs_hubgene_correlation_heatmap_significant.png"),
    file.path(output_dir, "immune_cell_vs_hubgene_correlation_heatmap_significant.svg"),
    width = 11,
    height = 6.5
  )

  cat("Saved filtered correlation heatmap.\n")
} else {
  cat("No significant correlations found (p < 0.05).\n")
}

cor_long <- as.data.frame(as.table(cor_matrix), stringsAsFactors = FALSE)
colnames(cor_long) <- c("Immune_Cell", "HubGene", "SpearmanR")

pval_long <- as.data.frame(as.table(pval_matrix), stringsAsFactors = FALSE)
colnames(pval_long) <- c("Immune_Cell", "HubGene", "P_value")

cor_stats <- cor_long %>%
  left_join(pval_long, by = c("Immune_Cell", "HubGene")) %>%
  arrange(P_value)

write.csv(
  cor_stats,
  file.path(output_dir, "immune_cell_vs_hubgene_correlation_statistics.csv"),
  row.names = FALSE
)

summary_stats <- data.frame(
  Metric = c("Hub genes requested",
             "Hub genes found",
             "Hub genes missing",
             "Immune cell types",
             "Total comparisons",
             "Significant (p < 0.05)",
             "Significant (p < 0.01)",
             "Significant (p < 0.001)",
             "Positive correlations",
             "Negative correlations",
             "Mean correlation (all)",
             "Mean correlation (significant)",
             "Max positive correlation",
             "Max negative correlation"),
  Value = c(
    length(hubgenes),
    length(present_expr_genes),
    length(missing_hubgenes),
    nrow(cor_matrix),
    length(cor_matrix),
    n_sig_05,
    n_sig_01,
    n_sig_001,
    sum(cor_matrix > 0, na.rm = TRUE),
    sum(cor_matrix < 0, na.rm = TRUE),
    round(mean(cor_matrix, na.rm = TRUE), 3),
    round(mean(cor_matrix[pval_matrix < 0.05], na.rm = TRUE), 3),
    round(max(cor_matrix, na.rm = TRUE), 3),
    round(min(cor_matrix, na.rm = TRUE), 3)
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_stats, file.path(output_dir, "immune_hubgene_correlation_summary.csv"), row.names = FALSE)

cat("\n=== Analysis Complete ===\n")
cat(paste0("All results saved to: ", output_dir, "\n"))
cat("\nOutput files:\n")
cat("  - hubgenes_in_dataset_status.csv\n")
cat("  - spearman_correlation_matrix_immune_vs_hubgenes.csv\n")
cat("  - spearman_pvalue_matrix_immune_vs_hubgenes.csv\n")
cat("  - immune_cell_vs_hubgene_correlation_heatmap.pdf/png\n")
cat("  - immune_cell_vs_hubgene_correlation_heatmap_significant.pdf/png\n")
cat("  - immune_cell_vs_hubgene_correlation_statistics.csv\n")
cat("  - immune_hubgene_correlation_summary.csv\n")

if (length(missing_hubgenes) > 0) {
  cat("\nGenes not in dataset (printed above) were excluded from analysis.\n")
}
