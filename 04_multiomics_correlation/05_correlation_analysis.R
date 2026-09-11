# ============================================================================
# Step 05: Spearman Correlation Analysis
# Differential taxa vs differential metabolites
#
# Plotting uses Hy4m/linkET source code:
# https://github.com/Hy4m/linkET
# The heatmap palette is kept consistent with step 06.
# ============================================================================

rm(list = ls())

required_packages <- c("dplyr", "ggplot2", "stringr", "tibble", "purrr", "magrittr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)
})

source_linket <- function() {
  if (requireNamespace("linkET", quietly = TRUE)) {
    suppressPackageStartupMessages(library(linkET))
    return("installed package")
  }
  stop("Missing required R package: linkET. Install Hy4m/linkET before running this script.")
}

save_plot_svg <- function(filename, plot, width, height) {
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = filename, width = width, height = height)
  } else {
    grDevices::svg(filename = filename, width = width, height = height, onefile = FALSE)
  }
  print(plot)
  grDevices::dev.off()
}

linket_source <- source_linket()

output_dir <- "05_correlation_analysis"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

short_label <- function(x, max_chars = 36) {
  x <- as.character(x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

cat("\n=== Step 05: Data Import ===\n")
cat("Using linkET from: ", linket_source, "\n", sep = "")

taxa_data <- read.csv(
  "04_differential_species_metabolites/differential_taxa_abundance_Model_vs_HQQD.csv",
  check.names = FALSE, row.names = 1, stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)

metabolite_data <- read.csv(
  "04_differential_species_metabolites/differential_metabolites_abundance_Model_vs_HQQD.csv",
  check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)

cat("Taxa loaded: ", nrow(taxa_data), "\n", sep = "")
cat("Metabolites loaded: ", nrow(metabolite_data), "\n", sep = "")

metabolite_name_col <- if ("Name" %in% colnames(metabolite_data)) "Name" else if ("name" %in% colnames(metabolite_data)) "name" else NA
if (is.na(metabolite_name_col)) stop("Metabolite table must contain Name or name column.")

metabolite_labels <- ifelse(
  is.na(metabolite_data[[metabolite_name_col]]) | metabolite_data[[metabolite_name_col]] == "",
  metabolite_data$ID,
  metabolite_data[[metabolite_name_col]]
)
metabolite_labels <- make.unique(short_label(metabolite_labels, 34))

taxa_matrix <- taxa_data %>%
  select(matches("^(M|Z|K)\\d+")) %>%
  as.matrix()

metabolite_matrix <- metabolite_data %>%
  select(matches("^(M|Z|K)\\d+")) %>%
  as.matrix()
rownames(metabolite_matrix) <- metabolite_labels

common_samples <- intersect(colnames(taxa_matrix), colnames(metabolite_matrix))
common_samples <- c(
  sort(common_samples[grepl("^K\\d+", common_samples)]),
  sort(common_samples[grepl("^M\\d+", common_samples)]),
  sort(common_samples[grepl("^Z\\d+", common_samples)])
)
if (length(common_samples) < 3) stop("Too few common samples for Spearman correlation.")

taxa_t <- as.data.frame(t(taxa_matrix[, common_samples, drop = FALSE]))
metabolite_t <- as.data.frame(t(metabolite_matrix[, common_samples, drop = FALSE]))

cat("Common samples: ", length(common_samples), "\n", sep = "")
cat("Final dimensions: ", ncol(taxa_t), " taxa x ", ncol(metabolite_t), " metabolites\n", sep = "")

cat("\n=== Spearman Correlation via linkET ===\n")

cor_obj <- correlate(taxa_t, metabolite_t, method = "spearman", use = "pairwise.complete.obs")
cor_matrix <- cor_obj$r
pval_matrix <- cor_obj$p

write.csv(cor_matrix, file.path(output_dir, "spearman_correlation_matrix.csv"), fileEncoding = "UTF-8")
write.csv(pval_matrix, file.path(output_dir, "spearman_pvalue_matrix.csv"), fileEncoding = "UTF-8")

sig_count <- sum(pval_matrix < 0.05, na.rm = TRUE)
cat("Significant correlations (p < 0.05): ", sig_count, "\n", sep = "")

cat("\n=== Plotting with linkET ===\n")

cor_plot <- cor_matrix
p_plot <- pval_matrix
if (nrow(cor_plot) > 1) {
  row_ord <- hclust(dist(cor_plot), method = "complete")$order
  cor_plot <- cor_plot[row_ord, , drop = FALSE]
  p_plot <- p_plot[row_ord, , drop = FALSE]
}
if (ncol(cor_plot) > 1) {
  col_ord <- hclust(dist(t(cor_plot)), method = "complete")$order
  cor_plot <- cor_plot[, col_ord, drop = FALSE]
  p_plot <- p_plot[, col_ord, drop = FALSE]
}

sig_df <- expand.grid(
  .rownames = rownames(p_plot),
  .colnames = colnames(p_plot),
  stringsAsFactors = FALSE
)
sig_df$p <- as.vector(p_plot)
sig_df$star <- ifelse(is.na(sig_df$p), "", ifelse(sig_df$p < 0.001, "***", ifelse(sig_df$p < 0.01, "**", ifelse(sig_df$p < 0.05, "*", ""))))
sig_df <- sig_df[sig_df$star != "", , drop = FALSE]

color_palette <- colorRampPalette(c(
  "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
  "#FFFFFF",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
))(100)

heatmap_plot <- qcorrplot(as_correlate(cor_plot, p = p_plot, is_corr = TRUE), is_corr = TRUE) +
  geom_tile(colour = "grey90", linewidth = 0.2, width = 1, height = 1) +
  geom_text(data = sig_df, aes(x = .colnames, y = .rownames, label = star), inherit.aes = FALSE, size = 3.2, color = "black", family = "sans") +
  scale_fill_gradientn(colours = color_palette, limits = c(-1, 1), name = "Spearman rho", guide = guide_colorbar(display = "rectangles")) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = "sans", base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  )

plot_width <- max(8.5, min(13, 6.5 + 0.18 * ncol(cor_plot)))
ggsave(file.path(output_dir, "spearman_correlation_heatmap.png"), heatmap_plot, width = plot_width, height = 6.5, dpi = 300)
ggsave(file.path(output_dir, "spearman_correlation_heatmap.pdf"), heatmap_plot, width = plot_width, height = 6.5)
save_plot_svg(file.path(output_dir, "spearman_correlation_heatmap.svg"), heatmap_plot, width = plot_width, height = 6.5)

cat("Step 05 complete. Outputs saved to: ", output_dir, "\n", sep = "")
