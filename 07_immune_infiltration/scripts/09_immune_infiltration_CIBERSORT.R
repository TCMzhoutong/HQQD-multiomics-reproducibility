# ============================================================================
# ============================================================================

required_pkgs <- c("ggplot2", "dplyr", "tidyr", "scales", "CIBERSORT")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
      ". Install them before running this script."
    )
  )
}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(scales))

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

batch_corrected_dir <- file.path(project_root, "03_batch_corrected")
immune_dir <- file.path(project_root, "06_immune_infiltration")

if (!dir.exists(immune_dir)) {
  dir.create(immune_dir, recursive = TRUE)
}

cat("=== Reading input data ===\n")

expr_file <- file.path(batch_corrected_dir, "expression_matrix_combat.csv")
meta_file <- file.path(batch_corrected_dir, "sample_metadata.csv")

if (!file.exists(expr_file)) {
  stop("Expression matrix not found: 03_batch_corrected/expression_matrix_combat.csv. Run scripts/04_batch_correction.R first.")
}

if (!file.exists(meta_file)) {
  stop("Sample metadata not found: 03_batch_corrected/sample_metadata.csv. Run scripts/04_batch_correction.R first.")
}

expr_mat <- read.csv(expr_file, row.names = 1, check.names = FALSE)
expr_mat <- as.matrix(expr_mat)
metadata <- read.csv(meta_file, stringsAsFactors = FALSE)

required_cols <- c("SampleID", "Batch", "Group")
if (!all(required_cols %in% colnames(metadata))) {
  stop("sample_metadata.csv is missing required columns: SampleID, Batch, Group")
}

common_samples <- intersect(colnames(expr_mat), metadata$SampleID)
if (length(common_samples) < 2) {
  stop("Too few shared samples between the expression matrix and metadata for CIBERSORT analysis")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]
metadata <- metadata[match(common_samples, metadata$SampleID), ]

cat("Expression-matrix dimensions:", dim(expr_mat), "\n")
cat("Sample groups:\n")
print(table(metadata$Group))

cat("\n=== Running CIBERSORT (LM22) ===\n")

sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
if (sig_matrix == "") {
  stop("The LM22 signature matrix was not found. Check the CIBERSORT installation.")
}

cibersort_res <- CIBERSORT::cibersort(
  sig_matrix = sig_matrix,
  mixture_file = expr_mat,
  perm = 100,
  QN = TRUE,
  maxSize = 500
)

cibersort_df <- as.data.frame(cibersort_res)
cibersort_df$SampleID <- rownames(cibersort_df)

cell_cols <- setdiff(colnames(cibersort_df), c("SampleID", "P-value", "Correlation", "RMSE"))

cibersort_with_meta <- metadata %>%
  inner_join(cibersort_df, by = "SampleID")

cat("CIBERSORT result dimensions:", dim(cibersort_with_meta), "\n")

cat("\n=== Saving result files ===\n")

write.csv(
  cibersort_with_meta,
  file.path(immune_dir, "cibersort_LM22_results_with_metadata.csv"),
  row.names = FALSE
)

immune_long <- cibersort_with_meta %>%
  select(SampleID, Group, Batch, all_of(cell_cols)) %>%
  pivot_longer(
    cols = all_of(cell_cols),
    names_to = "Cell_Type",
    values_to = "Proportion"
  )

write.csv(
  immune_long,
  file.path(immune_dir, "cibersort_LM22_results_long_format.csv"),
  row.names = FALSE
)

cat("\n=== Calculating between-group statistics ===\n")

group_levels <- unique(cibersort_with_meta$Group)
if (length(group_levels) != 2) {
  warning("Group does not contain exactly two levels; the Wilcoxon test will use the first two levels")
}

g1 <- group_levels[1]
g2 <- group_levels[2]

stats_df <- lapply(cell_cols, function(cell) {
  x <- cibersort_with_meta[cibersort_with_meta$Group == g1, cell, drop = TRUE]
  y <- cibersort_with_meta[cibersort_with_meta$Group == g2, cell, drop = TRUE]

  if (length(x) < 1 || length(y) < 1) {
    p_val <- NA_real_
  } else {
    p_val <- suppressWarnings(wilcox.test(x, y)$p.value)
  }

  data.frame(
    Cell_Type = cell,
    Group1 = g1,
    Group2 = g2,
    Mean_Group1 = mean(x, na.rm = TRUE),
    Mean_Group2 = mean(y, na.rm = TRUE),
    Median_Group1 = median(x, na.rm = TRUE),
    Median_Group2 = median(y, na.rm = TRUE),
    P_value = p_val,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

stats_df$P_adj_BH <- p.adjust(stats_df$P_value, method = "BH")
stats_df <- stats_df %>% arrange(P_adj_BH)

write.csv(
  stats_df,
  file.path(immune_dir, "cibersort_group_comparison_statistics.csv"),
  row.names = FALSE
)

cat("\n=== Plotting immune-infiltration results ===\n")

sample_order <- metadata %>%
  arrange(Group, Batch, SampleID) %>%
  pull(SampleID)

immune_long$SampleID <- factor(immune_long$SampleID, levels = sample_order)

lm22_palette <- c(
  "B cells naive"                = "#E41A1C",  # Set1 red
  "B cells memory"               = "#377EB8",  # Set1 blue
  "Plasma cells"                 = "#4DAF4A",  # Set1 green
  "T cells CD8"                  = "#984EA3",  # Set1 purple
  "T cells CD4 naive"            = "#FF7F00",  # Set1 orange
  "T cells CD4 memory resting"   = "#A65628",  # Set1 brown
  "T cells CD4 memory activated" = "#F781BF",  # Set1 pink
  "T cells follicular helper"    = "#999999",  # Set1 grey
  "T cells regulatory (Tregs)"   = "#1B9E77",  # Dark2 teal
  "T cells gamma delta"          = "#D95F02",  # Dark2 dark orange
  "NK cells resting"             = "#7570B3",  # Dark2 blue-purple
  "NK cells activated"           = "#E7298A",  # Dark2 magenta
  "Monocytes"                    = "#66A61E",  # Dark2 yellow-green
  "Macrophages M0"               = "#E6AB02",  # Dark2 gold
  "Macrophages M1"               = "#A6761D",  # Dark2 tan-brown
  "Macrophages M2"               = "#4E79A7",  # Tableau blue
  "Dendritic cells resting"      = "#F28E2B",  # Tableau orange
  "Dendritic cells activated"    = "#E15759",  # Tableau red
  "Mast cells resting"           = "#76B7B2",  # Tableau teal
  "Mast cells activated"         = "#59A14F",  # Tableau green
  "Eosinophils"                  = "#EDC948",  # Tableau yellow
  "Neutrophils"                  = "#B07AA1"   # Tableau purple
)

p_stacked <- ggplot(immune_long, aes(x = SampleID, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 0.85) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_fill_manual(values = lm22_palette) +
  labs(
    x = "Sample",
    y = "Estimated fraction",
    fill = "Cell type"
  ) +
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid    = element_blank(),
    legend.position = "right"
  )

boxplot_df <- immune_long %>%
  filter(Cell_Type %in% cell_cols)

group_colors <- c(
  "Bacterial_Pneumonia" = "#B2182B",
  "Control"             = "#2166AC"
)

global_y_top <- max(boxplot_df$Proportion, na.rm = TRUE) * 1.10

sig_map <- stats_df %>%
  filter(Cell_Type %in% cell_cols) %>%
  mutate(
    stars = dplyr::case_when(
      P_adj_BH < 0.001 ~ "***",
      P_adj_BH < 0.01  ~ "**",
      P_adj_BH < 0.05  ~ "*",
      TRUE              ~ ""
    ),
    y_pos = global_y_top
  ) %>%
  filter(stars != "")

p_box <- ggplot(boxplot_df, aes(x = Cell_Type, y = Proportion, fill = Group)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha         = 0.85,
    width         = 0.6,
    position      = position_dodge(0.72)
  ) +
  geom_jitter(
    aes(color = Group),
    position     = position_jitterdodge(dodge.width = 0.72, jitter.width = 0.1),
    alpha        = 0.65,
    size         = 1.2,
    show.legend  = FALSE
  ) +
  { if (nrow(sig_map) > 0)
      geom_text(
        data        = sig_map,
        aes(x = Cell_Type, y = y_pos, label = stars),
        inherit.aes = FALSE,
        size        = 4.5,
        vjust       = 0
      )
  } +
  scale_fill_manual(values  = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    x    = NULL,
    y    = "Estimated fraction",
    fill = NULL
  ) +
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.title      = element_blank(),
    legend.key.size   = unit(0.45, "cm"),
    axis.text.x       = element_text(angle = 35, hjust = 1, size = 9),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_blank()
  )

ggsave(
  file.path(immune_dir, "Figure_E_Immune_Composition_StackedBar.pdf"),
  p_stacked,
  width = 10,
  height = 4
)

ggsave(
  file.path(immune_dir, "Figure_E_Immune_Composition_StackedBar.png"),
  p_stacked,
  width = 10,
  height = 4,
  dpi = 300
)

ggsave(
  file.path(immune_dir, "Figure_E_Immune_Composition_StackedBar.svg"),
  p_stacked,
  width = 10,
  height = 4,
  device = "svg"
)

ggsave(
  file.path(immune_dir, "Figure_F_Immune_Cell_GroupComparison_Boxplot.pdf"),
  p_box,
  width = 11,
  height = 5
)

ggsave(
  file.path(immune_dir, "Figure_F_Immune_Cell_GroupComparison_Boxplot.png"),
  p_box,
  width = 11,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(immune_dir, "Figure_F_Immune_Cell_GroupComparison_Boxplot.svg"),
  p_box,
  width = 11,
  height = 5,
  device = "svg"
)

save(
  expr_mat, metadata, cibersort_res, cibersort_with_meta,
  immune_long, stats_df,
  file = file.path(immune_dir, "immune_infiltration_results.RData")
)

cat("\n")
cat("==========================================================\n")
cat("            CIBERSORT immune-infiltration analysis complete\n")
cat("==========================================================\n\n")
cat("Input data:\n")
cat("  - Samples:", ncol(expr_mat), "\n")
cat("  - Genes:", nrow(expr_mat), "\n\n")
cat("Output directory: 06_immune_infiltration/\n")
cat("  - cibersort_LM22_results_with_metadata.csv\n")
cat("  - cibersort_LM22_results_long_format.csv\n")
cat("  - cibersort_group_comparison_statistics.csv\n")
cat("  - Figure_E_Immune_Composition_StackedBar.pdf/png\n")
cat("  - Figure_F_Immune_Cell_GroupComparison_Boxplot.pdf/png\n")
cat("  - immune_infiltration_results.RData\n\n")
cat("This script runs independently of the DEG workflow.\n")
