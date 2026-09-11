# ============================================================================
# LEfSe Analysis for Quorum Sensing (QS) Class
# File: 09_LEfSe_QS_Class.R
# Description: Performs LEfSe analysis on QS classes and generates
#              LDA barplot + CSV results only.
# ============================================================================

rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readxl,
  microeco,
  magrittr,
  dplyr,
  tidyr,
  ggplot2,
  stringr,
  tibble
)

group_colors <- c(
  "Control" = "#3B4992",
  "Model" = "#EE0000",
  "LVX" = "#008B45",
  "HQQD" = "#FF8C00"
)

# ============================================================================
# Configuration Parameters
# ============================================================================

input_file <- "QS.class.relabundance.xlsx"

output_dir <- "09_QS_Analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ALPHA_VALUE <- 0.05
LDA_THRESHOLD <- 2.0

# ============================================================================
# Step 1: Load and Data Preparation
# ============================================================================

cat("\n=== Step 1: Loading and preparing data ===\n")

if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

raw_data <- read_excel(input_file)
cat("Original data dimensions:", dim(raw_data), "\n")

metadata_col <- "Quorum_Sensing_Class"
sample_cols <- setdiff(colnames(raw_data), metadata_col)

cat("Sample columns identified:", length(sample_cols), "\n")
print(sample_cols)

# 1. Feature metadata (for output merge)
feature_metadata <- raw_data %>%
  select(all_of(metadata_col)) %>%
  rename(Quorum_Sensing_Class = all_of(metadata_col))

# 2. OTU table
otu_table <- raw_data %>%
  select(all_of(c(metadata_col, sample_cols))) %>%
  tibble::column_to_rownames(metadata_col) %>%
  as.data.frame()

# Scale relative abundances to pseudo-counts (×10^6 → PPM)
otu_table <- otu_table * 1e6

# 3. Taxonomy table (pseudo-taxonomy for microeco)
tax_table <- raw_data %>%
  select(all_of(metadata_col)) %>%
  mutate(
    Kingdom = "QS_System",
    Phylum = "Unassigned",
    Class = "Unassigned",
    Order = "Unassigned",
    Family = "Unassigned",
    Genus = "Unassigned",
    Species = .[[1]]
  ) %>%
  tibble::column_to_rownames(metadata_col) %>%
  as.data.frame()

# 4. Sample metadata
sample_metadata <- data.frame(
  Sample = sample_cols,
  Group = case_when(
    grepl("^K", sample_cols) ~ "Control",
    grepl("^M", sample_cols) ~ "Model",
    grepl("^Y", sample_cols) ~ "LVX",
    grepl("^Z", sample_cols) ~ "HQQD",
    TRUE ~ "Unknown"
  ),
  stringsAsFactors = FALSE
)
rownames(sample_metadata) <- sample_metadata$Sample
sample_metadata$Group <- factor(sample_metadata$Group, levels = c("Control", "Model", "LVX", "HQQD"))

cat("\nSample Group Distribution:\n")
print(table(sample_metadata$Group))

# ============================================================================
# Step 2: Create Microtable Object
# ============================================================================

cat("\n=== Step 2: Creating microtable object ===\n")

dataset <- microtable$new(
  sample_table = sample_metadata,
  otu_table = otu_table,
  tax_table = tax_table
)

cat("Dataset created. Number of taxa:", nrow(dataset$tax_table), "\n")

# ============================================================================
# Step 3: Run LEfSe Analysis
# ============================================================================

cat("\n=== Step 3: Running LEfSe analysis ===\n")

t1 <- trans_diff$new(
  dataset = dataset,
  method = "lefse",
  group = "Group",
  alpha = ALPHA_VALUE,
  p_adjust_method = "fdr",
  taxa_level = "Species"
)

cat("LDA Threshold:", LDA_THRESHOLD, "\n")
cat("Alpha Value:", ALPHA_VALUE, "\n")

res_lefse <- t1$res_diff
cat("Total differential features found:", nrow(res_lefse), "\n")

# ============================================================================
# Step 4: Save Results as CSV
# ============================================================================

cat("\n=== Step 4: Processing and Saving Results ===\n")

if (nrow(res_lefse) > 0) {
  t1$res_diff$Taxa <- gsub(".*\\|", "", t1$res_diff$Taxa)
  res_lefse <- t1$res_diff

  final_results <- res_lefse %>%
    rename(Quorum_Sensing_Class = Taxa) %>%
    left_join(feature_metadata %>% distinct(Quorum_Sensing_Class, .keep_all = TRUE), by = "Quorum_Sensing_Class") %>%
    select(
      Quorum_Sensing_Class,
      Group,
      LDA,
      any_of(c("P.value", "P.adj", "P.unadj", "Significance"))
    ) %>%
    arrange(Group, desc(LDA))

  output_csv <- file.path(output_dir, "LEfSe_LDA_results.csv")
  write.csv(final_results, output_csv, row.names = FALSE)
  cat("Results table saved to:", output_csv, "\n")

  # ============================================================================
  # Step 5: Generate LDA Bar Plot
  # ============================================================================

  cat("\n=== Step 5: Generating LDA Bar Plot ===\n")

  full_res_diff <- t1$res_diff

  plot_data <- full_res_diff %>%
    filter(LDA >= LDA_THRESHOLD) %>%
    group_by(Group) %>%
    arrange(desc(LDA)) %>%
    slice_head(n = 10) %>%
    ungroup()

  t1$res_diff <- plot_data

  num_features_plot <- nrow(plot_data)
  cat("Number of features in plot:", num_features_plot, "\n")

  if (num_features_plot > 0) {
    plot_height <- max(6, num_features_plot * 0.25)
    cat("Plot height calculated:", round(plot_height, 2), "inches\n")

    p_bar <- t1$plot_diff_bar(
      use_number = 1:50,
      threshold = 0,
      group_order = c("Control", "Model", "LVX", "HQQD"),
      color_values = group_colors,
      width = 0.7
    )

    p_bar <- p_bar +
      ggtitle("Differential QS Classes (LEfSe)") +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        legend.title = element_blank()
      )

    output_plot_png <- file.path(output_dir, "LEfSe_LDA_barplot.png")
    ggsave(output_plot_png, p_bar, width = 10, height = plot_height, dpi = 300)
    cat("Bar plot saved to:", output_plot_png, "\n")

    output_plot_pdf <- file.path(output_dir, "LEfSe_LDA_barplot.pdf")
    ggsave(output_plot_pdf, p_bar, width = 10, height = plot_height)
    cat("Bar plot saved to:", output_plot_pdf, "\n")
  } else {
    cat("No features passed the LDA threshold for plotting.\n")
  }

  t1$res_diff <- full_res_diff

} else {
  cat("No significant features found with current Alpha threshold.\n")
}

cat("\nAnalysis Complete!\n")
