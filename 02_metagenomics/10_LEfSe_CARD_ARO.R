# ============================================================================
# LEfSe Analysis for CARD Antibiotic Resistance Genes (ARO)
# File: 10_LEfSe_CARD_ARO.R
# Description: Performs LEfSe analysis on CARD ARO and generates
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

input_file <- "CARD.ARO.relabundance.xlsx"

output_dir <- "10_CARD_Analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ALPHA_VALUE <- 0.05
LDA_THRESHOLD <- 2.0
TOP_N <- 5

# ============================================================================
# Step 1: Load and Data Preparation
# ============================================================================

cat("\n=== Step 1: Loading and preparing data ===\n")

if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

raw_data <- read_excel(input_file)
cat("Original data dimensions:", dim(raw_data), "\n")

metadata_cols <- c("ARO_accession", "ARO_name", "ARO_description", "Resistance", "Resistance Mechanism")
sample_cols <- setdiff(colnames(raw_data), metadata_cols)

cat("Sample columns identified:", length(sample_cols), "\n")
print(sample_cols)

# 1. Feature metadata
feature_metadata <- raw_data %>%
  select(all_of(metadata_cols)) %>%
  distinct(ARO_name, .keep_all = TRUE)

# 2. OTU table
if (any(duplicated(raw_data$ARO_name))) {
  warning("Duplicate ARO_name found! Aggregating by ARO_name.")
  otu_table <- raw_data %>%
    group_by(ARO_name) %>%
    summarise(across(all_of(sample_cols), sum), .groups = "drop") %>%
    tibble::column_to_rownames("ARO_name") %>%
    as.data.frame()
} else {
  otu_table <- raw_data %>%
    select(all_of(c("ARO_name", sample_cols))) %>%
    tibble::column_to_rownames("ARO_name") %>%
    as.data.frame()
}

# Scale relative abundances to pseudo-counts (×10^6 → PPM)
otu_table <- otu_table * 1e6

# 3. Taxonomy table (pseudo-taxonomy for microeco)
tax_table <- feature_metadata %>%
  mutate(
    Kingdom = "CARD_ARO",
    Phylum = "Unassigned",
    Class = "Unassigned",
    Order = "Unassigned",
    Family = "Unassigned",
    Genus = "Unassigned",
    Species = ARO_name
  ) %>%
  select(Species, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  as.data.frame()
rownames(tax_table) <- tax_table$Species

tax_table <- tax_table[rownames(otu_table), , drop = FALSE]

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
    rename(ARO_name = Taxa) %>%
    left_join(feature_metadata, by = "ARO_name") %>%
    select(
      ARO_name,
      ARO_accession,
      ARO_description,
      Resistance,
      `Resistance Mechanism`,
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
    slice_head(n = TOP_N) %>%
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
      ggtitle("Top 5 Differential CARD AROs (LEfSe)") +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        legend.title = element_blank()
      )

    output_plot_png <- file.path(output_dir, "LEfSe_LDA_barplot_top5.png")
    ggsave(output_plot_png, p_bar, width = 10, height = plot_height, dpi = 300)
    cat("Bar plot saved to:", output_plot_png, "\n")

    output_plot_pdf <- file.path(output_dir, "LEfSe_LDA_barplot_top5.pdf")
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
