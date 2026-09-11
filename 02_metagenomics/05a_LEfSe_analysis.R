# ============================================================================
# LEfSe Analysis (Linear discriminant analysis Effect Size)
# Author: Generated for metagenomics analysis
# Date: 2026-01-29
# Description: LEfSe analysis only - results saved for later plotting
#              This separates analysis from visualization for efficiency
# ============================================================================

# Clear environment
rm(list = ls())

# Load required libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readxl,      # For reading Excel files
  microeco,    # For LEfSe analysis
  magrittr,    # For pipe operations
  dplyr,       # For data manipulation
  tidyr,       # For data reshaping
  writexl      # For writing Excel files
)

# Define color palette (same as PCoA analysis)
group_colors <- c(
  "Control" = "#3B4992", 
  "Model" = "#EE0000", 
  "LVX" = "#008B45", 
  "HQQD" = "#FF8C00"
)

# ============================================================================
# Configuration Parameters
# ============================================================================

# Input file
input_file <- "taxonomy.all_level.relabundance.xlsx"

# Output directory
output_dir <- "05_LEfSe_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# LEfSe parameters
ALPHA_VALUE <- 0.05        # Significance level for Kruskal-Wallis test
LDA_THRESHOLD <- 2.0       # LDA score threshold (you can adjust to 3.0 or 4.0)
                           # Higher values = fewer but more robust features

# ============================================================================
# Step 1: Load and Prepare Data
# ============================================================================

cat("\n=== Step 1: Loading and preparing data ===\n")

# Read the Excel file
raw_data <- read_excel(input_file)
cat("Original data dimensions:", dim(raw_data), "\n")
cat("Column names:", paste(colnames(raw_data), collapse=", "), "\n\n")

# ============================================================================
# Step 2: Clean Taxonomic Data
# ============================================================================

cat("=== Step 2: Cleaning taxonomic data ===\n")

# Keep incomplete lower-rank annotations. Requiring every rank down to
# species to be classified removes valid higher-rank abundance and biases
# Phylum/Genus-level LEfSe results.
taxonomy_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
sample_cols <- setdiff(colnames(raw_data), taxonomy_cols)

cleaned_data <- raw_data %>%
  mutate(across(all_of(taxonomy_cols), ~ as.character(.x))) %>%
  filter(if_any(all_of(taxonomy_cols), ~ !is.na(.x) & trimws(.x) != ""))

cat("Data after cleaning:", dim(cleaned_data), "\n")
cat("Rows removed:", nrow(raw_data) - nrow(cleaned_data), "\n\n")

# Check if we still have data
if (nrow(cleaned_data) == 0) {
  stop("ERROR: All rows were filtered out! Please check your filtering criteria.")
}

# ============================================================================
# Step 3: Prepare Data for microeco
# ============================================================================

cat("=== Step 3: Preparing data for microeco ===\n")

# Separate taxonomy and abundance data
sample_cols <- setdiff(colnames(cleaned_data), taxonomy_cols)

cat("Number of samples:", length(sample_cols), "\n")
cat("Sample names:", paste(sample_cols, collapse=", "), "\n\n")

# Create OTU table (features × samples)
cleaned_data$OTU_ID <- paste0("OTU_", sprintf("%05d", 1:nrow(cleaned_data)))

otu_table <- cleaned_data %>%
  select(OTU_ID, all_of(sample_cols)) %>%
  tibble::column_to_rownames("OTU_ID") %>%
  as.data.frame()

# Create taxonomy table
tax_table <- cleaned_data %>%
  select(OTU_ID, all_of(taxonomy_cols)) %>%
  tibble::column_to_rownames("OTU_ID") %>%
  as.data.frame()

# Rename taxonomy columns to match microeco format
colnames(tax_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

cat("OTU table dimensions:", dim(otu_table), "\n")
cat("Taxonomy table dimensions:", dim(tax_table), "\n\n")

# ============================================================================
# Step 4: Create Sample Metadata
# ============================================================================

cat("=== Step 4: Creating sample metadata ===\n")

# Create sample metadata based on sample name patterns
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

# Convert Group to factor with proper ordering
sample_metadata$Group <- factor(sample_metadata$Group, 
                                levels = c("Control", "Model", "LVX", "HQQD"))

rownames(sample_metadata) <- sample_metadata$Sample

# Check for unknown groups
if (any(sample_metadata$Group == "Unknown")) {
  warning("Some samples could not be assigned to a group!")
  print(sample_metadata[sample_metadata$Group == "Unknown", ])
}

cat("Sample metadata:\n")
print(table(sample_metadata$Group))
cat("\n")

# ============================================================================
# Step 5: Create microtable Object
# ============================================================================

cat("=== Step 5: Creating microtable object ===\n")

# Create microtable object
dataset <- microtable$new(
  sample_table = sample_metadata,
  otu_table = otu_table,
  tax_table = tax_table
)

# Display basic information
cat("\nMicrotable object created successfully!\n")
cat("Number of samples:", ncol(dataset$otu_table), "\n")
cat("Number of OTUs:", nrow(dataset$otu_table), "\n\n")

# ============================================================================
# Step 6: Run LEfSe Analysis
# ============================================================================

cat("=== Step 6: Running LEfSe analysis ===\n")
cat("Parameters:\n")
cat("  Alpha:", ALPHA_VALUE, "\n")
cat("  LDA threshold:", LDA_THRESHOLD, "\n\n")

# Perform LEfSe analysis using trans_diff class
lefse_result <- trans_diff$new(
  dataset = dataset,
  method = "lefse",
  group = "Group",
  alpha = ALPHA_VALUE,
  lefse_subgroup = NULL,
  p_adjust_method = "fdr"
)

# Extract differential features and apply the configured LDA threshold to the
# primary result table. The unfiltered LEfSe output is saved separately below.
diff_features_all <- lefse_result$res_diff
diff_features <- diff_features_all %>%
  filter(LDA >= LDA_THRESHOLD)

lefse_result_filtered <- lefse_result
lefse_result_filtered$res_diff <- diff_features

cat("Total differential features found before LDA filtering:", nrow(diff_features_all), "\n")
cat("Differential features retained with LDA >=", LDA_THRESHOLD, ":", nrow(diff_features), "\n\n")

# Display top features for each group
cat("Top differential features by group:\n")
for (grp in levels(sample_metadata$Group)) {
  group_features <- diff_features %>%
    filter(Group == grp) %>%
    arrange(desc(LDA)) %>%
    head(10)
  
  if (nrow(group_features) > 0) {
    cat(paste0("\n", grp, " (n=", nrow(filter(diff_features, Group == grp)), "):\n"))
    cat("  Top features (by LDA score):\n")
    for (i in 1:min(5, nrow(group_features))) {
      cat(sprintf("    %d. %s (LDA=%.2f)\n", 
                  i, 
                  group_features$Taxa[i], 
                  group_features$LDA[i]))
    }
  }
}
cat("\n")

# ============================================================================
# Step 7: Save Results for Later Plotting
# ============================================================================

cat("=== Step 7: Saving results ===\n")

# Save LDA-filtered differential features to Excel
output_file_diff <- file.path(output_dir, "LEfSe_differential_features.xlsx")
write_xlsx(diff_features, output_file_diff)
cat("Differential features saved to:", output_file_diff, "\n")

# Save all LEfSe significant features before LDA filtering for auditability
output_file_diff_all <- file.path(output_dir, "LEfSe_differential_features_all_before_LDA_filter.xlsx")
write_xlsx(diff_features_all, output_file_diff_all)
cat("Unfiltered LEfSe features saved to:", output_file_diff_all, "\n")

# Save the entire LEfSe result object for plotting
output_file_rds <- file.path(output_dir, "LEfSe_result_object.rds")
saveRDS(lefse_result_filtered, output_file_rds)
cat("LDA-filtered LEfSe result object saved to:", output_file_rds, "\n")

# Save additional context for plotting
output_file_context <- file.path(output_dir, "LEfSe_analysis_context.rds")
saveRDS(list(
  group_colors = group_colors,
  sample_metadata = sample_metadata,
  LDA_THRESHOLD = LDA_THRESHOLD,
  ALPHA_VALUE = ALPHA_VALUE
), output_file_context)
cat("Analysis context saved to:", output_file_context, "\n\n")

# ============================================================================
# Step 8: Summary Statistics
# ============================================================================

cat("=== Step 8: Summary Statistics ===\n\n")

# Create summary by group
summary_stats <- diff_features %>%
  group_by(Group) %>%
  summarise(
    N_features = n(),
    Mean_LDA = mean(LDA, na.rm = TRUE),
    Median_LDA = median(LDA, na.rm = TRUE),
    Max_LDA = max(LDA, na.rm = TRUE),
    Min_LDA = min(LDA, na.rm = TRUE),
    .groups = "drop"
  )

cat("Summary by group:\n")
print(summary_stats)
cat("\n")

# Save summary statistics
output_file_summary <- file.path(output_dir, "LEfSe_summary_statistics.xlsx")
write_xlsx(summary_stats, output_file_summary)
cat("Summary statistics saved to:", output_file_summary, "\n\n")

# ============================================================================
# Analysis Complete
# ============================================================================

cat("\n", rep("=", 70), "\n", sep="")
cat("LEfSe analysis complete!\n")
cat("Results saved to:", output_dir, "\n")
cat("\nGenerated files:\n")
cat("  1. LEfSe_differential_features.xlsx - All differential features\n")
cat("  2. LEfSe_result_object.rds - LEfSe object for plotting\n")
cat("  3. LEfSe_analysis_context.rds - Analysis parameters\n")
cat("  4. LEfSe_summary_statistics.xlsx - Summary statistics by group\n")
cat("\nNext step:\n")
cat("  Run 05b_LEfSe_plot_cladogram.R to generate visualizations\n")
cat(rep("=", 70), "\n\n", sep="")
