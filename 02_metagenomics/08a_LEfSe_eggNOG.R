# ============================================================================
# LEfSe Analysis for eggNOG Class
# File: 08a_LEfSe_eggNOG_Class.R
# Description: Performs LEfSe analysis on eggNOG Classes and generates
#              barplots for top LDA features.
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
  ggplot2,     # For plotting
  stringr      # For string manipulation
)

# Define color palette (consistent with previous analyses)
group_colors <- c(
  "Control" = "#3B4992", 
  "Model" = "#EE0000", 
  "LVX" = "#008B45", 
  "HQQD" = "#FF8C00"
)

# ============================================================================
# Configuration Parameters
# ============================================================================

# Input file path
input_file <- "eggNOG.NOG.relabundance.xlsx"

# Output directory
output_dir <- "08_LEfSe_eggNOG_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# LEfSe parameters
ALPHA_VALUE <- 0.05        # Significance level for Kruskal-Wallis test
LDA_THRESHOLD <- 2.0       # LDA score threshold

# ============================================================================
# Step 1: Load and Data Preparation
# ============================================================================

cat("\n=== Step 1: Loading and preparing data ===\n")

# Check if file exists
if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

# Read the Excel file
raw_data <- read_excel(input_file)
cat("Original data dimensions:", dim(raw_data), "\n")

# Identify sample columns (assuming K, M, Y, Z prefixes based on requirements)
# Extract metadata columns first
metadata_cols <- c("eggNOG", "descrition", "eggNOG_Class", "Category")
all_cols <- colnames(raw_data)
sample_cols <- setdiff(all_cols, metadata_cols)

cat("Sample columns identified:", length(sample_cols), "\n")
print(sample_cols)

# Separate Expression Data and Feature Metadata
# Use eggNOG (e.g. COG0001) as unique feature identifier
raw_data$FeatureID <- raw_data$eggNOG

# 1. Prepare Feature Metadata (for later merging)
feature_metadata <- raw_data %>%
  select(FeatureID, all_of(metadata_cols))

# 2. Prepare OTU Table (Abundance Table)
# Row names: FeatureID, Columns: Samples
otu_table <- raw_data %>%
  select(FeatureID, all_of(sample_cols)) %>%
  tibble::column_to_rownames("FeatureID") %>%
  as.data.frame()

# Scale relative abundances to pseudo-counts (×10^6 → PPM)
# Required because MASS::lda() inside microeco LEfSe returns log10(effect_size),
# and near-zero relative abundance values produce log10(0) = -Inf LDA scores.
otu_table <- otu_table * 1e6

# 3. Prepare Taxonomy Table (Pseudo-taxonomy for microeco)
# We map 'Category' to 'Kingdom' and 'eggNOG' to 'Species'
tax_table <- feature_metadata %>%
  mutate(
    Kingdom = Category,
    Phylum = "Unassigned",
    Class = "Unassigned",
    Order = "Unassigned",
    Family = "Unassigned",
    Genus = "Unassigned",
    Species = eggNOG
  ) %>%
  select(FeatureID, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  tibble::column_to_rownames("FeatureID") %>%
  as.data.frame()

# 4. Prepare Sample Metadata
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

# Ensure Group levels are correct
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

# Initialize trans_diff object
# taxa_level = "Species" tells it to use the column 'Species' from tax_table which we mapped to eggNOG_Class
t1 <- trans_diff$new(
  dataset = dataset, 
  method = "lefse", 
  group = "Group", 
  alpha = ALPHA_VALUE, 
  p_adjust_method = "fdr", 
  taxa_level = "Species"
)

# Print comparison
cat("LDA Threshold:", LDA_THRESHOLD, "\n")
cat("Alpha Value:", ALPHA_VALUE, "\n")

# The results are stored in t1$res_diff
res_lefse <- t1$res_diff

cat("Total differential features found:", nrow(res_lefse), "\n")

# ============================================================================
# Step 4: Process and Save Results
# ============================================================================

cat("\n=== Step 4: Processing and Saving Results ===\n")

if (nrow(res_lefse) > 0) {
  # Fix Taxa names: microeco returns full lineage "Kingdom|Phylum|...|Species"
  # We only need the last part (Species level which corresponds to eggNOG_Class)
  t1$res_diff$Taxa <- gsub(".*\\|", "", t1$res_diff$Taxa)
  res_lefse <- t1$res_diff
  
  # [MATCHING 05a BEHAVIOR]
  # Do NOT filter by LDA threshold for the Excel table
  # We preserve all significant features (Alpha < 0.05) found by LEfSe
  res_lefse_export <- res_lefse
  
  cat("Total significant features (Alpha <", ALPHA_VALUE, "):", nrow(res_lefse_export), "\n")
  
  # Join with our feature_metadata to get Category info
  final_results <- res_lefse_export %>%
    rename(eggNOG = Taxa) %>%
    # Join with feature_metadata to get descrition, eggNOG_Class, Category
    left_join(feature_metadata %>% select(-FeatureID) %>% distinct(eggNOG, .keep_all = TRUE), by = "eggNOG") %>%
    select(
      eggNOG,
      descrition,
      eggNOG_Class,
      Category,
      Group,
      LDA,
      any_of(c("P.value", "P.adj", "P.unadj", "Significance"))
    ) %>%
    arrange(Group, desc(LDA))
  
  # Save to CSV
  output_csv <- file.path(output_dir, "LEfSe_LDA_results.csv")
  write.csv(final_results, output_csv, row.names = FALSE)
  cat("Results table saved to:", output_csv, "\n")
  
  # ============================================================================
  # Step 5: Generate Bar Plot
  # ============================================================================
  
  cat("\n=== Step 5: Generating Bar Plot ===\n")
  
  # 1. Back up full results
  full_res_diff <- t1$res_diff
  
  # 2. Filter top 10 per group based on LDA for plotting
  # [IMPORTANT] Apply LDA Threshold for plot
  plot_data <- full_res_diff %>%
    filter(LDA >= LDA_THRESHOLD) %>%
    group_by(Group) %>%
    arrange(desc(LDA)) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  # Update t1$res_diff for plotting function
  t1$res_diff <- plot_data
  
  num_features_plot <- nrow(t1$res_diff)
  cat("Number of features in plot (Top 5/group):", num_features_plot, "\n")
  
  if (num_features_plot > 0) {
    # Calculate dynamic height
    plot_height <- max(6, num_features_plot * 0.25)
    cat("Plot height calculated:", round(plot_height, 2), "inches\n")
    
    p_bar <- t1$plot_diff_bar(
      use_number = 1:50,       # Set high enough
      threshold = 0,           # Filter already applied
      group_order = c("Control", "Model", "LVX", "HQQD"),
      color_values = group_colors,
      width = 0.7
    )
    
    # Customize plot appearance
    p_bar <- p_bar + 
      ggtitle("Top 5 LEfSe Results (eggNOG NOG)") +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(), 
        legend.title = element_blank()
      )
    
    # Save plot (PNG and PDF)
    output_plot_png <- file.path(output_dir, "LEfSe_LDA_barplot_top5.png")
    ggsave(output_plot_png, p_bar, width = 10, height = plot_height, dpi = 300)
    cat("Bar plot saved to:", output_plot_png, "\n")
  
    output_plot_pdf <- file.path(output_dir, "LEfSe_LDA_barplot_top5.pdf")
    ggsave(output_plot_pdf, p_bar, width = 10, height = plot_height)
    cat("Bar plot saved to:", output_plot_pdf, "\n")
  } else {
    cat("No features passed the LDA threshold for plotting.\n")
  }
  
  # Restore full results
  t1$res_diff <- full_res_diff
  
} else {
  cat("No significant features found with current Alpha threshold.\n")
}

cat("\nAnalysis Complete!\n")
