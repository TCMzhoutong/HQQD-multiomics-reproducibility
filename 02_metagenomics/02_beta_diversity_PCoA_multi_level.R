# ============================================================================
# Beta Diversity PCoA Analysis (Multi-Taxonomic Levels)
# Author: Generated for metagenomics analysis
# Date: 2026-01-28
# Description: Principal Coordinates Analysis (PCoA) for beta diversity
#              at multiple taxonomic levels (Phylum, Genus, Species)
#              with PERMANOVA statistical testing
# ============================================================================

# Clear environment
rm(list = ls())

# Load required libraries
library(readxl)      # For reading Excel files
library(ggplot2)     # For visualization
library(dplyr)       # For data manipulation
library(vegan)       # For distance matrices and diversity analysis
library(ape)         # For PCoA
library(writexl)     # For writing Excel files
library(cowplot)     # For plot arrangement

# Define color palette
group_colors <- c("Control" = "#3B4992", 
                 "Model" = "#EE0000", 
                 "LVX" = "#008B45", 
                 "HQQD" = "#FF8C00")

# Define taxonomic levels to analyze
taxonomic_levels <- list(
  list(file = "taxonomy.phylum.relabundance.xlsx", 
       level = "Phylum", 
       output_dir = "02_beta_PCoA_Phylum"),
  list(file = "taxonomy.genus.relabundance.xlsx", 
       level = "Genus", 
       output_dir = "02_beta_PCoA_Genus"),
  list(file = "taxonomy.species.relabundance.xlsx", 
       level = "Species", 
       output_dir = "02_beta_PCoA_Species")
)

# ============================================================================
# Define Analysis Functions
# ============================================================================

# Function to extract PCoA results
extract_pcoa_results <- function(pcoa_obj, metadata, method_name) {
  # Extract coordinates
  coords <- as.data.frame(pcoa_obj$vectors)
  coords$Sample <- rownames(coords)
  
  # Add group information
  coords <- coords %>%
    left_join(metadata, by = "Sample")
  
  # Calculate variance explained
  variance_explained <- pcoa_obj$values$Relative_eig * 100  # Convert to percentage
  
  # Create variance data frame
  variance_df <- data.frame(
    PC = paste0("PC", 1:length(variance_explained)),
    Variance_Explained = variance_explained
  )
  
  cat(paste0("\n", method_name, ":\n"))
  cat("PC1 variance:", round(variance_df$Variance_Explained[1], 2), "%\n")
  cat("PC2 variance:", round(variance_df$Variance_Explained[2], 2), "%\n")
  cat("PC3 variance:", round(variance_df$Variance_Explained[3], 2), "%\n")
  
  return(list(
    coords = coords,
    variance = variance_df
  ))
}

# Function to create PCoA plot with PERMANOVA results
plot_pcoa <- function(coords, variance_df, method_name, taxonomic_level, permanova_result, pc_x = 1, pc_y = 2) {
  
  # Get PC column names
  pc_x_name <- paste0("Axis.", pc_x)
  pc_y_name <- paste0("Axis.", pc_y)
  
  # Get variance explained
  var_x <- variance_df$Variance_Explained[pc_x]
  var_y <- variance_df$Variance_Explained[pc_y]
  
  # Create axis labels with variance explained
  x_label <- sprintf("PC%d (%.2f%%)", pc_x, var_x)
  y_label <- sprintf("PC%d (%.2f%%)", pc_y, var_y)
  
  # Create title with taxonomic level
  plot_title <- paste0("PCoA - ", method_name, " (", taxonomic_level, " Level)")
  
  # Create PERMANOVA text annotation
  permanova_text <- sprintf("PERMANOVA: R² = %.3f, p = %.3f", 
                           permanova_result$R2, 
                           permanova_result$P_value)
  
  # Calculate centroids for each group
  centroids <- coords %>%
    group_by(Group) %>%
    summarise(
      centroid_x = mean(.data[[pc_x_name]], na.rm = TRUE),
      centroid_y = mean(.data[[pc_y_name]], na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Create data for connecting lines (samples to centroids)
  coords_with_centroids <- coords %>%
    left_join(centroids, by = "Group")
  
  # Create plot
  p <- ggplot(coords, aes(x = .data[[pc_x_name]], 
                          y = .data[[pc_y_name]], 
                          color = Group, 
                          fill = Group)) +
    # Add reference lines at 0 first (so they appear behind) - dashed lines
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    # Add 95% confidence ellipses with grey outline
    stat_ellipse(aes(fill = Group, color = NULL), 
                 geom = "polygon", 
                 alpha = 0.2, 
                 color = "grey60",
                 linewidth = 0.3,
                 level = 0.95,
                 show.legend = FALSE) +
    # Add connecting lines from samples to centroids - black dashed lines
    geom_segment(data = coords_with_centroids,
                 aes(x = .data[[pc_x_name]], 
                     y = .data[[pc_y_name]],
                     xend = centroid_x, 
                     yend = centroid_y),
                 color = "black",
                 linetype = "dashed",
                 alpha = 0.3,
                 linewidth = 0.3,
                 show.legend = FALSE) +
    # Add points (small, solid, no transparency)
    geom_point(size = 2, alpha = 1, shape = 16) +
    # Color scheme
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    # Format axis labels - show 0 instead of 0.0
    scale_x_continuous(labels = function(x) {
      sapply(x, function(val) {
        if (is.na(val)) return(NA)
        if (abs(val) < 1e-10) return("0")
        else return(as.character(val))
      })
    }) +
    scale_y_continuous(labels = function(x) {
      sapply(x, function(val) {
        if (is.na(val)) return(NA)
        if (abs(val) < 1e-10) return("0")
        else return(as.character(val))
      })
    }) +
    # Labels
    labs(title = plot_title,
         x = x_label,
         y = y_label) +
    # Add PERMANOVA annotation
    annotate("text", 
             x = -Inf, y = Inf, 
             label = permanova_text,
             hjust = -0.1, vjust = 1.5,
             size = 3,
             family = "sans",
             fontface = "italic") +
    # Theme
    theme_bw(base_size = 10, base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(size = 9, face = "bold"),
      legend.position = c(1.02, 0.8),
      legend.justification = c(0, 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      legend.spacing.x = unit(2, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(5, 50, 5, 5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      aspect.ratio = 1
    )
  
  return(p)
}

# Function to calculate group statistics
calculate_group_stats <- function(coords, method_name) {
  # Get all Axis columns
  axis_cols <- grep("^Axis\\.", colnames(coords), value = TRUE)
  
  # Calculate statistics for each group and axis
  stats_list <- list()
  
  for (axis_col in axis_cols[1:3]) {  # PC1, PC2, PC3
    group_stats <- coords %>%
      group_by(Group) %>%
      summarise(
        Mean = mean(.data[[axis_col]], na.rm = TRUE),
        SD = sd(.data[[axis_col]], na.rm = TRUE),
        SE = SD / sqrt(n()),
        .groups = 'drop'
      ) %>%
      mutate(
        PC = axis_col,
        Method = method_name
      )
    
    stats_list[[axis_col]] <- group_stats
  }
  
  return(bind_rows(stats_list))
}

# ============================================================================
# Main Analysis Loop for Each Taxonomic Level
# ============================================================================

for (tax_level in taxonomic_levels) {
  
  cat("\n")
  cat("============================================================================\n")
  cat(paste0("=== Analyzing ", tax_level$level, " Level ===\n"))
  cat("============================================================================\n")
  
  # Create output directory for this taxonomic level
  output_dir <- tax_level$output_dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(paste0("Created output directory: ", output_dir, "\n"))
  }
  
  # ==========================================================================
  # 1. Data Import and Preparation
  # ==========================================================================
  
  cat("\n=== Loading Abundance Data ===\n")
  
  # Read abundance data
  abundance_data <- read_excel(tax_level$file)
  
  cat("Data loaded successfully!\n")
  cat("Data dimensions:", dim(abundance_data), "\n")
  cat("First few column names:", paste(head(colnames(abundance_data), 10), collapse = ", "), "\n")
  
  # ==========================================================================
  # 2. Data Processing
  # ==========================================================================
  
  cat("\n=== Processing Data ===\n")
  
  # The first column is taxonomy names
  # Remaining columns are sample abundances
  taxonomy_col <- colnames(abundance_data)[1]
  sample_cols <- colnames(abundance_data)[-1]
  
  cat("Taxonomy column:", taxonomy_col, "\n")
  cat("Number of taxa:", nrow(abundance_data), "\n")
  cat("Total samples:", length(sample_cols), "\n")
  
  # Create sample metadata
  sample_metadata <- data.frame(
    Sample = sample_cols,
    Group = case_when(
      grepl("^K", sample_cols) ~ "Control",
      grepl("^M", sample_cols) ~ "Model",
      grepl("^Y", sample_cols) ~ "LVX",
      grepl("^Z", sample_cols) ~ "HQQD",
      TRUE ~ "Unknown"
    )
  )
  
  # Set factor levels for proper ordering
  sample_metadata$Group <- factor(sample_metadata$Group, 
                                   levels = c("Control", "Model", "LVX", "HQQD"))
  
  cat("Sample distribution:\n")
  print(table(sample_metadata$Group))
  
  # Prepare abundance matrix (samples as rows, taxa as columns)
  abundance_matrix <- abundance_data[, -1] %>%
    as.data.frame() %>%
    t()
  
  colnames(abundance_matrix) <- abundance_data[[taxonomy_col]]
  rownames(abundance_matrix) <- sample_cols
  
  # Remove taxa with zero abundance across all samples
  abundance_matrix <- abundance_matrix[, colSums(abundance_matrix, na.rm = TRUE) > 0]
  
  cat("After filtering zero-abundance taxa:", ncol(abundance_matrix), "taxa remain\n")
  
  # ==========================================================================
  # 3. Distance Matrix Calculation
  # ==========================================================================
  
  cat("\n=== Calculating Distance Matrices ===\n")
  
  # Binary Jaccard distance (presence/absence based)
  cat("Calculating Binary Jaccard distance...\n")
  dist_jaccard <- vegdist(abundance_matrix, method = "jaccard", binary = TRUE)
  
  # Bray-Curtis distance (abundance based)
  cat("Calculating Bray-Curtis distance...\n")
  dist_bray <- vegdist(abundance_matrix, method = "bray")
  
  cat("Distance matrices calculated successfully!\n")
  
  # ==========================================================================
  # 4. PCoA Analysis
  # ==========================================================================
  
  cat("\n=== Performing PCoA ===\n")
  
  # Perform PCoA
  cat("Binary Jaccard PCoA...\n")
  pcoa_jaccard <- pcoa(dist_jaccard)
  
  cat("Bray-Curtis PCoA...\n")
  pcoa_bray <- pcoa(dist_bray)
  
  cat("PCoA analysis completed!\n")
  
  # ==========================================================================
  # 5. Extract PCoA Results
  # ==========================================================================
  
  cat("\n=== Extracting PCoA Results ===\n")
  
  # Extract results
  jaccard_results <- extract_pcoa_results(pcoa_jaccard, sample_metadata, "Binary Jaccard")
  bray_results <- extract_pcoa_results(pcoa_bray, sample_metadata, "Bray-Curtis")
  
  # Save coordinate data
  write.csv(jaccard_results$coords, 
            file.path(output_dir, "PCoA_Binary_Jaccard_coordinates.csv"), 
            row.names = FALSE)
  
  write.csv(bray_results$coords, 
            file.path(output_dir, "PCoA_Bray_Curtis_coordinates.csv"), 
            row.names = FALSE)
  
  # Save variance explained
  write.csv(jaccard_results$variance, 
            file.path(output_dir, "PCoA_Binary_Jaccard_variance.csv"), 
            row.names = FALSE)
  
  write.csv(bray_results$variance, 
            file.path(output_dir, "PCoA_Bray_Curtis_variance.csv"), 
            row.names = FALSE)
  
  cat("Coordinate and variance data saved!\n")
  
  # ==========================================================================
  # 6. PERMANOVA Analysis
  # ==========================================================================
  
  cat("\n=== Performing PERMANOVA ===\n")
  
  # Binary Jaccard
  cat("Binary Jaccard PERMANOVA:\n")
  permanova_jaccard <- adonis2(dist_jaccard ~ Group, data = sample_metadata, permutations = 999)
  print(permanova_jaccard)
  
  # Bray-Curtis
  cat("\nBray-Curtis PERMANOVA:\n")
  permanova_bray <- adonis2(dist_bray ~ Group, data = sample_metadata, permutations = 999)
  print(permanova_bray)
  
  # Create PERMANOVA results data frame
  permanova_results <- data.frame(
    Method = c("Binary Jaccard", "Bray-Curtis"),
    R2 = c(permanova_jaccard$R2[1], permanova_bray$R2[1]),
    F_statistic = c(permanova_jaccard$F[1], permanova_bray$F[1]),
    P_value = c(permanova_jaccard$`Pr(>F)`[1], permanova_bray$`Pr(>F)`[1])
  )
  
  write.csv(permanova_results, 
            file.path(output_dir, "PERMANOVA_results.csv"), 
            row.names = FALSE)
  
  cat("\nPERMANOVA results saved!\n")
  
  # ==========================================================================
  # 7. Create PCoA Plots with PERMANOVA Results
  # ==========================================================================
  
  cat("\n=== Generating PCoA Plots ===\n")
  
  # Binary Jaccard PCoA
  p_jaccard <- plot_pcoa(jaccard_results$coords, 
                        jaccard_results$variance, 
                        "Binary Jaccard",
                        tax_level$level,
                        permanova_results[1, ],
                        pc_x = 1, pc_y = 2)
  
  # Bray-Curtis PCoA
  p_bray <- plot_pcoa(bray_results$coords, 
                     bray_results$variance, 
                     "Bray-Curtis",
                     tax_level$level,
                     permanova_results[2, ],
                     pc_x = 1, pc_y = 2)
  
  # Save individual plots
  ggsave(file.path(output_dir, "PCoA_Binary_Jaccard.png"), 
         plot = p_jaccard, 
         width = 6, 
         height = 5, 
         dpi = 300)
  
  ggsave(file.path(output_dir, "PCoA_Binary_Jaccard.pdf"), 
         plot = p_jaccard, 
         width = 6, 
         height = 5,
         device = cairo_pdf)
  
  ggsave(file.path(output_dir, "PCoA_Bray_Curtis.png"), 
         plot = p_bray, 
         width = 6, 
         height = 5, 
         dpi = 300)
  
  ggsave(file.path(output_dir, "PCoA_Bray_Curtis.pdf"), 
         plot = p_bray, 
         width = 6, 
         height = 5,
         device = cairo_pdf)
  
  cat("Individual PCoA plots saved!\n")
  
  # ==========================================================================
  # 8. Calculate Summary Statistics
  # ==========================================================================
  
  cat("\n=== Calculating Summary Statistics ===\n")
  
  # Calculate statistics for each method
  jaccard_stats <- calculate_group_stats(jaccard_results$coords, "Binary Jaccard")
  bray_stats <- calculate_group_stats(bray_results$coords, "Bray-Curtis")
  
  # Combine statistics
  all_stats <- bind_rows(jaccard_stats, bray_stats)
  
  # Save statistics
  write.csv(all_stats, 
            file.path(output_dir, "PCoA_group_statistics.csv"), 
            row.names = FALSE)
  
  cat("Statistics saved!\n")
  
  # ==========================================================================
  # 9. Create Variance Summary
  # ==========================================================================
  
  # Combine variance explained from both methods
  variance_summary <- bind_rows(
    jaccard_results$variance %>% mutate(Method = "Binary Jaccard"),
    bray_results$variance %>% mutate(Method = "Bray-Curtis")
  )
  
  # Save variance summary
  write.csv(variance_summary, 
            file.path(output_dir, "PCoA_variance_explained.csv"), 
            row.names = FALSE)
  
  # ==========================================================================
  # 10. Print Summary for this Taxonomic Level
  # ==========================================================================
  
  cat("\n============================================================================\n")
  cat(paste0("=== ", tax_level$level, " Level Analysis Complete ===\n"))
  cat("============================================================================\n\n")
  
  cat("=== Files Generated ===\n")
  cat("1. PCoA_Binary_Jaccard.png/pdf - Binary Jaccard PCoA plot\n")
  cat("2. PCoA_Bray_Curtis.png/pdf - Bray-Curtis PCoA plot\n")
  cat("3. PCoA_Binary_Jaccard_coordinates.csv - Sample coordinates\n")
  cat("4. PCoA_Bray_Curtis_coordinates.csv - Sample coordinates\n")
  cat("5. PCoA_Binary_Jaccard_variance.csv - Variance explained\n")
  cat("6. PCoA_Bray_Curtis_variance.csv - Variance explained\n")
  cat("7. PCoA_group_statistics.csv - Descriptive statistics by group\n")
  cat("8. PCoA_variance_explained.csv - Variance explained by each PC\n")
  cat("9. PERMANOVA_results.csv - Statistical test results\n")
  
  cat("\n=== Variance Explained ===\n")
  print(variance_summary)
  
  cat("\n=== PERMANOVA Results ===\n")
  print(permanova_results)
  
  cat("\n")
}

# ============================================================================
# Final Summary
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("=== All Taxonomic Level Analyses Complete ===\n")
cat("============================================================================\n\n")

cat("Analyzed taxonomic levels:\n")
for (tax_level in taxonomic_levels) {
  cat("  •", tax_level$level, "- Results in:", tax_level$output_dir, "\n")
}

cat("\n=== Key Features ===\n")
cat("1. Multiple taxonomic levels analyzed (Phylum, Genus, Species)\n")
cat("2. Two distance metrics: Binary Jaccard and Bray-Curtis\n")
cat("3. PERMANOVA statistical testing included\n")
cat("4. PERMANOVA results displayed on plots\n")
cat("5. Taxonomic level labeled in plot titles\n")
cat("6. 95% confidence ellipses for each group\n")
cat("7. Samples connected to group centroids\n")
cat("8. Academic publication-ready formatting\n")

cat("\n=== Interpretation Guide ===\n")
cat("• Binary Jaccard: Presence/absence based (unweighted)\n")
cat("• Bray-Curtis: Abundance based (weighted)\n")
cat("• PERMANOVA R²: Proportion of variance explained by grouping\n")
cat("• PERMANOVA p-value: Statistical significance (p < 0.05 = significant)\n")
cat("• If PERMANOVA is significant: Groups have distinct community compositions\n")
cat("• Compare patterns across taxonomic levels:\n")
cat("  - Phylum: Broad taxonomic differences\n")
cat("  - Genus: Intermediate resolution\n")
cat("  - Species: Finest resolution\n")

cat("\n=== Recommendations for Publication ===\n")
cat("1. Use combined plots (side-by-side comparison) for main figures\n")
cat("2. Report PERMANOVA results in figure captions or results section\n")
cat("3. Compare patterns across different taxonomic levels\n")
cat("4. Font: sans (standard for scientific publications)\n")
cat("5. Resolution: 300 DPI (PNG) + vector (PDF)\n")

cat("\n=== Methods Text Suggestion ===\n")
cat('  "Beta diversity was assessed at three taxonomic levels (phylum, genus,\n')
cat('   and species) using Principal Coordinates Analysis (PCoA) based on\n')
cat('   Binary Jaccard and Bray-Curtis distance matrices. PERMANOVA was used\n')
cat('   to test for significant differences among treatment groups (999\n')
cat('   permutations). All analyses were performed in R using the vegan\n')
cat('   package."\n')

cat("\n============================================================================\n")
cat("=== Analysis Complete ===\n")
cat("============================================================================\n\n")
