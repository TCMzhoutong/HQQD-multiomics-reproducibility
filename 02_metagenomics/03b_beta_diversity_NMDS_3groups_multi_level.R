# ============================================================================
# Beta Diversity NMDS Analysis (3 Groups, Multi-Taxonomic Levels)
# Author: Generated for metagenomics analysis
# Date: 2026-01-28
# Description: NMDS (Non-metric Multidimensional Scaling) for beta diversity
#              at multiple taxonomic levels (Phylum, Genus, Species)
#              for 3 groups (Control, Model, HQQD), with PERMANOVA and Stress
# ============================================================================

# Clear environment
rm(list = ls())

# Load required libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(vegan)
library(writexl)
library(cowplot)

# Define color palette (3 groups only)
group_colors <- c("Control" = "#3B4992", 
                 "Model" = "#EE0000", 
                 "HQQD" = "#FF8C00")

# Define taxonomic levels to analyze
taxonomic_levels <- list(
  list(file = "taxonomy.phylum.relabundance.xlsx", 
       level = "Phylum", 
       output_dir = "03b_beta_NMDS_3groups_Phylum"),
  list(file = "taxonomy.genus.relabundance.xlsx", 
       level = "Genus", 
       output_dir = "03b_beta_NMDS_3groups_Genus"),
  list(file = "taxonomy.species.relabundance.xlsx", 
       level = "Species", 
       output_dir = "03b_beta_NMDS_3groups_Species")
)

# NMDS plotting function (with PERMANOVA and Stress annotation)
plot_nmds <- function(nmds_coords, metadata, method_name, taxonomic_level, permanova_result, stress_value) {
  # Merge NMDS coordinates with metadata
  plot_data <- nmds_coords %>%
    left_join(metadata, by = "Sample")
  
  # Get column names for NMDS axes
  axis_cols <- colnames(plot_data)[1:2]
  nmds1_col <- axis_cols[1]
  nmds2_col <- axis_cols[2]
  
  # Calculate centroids
  centroids <- plot_data %>%
    group_by(Group) %>%
    summarise(
      centroid_x = mean(.data[[nmds1_col]], na.rm = TRUE),
      centroid_y = mean(.data[[nmds2_col]], na.rm = TRUE),
      .groups = 'drop'
    )
  plot_data <- plot_data %>% left_join(centroids, by = "Group")
  
  # PERMANOVA annotation (top-left, italic)
  permanova_text <- sprintf("PERMANOVA: R² = %.3f, p = %.3f", 
                           permanova_result$R2, 
                           permanova_result$P_value)
  
  # Stress annotation (bottom-right, regular font)
  stress_text <- sprintf("Stress = %.3f", stress_value)
  
  # Plot
  p <- ggplot(plot_data, aes(x = .data[[nmds1_col]], y = .data[[nmds2_col]], color = Group, fill = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    stat_ellipse(aes(fill = Group, color = NULL), 
                 geom = "polygon", 
                 alpha = 0.2, 
                 color = "grey60",
                 linewidth = 0.3,
                 level = 0.95,
                 show.legend = FALSE) +
    geom_segment(aes(xend = centroid_x, yend = centroid_y),
                 color = "black", linetype = "dashed", alpha = 0.3, linewidth = 0.3, show.legend = FALSE) +
    geom_point(size = 2, alpha = 1, shape = 16) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    labs(title = paste0("NMDS - ", method_name, " (", taxonomic_level, " Level)"),
         x = "NMDS1",
         y = "NMDS2") +
    # PERMANOVA annotation (top-left, italic)
    annotate("text", x = -Inf, y = Inf, label = permanova_text,
             hjust = -0.1, vjust = 1.5, size = 3, family = "sans", fontface = "italic") +
    # Stress annotation (bottom-right, regular font)
    annotate("text", x = Inf, y = -Inf, label = stress_text,
             hjust = 1.1, vjust = -0.5, size = 3, family = "sans", fontface = "plain") +
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

# Main analysis loop
for (tax_level in taxonomic_levels) {
  cat("\n============================================================================\n")
  cat(paste0("=== Analyzing ", tax_level$level, " Level (NMDS, 3 Groups) ===\n"))
  cat("============================================================================\n")
  output_dir <- tax_level$output_dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(paste0("Created output directory: ", output_dir, "\n"))
  }
  
  cat("\n=== Loading Abundance Data ===\n")
  abundance_data <- read_excel(tax_level$file)
  cat("Data loaded successfully!\n")
  
  taxonomy_col <- colnames(abundance_data)[1]
  sample_cols <- colnames(abundance_data)[-1]
  
  # Create sample metadata
  sample_metadata <- data.frame(
    Sample = sample_cols,
    Group = case_when(
      grepl("^K", sample_cols) ~ "Control",
      grepl("^M", sample_cols) ~ "Model",
      grepl("^Z", sample_cols) ~ "HQQD",
      TRUE ~ "Unknown"
    )
  )
  
  # Filter to keep only Control, Model, HQQD (exclude LVX)
  sample_metadata <- sample_metadata %>% filter(Group %in% c("Control", "Model", "HQQD"))
  sample_metadata$Group <- factor(sample_metadata$Group, levels = c("Control", "Model", "HQQD"))
  
  samples_to_keep <- sample_metadata$Sample
  
  cat("Sample distribution (after filtering to 3 groups):\n")
  print(table(sample_metadata$Group))
  
  # Prepare abundance matrix (only for 3 groups)
  abundance_matrix <- abundance_data[, samples_to_keep, drop = FALSE] %>%
    as.data.frame() %>%
    t()
  colnames(abundance_matrix) <- abundance_data[[taxonomy_col]]
  rownames(abundance_matrix) <- samples_to_keep
  
  # Remove taxa with zero abundance across all 3-group samples
  abundance_matrix <- abundance_matrix[, colSums(abundance_matrix, na.rm = TRUE) > 0]
  cat("After filtering zero-abundance taxa:", ncol(abundance_matrix), "taxa remain\n")
  
  # Distance matrices
  cat("\n=== Calculating Distance Matrices ===\n")
  dist_jaccard <- vegdist(abundance_matrix, method = "jaccard", binary = TRUE)
  dist_bray <- vegdist(abundance_matrix, method = "bray")
  cat("Distance matrices calculated!\n")
  
  # NMDS analysis
  cat("\n=== Performing NMDS ===\n")
  set.seed(123)
  cat("Binary Jaccard NMDS...\n")
  nmds_jaccard <- metaMDS(dist_jaccard, k = 2, trymax = 100)
  cat("Stress:", nmds_jaccard$stress, "\n")
  
  cat("Bray-Curtis NMDS...\n")
  nmds_bray <- metaMDS(dist_bray, k = 2, trymax = 100)
  cat("Stress:", nmds_bray$stress, "\n")
  
  # Extract NMDS coordinates
  nmds_jaccard_coords <- as.data.frame(nmds_jaccard$points)
  colnames(nmds_jaccard_coords) <- c("NMDS1", "NMDS2")
  nmds_jaccard_coords$Sample <- rownames(nmds_jaccard_coords)
  
  nmds_bray_coords <- as.data.frame(nmds_bray$points)
  colnames(nmds_bray_coords) <- c("NMDS1", "NMDS2")
  nmds_bray_coords$Sample <- rownames(nmds_bray_coords)
  
  # PERMANOVA
  cat("\n=== Performing PERMANOVA ===\n")
  permanova_jaccard <- adonis2(dist_jaccard ~ Group, data = sample_metadata, permutations = 999)
  permanova_bray <- adonis2(dist_bray ~ Group, data = sample_metadata, permutations = 999)
  
  permanova_results <- data.frame(
    Method = c("Binary Jaccard", "Bray-Curtis"),
    R2 = c(permanova_jaccard$R2[1], permanova_bray$R2[1]),
    F_statistic = c(permanova_jaccard$F[1], permanova_bray$F[1]),
    P_value = c(permanova_jaccard$`Pr(>F)`[1], permanova_bray$`Pr(>F)`[1]),
    Stress = c(nmds_jaccard$stress, nmds_bray$stress)
  )
  write.csv(permanova_results, file.path(output_dir, "PERMANOVA_and_Stress_3groups.csv"), row.names = FALSE)
  cat("PERMANOVA results saved!\n")
  
  # Plot with Stress values
  cat("\n=== Generating NMDS Plots ===\n")
  p_jaccard <- plot_nmds(nmds_jaccard_coords, sample_metadata, "Binary Jaccard", 
                         tax_level$level, permanova_results[1, ], nmds_jaccard$stress)
  p_bray <- plot_nmds(nmds_bray_coords, sample_metadata, "Bray-Curtis", 
                      tax_level$level, permanova_results[2, ], nmds_bray$stress)
  
  ggsave(file.path(output_dir, "NMDS_Binary_Jaccard_3Groups.png"), plot = p_jaccard, width = 6, height = 5, dpi = 300)
  ggsave(file.path(output_dir, "NMDS_Binary_Jaccard_3Groups.pdf"), plot = p_jaccard, width = 6, height = 5, device = cairo_pdf)
  ggsave(file.path(output_dir, "NMDS_Bray_Curtis_3Groups.png"), plot = p_bray, width = 6, height = 5, dpi = 300)
  ggsave(file.path(output_dir, "NMDS_Bray_Curtis_3Groups.pdf"), plot = p_bray, width = 6, height = 5, device = cairo_pdf)
  cat("Individual NMDS plots saved!\n")
  
  # Save NMDS coordinates
  write.csv(nmds_jaccard_coords, file.path(output_dir, "NMDS_Binary_Jaccard_coordinates_3groups.csv"), row.names = FALSE)
  write.csv(nmds_bray_coords, file.path(output_dir, "NMDS_Bray_Curtis_coordinates_3groups.csv"), row.names = FALSE)
  cat("NMDS coordinates saved!\n")
  
  # Print summary
  cat("\n=== Analysis Summary ===\n")
  cat("Binary Jaccard:\n")
  cat("  Stress:", round(nmds_jaccard$stress, 3), 
      ifelse(nmds_jaccard$stress < 0.05, "(Excellent)", 
      ifelse(nmds_jaccard$stress < 0.10, "(Good)",
      ifelse(nmds_jaccard$stress < 0.20, "(Acceptable)", "(Poor - consider reanalysis)"))), "\n")
  cat("  PERMANOVA R²:", round(permanova_results$R2[1], 3), "p =", round(permanova_results$P_value[1], 3), "\n")
  cat("\nBray-Curtis:\n")
  cat("  Stress:", round(nmds_bray$stress, 3),
      ifelse(nmds_bray$stress < 0.05, "(Excellent)", 
      ifelse(nmds_bray$stress < 0.10, "(Good)",
      ifelse(nmds_bray$stress < 0.20, "(Acceptable)", "(Poor - consider reanalysis)"))), "\n")
  cat("  PERMANOVA R²:", round(permanova_results$R2[2], 3), "p =", round(permanova_results$P_value[2], 3), "\n")
  
  cat("\n============================================================================\n")
  cat(paste0("=== ", tax_level$level, " Level (NMDS, 3 Groups) Analysis Complete ===\n"))
  cat("============================================================================\n\n")
}

cat("\n============================================================================\n")
cat("=== All 3-Group Taxonomic Level NMDS Analyses Complete ===\n")
cat("============================================================================\n\n")

cat("=== Key Points ===\n")
cat("1. Analyzed 3 groups only: Control, Model, HQQD (LVX excluded)\n")
cat("2. NMDS uses relative abundance data (correct for beta diversity)\n")
cat("3. Stress values are displayed on plots (bottom-right corner)\n")
cat("4. PERMANOVA tests group differences (top-left corner)\n")
cat("5. Stress interpretation:\n")
cat("   • < 0.05: Excellent representation\n")
cat("   • < 0.10: Good representation\n")
cat("   • < 0.20: Acceptable representation\n")
cat("   • > 0.20: Poor fit (consider adding dimensions or reanalysis)\n")
cat("\n6. Both metrics (Jaccard and Bray-Curtis) provide complementary views:\n")
cat("   • Binary Jaccard: Presence/absence (unweighted)\n")
cat("   • Bray-Curtis: Abundance-weighted\n")
cat("\n=== Comparison with 4-Group Analysis ===\n")
cat("• 4-group analysis: 03_beta_diversity_NMDS_multi_level.R\n")
cat("• 3-group analysis: 03b_beta_diversity_NMDS_3groups_multi_level.R\n")
cat("• Compare PERMANOVA R² and p-values to assess LVX group's effect\n")
cat("• Stress values may differ slightly due to different sample sizes\n")
cat("\n============================================================================\n\n")
