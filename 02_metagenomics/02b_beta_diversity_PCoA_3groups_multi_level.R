# ============================================================================
# Beta Diversity PCoA Analysis (3 Groups, Multi-Taxonomic Levels)
# Author: Generated for metagenomics analysis
# Date: 2026-01-28
# Description: Principal Coordinates Analysis (PCoA) for beta diversity
#              at multiple taxonomic levels (Phylum, Genus, Species)
#              for 3 groups (Control, Model, HQQD), with PERMANOVA results
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
                 "HQQD" = "#FF8C00")

# Define taxonomic levels to analyze
taxonomic_levels <- list(
  list(file = "taxonomy.phylum.relabundance.xlsx", 
       level = "Phylum", 
       output_dir = "02b_beta_PCoA_3groups_Phylum"),
  list(file = "taxonomy.genus.relabundance.xlsx", 
       level = "Genus", 
       output_dir = "02b_beta_PCoA_3groups_Genus"),
  list(file = "taxonomy.species.relabundance.xlsx", 
       level = "Species", 
       output_dir = "02b_beta_PCoA_3groups_Species")
)

# ============================================================================
# Define Analysis Functions
# ============================================================================

extract_pcoa_results <- function(pcoa_obj, metadata, method_name) {
  coords <- as.data.frame(pcoa_obj$vectors)
  coords$Sample <- rownames(coords)
  coords <- coords %>% left_join(metadata, by = "Sample")
  variance_explained <- pcoa_obj$values$Relative_eig * 100
  variance_df <- data.frame(
    PC = paste0("PC", 1:length(variance_explained)),
    Variance_Explained = variance_explained
  )
  return(list(coords = coords, variance = variance_df))
}

plot_pcoa <- function(coords, variance_df, method_name, taxonomic_level, permanova_result, pc_x = 1, pc_y = 2) {
  pc_x_name <- paste0("Axis.", pc_x)
  pc_y_name <- paste0("Axis.", pc_y)
  var_x <- variance_df$Variance_Explained[pc_x]
  var_y <- variance_df$Variance_Explained[pc_y]
  x_label <- sprintf("PC%d (%.2f%%)", pc_x, var_x)
  y_label <- sprintf("PC%d (%.2f%%)", pc_y, var_y)
  plot_title <- paste0("PCoA - ", method_name, " (", taxonomic_level, " Level)")
  permanova_text <- sprintf("PERMANOVA: R² = %.3f, p = %.3f", 
                           permanova_result$R2, 
                           permanova_result$P_value)
  centroids <- coords %>%
    group_by(Group) %>%
    summarise(
      centroid_x = mean(.data[[pc_x_name]], na.rm = TRUE),
      centroid_y = mean(.data[[pc_y_name]], na.rm = TRUE),
      .groups = 'drop'
    )
  coords_with_centroids <- coords %>% left_join(centroids, by = "Group")
  p <- ggplot(coords, aes(x = .data[[pc_x_name]], 
                          y = .data[[pc_y_name]], 
                          color = Group, 
                          fill = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    stat_ellipse(aes(fill = Group, color = NULL), 
                 geom = "polygon", 
                 alpha = 0.2, 
                 color = "grey60",
                 linewidth = 0.3,
                 level = 0.95,
                 show.legend = FALSE) +
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
    geom_point(size = 2, alpha = 1, shape = 16) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
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
    labs(title = plot_title,
         x = x_label,
         y = y_label) +
    annotate("text", 
             x = -Inf, y = Inf, 
             label = permanova_text,
             hjust = -0.1, vjust = 1.5,
             size = 3,
             family = "sans",
             fontface = "italic") +
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

calculate_group_stats <- function(coords, method_name) {
  axis_cols <- grep("^Axis\\.", colnames(coords), value = TRUE)
  stats_list <- list()
  for (axis_col in axis_cols[1:3]) {
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
  cat("\n============================================================================\n")
  cat(paste0("=== Analyzing ", tax_level$level, " Level (3 Groups) ===\n"))
  cat("============================================================================\n")
  output_dir <- tax_level$output_dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(paste0("Created output directory: ", output_dir, "\n"))
  }
  cat("\n=== Loading Abundance Data ===\n")
  abundance_data <- read_excel(tax_level$file)
  cat("Data loaded successfully!\n")
  cat("Data dimensions:", dim(abundance_data), "\n")
  taxonomy_col <- colnames(abundance_data)[1]
  sample_cols <- colnames(abundance_data)[-1]
  sample_metadata <- data.frame(
    Sample = sample_cols,
    Group = case_when(
      grepl("^K", sample_cols) ~ "Control",
      grepl("^M", sample_cols) ~ "Model",
      grepl("^Z", sample_cols) ~ "HQQD",
      TRUE ~ "Unknown"
    )
  )
  sample_metadata <- sample_metadata %>% filter(Group %in% c("Control", "Model", "HQQD"))
  sample_metadata$Group <- factor(sample_metadata$Group, levels = c("Control", "Model", "HQQD"))
  samples_to_keep <- sample_metadata$Sample
  abundance_matrix <- abundance_data[, samples_to_keep, drop = FALSE] %>%
    as.data.frame() %>%
    t()
  colnames(abundance_matrix) <- abundance_data[[taxonomy_col]]
  rownames(abundance_matrix) <- samples_to_keep
  abundance_matrix <- abundance_matrix[, colSums(abundance_matrix, na.rm = TRUE) > 0]
  cat("After filtering zero-abundance taxa:", ncol(abundance_matrix), "taxa remain\n")
  cat("Sample distribution:\n")
  print(table(sample_metadata$Group))
  cat("\n=== Calculating Distance Matrices ===\n")
  dist_jaccard <- vegdist(abundance_matrix, method = "jaccard", binary = TRUE)
  dist_bray <- vegdist(abundance_matrix, method = "bray")
  cat("Distance matrices calculated successfully!\n")
  cat("\n=== Performing PCoA ===\n")
  pcoa_jaccard <- pcoa(dist_jaccard)
  pcoa_bray <- pcoa(dist_bray)
  cat("PCoA analysis completed!\n")
  jaccard_results <- extract_pcoa_results(pcoa_jaccard, sample_metadata, "Binary Jaccard")
  bray_results <- extract_pcoa_results(pcoa_bray, sample_metadata, "Bray-Curtis")
  write.csv(jaccard_results$coords, file.path(output_dir, "PCoA_Binary_Jaccard_coordinates.csv"), row.names = FALSE)
  write.csv(bray_results$coords, file.path(output_dir, "PCoA_Bray_Curtis_coordinates.csv"), row.names = FALSE)
  write.csv(jaccard_results$variance, file.path(output_dir, "PCoA_Binary_Jaccard_variance.csv"), row.names = FALSE)
  write.csv(bray_results$variance, file.path(output_dir, "PCoA_Bray_Curtis_variance.csv"), row.names = FALSE)
  cat("Coordinate and variance data saved!\n")
  cat("\n=== Performing PERMANOVA ===\n")
  permanova_jaccard <- adonis2(dist_jaccard ~ Group, data = sample_metadata, permutations = 999)
  permanova_bray <- adonis2(dist_bray ~ Group, data = sample_metadata, permutations = 999)
  permanova_results <- data.frame(
    Method = c("Binary Jaccard", "Bray-Curtis"),
    R2 = c(permanova_jaccard$R2[1], permanova_bray$R2[1]),
    F_statistic = c(permanova_jaccard$F[1], permanova_bray$F[1]),
    P_value = c(permanova_jaccard$`Pr(>F)`[1], permanova_bray$`Pr(>F)`[1])
  )
  write.csv(permanova_results, file.path(output_dir, "PERMANOVA_results.csv"), row.names = FALSE)
  cat("PERMANOVA results saved!\n")
  cat("\n=== Generating PCoA Plots ===\n")
  p_jaccard <- plot_pcoa(jaccard_results$coords, jaccard_results$variance, "Binary Jaccard", tax_level$level, permanova_results[1, ], pc_x = 1, pc_y = 2)
  p_bray <- plot_pcoa(bray_results$coords, bray_results$variance, "Bray-Curtis", tax_level$level, permanova_results[2, ], pc_x = 1, pc_y = 2)
  ggsave(file.path(output_dir, "PCoA_Binary_Jaccard_3Groups.png"), plot = p_jaccard, width = 6, height = 5, dpi = 300)
  ggsave(file.path(output_dir, "PCoA_Binary_Jaccard_3Groups.pdf"), plot = p_jaccard, width = 6, height = 5, device = cairo_pdf)
  ggsave(file.path(output_dir, "PCoA_Bray_Curtis_3Groups.png"), plot = p_bray, width = 6, height = 5, dpi = 300)
  ggsave(file.path(output_dir, "PCoA_Bray_Curtis_3Groups.pdf"), plot = p_bray, width = 6, height = 5, device = cairo_pdf)
  cat("Individual PCoA plots saved!\n")
  
  cat("\n=== Calculating Summary Statistics ===\n")
  jaccard_stats <- calculate_group_stats(jaccard_results$coords, "Binary Jaccard")
  bray_stats <- calculate_group_stats(bray_results$coords, "Bray-Curtis")
  all_stats <- bind_rows(jaccard_stats, bray_stats)
  write.csv(all_stats, file.path(output_dir, "PCoA_group_statistics_3groups.csv"), row.names = FALSE)
  cat("Statistics saved!\n")
  variance_summary <- bind_rows(
    jaccard_results$variance %>% mutate(Method = "Binary Jaccard"),
    bray_results$variance %>% mutate(Method = "Bray-Curtis")
  )
  write.csv(variance_summary, file.path(output_dir, "PCoA_variance_explained_3groups.csv"), row.names = FALSE)
  cat("\n============================================================================\n")
  cat(paste0("=== ", tax_level$level, " Level (3 Groups) Analysis Complete ===\n"))
  cat("============================================================================\n\n")
}
cat("\n============================================================================\n")
cat("=== All 3-Group Taxonomic Level Analyses Complete ===\n")
cat("============================================================================\n\n")
