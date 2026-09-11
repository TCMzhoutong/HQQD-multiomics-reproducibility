# ============================================================================
# Taxonomy Composition Stacked Bar Charts
# Author: Generated for metagenomics analysis
# Date: 2026-01-28
# Description: Stacked bar charts showing taxonomic composition at different
#              levels (Phylum, Genus, Species) for individual samples and groups
#              with top 10 taxa, Others, and connecting lines between bars
# ============================================================================

# Clear environment
rm(list = ls())

# Load required libraries
library(readxl)      # For reading Excel files
library(ggplot2)     # For visualization
library(dplyr)       # For data manipulation
library(tidyr)       # For data reshaping
library(RColorBrewer)  # For color palettes
library(ggalluvial)  # For alluvial/flow diagrams
library(cowplot)     # For fixed-width legend layout

# Define Morandi color palette (soft, low-saturation colors for scientific publications)
# Expanded palette with 20 distinct colors to ensure no duplication
morandi_colors <- c(
  "#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2", "#BEB8DC",
  "#A8E6CF", "#FFD3B6", "#FFAAA5", "#B0C4DE", "#D4A5A5",
  "#C9E4CA", "#F7C6C7", "#B5EAD7", "#C7CEEA", "#FFDAC1",
  "#E0BBE4", "#957DAD", "#D291BC", "#FEC8D8", "#C3B1E1"
)

# Define sample group mapping
sample_groups <- data.frame(
  Sample = c("C1", "C2", "C3", "C4", "C5", "C6",
             "M1", "M2", "M3", "M4", "M5", "M6",
             "L1", "L2", "L3", "L4", "L5", "L6",
             "H1", "H2", "H3", "H4", "H5", "H6"),
  Group = rep(c("Control", "Model", "LVX", "HQQD"), each = 6),
  stringsAsFactors = FALSE
)

# Define taxonomic levels to analyze
taxonomic_levels <- list(
  list(file = "taxonomy.phylum.relabundance.xlsx", 
       level = "Phylum", 
       output_dir = "04_composition_Phylum"),
  list(file = "taxonomy.genus.relabundance.xlsx", 
       level = "Genus", 
       output_dir = "04_composition_Genus"),
  list(file = "taxonomy.species.relabundance.xlsx", 
       level = "Species", 
       output_dir = "04_composition_Species")
)

# ============================================================================
# Define Helper Functions
# ============================================================================

# Function to rename samples from original naming to standard naming
rename_samples <- function(sample_name) {
  # K1-K6 -> C1-C6 (Control)
  # Y1-Y6 -> L1-L6 (LVX)
  # Z1-Z6 -> H1-H6 (HQQD)
  # M1-M6 remains M1-M6 (Model)
  sample_name <- gsub("^K(\\d+)$", "C\\1", sample_name)
  sample_name <- gsub("^Y(\\d+)$", "L\\1", sample_name)
  sample_name <- gsub("^Z(\\d+)$", "H\\1", sample_name)
  return(sample_name)
}

# Function to prepare data for stacked bar chart
# SCI Standard: Use mean abundance across all samples to rank taxa
prepare_taxonomy_data <- function(data, level_name) {
  # Convert to long format
  data_long <- data %>%
    pivot_longer(cols = -1, names_to = "Sample", values_to = "Abundance") %>%
    rename(Taxon = 1)
  
  # Rename samples
  data_long$Sample <- sapply(data_long$Sample, rename_samples)
  
  # Add group information
  data_long <- data_long %>%
    left_join(sample_groups, by = "Sample")
  
  # Identify Unassigned and Unclassified taxa separately
  unassigned_pattern <- "Unassigned|unassigned|Unknown|unknown"
  unclassified_pattern <- "Unclassified|unclassified|uncultured|Uncultured"
  
  # Separate into Unassigned and Unclassified
  unassigned_taxa <- data_long %>%
    filter(grepl(unassigned_pattern, Taxon)) %>%
    pull(Taxon) %>%
    unique()
  
  unclassified_taxa <- data_long %>%
    filter(grepl(unclassified_pattern, Taxon) & !Taxon %in% unassigned_taxa) %>%
    pull(Taxon) %>%
    unique()
  
  all_special_taxa <- c(unassigned_taxa, unclassified_taxa)
  
  # Calculate mean abundance for each taxon
  # EXCLUDING Unassigned/Unclassified from Top10 calculation
  taxon_means <- data_long %>%
    filter(!Taxon %in% all_special_taxa) %>%
    group_by(Taxon) %>%
    summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
    arrange(desc(MeanAbundance))
  
  # Get top 10 taxa (from normal taxa only)
  top10_taxa <- taxon_means$Taxon[1:min(10, nrow(taxon_means))]
  
  cat(paste0("  Identified ", length(unassigned_taxa), " Unassigned taxa\n"))
  cat(paste0("  Identified ", length(unclassified_taxa), " Unclassified taxa\n"))
  cat(paste0("  Selected Top 10 taxa (excluding special categories)\n"))
  
  # Classify taxa with clear hierarchy
  # Order: Top1-10, Others, Unassigned, Unclassified
  data_long <- data_long %>%
    mutate(TaxonCategory = case_when(
      Taxon %in% top10_taxa ~ Taxon,
      Taxon %in% unassigned_taxa ~ "Unassigned",
      Taxon %in% unclassified_taxa ~ "Unclassified",
      TRUE ~ "Others"
    ))
  
  return(list(data = data_long, top10 = top10_taxa))
}

# Function to create sample-level stacked bar chart
plot_sample_composition <- function(data_long, top10_taxa, level_name, output_dir) {
  cat(paste0("Creating sample-level alluvial chart for ", level_name, "...\n"))
  
  # Aggregate data
  plot_data <- data_long %>%
    group_by(Sample, Group, TaxonCategory) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  # Check if abundances sum to 100% (or close to 1.0)
  abundance_check <- plot_data %>%
    group_by(Sample) %>%
    summarise(Total = sum(Abundance)) %>%
    filter(abs(Total - 1.0) > 0.01)
  
  if (nrow(abundance_check) > 0) {
    cat("  Warning: Some samples do not sum to 100%:\n")
    print(abundance_check)
  }
  
  # Set consistent factor order: Top1-10 (ascending), Others, Unassigned, Unclassified
  # Bottom to top in stacked bar: Top1, Top2, ..., Top10, Others, Unassigned, Unclassified
  # Reverse the order for correct stacking direction
  taxon_order <- rev(c(top10_taxa, "Others", "Unassigned", "Unclassified"))
  plot_data$TaxonCategory <- factor(plot_data$TaxonCategory, levels = taxon_order)
  
  # Order samples by group
  sample_order <- sample_groups$Sample
  plot_data$Sample <- factor(plot_data$Sample, levels = sample_order)
  
  # Assign colors - consistent order (reverse back for color assignment)
  original_order <- c(top10_taxa, "Others", "Unassigned", "Unclassified")
  n_taxa <- length(original_order)
  if (n_taxa <= length(morandi_colors)) {
    taxa_colors <- setNames(morandi_colors[1:n_taxa], original_order)
  } else {
    taxa_colors <- setNames(colorRampPalette(morandi_colors)(n_taxa), original_order)
  }
  # Make Others, Unassigned, and Unclassified grey colors
  if ("Others" %in% names(taxa_colors)) taxa_colors["Others"] <- "#999999"
  if ("Unassigned" %in% names(taxa_colors)) taxa_colors["Unassigned"] <- "#BBBBBB"
  if ("Unclassified" %in% names(taxa_colors)) taxa_colors["Unclassified"] <- "#DDDDDD"
  
  # Create alluvial plot
  p <- ggplot(plot_data,
              aes(x = Sample, y = Abundance * 100, 
                  alluvium = TaxonCategory, stratum = TaxonCategory, fill = TaxonCategory)) +
    geom_alluvium(alpha = 0.7, decreasing = NA) +
    geom_stratum(width = 0.8, decreasing = NA, color = "#000000", linewidth = 0.2) +
    scale_fill_manual(values = taxa_colors, name = NULL) +
    labs(title = paste0("Taxonomic Composition at ", level_name, " Level (By Sample)"),
         x = NULL,
         y = "Relative Abundance (%)") +
    theme_bw(base_size = 10, base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
      axis.text.y = element_text(size = 9, face = "bold"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      legend.spacing.y = unit(0.1, "cm"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 1))
  
  # Save plot
  ggsave(file.path(output_dir, paste0(level_name, "_composition_by_sample.png")), 
         plot = p, width = 14, height = 6, dpi = 300)
  ggsave(file.path(output_dir, paste0(level_name, "_composition_by_sample.pdf")), 
         plot = p, width = 14, height = 6, device = cairo_pdf)
  
  cat(paste0("Sample-level plot saved to ", output_dir, "\n"))
  
  return(p)
}

# Function to create group-level stacked bar chart
plot_group_composition <- function(data_long, top10_taxa, level_name, output_dir) {
  cat(paste0("Creating group-level alluvial chart for ", level_name, "...\n"))
  
  # Aggregate data by group - calculate mean for each taxon within each group
  plot_data <- data_long %>%
    group_by(Sample, Group, TaxonCategory) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Group, TaxonCategory) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")
  
  # Check if abundances sum to 100% (or close to 1.0)
  abundance_check <- plot_data %>%
    group_by(Group) %>%
    summarise(Total = sum(Abundance)) %>%
    mutate(Percent = Total * 100)
  
  cat("  Group abundance totals:\n")
  print(abundance_check)
  
  # Set consistent factor order: Top1-10 (ascending), Others, Unassigned, Unclassified
  # Bottom to top in stacked bar: Top1, Top2, ..., Top10, Others, Unassigned, Unclassified
  # Reverse the order for correct stacking direction
  taxon_order <- rev(c(top10_taxa, "Others", "Unassigned", "Unclassified"))
  plot_data$TaxonCategory <- factor(plot_data$TaxonCategory, levels = taxon_order)
  
  # Order groups
  group_order <- c("Control", "Model", "LVX", "HQQD")
  plot_data$Group <- factor(plot_data$Group, levels = group_order)
  
  # Assign colors - consistent order (reverse back for color assignment)
  original_order <- c(top10_taxa, "Others", "Unassigned", "Unclassified")
  n_taxa <- length(original_order)
  if (n_taxa <= length(morandi_colors)) {
    taxa_colors <- setNames(morandi_colors[1:n_taxa], original_order)
  } else {
    taxa_colors <- setNames(colorRampPalette(morandi_colors)(n_taxa), original_order)
  }
  # Make Others, Unassigned, and Unclassified grey colors
  if ("Others" %in% names(taxa_colors)) taxa_colors["Others"] <- "#999999"
  if ("Unassigned" %in% names(taxa_colors)) taxa_colors["Unassigned"] <- "#BBBBBB"
  if ("Unclassified" %in% names(taxa_colors)) taxa_colors["Unclassified"] <- "#DDDDDD"
  
  wrap_legend_label <- function(x) {
    vapply(x, function(label) paste(strwrap(label, width = 18), collapse = "\n"), character(1))
  }

  # Create alluvial plot
  p <- ggplot(plot_data,
              aes(x = Group, y = Abundance * 100, 
                  alluvium = TaxonCategory, stratum = TaxonCategory, fill = TaxonCategory)) +
    geom_alluvium(alpha = 0.7, decreasing = NA) +
    geom_stratum(width = 0.55, decreasing = NA, color = "#000000", linewidth = 0.2) +
    scale_fill_manual(values = taxa_colors, name = NULL, labels = wrap_legend_label) +
    labs(title = paste0("Taxonomic Composition at ", level_name, " Level"),
         x = NULL,
         y = "Relative Abundance (%)") +
    theme_bw(base_size = 11, base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      legend.key.size = unit(0.5, "cm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      legend.spacing.y = unit(0.1, "cm"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 1))

  p_main <- p +
    theme(
      legend.position = "none",
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )

  legend_levels <- levels(plot_data$TaxonCategory)
  legend_data <- data.frame(
    TaxonCategory = factor(legend_levels, levels = legend_levels),
    Label = wrap_legend_label(legend_levels),
    y = seq_along(legend_levels),
    stringsAsFactors = FALSE
  )
  legend_data$fill <- taxa_colors[as.character(legend_data$TaxonCategory)]

  p_legend <- ggplot(legend_data, aes(y = y)) +
    geom_tile(aes(x = 0.08, fill = TaxonCategory), width = 0.16, height = 0.9,
              color = "black", linewidth = 0.2) +
    geom_text(aes(x = 0.22, label = Label), hjust = 0, vjust = 0.5,
              size = 11 / ggplot2::.pt, lineheight = 0.9, family = "sans") +
    scale_fill_manual(values = taxa_colors, guide = "none") +
    scale_y_reverse(limits = c(length(legend_levels) + 0.5, 0.5)) +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    theme_void(base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(0, 0, 0, 0)
    )

  p_final <- cowplot::plot_grid(
    p_main,
    p_legend,
    ncol = 2,
    rel_widths = c(0.64, 0.36),
    align = "h",
    axis = "tb"
  ) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  # Save plot
  ggsave(file.path(output_dir, paste0(level_name, "_composition_by_group.png")), 
         plot = p_final, width = 6.4, height = 5.2, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0(level_name, "_composition_by_group.pdf")), 
         plot = p_final, width = 6.4, height = 5.2, device = cairo_pdf, bg = "white")
  ggsave(file.path(output_dir, paste0(level_name, "_composition_by_group.svg")),
         plot = p_final, width = 6.4, height = 5.2, bg = "white")
  
  cat(paste0("Group-level plot saved to ", output_dir, "\n"))
  
  return(p_final)
}

# ============================================================================
# Main Analysis Loop
# ============================================================================

for (tax_level in taxonomic_levels) {
  cat("\n============================================================================\n")
  cat(paste0("=== Analyzing ", tax_level$level, " Level (Composition Bar Charts) ===\n"))
  cat("============================================================================\n")
  
  output_dir <- tax_level$output_dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(paste0("Created output directory: ", output_dir, "\n"))
  }
  
  cat("\n=== Loading Abundance Data ===\n")
  data <- read_excel(tax_level$file)
  cat(paste0("Loaded data: ", nrow(data), " taxa, ", ncol(data) - 1, " samples\n"))
  
  cat("\n=== Preparing Data ===\n")
  prepared <- prepare_taxonomy_data(data, tax_level$level)
  data_long <- prepared$data
  top10_taxa <- prepared$top10
  
  cat(paste0("Top 10 taxa identified:\n"))
  for (i in 1:length(top10_taxa)) {
    cat(paste0("  ", i, ". ", top10_taxa[i], "\n"))
  }
  
  cat("\n=== Creating Plots ===\n")
  # Sample-level plot
  p_sample <- plot_sample_composition(data_long, top10_taxa, tax_level$level, output_dir)
  
  # Group-level plot
  p_group <- plot_group_composition(data_long, top10_taxa, tax_level$level, output_dir)
  
  cat(paste0("\n", tax_level$level, " level analysis complete!\n"))
}

cat("\n============================================================================\n")
cat("=== All Taxonomic Composition Analyses Complete! ===\n")
cat("============================================================================\n")
cat("\nOutput directories:\n")
for (tax_level in taxonomic_levels) {
  cat(paste0("- ", tax_level$output_dir, "/\n"))
  cat(paste0("  - ", tax_level$level, "_composition_by_sample.png/pdf\n"))
  cat(paste0("  - ", tax_level$level, "_composition_by_group.png/pdf\n"))
}
