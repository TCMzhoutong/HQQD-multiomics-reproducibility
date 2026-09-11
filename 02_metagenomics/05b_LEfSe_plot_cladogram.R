# ============================================================================
# LEfSe Visualization (Cladogram and Bar plots)
# Author: Generated for metagenomics analysis
# Date: 2026-01-29
# Description: Load LEfSe results and generate visualizations
#              Standard Python LEfSe (Segata style) Cladogram
# ============================================================================

# Clear environment
rm(list = ls())

# Load required libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  microeco,    # For LEfSe result plotting methods
  phyloseq,    # For cladogram helper functions used by microeco
  ggplot2,     # For plotting
  dplyr,       # For data manipulation
  ggtree,      # For tree visualization
  lifecycle,   # For compatibility helpers used by plotting methods
  writexl      # For writing Excel files
)

# ============================================================================
# Configuration Parameters
# ============================================================================

# Input/Output directory
input_dir <- "05_LEfSe_results"
output_dir <- input_dir

# Check if analysis results exist
if (!dir.exists(input_dir)) {
  stop("ERROR: LEfSe results directory not found!\n",
       "Please run 05a_LEfSe_analysis.R first!")
}

# Cladogram parameters
USE_TAXA_NUM <- 200        # Number of taxa to use in cladogram
USE_FEATURE_NUM <- 50      # Number of differential features to show
CLADE_LABEL_LEVEL <- 99    # Set very high to force ALL labels to use alphabetic mapping
                            # All labels will show as a, b, c... with full names in legend

# ============================================================================
# Step 1: Load Saved Results
# ============================================================================

cat("\n=== Step 1: Loading saved LEfSe results ===\n")

# Load LEfSe result object
lefse_result_file <- file.path(input_dir, "LEfSe_result_object.rds")
if (!file.exists(lefse_result_file)) {
  stop("ERROR: LEfSe result object not found!\n",
       "Please run 05a_LEfSe_analysis.R first!")
}

lefse_result <- readRDS(lefse_result_file)
cat("LEfSe result object loaded successfully\n")

# Load analysis context
context_file <- file.path(input_dir, "LEfSe_analysis_context.rds")
if (!file.exists(context_file)) {
  stop("ERROR: Analysis context not found!\n",
       "Please run 05a_LEfSe_analysis.R first!")
}

context <- readRDS(context_file)
group_colors <- context$group_colors
sample_metadata <- context$sample_metadata
LDA_THRESHOLD <- context$LDA_THRESHOLD

cat("Analysis context loaded successfully\n")
cat("LDA threshold:", LDA_THRESHOLD, "\n")
cat("Number of groups:", length(levels(sample_metadata$Group)), "\n\n")

# Standardize taxonomy paths used by cladogram plotting. The LEfSe statistics
# can include paths such as "Unassigned|Unassigned" in the cached abundance
# table; plot_diff_cladogram requires every path element to carry a rank prefix.
standardize_taxa_path <- function(x) {
  ranks <- c("k", "p", "c", "o", "f", "g", "s")
  
  vapply(as.character(x), function(path) {
    parts <- strsplit(path, "\\|")[[1]]
    parts <- trimws(parts)
    
    fixed <- mapply(function(part, rank) {
      if (is.na(part) || part == "") {
        return(paste0(rank, "__Unclassified"))
      }
      if (grepl("^[kpcofgs]__", part)) {
        return(part)
      }
      paste0(rank, "__", part)
    }, parts, ranks[seq_along(parts)], USE.NAMES = FALSE)
    
    paste(fixed, collapse = "|")
  }, character(1), USE.NAMES = FALSE)
}

if (!is.null(lefse_result$abund_table)) {
  rownames(lefse_result$abund_table) <- standardize_taxa_path(rownames(lefse_result$abund_table))
}

if (!is.null(lefse_result$res_abund) && "Taxa" %in% colnames(lefse_result$res_abund)) {
  lefse_result$res_abund$Taxa <- standardize_taxa_path(lefse_result$res_abund$Taxa)
}

if (!is.null(lefse_result$res_diff) && "Taxa" %in% colnames(lefse_result$res_diff)) {
  lefse_result$res_diff$Taxa <- standardize_taxa_path(lefse_result$res_diff$Taxa)
  rownames(lefse_result$res_diff) <- lefse_result$res_diff$Taxa
}

# Exclude unresolved taxonomic labels from visualizations. The full LEfSe
# result files remain unchanged; this only controls plotted features.
is_unresolved_taxon <- function(x) {
  unresolved_pattern <- "(^|[|_[:space:]])(unclassified|unassigned|unknown)([|_[:space:]]|$)"
  grepl(unresolved_pattern, as.character(x), ignore.case = TRUE)
}

if (!is.null(lefse_result$abund_table)) {
  n_before_abund_filter <- nrow(lefse_result$abund_table)
  lefse_result$abund_table <- lefse_result$abund_table[!is_unresolved_taxon(rownames(lefse_result$abund_table)), , drop = FALSE]
  cat(
    "Removed unresolved taxa from cladogram background:",
    n_before_abund_filter - nrow(lefse_result$abund_table),
    "\n"
  )
}

if (!is.null(lefse_result$res_abund) && "Taxa" %in% colnames(lefse_result$res_abund)) {
  lefse_result$res_abund <- lefse_result$res_abund %>%
    filter(!is_unresolved_taxon(Taxa))
}

if (!is.null(lefse_result$res_diff) && "Taxa" %in% colnames(lefse_result$res_diff)) {
  n_before_unresolved_filter <- nrow(lefse_result$res_diff)
  lefse_result$res_diff <- lefse_result$res_diff %>%
    filter(!is_unresolved_taxon(Taxa))
  rownames(lefse_result$res_diff) <- lefse_result$res_diff$Taxa
  cat(
    "Removed unresolved taxa from visualization:",
    n_before_unresolved_filter - nrow(lefse_result$res_diff),
    "\n"
  )
}

# Extract differential features
diff_features <- lefse_result$res_diff
cat("Total differential features:", nrow(diff_features), "\n\n")

save_plot_svg <- function(filename, plot, width, height) {
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = filename, width = width, height = height)
  } else {
    grDevices::svg(filename = filename, width = width, height = height, onefile = FALSE)
  }
  print(plot)
  grDevices::dev.off()
}

# ============================================================================
# Step 2: Generate Standard LEfSe Cladogram (Segata Style)
# ============================================================================

cat("=== Step 2: Generating LEfSe Cladogram ===\n")
cat("Style: Standard Python LEfSe (Segata style)\n")
cat("Layout: Circular hierarchical tree with sector highlighting\n")
cat("Parameters:\n")
cat("  use_taxa_num:", USE_TAXA_NUM, "\n")
cat("  use_feature_num:", USE_FEATURE_NUM, "\n")
cat("  clade_label_level:", CLADE_LABEL_LEVEL, "\n\n")

# Set group order to match color palette
group_order <- c("Control", "Model", "LVX", "HQQD")

# Generate cladogram with automatic labeling
tryCatch({
  
  # Plot 1: Full automatic cladogram with alphabetic mapping
  cat("Generating full cladogram with alphabetic mapping...\n")
  cat("Note: Outer circle shows letters (a, b, c...), full names in legend\n\n")
  
  p1 <- lefse_result$plot_diff_cladogram(
    use_taxa_num = USE_TAXA_NUM,
    use_feature_num = USE_FEATURE_NUM,
    clade_label_level = CLADE_LABEL_LEVEL,
    group_order = group_order,
    alpha = 0.2,                         # Shading transparency for sector highlighting
    select_show_labels = NULL,           # NULL = automatic alphabetic mapping
    color = group_colors,                # Use our color palette
    branch_size = 0.2,                   # Branch line size
    node_size_scale = 1,                 # Node size scaling
    node_size_offset = 1,                # Node size offset
    annotation_shape = 22,               # Legend shape (square)
    annotation_shape_size = 5,
    clade_label_size = 2,
    clade_label_size_add = 3,
    clade_label_size_log = exp(1)        # Keep default logarithmic scaling
  ) + 
  theme(
    legend.text = element_text(size = 11.5)
  )
  # Save full cladogram with larger size to show legend properly
  output_file_p1 <- file.path(output_dir, "LEfSe_cladogram_full.pdf")
  ggsave(
    filename = output_file_p1,
    plot = p1,
    width = 16,                          # Increased width for legend space
    height = 12,
    units = "in",
    dpi = 300
  )
  cat("Full cladogram saved to:", output_file_p1, "\n\n")
  
  # Also save as PNG with larger size
  output_file_p1_png <- file.path(output_dir, "LEfSe_cladogram_full.png")
  ggsave(
    filename = output_file_p1_png,
    plot = p1,
    width = 16,                          # Increased width for legend space
    height = 12,
    units = "in",
    dpi = 300
  )
  cat("Full cladogram (PNG) saved to:", output_file_p1_png, "\n\n")
  
  output_file_p1_svg <- file.path(output_dir, "LEfSe_cladogram_full.svg")
  save_plot_svg(output_file_p1_svg, p1, width = 16, height = 12)
  cat("Full cladogram (SVG) saved to:", output_file_p1_svg, "\n\n")
  
}, error = function(e) {
  cat("ERROR in cladogram generation:", conditionMessage(e), "\n")
  cat("This might be due to insufficient differential features.\n")
  cat("Try reducing USE_FEATURE_NUM or USE_TAXA_NUM parameters.\n\n")
})

# ============================================================================
# Step 3: Generate Cladogram with Selected Labels
# ============================================================================

cat("=== Step 3: Generating cladogram with selected labels ===\n")
cat("This version uses strategic label selection to minimize overlap\n\n")

tryCatch({
  
  # Get the top differential features and select diverse taxonomic levels
  top_features <- diff_features %>%
    arrange(desc(LDA)) %>%
    head(50)
  
  # Extract taxa names and filter by the terminal taxonomic rank for diversity
  selected_phyla <- grep("(^|\\|)p__[^|]+$", top_features$Taxa, value = TRUE)
  selected_classes <- grep("(^|\\|)c__[^|]+$", top_features$Taxa, value = TRUE)
  selected_orders <- grep("(^|\\|)o__[^|]+$", top_features$Taxa, value = TRUE)
  selected_families <- grep("(^|\\|)f__[^|]+$", top_features$Taxa, value = TRUE)
  
  # Take top representatives from each level
  use_labels <- c(
    head(selected_phyla, 5),
    head(selected_classes, 5),
    head(selected_orders, 5),
    head(selected_families, 5)
  )
  
  # Remove duplicates
  use_labels <- unique(use_labels)
  
  if (length(use_labels) > 0) {
    cat("Selected labels for display:\n")
    cat(paste(use_labels, collapse = "\n"), "\n\n")
    
    cat("Generating cladogram with selected labels...\n")
    
    p2 <- lefse_result$plot_diff_cladogram(
      use_taxa_num = USE_TAXA_NUM,
      use_feature_num = USE_FEATURE_NUM,
      select_show_labels = use_labels,
      group_order = group_order,
      alpha = 0.2,
      color = group_colors,
      branch_size = 0.2,
      node_size_scale = 1,
      node_size_offset = 1,
      annotation_shape = 22,
      annotation_shape_size = 5
    )
    
    # Save selected labels cladogram
    output_file_p2 <- file.path(output_dir, "LEfSe_cladogram_selected_labels.pdf")
    ggsave(
      filename = output_file_p2,
      plot = p2,
      width = 12,
      height = 12,
      units = "in",
      dpi = 300
    )
    cat("Selected labels cladogram saved to:", output_file_p2, "\n\n")
    
    # Also save as PNG
    output_file_p2_png <- file.path(output_dir, "LEfSe_cladogram_selected_labels.png")
    ggsave(
      filename = output_file_p2_png,
      plot = p2,
      width = 12,
      height = 12,
      units = "in",
      dpi = 300
    )
    cat("Selected labels cladogram (PNG) saved to:", output_file_p2_png, "\n\n")
    
    output_file_p2_svg <- file.path(output_dir, "LEfSe_cladogram_selected_labels.svg")
    save_plot_svg(output_file_p2_svg, p2, width = 12, height = 12)
    cat("Selected labels cladogram (SVG) saved to:", output_file_p2_svg, "\n\n")
  } else {
    cat("No labels selected. Skipping this plot.\n\n")
  }
  
}, error = function(e) {
  cat("ERROR in selected labels cladogram:", conditionMessage(e), "\n\n")
})

# ============================================================================
# Step 4: Generate LDA Score Bar Plot (Top 10 per group)
# ============================================================================

cat("=== Step 4: Generating LDA score bar plot (Top 10 per group) ===\n")

tryCatch({
  
  # Filter to get top 10 features per group
  top10_per_group <- diff_features %>%
    group_by(Group) %>%
    arrange(desc(LDA)) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  n_features_top10 <- nrow(top10_per_group)
  cat("Total features after top 10 per group filtering:", n_features_top10, "\n")
  
  # Calculate appropriate height
  plot_height <- max(8, n_features_top10 * 0.2)
  cat("Plot height:", round(plot_height, 1), "inches\n\n")
  
  # Temporarily update the result object with filtered data
  lefse_result_temp <- lefse_result
  lefse_result_temp$res_diff <- top10_per_group
  
  # Generate bar plot of LDA scores (top 10 per group)
  p3 <- lefse_result_temp$plot_diff_bar(
    threshold = LDA_THRESHOLD,
    group_order = group_order,
    color_values = group_colors,
    width = 0.7
  ) + theme(axis.title.x = element_blank(), legend.title = element_blank())
  
  # Save bar plot
  output_file_p3 <- file.path(output_dir, "LEfSe_LDA_barplot_top10.pdf")
  ggsave(
    filename = output_file_p3,
    plot = p3,
    width = 10,
    height = plot_height,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
  cat("LDA bar plot (top 10 per group) saved to:", output_file_p3, "\n")
  
  # Also save as PNG
  output_file_p3_png <- file.path(output_dir, "LEfSe_LDA_barplot_top10.png")
  ggsave(
    filename = output_file_p3_png,
    plot = p3,
    width = 10,
    height = plot_height,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
  cat("LDA bar plot (top 10 per group, PNG) saved to:", output_file_p3_png, "\n\n")
  
}, error = function(e) {
  cat("ERROR in bar plot generation:", conditionMessage(e), "\n\n")
})

# ============================================================================
# Step 5: Generate Taxonomic Level-Specific LDA Plots (Top 5 per group)
# ============================================================================

cat("=== Step 5: Generating taxonomic level-specific LDA plots ===\n")
cat("Plotting top 5 features per group for each taxonomic level\n\n")

# Define taxonomic levels and their prefixes
# Note: Use \\|p__ to match phylum in full taxonomic path (e.g., k__Bacteria|p__Firmicutes)
# Or use p__[^|]+$ to match features ending at phylum level
taxonomic_levels <- list(
  list(name = "Phylum", prefix = "\\|p__[^|]+$", label = "Phylum"),
  list(name = "Genus", prefix = "\\|g__[^|]+$", label = "Genus"),
  list(name = "Species", prefix = "\\|s__[^|]+$", label = "Species")
)

for (tax_level in taxonomic_levels) {
  
  tryCatch({
    
    cat("Processing", tax_level$name, "level...\n")
    
    # Filter features for this taxonomic level
    level_features <- diff_features %>%
      filter(grepl(tax_level$prefix, Taxa))
    
    if (nrow(level_features) == 0) {
      cat("  No features found at", tax_level$name, "level. Skipping.\n\n")
      next
    }
    
    cat("  Total features at", tax_level$name, "level:", nrow(level_features), "\n")
    
    # Get top 5 per group
    top5_per_group <- level_features %>%
      group_by(Group) %>%
      arrange(desc(LDA)) %>%
      slice_head(n = 5) %>%
      ungroup()
    
    n_features_top5 <- nrow(top5_per_group)
    cat("  Features after top 5 per group filtering:", n_features_top5, "\n")
    
    if (n_features_top5 == 0) {
      cat("  No features to plot. Skipping.\n\n")
      next
    }
    
    # Calculate plot height
    plot_height <- max(6, n_features_top5 * 0.25)
    
    # Temporarily update the result object
    lefse_result_temp <- lefse_result
    lefse_result_temp$res_diff <- top5_per_group
    
    # Generate bar plot
    p_level <- lefse_result_temp$plot_diff_bar(
      threshold = LDA_THRESHOLD,
      group_order = group_order,
      color_values = group_colors,
      width = 0.7
    ) + theme(axis.title.x = element_blank(), legend.title = element_blank())
    
    # Save plots
    output_file_pdf <- file.path(output_dir, paste0("LEfSe_LDA_barplot_", tax_level$name, "_top5.pdf"))
    ggsave(
      filename = output_file_pdf,
      plot = p_level,
      width = 10,
      height = plot_height,
      units = "in",
      dpi = 300,
      limitsize = FALSE
    )
    cat("  PDF saved to:", output_file_pdf, "\n")
    
    output_file_png <- file.path(output_dir, paste0("LEfSe_LDA_barplot_", tax_level$name, "_top5.png"))
    ggsave(
      filename = output_file_png,
      plot = p_level,
      width = 10,
      height = plot_height,
      units = "in",
      dpi = 300,
      limitsize = FALSE
    )
    cat("  PNG saved to:", output_file_png, "\n\n")
    
  }, error = function(e) {
    cat("  ERROR in", tax_level$name, "level plot:", conditionMessage(e), "\n\n")
  })
}

# ============================================================================
# Visualization Complete
# ============================================================================

cat("\n", rep("=", 70), "\n", sep="")
cat("LEfSe visualization complete!\n")
cat("Output files saved to:", output_dir, "\n")
cat("\nGenerated plots:\n")
cat("  1. LEfSe_cladogram_full.pdf/png - Full cladogram\n")
cat("  2. LEfSe_cladogram_selected_labels.pdf/png - Cladogram with selected labels\n")
cat("  3. LEfSe_LDA_barplot_top10.pdf/png - LDA bar plot (top 10 per group)\n")
cat("  4. LEfSe_LDA_barplot_Phylum_top5.pdf/png - Phylum level (top 5 per group)\n")
cat("  5. LEfSe_LDA_barplot_Genus_top5.pdf/png - Genus level (top 5 per group)\n")
cat("  6. LEfSe_LDA_barplot_Species_top5.pdf/png - Species level (top 5 per group)\n")
cat("\nInterpretation notes:\n")
cat("  - Circular layout represents phylogenetic relationships\n")
cat("  - Sector shading highlights enriched clades\n")
cat("  - Node size reflects relative abundance\n")
cat("  - Letters in outer circle map to full taxonomic names in legend\n")
cat("  - LDA scores indicate strength of differential abundance\n")
cat("  - Top N filtering shows most discriminative features per group\n")
cat("  - Color coding: Control=", group_colors["Control"], 
    ", Model=", group_colors["Model"],
    ", LVX=", group_colors["LVX"],
    ", HQQD=", group_colors["HQQD"], "\n")
cat(rep("=", 70), "\n\n", sep="")
