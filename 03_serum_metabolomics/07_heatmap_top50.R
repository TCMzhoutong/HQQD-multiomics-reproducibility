################################################################################
# Top 50 Differential Metabolites Heatmap with ComplexHeatmap
# Author: Data Analysis Script
# Date: 2026-02-04
# Description: Generate combined heatmaps for top 50 differential metabolites
################################################################################

# Load required packages
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)
library(stringr)
library(grid)

# Create output directory
out_dir <- "07_heatmap_top50"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Define color palettes (Academic/NPG style)
# Note: Color ranges will be dynamically determined for each comparison

# Sequential palette for -log(p) (White-Purple, academic style)
col_sequential <- colorRamp2(c(0, 2, 4), c("#F5F5F5", "#8856A7", "#4D004B"))

# Group colors (Blue and Red for two-group comparison)
group_colors <- c("Control" = "#3C5488", 
                  "Model" = "#DC0000", 
                  "LVX" = "#3C5488", 
                  "HQQD" = "#DC0000")

################################################################################
# 1. Read Annotation Data
################################################################################

cat("=== Loading Annotation Data ===\n")

# Read HMDB annotations
anno_neg <- read_csv("meta_hmdb_anno_neg.csv", show_col_types = FALSE)
anno_pos <- read_csv("meta_hmdb_anno_pos.csv", show_col_types = FALSE)

# Standardize column names
colnames(anno_neg) <- gsub("^#", "", colnames(anno_neg))
colnames(anno_pos) <- gsub("^#", "", colnames(anno_pos))

# Combine annotations
anno_all <- bind_rows(anno_neg, anno_pos)

# Create lookup table: ID -> super_class
anno_lookup <- anno_all %>%
  select(ID, super_class = `super_class(HMDB)`) %>%
  distinct()

cat(sprintf("  Loaded %d annotations\n", nrow(anno_lookup)))

################################################################################
# 2. Read Expression Data
################################################################################

cat("=== Loading Expression Data ===\n")

# Read raw expression data
exp_neg <- read_csv("metabolites_exp_neg.csv", show_col_types = FALSE)
exp_pos <- read_csv("metabolites_exp_pos.csv", show_col_types = FALSE)

# Combine expression data
exp_all <- bind_rows(exp_neg, exp_pos)

# Rename sample columns (K->Control, M->Model, Y->LVX, Z->HQQD)
colnames(exp_all) <- gsub("^K(\\d+)$", "Control\\1", colnames(exp_all))
colnames(exp_all) <- gsub("^M(\\d+)$", "Model\\1", colnames(exp_all))
colnames(exp_all) <- gsub("^Y(\\d+)$", "LVX\\1", colnames(exp_all))
colnames(exp_all) <- gsub("^Z(\\d+)$", "HQQD\\1", colnames(exp_all))

# Extract sample columns (exclude ID, name, QC columns)
sample_cols <- grep("^(Control|Model|LVX|HQQD)\\d+$", colnames(exp_all), value = TRUE)

cat(sprintf("  Loaded expression data for %d metabolites\n", nrow(exp_all)))
cat(sprintf("  Found %d sample columns\n", length(sample_cols)))

################################################################################
# 3. Process Each Comparison
################################################################################

# Get differential metabolite files from revised statistics
diff_dir <- "06_merge_diff_revised"
diff_files <- list.files(diff_dir, pattern = "^06_merged_.*_diff\\.csv$", full.names = FALSE)

if (length(diff_files) == 0) {
  stop("No differential metabolite files found in ", diff_dir, "/")
}

cat(sprintf("\n=== Processing %d Comparisons ===\n\n", length(diff_files)))

for (f in diff_files) {
  
  # Extract comparison name
  comp_name <- str_remove(f, "^06_merged_")
  comp_name <- str_remove(comp_name, "_diff\\.csv$")
  
  cat(sprintf("Processing: %s\n", comp_name))
  
  # Read differential metabolites
  diff_data <- read_csv(file.path(diff_dir, f), show_col_types = FALSE)
  
  if (nrow(diff_data) == 0) {
    cat("  Skipping: No differential metabolites found\n\n")
    next
  }
  
  # Extract group names from comparison (e.g., "Control_vs_HQQD" -> c("Control", "HQQD"))
  groups <- str_split(comp_name, "_vs_")[[1]]
  group1 <- groups[1]
  group2 <- groups[2]
  
  cat(sprintf("  Groups: %s vs %s\n", group1, group2))
  
  # Sort by absolute Log2FC and select Top 50
  diff_data <- diff_data %>%
    mutate(abs_Log2FC = abs(Log2FC)) %>%
    arrange(desc(abs_Log2FC)) %>%
    head(50)
  
  # Then sort by original Log2FC value for display order (descending: high to low)
  diff_data <- diff_data %>%
    arrange(desc(Log2FC))
  
  cat(sprintf("  Selected top %d metabolites\n", nrow(diff_data)))
  
  # Get IDs
  top_ids <- diff_data$ID
  
  # Calculate -log10(p)
  diff_data <- diff_data %>%
    mutate(neg_log_p = -log10(pmax(p_value, 1e-10))) # Avoid log(0)
  
  # Determine color range for log2(FC) based on actual data
  fc_range <- range(diff_data$Log2FC, na.rm = TRUE)
  fc_max <- max(abs(fc_range))
  fc_max <- ceiling(fc_max)  # Round up to integer
  
  # Create color function for log2(FC)
  col_fc <- colorRamp2(c(-fc_max, 0, fc_max), 
                       c("#3C5488", "white", "#DC0000"))
  
  cat(sprintf("  log2(FC) range: [%.2f, %.2f], color range: [-%d, %d]\n", 
              fc_range[1], fc_range[2], fc_max, fc_max))
  
  # Cap extreme values for visualization
  diff_data$Log2FC_capped <- pmax(pmin(diff_data$Log2FC, fc_max), -fc_max)
  diff_data$neg_log_p_capped <- pmin(diff_data$neg_log_p, 4)
  
  # Join with ontology
  diff_data <- diff_data %>%
    left_join(anno_lookup, by = "ID")
  
  # Handle missing super_class
  diff_data$super_class[is.na(diff_data$super_class) | diff_data$super_class == ""] <- "Unknown"
  
  # Create row labels with truncation
  # Format: "Name... (ID)" only if name is too long; otherwise just "Name"
  max_name_length <- 40
  diff_data$row_label <- sapply(1:nrow(diff_data), function(i) {
    name <- diff_data$Name[i]
    id <- diff_data$ID[i]
    if (is.na(name) || name == "") {
      return(id)
    }
    if (nchar(name) > max_name_length) {
      # Truncate and add ID for long names
      name_truncated <- paste0(substr(name, 1, max_name_length), "...")
      return(paste0(name_truncated, " (", id, ")"))
    }
    # Short names: no ID needed
    return(name)
  })
  
  # Extract expression matrix for these metabolites
  exp_subset <- exp_all %>%
    filter(ID %in% top_ids) %>%
    select(ID, all_of(sample_cols))
  
  # Match order with diff_data
  exp_subset <- exp_subset[match(top_ids, exp_subset$ID), ]
  
  # Create expression matrix (rows = metabolites, columns = samples)
  expr_mat <- as.matrix(exp_subset[, -1])
  rownames(expr_mat) <- diff_data$row_label
  
  # Filter columns to only include samples from the two groups
  sample_pattern <- paste0("^(", group1, "|", group2, ")\\d+$")
  selected_samples <- grep(sample_pattern, colnames(expr_mat), value = TRUE)
  
  if (length(selected_samples) == 0) {
    cat("  Warning: No matching samples found. Skipping.\n\n")
    next
  }
  
  expr_mat <- expr_mat[, selected_samples, drop = FALSE]
  
  # Log transformation and Z-score normalization
  expr_mat_log <- log2(expr_mat + 1)
  
  # Z-score normalization (row-wise)
  expr_mat_scaled <- t(scale(t(expr_mat_log)))
  
  # Replace NA with 0 (in case of zero variance metabolites)
  expr_mat_scaled[is.na(expr_mat_scaled)] <- 0
  
  # Determine color range for Scaled Intensity based on actual data
  # Use 1st and 99th percentile to handle outliers
  intensity_range <- quantile(expr_mat_scaled, probs = c(0.01, 0.99), na.rm = TRUE)
  intensity_max <- max(abs(intensity_range))
  intensity_max <- ceiling(intensity_max)  # Round up to integer
  
  # Cap extreme values for visualization
  expr_mat_scaled[expr_mat_scaled > intensity_max] <- intensity_max
  expr_mat_scaled[expr_mat_scaled < -intensity_max] <- -intensity_max
  
  # Create color function for Scaled Intensity
  col_intensity <- colorRamp2(c(-intensity_max, 0, intensity_max), 
                              c("#3C5488", "white", "#DC0000"))
  
  cat(sprintf("  Scaled Intensity range: [-%d, %d]\n", intensity_max, intensity_max))
  
  cat(sprintf("  Matrix dimensions: %d metabolites x %d samples\n", 
              nrow(expr_mat_scaled), ncol(expr_mat_scaled)))
  
  ############################################################################
  # 4. Create Sample Group Annotation
  ############################################################################
  
  # Determine group for each sample
  sample_groups <- sapply(colnames(expr_mat_scaled), function(s) {
    if (grepl(paste0("^", group1), s)) return(group1)
    if (grepl(paste0("^", group2), s)) return(group2)
    return("Other")
  })
  
  # Define colors for the two groups (Blue vs Red)
  pair_colors <- c("#3C5488", "#DC0000")
  names(pair_colors) <- c(group1, group2)
  
  # Create column annotation (Top)
  col_anno <- HeatmapAnnotation(
    Group = sample_groups,
    col = list(Group = pair_colors),
    annotation_name_side = "right",  # Label on right side
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = 9),  # Smaller than default but > sample names
    show_legend = TRUE
  )
  
  ############################################################################
  # 5. Create Row Annotations (Right side)
  ############################################################################
  
  # Ontology colors (discrete) - Academic color palette
  unique_classes <- unique(diff_data$super_class)
  n_classes <- length(unique_classes)
  
  # Use NPG-style academic color palette for ontology
  npg_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", 
                  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                  "#7E6148", "#B09C85", "#00A1D5", "#79AF97")
  
  if (n_classes <= length(npg_colors)) {
    ontology_colors <- setNames(npg_colors[1:n_classes], unique_classes)
  } else {
    # For more classes, extend with additional colors
    ontology_colors <- setNames(
      c(npg_colors, rainbow(n_classes - length(npg_colors)))[1:n_classes],
      unique_classes
    )
  }
  
  # Row annotation data frame
  row_anno_df <- data.frame(
    Ontology = diff_data$super_class,
    Log2FC = diff_data$Log2FC_capped,
    neg_log_p = diff_data$neg_log_p_capped,
    VIP = as.numeric(diff_data$VIP),
    row.names = diff_data$row_label
  )
  vip_available_for_heatmap <- any(is.finite(row_anno_df$VIP))
  
  # Create right row annotations - Split into two parts for dendrogram positioning
  # First annotation: Ontology only (will be adjacent to dendrogram)
  row_anno_1 <- rowAnnotation(
    Ontology = row_anno_df$Ontology,
    col = list(Ontology = ontology_colors),
    annotation_name_rot = 90,
    annotation_name_side = "bottom",
    annotation_name_gp = gpar(fontsize = 9),
    show_legend = TRUE,
    width = unit(8, "mm")
  )
  
  # Second annotation: log2(FC), -log(p), and VIP only when a validated VIP exists.
  if (vip_available_for_heatmap) {
    row_anno_2 <- rowAnnotation(
      `log2(FC)` = anno_simple(
        row_anno_df$Log2FC,
        col = col_fc,
        pch = NA,
        width = unit(8, "mm")
      ),
      
      `-log(p)` = anno_simple(
        row_anno_df$neg_log_p,
        col = col_sequential,
        pch = NA,
        width = unit(8, "mm")
      ),
      
      VIP = anno_barplot(
        row_anno_df$VIP,
        bar_width = 0.6,
        gp = gpar(fill = "#2166AC", col = NA),
        axis = TRUE,
        axis_param = list(
          side = "bottom", 
          labels_rot = 0,
          gp = gpar(fontsize = 7)
        ),
        border = FALSE,
        width = unit(20, "mm")
      ),
      
      annotation_name_rot = c(90, 90, 0),
      annotation_name_side = "bottom",
      annotation_name_gp = gpar(fontsize = 9),
      show_legend = FALSE,
      gap = unit(2, "mm")
    )
  } else {
    row_anno_2 <- rowAnnotation(
      `log2(FC)` = anno_simple(
        row_anno_df$Log2FC,
        col = col_fc,
        pch = NA,
        width = unit(8, "mm")
      ),
      
      `-log(p)` = anno_simple(
        row_anno_df$neg_log_p,
        col = col_sequential,
        pch = NA,
        width = unit(8, "mm")
      ),
      
      annotation_name_rot = c(90, 90),
      annotation_name_side = "bottom",
      annotation_name_gp = gpar(fontsize = 9),
      show_legend = FALSE,
      gap = unit(2, "mm")
    )
  }
  
  ############################################################################
  # 6. Create Main Heatmap
  ############################################################################
  
  # Calculate cell size for square cells
  n_samples <- ncol(expr_mat_scaled)
  n_metabolites <- nrow(expr_mat_scaled)
  cell_size <- unit(5, "mm")  # Square cell size
  
  ht <- Heatmap(
    expr_mat_scaled,
    name = "Scaled Intensity",
    col = col_intensity,  # Dynamic color range
    
    # Column settings
    top_annotation = col_anno,
    column_names_rot = 90,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 8),
    cluster_columns = FALSE,
    column_split = sample_groups,
    column_title = NULL,
    
    # Row settings
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    cluster_rows = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    row_dend_side = "right",
    row_dend_width = unit(15, "mm"),
    show_row_dend = TRUE,
    
    # Right annotations - only Ontology (dendrogram will be between heatmap and this)
    right_annotation = row_anno_1,
    
    # Dimensions - Square cells
    width = n_samples * cell_size,
    height = n_metabolites * cell_size,
    
    # Heatmap body
    rect_gp = gpar(col = "black", lwd = 0.3),
    
    # Legend
    heatmap_legend_param = list(
      title = "Scaled Intensity",
      direction = "vertical",
      legend_height = unit(40, "mm")
    )
  ) + row_anno_2  # Add second annotation after Ontology
  
  ############################################################################
  # 7. Create Custom Legend for Annotations
  ############################################################################
  
  # Legend for log2(FC) - use dynamic range
  lgd_fc <- Legend(
    title = "log2(FC)",
    col_fun = col_fc,
    direction = "vertical",
    legend_height = unit(30, "mm")
  )
  
  # Legend for -log(p)
  lgd_p <- Legend(
    title = "-log(p)",
    col_fun = col_sequential,
    direction = "vertical",
    legend_height = unit(30, "mm")
  )
  
  ############################################################################
  # 8. Draw and Save
  ############################################################################
  
  # Combine all legends
  lgd_list <- list(lgd_fc, lgd_p)
  
  # Calculate appropriate output dimensions
  # Base dimensions + extra space for annotations, labels, legends
  heatmap_width_mm <- n_samples * 5 + 150  # cell width + row labels + right annotations + legends
  heatmap_height_mm <- n_metabolites * 5 + 80  # cell height + column labels + top annotation
  
  # Convert to inches (1 inch = 25.4 mm)
  pdf_width <- max(heatmap_width_mm / 25.4, 14)
  pdf_height <- max(heatmap_height_mm / 25.4, 12)
  
  # Save as PDF
  pdf_file <- file.path(out_dir, paste0("07_heatmap_", comp_name, ".pdf"))
  pdf(pdf_file, width = pdf_width, height = pdf_height)
  
  draw(ht, 
       annotation_legend_list = lgd_list,
       heatmap_legend_side = "right",
       annotation_legend_side = "right",
       merge_legend = TRUE,
       padding = unit(c(10, 10, 10, 10), "mm"))
  
  dev.off()
  
  # Save as PNG (higher resolution for large heatmaps)
  png_file <- file.path(out_dir, paste0("07_heatmap_", comp_name, ".png"))
  png(png_file, width = pdf_width * 300, height = pdf_height * 300, res = 300)
  
  draw(ht, 
       annotation_legend_list = lgd_list,
       heatmap_legend_side = "right",
       annotation_legend_side = "right",
       merge_legend = TRUE,
       padding = unit(c(10, 10, 10, 10), "mm"))
  
  dev.off()
  
  cat(sprintf("  Saved: %s\n", pdf_file))
  cat(sprintf("  Saved: %s\n\n", png_file))
}

cat("=== All Heatmaps Complete ===\n")
cat(sprintf("Results saved to: %s/\n", out_dir))
