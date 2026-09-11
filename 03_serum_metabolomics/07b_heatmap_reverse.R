################################################################################
# Reverse Metabolites Heatmap (Control_vs_Model & Model_vs_HQQD)
# Description: Old ComplexHeatmap drawing style, using revised reverse candidates.
################################################################################

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)
library(grid)

open_svg_device <- function(file, width, height) {
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = file, width = width, height = height)
  } else {
    svg(filename = file, width = width, height = height, onefile = FALSE)
  }
}

out_dir <- "07b_heatmap_reverse"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

col_sequential <- colorRamp2(c(0, 2, 4), c("#F5F5F5", "#8856A7", "#4D004B"))

cat("=== Loading Annotation Data ===\n")

anno_neg <- read_csv("meta_hmdb_anno_neg.csv", show_col_types = FALSE)
anno_pos <- read_csv("meta_hmdb_anno_pos.csv", show_col_types = FALSE)

colnames(anno_neg) <- gsub("^#", "", colnames(anno_neg))
colnames(anno_pos) <- gsub("^#", "", colnames(anno_pos))

anno_all <- bind_rows(anno_neg, anno_pos)

anno_lookup <- anno_all %>%
  select(ID, super_class = `super_class(HMDB)`) %>%
  distinct()

cat(sprintf("  Loaded %d annotations\n", nrow(anno_lookup)))

cat("=== Loading Expression Data ===\n")

exp_neg <- read_csv("metabolites_exp_neg.csv", show_col_types = FALSE)
exp_pos <- read_csv("metabolites_exp_pos.csv", show_col_types = FALSE)

exp_all <- bind_rows(exp_neg, exp_pos)

colnames(exp_all) <- gsub("^K(\\d+)$", "Control\\1", colnames(exp_all))
colnames(exp_all) <- gsub("^M(\\d+)$", "Model\\1", colnames(exp_all))
colnames(exp_all) <- gsub("^Y(\\d+)$", "LVX\\1", colnames(exp_all))
colnames(exp_all) <- gsub("^Z(\\d+)$", "HQQD\\1", colnames(exp_all))

sample_cols <- grep("^(Control|Model|LVX|HQQD)\\d+$", colnames(exp_all), value = TRUE)

cat(sprintf("  Loaded expression data for %d metabolites\n", nrow(exp_all)))
cat(sprintf("  Found %d sample columns\n", length(sample_cols)))

reverse_file <- "06_HQQD_reverse_revised/HQQD_reverse_candidates.csv"
if (!file.exists(reverse_file)) {
  stop("Reverse file not found: ", reverse_file)
}

reverse_data <- read_csv(reverse_file, show_col_types = FALSE)
if (nrow(reverse_data) == 0) {
  stop("HQQD_reverse_candidates.csv is empty. Cannot plot heatmaps.")
}

reverse_rank <- reverse_data %>%
  mutate(
    CM_Log2FC = as.numeric(CM_Log2FC),
    abs_CM_Log2FC = abs(CM_Log2FC)
  ) %>%
  arrange(desc(abs_CM_Log2FC), desc(CM_Log2FC)) %>%
  head(50)

shared_keys <- paste(reverse_rank$ID, reverse_rank$Ion_Mode, sep = "||")
cat(sprintf("=== Using %d shared metabolites (same row order in both heatmaps) ===\n", length(shared_keys)))

plot_configs <- list(
  list(comp_name = "Control_vs_Model", group1 = "Control", group2 = "Model", prefix = "CM"),
  list(comp_name = "Model_vs_HQQD", group1 = "Model", group2 = "HQQD", prefix = "MH")
)

for (cfg in plot_configs) {
  comp_name <- cfg$comp_name
  group1 <- cfg$group1
  group2 <- cfg$group2
  prefix <- cfg$prefix

  cat(sprintf("\nProcessing: %s\n", comp_name))

  diff_data <- reverse_data %>%
    mutate(
      key = paste(ID, Ion_Mode, sep = "||"),
      Name = ifelse(!is.na(CM_Name) & CM_Name != "", CM_Name,
                    ifelse(!is.na(MH_Name) & MH_Name != "", MH_Name, ID)),
      Log2FC = as.numeric(.data[[paste0(prefix, "_Log2FC")]]),
      p_value = as.numeric(.data[[paste0(prefix, "_p_value")]]),
      VIP = as.numeric(.data[[paste0(prefix, "_VIP")]])
    ) %>%
    filter(key %in% shared_keys) %>%
    mutate(key = factor(key, levels = shared_keys)) %>%
    arrange(key)

  if (nrow(diff_data) == 0) {
    cat("  Skipping: No metabolites after shared-order filtering\n")
    next
  }

  cat(sprintf("  Groups: %s vs %s\n", group1, group2))
  cat(sprintf("  Metabolites used: %d\n", nrow(diff_data)))

  diff_data <- diff_data %>%
    mutate(neg_log_p = -log10(pmax(p_value, 1e-10)))

  fc_range <- range(diff_data$Log2FC, na.rm = TRUE)
  fc_max <- max(abs(fc_range))
  fc_max <- ceiling(fc_max)
  if (!is.finite(fc_max) || fc_max == 0) fc_max <- 1

  col_fc <- colorRamp2(c(-fc_max, 0, fc_max), c("#3C5488", "white", "#DC0000"))

  diff_data$Log2FC_capped <- pmax(pmin(diff_data$Log2FC, fc_max), -fc_max)
  diff_data$neg_log_p_capped <- pmin(diff_data$neg_log_p, 4)

  diff_data <- diff_data %>% left_join(anno_lookup, by = "ID")
  diff_data$super_class[is.na(diff_data$super_class) | diff_data$super_class == ""] <- "Unknown"

  max_name_length <- 40
  diff_data$row_label <- sapply(seq_len(nrow(diff_data)), function(i) {
    name <- diff_data$Name[i]
    id <- diff_data$ID[i]
    if (is.na(name) || name == "") return(id)
    if (nchar(name) > max_name_length) {
      name_truncated <- paste0(substr(name, 1, max_name_length), "...")
      return(paste0(name_truncated, " (", id, ")"))
    }
    return(name)
  })
  diff_data$row_label <- make.unique(diff_data$row_label)

  exp_subset <- exp_all %>%
    filter(ID %in% diff_data$ID) %>%
    select(ID, all_of(sample_cols))

  exp_subset <- exp_subset[match(diff_data$ID, exp_subset$ID), ]
  expr_mat <- as.matrix(exp_subset[, -1, drop = FALSE])
  rownames(expr_mat) <- diff_data$row_label

  sample_pattern <- paste0("^(", group1, "|", group2, ")\\d+$")
  selected_samples <- grep(sample_pattern, colnames(expr_mat), value = TRUE)
  if (length(selected_samples) == 0) {
    cat("  Warning: No matching samples found. Skipping.\n")
    next
  }
  expr_mat <- expr_mat[, selected_samples, drop = FALSE]

  expr_mat_log <- log2(expr_mat + 1)
  expr_mat_scaled <- t(scale(t(expr_mat_log)))
  expr_mat_scaled[is.na(expr_mat_scaled)] <- 0

  intensity_range <- quantile(expr_mat_scaled, probs = c(0.01, 0.99), na.rm = TRUE)
  intensity_max <- max(abs(intensity_range))
  intensity_max <- ceiling(intensity_max)
  if (!is.finite(intensity_max) || intensity_max == 0) intensity_max <- 1

  expr_mat_scaled[expr_mat_scaled > intensity_max] <- intensity_max
  expr_mat_scaled[expr_mat_scaled < -intensity_max] <- -intensity_max

  col_intensity <- colorRamp2(c(-intensity_max, 0, intensity_max),
                              c("#3C5488", "white", "#DC0000"))

  cat(sprintf("  Matrix dimensions: %d metabolites x %d samples\n",
              nrow(expr_mat_scaled), ncol(expr_mat_scaled)))

  sample_groups <- sapply(colnames(expr_mat_scaled), function(s) {
    if (grepl(paste0("^", group1), s)) return(group1)
    if (grepl(paste0("^", group2), s)) return(group2)
    return("Other")
  })

  pair_colors <- c("#3C5488", "#DC0000")
  names(pair_colors) <- c(group1, group2)

  col_anno <- HeatmapAnnotation(
    Group = sample_groups,
    col = list(Group = pair_colors),
    annotation_name_side = "right",
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    show_legend = TRUE,
    annotation_legend_param = list(
      Group = list(
        title_gp = gpar(fontface = "bold"),
        labels_gp = gpar()
      )
    )
  )

  unique_classes <- unique(diff_data$super_class)
  n_classes <- length(unique_classes)
  npg_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                  "#7E6148", "#B09C85", "#00A1D5", "#79AF97")

  if (n_classes <= length(npg_colors)) {
    ontology_colors <- setNames(npg_colors[1:n_classes], unique_classes)
  } else {
    ontology_colors <- setNames(
      c(npg_colors, rainbow(n_classes - length(npg_colors)))[1:n_classes],
      unique_classes
    )
  }

  row_anno_df <- data.frame(
    Ontology = diff_data$super_class,
    Log2FC = diff_data$Log2FC_capped,
    neg_log_p = diff_data$neg_log_p_capped,
    VIP = as.numeric(diff_data$VIP),
    row.names = diff_data$row_label
  )
  vip_available_for_heatmap <- any(is.finite(row_anno_df$VIP))

  row_anno_1 <- rowAnnotation(
    Ontology = row_anno_df$Ontology,
    col = list(Ontology = ontology_colors),
    annotation_name_rot = 90,
    annotation_name_side = "bottom",
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    show_legend = TRUE,
    annotation_legend_param = list(
      Ontology = list(
        title_gp = gpar(fontface = "bold"),
        labels_gp = gpar()
      )
    ),
    width = unit(8, "mm")
  )

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
      annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
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
      annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
      show_legend = FALSE,
      gap = unit(2, "mm")
    )
  }

  n_samples <- ncol(expr_mat_scaled)
  n_metabolites <- nrow(expr_mat_scaled)
  cell_size <- unit(5, "mm")

  ht <- Heatmap(
    expr_mat_scaled,
    name = "Scaled Intensity",
    col = col_intensity,
    top_annotation = col_anno,
    column_names_rot = 90,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 8),
    cluster_columns = FALSE,
    column_split = sample_groups,
    column_title = NULL,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    right_annotation = row_anno_1,
    width = n_samples * cell_size,
    height = n_metabolites * cell_size,
    rect_gp = gpar(col = "black", lwd = 0.3),
    heatmap_legend_param = list(
      title = "Scaled Intensity",
      direction = "vertical",
      legend_height = unit(24, "mm"),
      title_gp = gpar(fontface = "bold"),
      labels_gp = gpar()
    )
  ) + row_anno_2

  lgd_fc <- Legend(
    title = "log2(FC)",
    col_fun = col_fc,
    direction = "vertical",
    legend_height = unit(20, "mm"),
    title_gp = gpar(fontface = "bold"),
    labels_gp = gpar()
  )

  lgd_p <- Legend(
    title = "-log(p)",
    col_fun = col_sequential,
    direction = "vertical",
    legend_height = unit(20, "mm"),
    title_gp = gpar(fontface = "bold"),
    labels_gp = gpar()
  )

  lgd_list <- list(lgd_fc, lgd_p)

  heatmap_width_mm <- n_samples * 5 + 150
  heatmap_height_mm <- n_metabolites * 5 + 80
  pdf_width <- max(heatmap_width_mm / 25.4, 14)
  # Scale the canvas to the retained row count after curated exclusions.
  pdf_height <- max(heatmap_height_mm / 25.4, 8.5)

  pdf_file <- file.path(out_dir, paste0("07b_heatmap_", comp_name, "_reverse.pdf"))
  pdf(pdf_file, width = pdf_width, height = pdf_height)
  draw(ht,
       annotation_legend_list = lgd_list,
       heatmap_legend_side = "right",
       annotation_legend_side = "right",
       merge_legend = TRUE,
       padding = unit(c(10, 10, 10, 10), "mm"))
  dev.off()

  png_file <- file.path(out_dir, paste0("07b_heatmap_", comp_name, "_reverse.png"))
  png(png_file, width = pdf_width * 300, height = pdf_height * 300, res = 300)
  draw(ht,
       annotation_legend_list = lgd_list,
       heatmap_legend_side = "right",
       annotation_legend_side = "right",
       merge_legend = TRUE,
       padding = unit(c(10, 10, 10, 10), "mm"))
  dev.off()

  svg_file <- file.path(out_dir, paste0("07b_heatmap_", comp_name, "_reverse.svg"))
  open_svg_device(svg_file, width = pdf_width, height = pdf_height)
  draw(ht,
       annotation_legend_list = lgd_list,
       heatmap_legend_side = "right",
       annotation_legend_side = "right",
       merge_legend = TRUE,
       padding = unit(c(10, 10, 10, 10), "mm"))
  dev.off()

  cat(sprintf("  Saved: %s\n", pdf_file))
  cat(sprintf("  Saved: %s\n", png_file))
  cat(sprintf("  Saved: %s\n", svg_file))
}

cat("\n=== 07b Heatmaps Complete ===\n")
cat(sprintf("Results saved to: %s/\n", out_dir))
