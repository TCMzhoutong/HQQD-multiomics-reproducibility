suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

if (!exists("cfg")) {
  source("config/analysis_config_gse252663.R")
  cfg <- get_gse252663_config()
}

if (!exists("combined")) {
  combined <- readRDS(file.path(cfg$intermediate_dir, "GSE252663_annotated_aucell_combined.rds"))
}
if (!exists("t_cells")) {
  t_cells <- readRDS(file.path(cfg$intermediate_dir, "GSE252663_tcells_subclustered.rds"))
}

dir.create(cfg$figure_dir, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(cfg$figure_dir, full.names = TRUE), recursive = TRUE, force = TRUE)

save_plot_pdf_png <- function(plot_obj, file_stem, out_dir, width = 8, height = 6, dpi = 320) {
  pdf_file <- file.path(out_dir, paste0(file_stem, ".pdf"))
  png_file <- file.path(out_dir, paste0(file_stem, ".png"))
  svg_file <- file.path(out_dir, paste0(file_stem, ".svg"))
  ggsave(pdf_file, plot = plot_obj, width = width, height = height, device = cairo_pdf)
  ggsave(png_file, plot = plot_obj, width = width, height = height, dpi = dpi)
  ggsave(svg_file, plot = plot_obj, width = width, height = height, device = grDevices::svg)
  invisible(list(pdf = pdf_file, png = png_file, svg = svg_file))
}

read_core_targets <- function(path) {
  hub <- read.csv(path, stringsAsFactors = FALSE)
  unique(toupper(trimws(as.character(hub$name))))
}

p_to_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

pairwise_tests_by_label <- function(df, label_col, score_col) {
  bind_rows(lapply(unique(as.character(df[[label_col]])), function(label_value) {
    sub_df <- df[df[[label_col]] == label_value, , drop = FALSE]
    sub_df <- sub_df[!is.na(sub_df[[score_col]]) & !is.na(sub_df$Infection_Status), , drop = FALSE]
    if (length(unique(sub_df$Infection_Status)) < 2) {
      return(NULL)
    }
    p_value <- tryCatch(
      wilcox.test(sub_df[[score_col]] ~ droplevels(sub_df$Infection_Status))$p.value,
      error = function(e) NA_real_
    )
    data.frame(
      label = label_value,
      group1 = "Control",
      group2 = "KP",
      n_group1 = sum(sub_df$Infection_Status == "Control"),
      n_group2 = sum(sub_df$Infection_Status == "KP"),
      median_group1 = median(sub_df[[score_col]][sub_df$Infection_Status == "Control"], na.rm = TRUE),
      median_group2 = median(sub_df[[score_col]][sub_df$Infection_Status == "KP"], na.rm = TRUE),
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  })) %>%
    mutate(
      p_adj_BH_within_label = p.adjust(p_value, method = "BH"),
      p_adj_BH_global = p.adjust(p_value, method = "BH"),
      p_star = p_to_star(p_adj_BH_within_label)
    )
}

build_two_group_brackets <- function(df, label_col, score_col, pairwise_df) {
  if (nrow(pairwise_df) == 0) {
    return(data.frame())
  }
  label_levels <- levels(df[[label_col]])
  offsets <- c(Control = -0.18, KP = 0.18)
  score_range <- range(df[[score_col]], na.rm = TRUE)
  y_step <- diff(score_range) * 0.055
  if (!is.finite(y_step) || y_step <= 0) y_step <- 0.02
  y_tip <- y_step * 0.23
  y_df <- df %>%
    group_by(.data[[label_col]]) %>%
    summarise(y_base = max(.data[[score_col]], na.rm = TRUE) + y_step, .groups = "drop")

  pairwise_df %>%
    mutate(
      star = p_to_star(p_adj_BH_within_label),
      label_index = match(label, label_levels),
      x = label_index + unname(offsets[group1]),
      xend = label_index + unname(offsets[group2])
    ) %>%
    filter(star != "", !is.na(label_index), !is.na(x), !is.na(xend)) %>%
    left_join(y_df, by = setNames(label_col, "label")) %>%
    mutate(
      y = y_base,
      y_tip = y - y_tip,
      xmid = (x + xend) / 2
    )
}

score_gene_sets_average <- function(seu, gene_sets) {
  mat <- GetAssayData(seu, assay = "RNA", layer = "data")
  rownames(mat) <- toupper(rownames(mat))
  matched <- lapply(gene_sets, function(genes) intersect(unique(toupper(genes)), rownames(mat)))
  score_df <- data.frame(cell = colnames(seu), stringsAsFactors = FALSE)
  for (nm in names(matched)) {
    score_df[[nm]] <- if (length(matched[[nm]]) == 0) {
      NA_real_
    } else {
      Matrix::colMeans(mat[matched[[nm]], , drop = FALSE])
    }
  }
  rownames(score_df) <- score_df$cell
  list(scores = score_df, matched = matched)
}

score_name <- cfg$aucell_score_name
core_targets <- read_core_targets(cfg$hub_gene_file)
th17_sig <- c("Rorc", "Ccr6", "Il17a", "Il17f", "Il23r", "Stat3", "Rora", "Batf", "Maf")
treg_sig <- c("Foxp3", "Il2ra", "Ctla4", "Ikzf2", "Tigit", "Tnfrsf18", "Il10", "Lrrc32")

combined <- JoinLayers(combined, assay = "RNA")
t_cells <- JoinLayers(t_cells, assay = "RNA")

combined$Infection_Status <- factor(as.character(combined$Infection_Status), levels = cfg$group_order)
t_cells$Infection_Status <- factor(as.character(t_cells$Infection_Status), levels = cfg$group_order)

celltype_order <- c(
  "T_cells", "B_cells", "NK_cells", "Monocyte_Macrophage",
  "Dendritic", "Neutrophils", "Epithelial", "Endothelial",
  "Fibroblast", "Plasma", "Unknown"
)
celltype_levels <- intersect(celltype_order, unique(as.character(combined$CellType)))
extra_celltypes <- setdiff(unique(as.character(combined$CellType)), celltype_levels)
combined$CellType_paper <- factor(as.character(combined$CellType), levels = c(celltype_levels, sort(extra_celltypes)))

if (!"T_subtype" %in% colnames(t_cells@meta.data)) {
  t_cells$T_subtype <- "T_cells"
}

# Keep the same T-subtype label system used in the GSE263380 publication figures.
# The dataset is mouse lung tissue, so labels are assigned with mouse ortholog
# marker/signature scores rather than the human STCAT model.
reference_t_levels <- c(
  "CD4 Tn", "CD4 Tcm", "CD4 Tc", "CD8 Temra", "CD8 Tem", "Tgd",
  "CD8 Tn", "CD8 Tc", "CD4 Treg", "CD4 Tem", "CD8 Tcm", "MAIT",
  "CD4 Temra", "CD4 Th17", "CD4 Th1"
)
reference_t_markers <- list(
  "CD4 Tn" = c("Cd4", "Ccr7", "Sell", "Lef1", "Tcf7", "Il7r"),
  "CD4 Tcm" = c("Cd4", "Ccr7", "Il7r", "Tcf7", "Cd44"),
  "CD4 Tc" = c("Cd4", "Nkg7", "Gzmb", "Prf1", "Ccl5", "Ctsw"),
  "CD8 Temra" = c("Cd8a", "Cd8b1", "Klrg1", "Cx3cr1", "Gzmb", "Prf1"),
  "CD8 Tem" = c("Cd8a", "Cd8b1", "Gzmk", "Ccl5", "Nkg7", "Cxcr3"),
  "Tgd" = c("Trdc", "Trgc1", "Trgc2", "Cd3d", "Cd3e"),
  "CD8 Tn" = c("Cd8a", "Cd8b1", "Ccr7", "Sell", "Lef1", "Tcf7"),
  "CD8 Tc" = c("Cd8a", "Cd8b1", "Nkg7", "Gzmb", "Prf1", "Ctsw"),
  "CD4 Treg" = c("Cd4", "Foxp3", "Il2ra", "Ctla4", "Ikzf2", "Tigit", "Tnfrsf18"),
  "CD4 Tem" = c("Cd4", "Il7r", "Gzmk", "Cxcr3", "Ccl5", "Cd44"),
  "CD8 Tcm" = c("Cd8a", "Cd8b1", "Ccr7", "Il7r", "Tcf7", "Cd44"),
  "MAIT" = c("Zbtb16", "Rorc", "Slamf6", "Il18r1", "Klrb1c", "Slc4a10"),
  "CD4 Temra" = c("Cd4", "Klrg1", "Cx3cr1", "Gzmb", "Prf1", "Ccl5"),
  "CD4 Th17" = c("Cd4", "Rorc", "Ccr6", "Il17a", "Il17f", "Il23r", "Maf", "Batf"),
  "CD4 Th1" = c("Cd4", "Tbx21", "Ifng", "Cxcr3", "Stat4", "Il12rb2")
)

reference_t_score_res <- score_gene_sets_average(t_cells, reference_t_markers)
reference_t_score_mat <- as.matrix(reference_t_score_res$scores[, reference_t_levels, drop = FALSE])
reference_t_score_z <- scale(reference_t_score_mat)
reference_t_score_z[is.na(reference_t_score_z)] <- -Inf
t_display <- reference_t_levels[max.col(reference_t_score_z, ties.method = "first")]
t_display[!is.finite(rowSums(reference_t_score_z))] <- "Other / rare T"
t_cells$T_major_paper <- factor(t_display, levels = c(reference_t_levels, "Other / rare T"))

reference_t_matched_tbl <- data.frame(
  T_subtype = names(reference_t_score_res$matched),
  n_matched = vapply(reference_t_score_res$matched, length, integer(1)),
  matched_genes = vapply(reference_t_score_res$matched, paste, collapse = ";", FUN.VALUE = character(1)),
  stringsAsFactors = FALSE
)
write.csv(reference_t_matched_tbl, file.path(cfg$figure_dir, "01A_Table_Tsubtype_reference_marker_matched_genes.csv"), row.names = FALSE)

sig_scores_t <- score_gene_sets_average(t_cells, list(Th17 = th17_sig, Treg = treg_sig))
for (nm in setdiff(colnames(sig_scores_t$scores), "cell")) {
  t_cells[[paste0(nm, "_score")]] <- sig_scores_t$scores[colnames(t_cells), nm]
}
sig_scores_combined <- score_gene_sets_average(combined, list(Th17 = th17_sig, Treg = treg_sig))
for (nm in setdiff(colnames(sig_scores_combined$scores), "cell")) {
  combined[[paste0(nm, "_score")]] <- sig_scores_combined$scores[colnames(combined), nm]
}

matched_tbl <- data.frame(
  signature = c("CoreTargets", "Th17", "Treg"),
  n_matched = c(
    sum(core_targets %in% toupper(rownames(combined))),
    length(sig_scores_t$matched$Th17),
    length(sig_scores_t$matched$Treg)
  ),
  matched_genes = c(
    paste(intersect(core_targets, toupper(rownames(combined))), collapse = ";"),
    paste(sig_scores_t$matched$Th17, collapse = ";"),
    paste(sig_scores_t$matched$Treg, collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write.csv(matched_tbl, file.path(cfg$figure_dir, "01_Table_signature_matched_genes.csv"), row.names = FALSE)

global_df <- combined@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(!is.na(CellType_paper), !is.na(Infection_Status), !is.na(.data[[score_name]])) %>%
  mutate(CellType_paper = factor(as.character(CellType_paper), levels = levels(combined$CellType_paper)))
t_df <- t_cells@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(!is.na(T_major_paper), !is.na(Infection_Status), !is.na(.data[[score_name]])) %>%
  mutate(T_major_paper = factor(as.character(T_major_paper), levels = levels(t_cells$T_major_paper)))

write.csv(
  global_df %>% count(CellType_paper, Infection_Status, name = "n_cells") %>% arrange(CellType_paper, Infection_Status),
  file.path(cfg$figure_dir, "02_Table_Global_celltype_counts_by_InfectionStatus.csv"),
  row.names = FALSE
)
write.csv(
  t_df %>% count(T_major_paper, Infection_Status, name = "n_cells") %>% arrange(T_major_paper, Infection_Status),
  file.path(cfg$figure_dir, "03_Table_Tsubtype_counts_for_paper_landscape.csv"),
  row.names = FALSE
)
write.csv(
  global_df %>%
    group_by(CellType_paper, Infection_Status) %>%
    summarise(n_cells = n(),
              CoreTargets_AUC_median = median(.data[[score_name]], na.rm = TRUE),
              CoreTargets_AUC_mean = mean(.data[[score_name]], na.rm = TRUE),
              .groups = "drop"),
  file.path(cfg$figure_dir, "04_Table_Global_celltype_scores_by_InfectionStatus.csv"),
  row.names = FALSE
)
write.csv(
  t_df %>%
    group_by(T_major_paper, Infection_Status) %>%
    summarise(n_cells = n(),
              CoreTargets_AUC_median = median(.data[[score_name]], na.rm = TRUE),
              CoreTargets_AUC_mean = mean(.data[[score_name]], na.rm = TRUE),
              Th17_score_median = median(Th17_score, na.rm = TRUE),
              Treg_score_median = median(Treg_score, na.rm = TRUE),
              .groups = "drop"),
  file.path(cfg$figure_dir, "05_Table_Tsubtype_scores_by_InfectionStatus.csv"),
  row.names = FALSE
)

global_pairwise <- pairwise_tests_by_label(global_df, "CellType_paper", score_name)
t_pairwise <- pairwise_tests_by_label(t_df, "T_major_paper", score_name)
write.csv(global_pairwise, file.path(cfg$figure_dir, "06A_Table_Global_celltype_pairwise_Wilcoxon_by_InfectionStatus.csv"), row.names = FALSE)
write.csv(t_pairwise, file.path(cfg$figure_dir, "06B_Table_Tsubtype_pairwise_Wilcoxon_by_InfectionStatus.csv"), row.names = FALSE)
global_brackets <- build_two_group_brackets(global_df, "CellType_paper", score_name, global_pairwise)
t_brackets <- build_two_group_brackets(t_df, "T_major_paper", score_name, t_pairwise)

base_dim_theme <- theme_classic(base_size = 10, base_family = "sans") +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.65, "cm"),
    legend.key.height = unit(0.72, "cm"),
    legend.spacing.y = unit(0.35, "cm"),
    plot.margin = margin(8, 8, 8, 8, "pt")
  )

p_global_umap <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "CellType_paper",
  label = TRUE,
  label.size = 5,
  repel = TRUE,
  pt.size = 0.16,
  shuffle = TRUE
) +
  labs(title = "Global cell types", color = "Cell type") +
  coord_fixed() +
  base_dim_theme
save_plot_pdf_png(p_global_umap, "07_Fig1_Global_UMAP_CellType_Labelled", cfg$figure_dir, width = 10, height = 7)

p_t_umap <- DimPlot(
  t_cells,
  reduction = "umap",
  group.by = "T_major_paper",
  label = TRUE,
  label.size = 5,
  repel = TRUE,
  pt.size = 0.3,
  shuffle = TRUE
) +
  labs(title = "T-cell subtypes", color = "T subtype") +
  coord_fixed() +
  base_dim_theme
save_plot_pdf_png(p_t_umap, "08_Fig2_Tcell_UMAP_MajorSubtype_Labelled", cfg$figure_dir, width = 10, height = 7)

feature_colors <- c("#00E5FF", "#8A2BE2", "#4B0082")

make_split_feature_plot <- function(seu, split_by, point_size) {
  split_labels <- levels(droplevels(factor(seu@meta.data[[split_by]])))
  plots <- FeaturePlot(
    seu,
    features = score_name,
    reduction = "umap",
    split.by = split_by,
    keep.scale = "all",
    min.cutoff = "q5",
    max.cutoff = "q95",
    pt.size = point_size,
    raster = FALSE,
    combine = FALSE
  )
  plots <- lapply(seq_along(plots), function(i) {
    p <- plots[[i]] +
      labs(title = split_labels[[i]], color = score_name) +
      scale_color_gradientn(
        colors = feature_colors,
        name = score_name,
        guide = guide_colorbar(display = "rectangles")
      ) +
      coord_fixed() +
      theme_classic(base_size = 9, base_family = "sans") +
      theme(
        aspect.ratio = 1,
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank()
      )
    if (i > 1) {
      p <- p +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()
        )
    }
    p
  })
  wrap_plots(plots, nrow = 1, guides = "collect") &
    theme(legend.position = "right")
}

p_feature_global <- make_split_feature_plot(combined, "Infection_Status", 0.07)
save_plot_pdf_png(p_feature_global, "09_Fig3_Global_CoreTargets_AUCell_FeaturePlot_by_InfectionStatus",
                  cfg$figure_dir, width = 10.5, height = 4.4)

p_feature_t <- make_split_feature_plot(t_cells, "Infection_Status", 0.09)
save_plot_pdf_png(p_feature_t, "10_Fig4_Tcell_CoreTargets_AUCell_FeaturePlot_by_InfectionStatus",
                  cfg$figure_dir, width = 10.5, height = 4.4)

group_palette <- c(Control = "#2166AC", KP = "#C51B2D")
box_theme <- theme_bw(base_size = 10, base_family = "sans") +
  theme(
    panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11)
  )

make_consistent_style_box <- function(df, x_var, title_text, bracket_df = NULL) {
  p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[score_name]], fill = Infection_Status, color = Infection_Status)) +
    geom_boxplot(position = position_dodge(width = 0.72), width = 0.55,
                 outlier.shape = NA, alpha = 0.82, color = "#2F2F2F") +
    geom_point(position = position_jitterdodge(jitter.width = 0.13, dodge.width = 0.72),
               size = 0.16, alpha = 0.18, show.legend = FALSE) +
    scale_fill_manual(values = group_palette, drop = FALSE) +
    scale_color_manual(values = group_palette, drop = FALSE) +
    labs(title = title_text, x = NULL, y = "Core target AUCell score") +
    coord_cartesian(clip = "off") +
    box_theme
  if (!is.null(bracket_df) && nrow(bracket_df) > 0) {
    p <- p +
      geom_segment(data = bracket_df, aes(x = x, xend = xend, y = y, yend = y),
                   inherit.aes = FALSE, color = "black", linewidth = 0.28) +
      geom_segment(data = bracket_df, aes(x = x, xend = x, y = y, yend = y_tip),
                   inherit.aes = FALSE, color = "black", linewidth = 0.28) +
      geom_segment(data = bracket_df, aes(x = xend, xend = xend, y = y, yend = y_tip),
                   inherit.aes = FALSE, color = "black", linewidth = 0.28) +
      geom_text(data = bracket_df, aes(x = xmid, y = y + 0.004, label = star),
                inherit.aes = FALSE, size = 3.0, fontface = "bold")
  }
  p
}

p_box_global <- make_consistent_style_box(global_df, "CellType_paper", "Global cell types", global_brackets)
p_box_t <- make_consistent_style_box(t_df, "T_major_paper", "T-cell subtypes", t_brackets) +
  labs(y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
box_y_limit <- range(c(global_df[[score_name]], t_df[[score_name]],
                       global_brackets$y, t_brackets$y), na.rm = TRUE)
box_y_limit[1] <- min(0, box_y_limit[1])
p_box_global <- p_box_global + coord_cartesian(ylim = box_y_limit, clip = "off")
p_box_t <- p_box_t + coord_cartesian(ylim = box_y_limit, clip = "off")
box_widths <- c(
  length(levels(droplevels(global_df$CellType_paper))),
  length(levels(droplevels(t_df$T_major_paper)))
)
p_box_combined <- p_box_global + p_box_t +
  plot_layout(widths = box_widths, guides = "collect") &
  theme(legend.position = "right")
save_plot_pdf_png(p_box_combined, "11_Fig5_CoreTargets_AUC_Global_and_Tcell_by_InfectionStatus",
                  cfg$figure_dir, width = 16, height = 5.6)

core_features_global <- rownames(combined)[toupper(rownames(combined)) %in% core_targets]
core_features_t <- rownames(t_cells)[toupper(rownames(t_cells)) %in% core_targets]
paper_dot_gene_order <- c("PTGS2", "PTPRC", "PPARG", "MAPK14", "SIRT1", "CASP1", "AKT1", "BCL2", "MMP9")
dot_feature_map <- setNames(intersect(core_features_global, core_features_t),
                            toupper(intersect(core_features_global, core_features_t)))
dot_features <- unname(dot_feature_map[paper_dot_gene_order[paper_dot_gene_order %in% names(dot_feature_map)]])
if (length(dot_features) > 0) {
  p_dot_global <- DotPlot(
    combined,
    features = dot_features,
    group.by = "CellType_paper",
    dot.scale = 5,
    col.min = -2,
    col.max = 2
  ) +
    RotatedAxis() +
    scale_x_discrete(labels = toupper) +
    labs(title = "Global cell types", x = NULL, y = NULL) +
    guides(color = guide_colorbar(display = "rectangles")) +
    theme_bw(base_size = 12, base_family = "sans") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12),
          legend.position = "none")

  p_dot_t <- DotPlot(
    t_cells,
    features = dot_features,
    group.by = "T_major_paper",
    dot.scale = 5,
    col.min = -2,
    col.max = 2
  ) +
    RotatedAxis() +
    scale_x_discrete(labels = toupper) +
    labs(title = "T-cell subtypes", x = NULL, y = NULL) +
    guides(color = guide_colorbar(display = "rectangles")) +
    theme_bw(base_size = 12, base_family = "sans") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.55, "cm"))

  p_dot_combined <- p_dot_global + p_dot_t +
    plot_layout(widths = c(1, 1))
  save_plot_pdf_png(p_dot_combined, "12_Fig6_CoreTarget_Gene_DotPlot_Global_and_Tcell",
                    cfg$figure_dir, width = 15, height = 6.4)
}

th_ref_df <- t_df %>%
  mutate(T_reference = ifelse(grepl("Th17", as.character(T_major_paper)), "CD4 Th17",
                              ifelse(grepl("Treg", as.character(T_major_paper)), "CD4 Treg", NA_character_))) %>%
  filter(!is.na(T_reference))
write.csv(
  th_ref_df %>%
    group_by(T_reference, Infection_Status) %>%
    summarise(n_cells = n(),
              Th17_score_median = median(Th17_score, na.rm = TRUE),
              Treg_score_median = median(Treg_score, na.rm = TRUE),
              CoreTargets_AUC_median = median(.data[[score_name]], na.rm = TRUE),
              .groups = "drop"),
  file.path(cfg$figure_dir, "06C_Table_Tcell_Th17_Treg_signature_validation.csv"),
  row.names = FALSE
)

if (nrow(th_ref_df) > 0) {
  p_th17 <- ggplot(th_ref_df, aes(T_reference, Th17_score, fill = Infection_Status, color = Infection_Status)) +
    geom_boxplot(position = position_dodge(width = 0.72), width = 0.55,
                 outlier.shape = NA, alpha = 0.82, color = "#2F2F2F") +
    geom_point(position = position_jitterdodge(jitter.width = 0.13, dodge.width = 0.72),
               size = 0.35, alpha = 0.28, show.legend = FALSE) +
    scale_fill_manual(values = group_palette, drop = FALSE) +
    scale_color_manual(values = group_palette, drop = FALSE) +
    labs(x = NULL, y = "Th17 signature score") +
    box_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  save_plot_pdf_png(p_th17, "13_Fig7_CD4Th17_Signature_TregReference_by_InfectionStatus",
                    cfg$figure_dir, width = 7, height = 4.5)
}

validation_features <- rownames(t_cells)[toupper(rownames(t_cells)) %in% c(toupper(th17_sig), toupper(treg_sig))]
if (length(validation_features) > 0) {
  p_validation_dot <- DotPlot(
    t_cells,
    features = validation_features,
    group.by = "T_major_paper",
    dot.scale = 5,
    col.min = -2,
    col.max = 2
  ) +
    RotatedAxis() +
    labs(title = "Marker expression", x = NULL, y = NULL) +
    guides(color = guide_colorbar(display = "rectangles")) +
    theme_bw(base_size = 10, base_family = "sans") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  sig_long <- t_df %>%
    select(cell, T_major_paper, Th17_score, Treg_score) %>%
    pivot_longer(cols = c(Th17_score, Treg_score), names_to = "Signature", values_to = "Score") %>%
    mutate(Signature = recode(Signature, Th17_score = "Th17 signature", Treg_score = "Treg signature"))

  p_validation_sig <- ggplot(sig_long, aes(T_major_paper, Score, fill = T_major_paper)) +
    geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.85, color = "#2F2F2F") +
    facet_wrap(~Signature, ncol = 1, scales = "free_y") +
    labs(title = "Signature score", x = NULL, y = "Average log-normalized expression") +
    theme_bw(base_size = 10, base_family = "sans") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "none")

  p_validation <- p_validation_dot + p_validation_sig +
    plot_layout(widths = c(1.4, 1))
  save_plot_pdf_png(p_validation, "14_Fig8_Tcell_Th17_Treg_marker_signature_validation",
                    cfg$figure_dir, width = 14, height = 6.5)
}

report_lines <- c(
  "# GSE252663 AUCell Result Report",
  "",
  "## Dataset and Grouping",
  "",
  "- Dataset: GSE252663, mouse lung scRNA-seq after Klebsiella pneumoniae infection.",
  "- Samples used: KO-KP1, KO-KP2, WT-KP1, WT-KP2, WT-Control, and KO-Control.",
  "- Primary comparison: Control versus KP infection.",
  "- Genotype is retained in metadata while the main analysis uses the primary disease/infection grouping.",
  "",
  "## Output Style",
  "",
  "- Outputs include global and T-cell UMAPs, feature plots, score comparisons, dot plots, and Th17/Treg marker-signature validation.",
  "- Human Monaco/STCAT annotation was not reused because this is mouse lung tissue; mouse marker/signature annotation is used instead.",
  "- T-cell subtype labels follow the previous GSE263380 label system; labels are assigned using mouse ortholog marker/signature scores.",
  "",
  "## Interpretation Notes",
  "",
  "- All nine core targets were matched in the mouse expression matrix.",
  "- T-cell subtypes are shown with the same old-style labels used in the previous result folder.",
  "- Th17/Treg labels should still be interpreted through the marker/signature validation figure, because this is mouse lung tissue and not the original human PBMC STCAT setting."
)
writeLines(report_lines, file.path(cfg$figure_dir, "00_Result_report.md"))

message("[GSE252663 Step6] Publication-style outputs saved in: ", cfg$figure_dir)
