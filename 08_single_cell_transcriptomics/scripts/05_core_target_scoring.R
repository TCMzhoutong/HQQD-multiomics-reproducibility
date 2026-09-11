suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
})

source("scripts/utils.R")

message("[Step5] Loading AUCell objects for paper figures...")
combined <- readRDS(file.path(cfg$step3_dir, "step3_annotated_with_aucell.rds"))
t_cells <- readRDS(file.path(cfg$step3_dir, "step3_tcells_with_aucell.rds"))
combined <- JoinLayers(combined, assay = "RNA")
t_cells <- JoinLayers(t_cells, assay = "RNA")

out_dir <- file.path(cfg$results_dir, "source_tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(out_dir, full.names = TRUE), recursive = TRUE, force = TRUE)

score_name <- cfg$aucell_score_name
group_order <- cfg$comparison$group_order
exclude_t <- cfg$comparison$t_exclude_labels
group_palette <- c(
  Moderate = "#C51B2D",
  Severe = "#E69F00",
  Convalescent = "#2166AC"
)

p_to_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

ensure_disease_state <- function(seu) {
  if (!"Disease_State" %in% colnames(seu@meta.data)) {
    seu$Disease_State <- ifelse(
      as.character(seu$Phase) == "Convalescent",
      "Convalescent",
      as.character(seu$Severity)
    )
  }
  seu$Disease_State <- factor(as.character(seu$Disease_State), levels = group_order)
  seu
}

combined <- ensure_disease_state(combined)
t_cells <- ensure_disease_state(t_cells)

if (!(score_name %in% colnames(combined@meta.data))) {
  stop("Core-target score not found in combined metadata: ", score_name)
}
if (!(score_name %in% colnames(t_cells@meta.data))) {
  stop("Core-target score not found in T-cell metadata: ", score_name)
}

core_target_sig <- read_hub_genes(cfg)
th17_sig <- c("RORC", "CCR6", "IL17A", "IL17F", "KLRB1", "IL23R", "STAT3", "RORA", "BATF", "MAF")
treg_sig <- c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "TNFRSF18", "CCR8", "ENTPD1", "IL10", "LRRC32")

score_gene_sets_auc <- function(seu, gene_sets) {
  expr <- get_matrix_layer(seu, assay = "RNA", layer = "data")
  rownames(expr) <- toupper(rownames(expr))
  gene_sets <- lapply(gene_sets, function(x) intersect(unique(toupper(x)), rownames(expr)))
  matched <- data.frame(
    signature = names(gene_sets),
    n_matched = vapply(gene_sets, length, integer(1)),
    matched_genes = vapply(gene_sets, paste, collapse = ";", FUN.VALUE = character(1)),
    stringsAsFactors = FALSE
  )
  rankings <- AUCell_buildRankings(expr, plotStats = FALSE, verbose = FALSE)
  auc <- AUCell_calcAUC(gene_sets, rankings, aucMaxRank = ceiling(0.05 * nrow(rankings)), verbose = FALSE)
  score_mat <- as.data.frame(t(as.matrix(getAUC(auc))))
  colnames(score_mat) <- paste0(colnames(score_mat), "_AUC")
  score_mat$cell <- rownames(score_mat)
  list(scores = score_mat, matched = matched)
}

sig_res <- score_gene_sets_auc(
  t_cells,
  list(Th17 = th17_sig, Treg = treg_sig, CoreTargets = core_target_sig)
)
write.csv(sig_res$matched, file.path(out_dir, "01_Table_signature_matched_genes.csv"), row.names = FALSE)

score_df <- sig_res$scores
rownames(score_df) <- score_df$cell
for (nm in setdiff(colnames(score_df), "cell")) {
  t_cells[[nm]] <- score_df[colnames(t_cells), nm]
}

valid_t <- t_cells@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(!is.na(T_subtype), !T_subtype %in% exclude_t)

major_t_counts <- valid_t %>%
  count(T_subtype, name = "n_cells") %>%
  arrange(desc(n_cells))
force_keep <- c("CD4 Th17", "CD4 Treg", "CD4 Th1", "MAIT", "Tgd")
major_t_labels <- major_t_counts %>%
  filter(n_cells >= 250 | T_subtype %in% force_keep) %>%
  pull(T_subtype) %>%
  as.character()
major_t_levels <- as.character(major_t_counts$T_subtype[major_t_counts$T_subtype %in% major_t_labels])

t_cells$T_major_paper <- ifelse(
  as.character(t_cells$T_subtype) %in% major_t_labels,
  as.character(t_cells$T_subtype),
  "Other / rare T"
)
t_cells$T_major_paper <- factor(t_cells$T_major_paper, levels = c(major_t_levels, "Other / rare T"))

celltype_order <- c(
  "T_cells", "Myeloid", "B_cells", "NK_cells", "Dendritic",
  "Neutrophils", "Plasma", "Epithelial", "Basophils", "Progenitors", "Unknown"
)
combined$CellType_paper <- factor(as.character(combined$CellType),
                                  levels = intersect(celltype_order, unique(as.character(combined$CellType))))

global_df <- combined@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(!is.na(CellType_paper), !is.na(Disease_State), !is.na(.data[[score_name]])) %>%
  mutate(Disease_State = factor(as.character(Disease_State), levels = group_order))

major_t_df <- t_cells@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(T_major_paper != "Other / rare T",
         !is.na(Disease_State),
         !is.na(.data[[score_name]])) %>%
  mutate(
    T_major_paper = factor(as.character(T_major_paper), levels = major_t_levels),
    Disease_State = factor(as.character(Disease_State), levels = group_order)
  )

write.csv(
  global_df %>%
    count(CellType_paper, Disease_State, name = "n_cells") %>%
    arrange(CellType_paper, Disease_State),
  file.path(out_dir, "02_Table_Global_celltype_counts_by_DiseaseState.csv"),
  row.names = FALSE
)

write.csv(
  major_t_counts %>% mutate(in_paper_landscape = T_subtype %in% major_t_labels),
  file.path(out_dir, "03_Table_Tsubtype_counts_for_paper_landscape.csv"),
  row.names = FALSE
)

write.csv(
  global_df %>%
    group_by(CellType_paper, Disease_State) %>%
    summarise(
      n_cells = n(),
      CoreTargets_AUC_median = median(.data[[score_name]], na.rm = TRUE),
      CoreTargets_AUC_mean = mean(.data[[score_name]], na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(out_dir, "04_Table_Global_celltype_scores_by_DiseaseState.csv"),
  row.names = FALSE
)

write.csv(
  major_t_df %>%
    group_by(T_major_paper, Disease_State) %>%
    summarise(
      n_cells = n(),
      CoreTargets_AUC_median = median(.data[[score_name]], na.rm = TRUE),
      CoreTargets_AUC_mean = mean(.data[[score_name]], na.rm = TRUE),
      Th17_AUC_median = median(Th17_AUC, na.rm = TRUE),
      Treg_AUC_median = median(Treg_AUC, na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(out_dir, "05_Table_Tsubtype_scores_by_DiseaseState.csv"),
  row.names = FALSE
)

pairwise_tests_by_label <- function(df, label_col) {
  label_values <- unique(as.character(df[[label_col]]))
  test_rows <- lapply(label_values, function(label_value) {
    sub_df <- df[df[[label_col]] == label_value, , drop = FALSE]
    groups_present <- group_order[group_order %in% as.character(unique(sub_df$Disease_State))]
    if (length(groups_present) < 2) {
      return(NULL)
    }
    pairs <- combn(groups_present, 2, simplify = FALSE)
    do.call(rbind, lapply(pairs, function(pair) {
      pair_df <- sub_df[as.character(sub_df$Disease_State) %in% pair, , drop = FALSE]
      p_value <- tryCatch(
        wilcox.test(pair_df[[score_name]] ~ droplevels(pair_df$Disease_State))$p.value,
        error = function(e) NA_real_
      )
      data.frame(
        label = label_value,
        group1 = pair[[1]],
        group2 = pair[[2]],
        n_group1 = sum(as.character(pair_df$Disease_State) == pair[[1]]),
        n_group2 = sum(as.character(pair_df$Disease_State) == pair[[2]]),
        median_group1 = median(pair_df[[score_name]][as.character(pair_df$Disease_State) == pair[[1]]], na.rm = TRUE),
        median_group2 = median(pair_df[[score_name]][as.character(pair_df$Disease_State) == pair[[2]]], na.rm = TRUE),
        p_value = p_value,
        stringsAsFactors = FALSE
      )
    }))
  })
  res <- bind_rows(test_rows)
  if (nrow(res) == 0) {
    return(res)
  }
  res %>%
    group_by(label) %>%
    mutate(p_adj_BH_within_label = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(
      p_adj_BH_global = p.adjust(p_value, method = "BH"),
      p_star = p_to_star(p_adj_BH_within_label)
    )
}

build_pairwise_brackets <- function(df, label_col, pairwise_df) {
  if (nrow(pairwise_df) == 0) {
    return(data.frame())
  }
  label_levels <- levels(df[[label_col]])
  group_offsets <- c(Moderate = -0.24, Severe = 0, Convalescent = 0.24)
  score_range <- range(df[[score_name]], na.rm = TRUE)
  y_step <- diff(score_range) * 0.055
  if (!is.finite(y_step) || y_step <= 0) {
    y_step <- max(df[[score_name]], na.rm = TRUE) * 0.08
  }
  y_tip <- y_step * 0.23
  y_df <- df %>%
    group_by(.data[[label_col]]) %>%
    summarise(y_base = max(.data[[score_name]], na.rm = TRUE) + y_step, .groups = "drop")

  pairwise_df %>%
    mutate(
      star = p_to_star(p_adj_BH_within_label),
      label_index = match(label, label_levels),
      x = label_index + unname(group_offsets[group1]),
      xend = label_index + unname(group_offsets[group2])
    ) %>%
    filter(star != "", !is.na(label_index), !is.na(x), !is.na(xend)) %>%
    arrange(label, p_adj_BH_within_label) %>%
    group_by(label) %>%
    mutate(rank_in_label = row_number()) %>%
    ungroup() %>%
    left_join(y_df, by = setNames(label_col, "label")) %>%
    mutate(
      y = y_base + (rank_in_label - 1) * y_step,
      y_tip = y - y_tip,
      xmid = (x + xend) / 2
    )
}

global_pairwise <- pairwise_tests_by_label(global_df, "CellType_paper")
t_pairwise <- pairwise_tests_by_label(major_t_df, "T_major_paper")
write.csv(
  global_pairwise,
  file.path(out_dir, "06A_Table_Global_celltype_pairwise_Wilcoxon_by_DiseaseState.csv"),
  row.names = FALSE
)
write.csv(
  t_pairwise,
  file.path(out_dir, "06B_Table_Tsubtype_pairwise_Wilcoxon_by_DiseaseState.csv"),
  row.names = FALSE
)
global_bracket_df <- build_pairwise_brackets(global_df, "CellType_paper", global_pairwise)
t_bracket_df <- build_pairwise_brackets(major_t_df, "T_major_paper", t_pairwise)

th17_df <- t_cells@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(T_subtype %in% c("CD4 Th17", "CD4 Treg"),
         !is.na(Disease_State),
         !is.na(.data[[score_name]])) %>%
  mutate(
    T_subtype = factor(T_subtype, levels = c("CD4 Th17", "CD4 Treg")),
    Disease_State = factor(as.character(Disease_State), levels = group_order)
  )

write.csv(
  valid_t %>%
    group_by(Disease_State = ifelse(as.character(Phase) == "Convalescent", "Convalescent", as.character(Severity))) %>%
    summarise(
      valid_t_cells = n(),
      cd4_like_cells = sum(grepl("^CD4", T_subtype)),
      cd4_th17_cells = sum(T_subtype == "CD4 Th17"),
      cd4_treg_cells = sum(T_subtype == "CD4 Treg"),
      th17_pct_valid_t = cd4_th17_cells / valid_t_cells * 100,
      th17_pct_cd4_like = cd4_th17_cells / cd4_like_cells * 100,
      treg_pct_valid_t = cd4_treg_cells / valid_t_cells * 100,
      treg_pct_cd4_like = cd4_treg_cells / cd4_like_cells * 100,
      .groups = "drop"
    ),
  file.path(out_dir, "06_Table_CD4_Th17_frequency_Treg_reference_by_DiseaseState.csv"),
  row.names = FALSE
)

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

# Figure 1: global parent cell-type landscape.
p_global_umap <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "CellType_paper",
  label = TRUE,
  label.size = 5,
  repel = TRUE,
  pt.size = 0.2,
  shuffle = TRUE
) +
  labs(title = "Global cell types", color = "Cell type") +
  coord_fixed() +
  base_dim_theme
save_plot_pdf_png(p_global_umap, "07_Fig1_Global_UMAP_CellType_Labelled",
                  out_dir, width = 10, height = 7)

# Figure 2: labelled T-cell subtype landscape.
p_t_umap <- DimPlot(
  t_cells,
  reduction = "umap",
  group.by = "T_major_paper",
  label = TRUE,
  label.size = 5,
  repel = TRUE,
  pt.size = 0.25,
  shuffle = TRUE
) +
  labs(title = "T-cell subtypes", color = "T subtype") +
  coord_fixed() +
  base_dim_theme
save_plot_pdf_png(p_t_umap, "08_Fig2_Tcell_UMAP_MajorSubtype_Labelled",
                  out_dir, width = 10, height = 7)

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

p_feature_global <- make_split_feature_plot(combined, "Disease_State", 0.08)
save_plot_pdf_png(p_feature_global, "09_Fig3_Global_CoreTargets_AUCell_FeaturePlot_by_DiseaseState",
                  out_dir, width = 13, height = 4.2)

p_feature_t <- make_split_feature_plot(t_cells, "Disease_State", 0.1)
save_plot_pdf_png(p_feature_t, "10_Fig4_Tcell_CoreTargets_AUCell_FeaturePlot_by_DiseaseState",
                  out_dir, width = 13, height = 4.2)

box_theme <- theme_bw(base_size = 10, base_family = "sans") +
  theme(
    panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.55, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11)
  )

make_old_style_box <- function(df, x_var, title_text, bracket_df = NULL) {
  p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[score_name]], fill = Disease_State, color = Disease_State)) +
    geom_boxplot(position = position_dodge(width = 0.72), width = 0.55,
                 outlier.shape = NA, alpha = 0.82, color = "#2F2F2F") +
    geom_point(position = position_jitterdodge(jitter.width = 0.14, dodge.width = 0.72),
               size = 0.18, alpha = 0.22, show.legend = FALSE) +
    scale_fill_manual(values = group_palette, drop = FALSE) +
    scale_color_manual(values = group_palette, drop = FALSE) +
    labs(title = title_text, x = NULL, y = "Core target AUCell score") +
    coord_cartesian(clip = "off") +
    box_theme
  if (!is.null(bracket_df) && nrow(bracket_df) > 0) {
    p <- p +
      geom_segment(
        data = bracket_df,
        aes(x = x, xend = xend, y = y, yend = y),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.28
      ) +
      geom_segment(
        data = bracket_df,
        aes(x = x, xend = x, y = y, yend = y_tip),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.28
      ) +
      geom_segment(
        data = bracket_df,
        aes(x = xend, xend = xend, y = y, yend = y_tip),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.28
      ) +
      geom_text(
        data = bracket_df,
        aes(x = xmid, y = y + 0.004, label = star),
        inherit.aes = FALSE,
        size = 3.0,
        fontface = "bold"
      )
  }
  p
}

p_box_global <- make_old_style_box(global_df, "CellType_paper", "Global cell types", global_bracket_df)
p_box_t <- make_old_style_box(major_t_df, "T_major_paper", "T-cell subtypes", t_bracket_df) +
  labs(y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
box_y_limit <- range(c(global_df[[score_name]], major_t_df[[score_name]],
                       global_bracket_df$y, t_bracket_df$y), na.rm = TRUE)
box_y_limit[1] <- min(0, box_y_limit[1])
p_box_global <- p_box_global +
  coord_cartesian(ylim = box_y_limit, clip = "off") +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf,
           color = "white", linewidth = 1.1) +
  theme(legend.position = "none")
p_box_t <- p_box_t +
  coord_cartesian(ylim = box_y_limit, clip = "off") +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf,
           color = "white", linewidth = 1.1) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0.985, 0.86),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = ggplot2::alpha("white", 0.88), color = "#BDBDBD", linewidth = 0.25),
    legend.box.background = element_blank()
  )
box_widths <- c(
  length(levels(droplevels(global_df$CellType_paper))),
  length(levels(droplevels(major_t_df$T_major_paper)))
)
p_box_combined <- p_box_global + p_box_t +
  plot_layout(widths = box_widths, guides = "keep")
save_plot_pdf_png(p_box_combined, "11_Fig5_CoreTargets_AUC_Global_and_Tcell_by_DiseaseState",
                  out_dir, width = 16, height = 5.6)

core_target_features_global <- rownames(combined)[toupper(rownames(combined)) %in% core_target_sig]
core_target_features_t <- rownames(t_cells)[toupper(rownames(t_cells)) %in% core_target_sig]
dot_features <- intersect(core_target_features_global, core_target_features_t)

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
    labs(title = "Global cell types", x = NULL, y = NULL) +
    guides(color = guide_colorbar(display = "rectangles")) +
    theme_bw(base_size = 12, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.position = "none"
    )

  p_dot_t <- DotPlot(
    subset(t_cells, subset = T_major_paper != "Other / rare T"),
    features = dot_features,
    group.by = "T_major_paper",
    dot.scale = 5,
    col.min = -2,
    col.max = 2
  ) +
    RotatedAxis() +
    labs(title = "T-cell subtypes", x = NULL, y = NULL) +
    guides(color = guide_colorbar(display = "rectangles")) +
    theme_bw(base_size = 12, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.key.size = unit(0.55, "cm")
    )

  p_dot_combined <- p_dot_global + p_dot_t +
    plot_layout(widths = c(1, 1))
  save_plot_pdf_png(p_dot_combined, "12_Fig6_CoreTarget_Gene_DotPlot_Global_and_Tcell",
                    out_dir, width = 15, height = 6.4)
}

if (nrow(th17_df) > 0) {
  p_th17 <- ggplot(th17_df, aes(T_subtype, Th17_AUC, fill = Disease_State, color = Disease_State)) +
    geom_boxplot(position = position_dodge(width = 0.72), width = 0.55,
                 outlier.shape = NA, alpha = 0.82, color = "#2F2F2F") +
    geom_point(position = position_jitterdodge(jitter.width = 0.14, dodge.width = 0.72),
               size = 0.35, alpha = 0.28, show.legend = FALSE) +
    scale_fill_manual(values = group_palette, drop = FALSE) +
    scale_color_manual(values = group_palette, drop = FALSE) +
    labs(x = NULL, y = "Th17 AUCell score", fill = NULL) +
    box_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  save_plot_pdf_png(p_th17, "13_Fig7_CD4Th17_Signature_TregReference_by_DiseaseState",
                    out_dir, width = 7, height = 4.5)
}

cd4_validation_levels <- c(
  "CD4 Th17", "CD4 Treg", "CD4 Th1", "CD4 Tn", "CD4 Tcm",
  "CD4 Tem", "CD4 Tc", "CD4 Temra", "Other CD4 T"
)
t_cells$CD4_validation_group <- ifelse(
  grepl("^CD4", as.character(t_cells$T_subtype)) &
    as.character(t_cells$T_subtype) %in% cd4_validation_levels,
  as.character(t_cells$T_subtype),
  ifelse(grepl("^CD4", as.character(t_cells$T_subtype)), "Other CD4 T", NA_character_)
)
t_cells$CD4_validation_group <- factor(t_cells$CD4_validation_group, levels = cd4_validation_levels)

cd4_validation_df <- t_cells@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(!is.na(CD4_validation_group),
         !is.na(Th17_AUC),
         !is.na(Treg_AUC)) %>%
  mutate(CD4_validation_group = factor(as.character(CD4_validation_group), levels = cd4_validation_levels))

compare_signature_target <- function(df, signature_col, target_label) {
  comparison_labels <- setdiff(levels(df$CD4_validation_group), target_label)
  res <- lapply(comparison_labels, function(other_label) {
    pair_df <- df %>%
      filter(CD4_validation_group %in% c(target_label, other_label)) %>%
      mutate(group2 = factor(as.character(CD4_validation_group), levels = c(target_label, other_label)))
    if (n_distinct(pair_df$group2) < 2) {
      return(NULL)
    }
    p_value <- tryCatch(
      wilcox.test(pair_df[[signature_col]] ~ pair_df$group2)$p.value,
      error = function(e) NA_real_
    )
    data.frame(
      signature = signature_col,
      target_label = target_label,
      reference_label = other_label,
      n_target = sum(pair_df$group2 == target_label),
      n_reference = sum(pair_df$group2 == other_label),
      median_target = median(pair_df[[signature_col]][pair_df$group2 == target_label], na.rm = TRUE),
      median_reference = median(pair_df[[signature_col]][pair_df$group2 == other_label], na.rm = TRUE),
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(res) %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      p_star = p_to_star(p_adj_BH)
    )
}

validation_stats <- bind_rows(
  compare_signature_target(cd4_validation_df, "Th17_AUC", "CD4 Th17"),
  compare_signature_target(cd4_validation_df, "Treg_AUC", "CD4 Treg")
)
write.csv(
  validation_stats,
  file.path(out_dir, "06C_Table_CD4_Th17_Treg_signature_validation.csv"),
  row.names = FALSE
)

if (nrow(cd4_validation_df) > 0) {
  validation_features <- rownames(t_cells)[toupper(rownames(t_cells)) %in% c(th17_sig, treg_sig)]
  t_cells_cd4_validation <- subset(t_cells, subset = !is.na(CD4_validation_group))
  if (length(validation_features) > 0) {
    p_validation_dot <- DotPlot(
      t_cells_cd4_validation,
      features = validation_features,
      group.by = "CD4_validation_group",
      dot.scale = 5,
      col.min = -2,
      col.max = 2
    ) +
      RotatedAxis() +
      labs(title = "Marker expression", x = NULL, y = NULL) +
      guides(color = guide_colorbar(display = "rectangles")) +
      theme_bw(base_size = 10, base_family = "sans") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
      )

    cd4_sig_long <- cd4_validation_df %>%
      select(cell, CD4_validation_group, Th17_AUC, Treg_AUC) %>%
      pivot_longer(cols = c(Th17_AUC, Treg_AUC), names_to = "Signature", values_to = "AUC") %>%
      mutate(Signature = recode(Signature, Th17_AUC = "Th17 signature", Treg_AUC = "Treg signature"))

    p_validation_sig <- ggplot(cd4_sig_long, aes(CD4_validation_group, AUC, fill = CD4_validation_group)) +
      geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.85, color = "#2F2F2F") +
      facet_wrap(~Signature, ncol = 1, scales = "free_y") +
      labs(title = "Signature score", x = NULL, y = "AUCell score") +
      theme_bw(base_size = 10, base_family = "sans") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none"
      )

    p_validation <- p_validation_dot + p_validation_sig +
      plot_layout(widths = c(1.4, 1))
    save_plot_pdf_png(p_validation, "14_Fig8_CD4_Th17_Treg_marker_signature_validation",
                      out_dir, width = 14, height = 6.5)
  }
}

report_lines <- c(
  "# Core Target AUCell Result Report",
  "",
  "## Dataset and Grouping",
  "",
  "- Dataset: GSE263380, five PBMC scRNA-seq samples from CRKP pneumonia patients.",
  "- Samples used: two Moderate acute, one Severe acute, and two Convalescent samples.",
  "- Main grouping used for paper figures: Moderate, Severe, Convalescent.",
  "- Caveat: Severe acute contains one sample, so group-level statistical claims should be conservative.",
  "",
  "## Published-Figure Folder",
  "",
  "- 07_Fig1: global parent cell-type UMAP.",
  "- 08_Fig2: T-cell subtype UMAP.",
  "- 09_Fig3 and 10_Fig4: Global and T-cell AUCell FeaturePlots by disease state.",
  "- 11_Fig5: Global and T-cell AUCell boxplots in the old two-panel style.",
  "- 06A and 06B: pairwise Wilcoxon tests for the three disease-state groups; brackets in 11_Fig5 mark significant within-label pairwise BH-adjusted comparisons.",
  "- 12_Fig6: Global and T-cell single-gene DotPlots for the nine core targets.",
  "- 13_Fig7: CD4 Th17 signature with CD4 Treg as reference.",
  "- 14_Fig8 and 06C: CD4 Th17/Treg marker and signature validation after STCAT annotation.",
  "",
  "## Key Findings",
  "",
  "- All nine core targets were matched: AKT1, BCL2, CASP1, MAPK14, MMP9, PPARG, PTGS2, PTPRC, and SIRT1.",
  "- CD4 Th17 cells were detectable but sparse, so Th17 is suitable as a mechanistic anchor but not as a strong standalone population-level claim.",
  "- The single-cell result should not be framed as evidence for restoration of the Th17/Treg ratio.",
  "- PTPRC is broadly expressed across leukocytes; CASP1 and BCL2 show clearer T-cell expression than PTGS2, MMP9, or PPARG.",
  "- PTGS2, MMP9, and PPARG may be more informative in non-T immune compartments than in T-cell-centered analysis."
)
writeLines(report_lines, file.path(out_dir, "00_Result_report.md"))

message("[Step5] Completed. Paper-ready outputs saved in ", out_dir)
