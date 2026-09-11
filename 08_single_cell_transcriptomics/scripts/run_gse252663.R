suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

source("config/analysis_config_gse252663.R")

cfg <- get_gse252663_config()
dir.create(cfg$results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$intermediate_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$figure_dir, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(cfg$figure_dir, full.names = TRUE), recursive = TRUE, force = TRUE)

set.seed(252663)

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
  if (!"name" %in% colnames(hub)) {
    stop("Column 'name' was not found in core target file: ", path)
  }
  unique(toupper(trimws(as.character(hub$name))))
}

read_one_sample <- function(smp) {
  counts <- ReadMtx(
    mtx = file.path(cfg$raw_dir, paste0(smp$file_prefix, "_matrix.mtx.gz")),
    cells = file.path(cfg$raw_dir, paste0(smp$file_prefix, "_barcodes.tsv.gz")),
    features = file.path(cfg$raw_dir, paste0(smp$file_prefix, "_features.tsv.gz")),
    feature.column = 2,
    unique.features = TRUE
  )
  seu <- CreateSeuratObject(
    counts = counts,
    project = smp$sample_name,
    min.cells = 3,
    min.features = 100
  )
  seu$Sample <- smp$sample_name
  seu$GEO_Accession <- smp$geo_accession
  seu$Genotype <- smp$genotype
  seu$Infection_Status <- smp$infection_status
  seu$Condition <- smp$condition
  seu$Replicate <- smp$replicate
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
  seu
}

get_data_matrix_upper <- function(seu) {
  mat <- GetAssayData(seu, assay = "RNA", layer = "data")
  rownames(mat) <- toupper(rownames(mat))
  mat
}

score_gene_sets <- function(seu, gene_sets, prefix) {
  mat <- get_data_matrix_upper(seu)
  matched <- lapply(gene_sets, function(genes) intersect(unique(toupper(genes)), rownames(mat)))
  score_df <- data.frame(cell = colnames(seu), stringsAsFactors = FALSE)
  for (nm in names(matched)) {
    genes <- matched[[nm]]
    score_name <- paste0(prefix, nm)
    if (length(genes) == 0) {
      score_df[[score_name]] <- NA_real_
    } else {
      score_df[[score_name]] <- Matrix::colMeans(mat[genes, , drop = FALSE])
    }
  }
  rownames(score_df) <- score_df$cell
  list(scores = score_df, matched = matched)
}

assign_cluster_by_signature <- function(seu, cluster_col, gene_sets, prefix, unknown_label = "Unknown") {
  res <- score_gene_sets(seu, gene_sets, prefix = prefix)
  for (nm in setdiff(colnames(res$scores), "cell")) {
    seu[[nm]] <- res$scores[colnames(seu), nm]
  }
  score_cols <- paste0(prefix, names(gene_sets))
  score_cols <- score_cols[score_cols %in% colnames(seu@meta.data)]
  cluster_scores <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    group_by(.data[[cluster_col]]) %>%
    summarise(across(all_of(score_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  cluster_label <- apply(cluster_scores[, score_cols, drop = FALSE], 1, function(x) {
    if (all(is.na(x)) || max(x, na.rm = TRUE) <= 0) {
      return(unknown_label)
    }
    gsub(paste0("^", prefix), "", names(which.max(x)))
  })
  cluster_map <- setNames(cluster_label, as.character(cluster_scores[[cluster_col]]))
  list(
    object = seu,
    labels = unname(cluster_map[as.character(seu@meta.data[[cluster_col]])]),
    cluster_scores = cluster_scores,
    matched = res$matched
  )
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

pairwise_infection_tests <- function(df, label_col, score_col) {
  bind_rows(lapply(unique(as.character(df[[label_col]])), function(label_value) {
    sub_df <- df[df[[label_col]] == label_value, , drop = FALSE]
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
      n_control = sum(sub_df$Infection_Status == "Control"),
      n_kp = sum(sub_df$Infection_Status == "KP"),
      median_control = median(sub_df[[score_col]][sub_df$Infection_Status == "Control"], na.rm = TRUE),
      median_kp = median(sub_df[[score_col]][sub_df$Infection_Status == "KP"], na.rm = TRUE),
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  })) %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      p_star = p_to_star(p_adj_BH)
    )
}

message("[GSE252663] Loading 10X matrices...")
if (!dir.exists(cfg$raw_dir)) {
  if (!file.exists(cfg$tar_file)) {
    stop("Raw directory and tar file were not found. Expected: ", cfg$raw_dir)
  }
  dir.create(cfg$raw_dir, showWarnings = FALSE, recursive = TRUE)
  untar(cfg$tar_file, exdir = cfg$raw_dir)
}

seu_list <- lapply(seq_len(nrow(cfg$sample_info)), function(i) read_one_sample(cfg$sample_info[i, ]))
combined <- Reduce(function(x, y) merge(x, y, add.cell.ids = c(unique(x$Sample)[1], unique(y$Sample)[1])), seu_list)
combined$Infection_Status <- factor(as.character(combined$Infection_Status), levels = cfg$group_order)
combined$Condition <- factor(as.character(combined$Condition), levels = cfg$condition_order)

message("[GSE252663] QC and integration...")
p_qc_pre <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    group.by = "Sample", pt.size = 0) +
  plot_annotation(title = "GSE252663 QC before filtering")
save_plot_pdf_png(p_qc_pre, "01_QC_before_filter", cfg$intermediate_dir, width = 13, height = 5)

qc <- cfg$qc
combined <- subset(
  combined,
  subset = nFeature_RNA >= qc$min_features &
    nFeature_RNA <= qc$max_features &
    percent.mt <= qc$max_percent_mt
)

combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- JoinLayers(combined, assay = "RNA")
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = cfg$seurat$nfeatures)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = cfg$seurat$npcs, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = cfg$seurat$dims_use)
combined <- FindNeighbors(combined, reduction = "pca", dims = cfg$seurat$dims_use)
combined <- FindClusters(combined, resolution = cfg$seurat$cluster_resolution)

p_qc_post <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     group.by = "Sample", pt.size = 0) +
  plot_annotation(title = "GSE252663 QC after filtering")
save_plot_pdf_png(p_qc_post, "02_QC_after_filter", cfg$intermediate_dir, width = 13, height = 5)

message("[GSE252663] Marker/signature annotation for mouse lung cells...")
major_markers <- list(
  T_cells = c("Cd3d", "Cd3e", "Cd3g", "Trac", "Cd247"),
  B_cells = c("Ms4a1", "Cd79a", "Cd79b", "Cd19", "Bank1"),
  Plasma = c("Jchain", "Mzb1", "Xbp1", "Sdc1"),
  NK_cells = c("Nkg7", "Klrb1c", "Gzma", "Gzmb", "Prf1"),
  Monocyte_Macrophage = c("Lyz2", "Csf1r", "Adgre1", "Fcgr3", "Lgals3", "C1qa"),
  Neutrophils = c("S100a8", "S100a9", "Retnlg", "Cxcr2", "Mmp8", "Elane"),
  Dendritic = c("Itgax", "Flt3", "Clec9a", "Cd209a", "Xcr1"),
  Epithelial = c("Epcam", "Krt8", "Krt18", "Sftpc", "Scgb1a1"),
  Endothelial = c("Pecam1", "Kdr", "Cdh5", "Vwf", "Esam"),
  Fibroblast = c("Col1a1", "Col1a2", "Dcn", "Lum", "Pdgfra")
)

major_res <- assign_cluster_by_signature(combined, "seurat_clusters", major_markers, "MajorScore_")
combined <- major_res$object
combined$CellType <- factor(major_res$labels)
write.csv(major_res$cluster_scores, file.path(cfg$intermediate_dir, "major_celltype_cluster_signature_scores.csv"), row.names = FALSE)

core_targets <- read_core_targets(cfg$hub_gene_file)
th17_sig <- c("Rorc", "Ccr6", "Il17a", "Il17f", "Il23r", "Stat3", "Rora", "Batf", "Maf")
treg_sig <- c("Foxp3", "Il2ra", "Ctla4", "Ikzf2", "Tigit", "Tnfrsf18", "Il10", "Lrrc32")

message("[GSE252663] AUCell scoring...")
expr_mat <- get_data_matrix_upper(combined)
core_matched <- intersect(core_targets, rownames(expr_mat))
if (length(core_matched) < 5) {
  stop("Too few core targets matched in mouse expression matrix: ", paste(core_matched, collapse = ", "))
}
rankings <- AUCell_buildRankings(expr_mat, plotStats = FALSE, verbose = FALSE)
auc <- AUCell_calcAUC(list(CoreTargets = core_matched), rankings,
                      aucMaxRank = ceiling(0.05 * nrow(rankings)), verbose = FALSE)
auc_score <- as.numeric(getAUC(auc)[1, ])
names(auc_score) <- colnames(expr_mat)
combined[[cfg$aucell_score_name]] <- auc_score[colnames(combined)]

sig_scores <- score_gene_sets(combined, list(Th17 = th17_sig, Treg = treg_sig), prefix = "")
for (nm in setdiff(colnames(sig_scores$scores), "cell")) {
  combined[[nm]] <- sig_scores$scores[colnames(combined), nm]
}

message("[GSE252663] T-cell subclustering...")
t_cells <- subset(combined, subset = CellType == "T_cells")
if (ncol(t_cells) >= 100) {
  t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 2000)
  t_cells[["RNA"]]$scale.data <- NULL
  t_cells <- ScaleData(t_cells)
  t_cells <- RunPCA(t_cells, npcs = 25, verbose = FALSE)
  t_cells <- RunUMAP(t_cells, reduction = "pca", dims = 1:20)
  t_cells <- FindNeighbors(t_cells, reduction = "pca", dims = 1:20)
  t_cells <- FindClusters(t_cells, resolution = cfg$seurat$tcell_cluster_resolution)
  t_cell_markers <- list(
    CD4_T = c("Cd4", "Il7r", "Ccr7", "Lef1"),
    CD8_T = c("Cd8a", "Cd8b1", "Gzmk", "Nkg7"),
    Treg = c("Foxp3", "Il2ra", "Ctla4", "Ikzf2"),
    Th17_like = c("Rorc", "Il17a", "Il17f", "Ccr6", "Il23r"),
    GdT = c("Trdc", "Trgc1", "Trgc2"),
    Proliferating_T = c("Mki67", "Top2a", "Stmn1"),
    Cytotoxic_T_NKlike = c("Nkg7", "Gzmb", "Prf1", "Klrb1c")
  )
  t_res <- assign_cluster_by_signature(t_cells, "seurat_clusters", t_cell_markers, "TScore_", unknown_label = "T_unknown")
  t_cells <- t_res$object
  t_cells$T_subtype <- factor(t_res$labels)
  write.csv(t_res$cluster_scores, file.path(cfg$intermediate_dir, "tcell_cluster_signature_scores.csv"), row.names = FALSE)
} else {
  warning("Too few T cells for subclustering; skipping T-cell figures.")
  t_cells$T_subtype <- "T_cells"
}

saveRDS(combined, file.path(cfg$intermediate_dir, "GSE252663_annotated_aucell_combined.rds"))
saveRDS(t_cells, file.path(cfg$intermediate_dir, "GSE252663_tcells_subclustered.rds"))

source("scripts/06_gse252663_core_target_scoring.R")
message("[GSE252663] Completed. Outputs saved in: ", cfg$figure_dir)
quit(save = "no", status = 0, runLast = FALSE)

message("[GSE252663] Publication tables and figures...")
write.csv(cfg$sample_info, file.path(cfg$figure_dir, "01_Table_sample_metadata.csv"), row.names = FALSE)
matched_tbl <- data.frame(
  signature = c("CoreTargets", "Th17", "Treg"),
  n_matched = c(length(core_matched), length(sig_scores$matched$Th17), length(sig_scores$matched$Treg)),
  matched_genes = c(
    paste(core_matched, collapse = ";"),
    paste(sig_scores$matched$Th17, collapse = ";"),
    paste(sig_scores$matched$Treg, collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write.csv(matched_tbl, file.path(cfg$figure_dir, "02_Table_signature_matched_genes.csv"), row.names = FALSE)

write.csv(
  as.data.frame(table(CellType = combined$CellType, Infection_Status = combined$Infection_Status, Condition = combined$Condition)),
  file.path(cfg$figure_dir, "03_Table_celltype_counts_by_condition.csv"),
  row.names = FALSE
)
if (ncol(t_cells) > 0) {
  write.csv(
    as.data.frame(table(T_subtype = t_cells$T_subtype, Infection_Status = t_cells$Infection_Status, Condition = t_cells$Condition)),
    file.path(cfg$figure_dir, "04_Table_tcell_subtype_counts_by_condition.csv"),
    row.names = FALSE
  )
}

score_name <- cfg$aucell_score_name
combined_meta <- combined@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    Infection_Status = factor(as.character(Infection_Status), levels = cfg$group_order),
    Condition = factor(as.character(Condition), levels = cfg$condition_order)
  )

celltype_order <- combined_meta %>%
  count(CellType, sort = TRUE) %>%
  pull(CellType) %>%
  as.character()
combined$CellType_paper <- factor(as.character(combined$CellType), levels = celltype_order)
combined_meta$CellType_paper <- factor(as.character(combined_meta$CellType), levels = celltype_order)

write.csv(
  combined_meta %>%
    group_by(CellType_paper, Infection_Status) %>%
    summarise(
      n_cells = n(),
      CoreTargets_AUC_median = median(.data[[score_name]], na.rm = TRUE),
      CoreTargets_AUC_mean = mean(.data[[score_name]], na.rm = TRUE),
      Th17_score_median = median(Th17, na.rm = TRUE),
      Treg_score_median = median(Treg, na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(cfg$figure_dir, "05_Table_coretarget_auc_by_celltype_infection.csv"),
  row.names = FALSE
)
global_pairwise <- pairwise_infection_tests(combined_meta, "CellType_paper", score_name)
write.csv(global_pairwise, file.path(cfg$figure_dir, "06_Table_celltype_KP_vs_Control_Wilcoxon.csv"), row.names = FALSE)

base_dim_theme <- theme_classic(base_size = 10, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 7)
  )

p_global <- DimPlot(combined, reduction = "umap", group.by = "CellType_paper",
                    label = TRUE, repel = TRUE, pt.size = 0.15, shuffle = TRUE) +
  labs(title = "GSE252663 mouse lung cell types", color = "Cell type") +
  base_dim_theme
save_plot_pdf_png(p_global, "07_Fig1_Global_UMAP_CellType_Labelled", cfg$figure_dir, width = 9, height = 7)

p_condition <- DimPlot(combined, reduction = "umap", group.by = "Condition",
                       pt.size = 0.15, shuffle = TRUE) +
  labs(title = "Sample condition", color = "Condition") +
  base_dim_theme
save_plot_pdf_png(p_condition, "08_Fig2_Global_UMAP_Condition", cfg$figure_dir, width = 8, height = 6)

feature_colors <- c("#D3D3D3", "#C6DBEF", "#6BAED6", "#2171B5", "#084594")
p_feature <- FeaturePlot(
  combined,
  features = score_name,
  reduction = "umap",
  split.by = "Infection_Status",
  keep.scale = "all",
  min.cutoff = "q5",
  max.cutoff = "q95",
  pt.size = 0.12
) &
  scale_color_gradientn(colors = feature_colors) &
  theme_classic(base_size = 9, base_family = "sans") &
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
save_plot_pdf_png(p_feature, "09_Fig3_CoreTargets_AUCell_FeaturePlot_by_InfectionStatus",
                  cfg$figure_dir, width = 10, height = 4.5)

box_theme <- theme_bw(base_size = 10, base_family = "sans") +
  theme(
    panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11)
  )
group_palette <- c(Control = "#2166AC", KP = "#C51B2D")

p_box <- ggplot(combined_meta, aes(CellType_paper, .data[[score_name]], fill = Infection_Status, color = Infection_Status)) +
  geom_boxplot(position = position_dodge(width = 0.72), width = 0.55,
               outlier.shape = NA, alpha = 0.82, color = "#2F2F2F") +
  geom_point(position = position_jitterdodge(jitter.width = 0.13, dodge.width = 0.72),
             size = 0.12, alpha = 0.16, show.legend = FALSE) +
  scale_fill_manual(values = group_palette, drop = FALSE) +
  scale_color_manual(values = group_palette, drop = FALSE) +
  labs(title = "Core target AUCell score by lung cell type", x = NULL, y = "Core target AUCell score") +
  box_theme
save_plot_pdf_png(p_box, "10_Fig4_CoreTargets_AUC_by_CellType_InfectionStatus",
                  cfg$figure_dir, width = 12, height = 5.5)

core_target_features <- rownames(combined)[toupper(rownames(combined)) %in% core_targets]
if (length(core_target_features) > 0) {
  p_dot <- DotPlot(
    combined,
    features = core_target_features,
    group.by = "CellType_paper",
    dot.scale = 5,
    col.min = -2,
    col.max = 2
  ) +
    RotatedAxis() +
    labs(title = "Core target gene expression by lung cell type", x = NULL, y = NULL) +
    theme_bw(base_size = 10, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )
  save_plot_pdf_png(p_dot, "11_Fig5_CoreTarget_Gene_DotPlot_by_CellType",
                    cfg$figure_dir, width = 10, height = 6)
}

if (ncol(t_cells) >= 100) {
  t_meta <- t_cells@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(
      Infection_Status = factor(as.character(Infection_Status), levels = cfg$group_order),
      Condition = factor(as.character(Condition), levels = cfg$condition_order)
    )
  t_order <- t_meta %>% count(T_subtype, sort = TRUE) %>% pull(T_subtype) %>% as.character()
  t_cells$T_subtype_paper <- factor(as.character(t_cells$T_subtype), levels = t_order)
  t_meta$T_subtype_paper <- factor(as.character(t_meta$T_subtype), levels = t_order)

  write.csv(
    t_meta %>%
      group_by(T_subtype_paper, Infection_Status) %>%
      summarise(
        n_cells = n(),
        CoreTargets_AUC_median = median(.data[[score_name]], na.rm = TRUE),
        Th17_score_median = median(Th17, na.rm = TRUE),
        Treg_score_median = median(Treg, na.rm = TRUE),
        .groups = "drop"
      ),
    file.path(cfg$figure_dir, "12_Table_tcell_coretarget_auc_by_infection.csv"),
    row.names = FALSE
  )
  t_pairwise <- pairwise_infection_tests(t_meta, "T_subtype_paper", score_name)
  write.csv(t_pairwise, file.path(cfg$figure_dir, "13_Table_tcell_KP_vs_Control_Wilcoxon.csv"), row.names = FALSE)

  p_t_umap <- DimPlot(t_cells, reduction = "umap", group.by = "T_subtype_paper",
                      label = TRUE, repel = TRUE, pt.size = 0.35, shuffle = TRUE) +
    labs(title = "T-cell subtypes", color = "T subtype") +
    base_dim_theme
  save_plot_pdf_png(p_t_umap, "14_Fig6_Tcell_UMAP_Subtype_Labelled", cfg$figure_dir, width = 9, height = 7)

  p_t_box <- ggplot(t_meta, aes(T_subtype_paper, .data[[score_name]], fill = Infection_Status, color = Infection_Status)) +
    geom_boxplot(position = position_dodge(width = 0.72), width = 0.55,
                 outlier.shape = NA, alpha = 0.82, color = "#2F2F2F") +
    geom_point(position = position_jitterdodge(jitter.width = 0.13, dodge.width = 0.72),
               size = 0.22, alpha = 0.22, show.legend = FALSE) +
    scale_fill_manual(values = group_palette, drop = FALSE) +
    scale_color_manual(values = group_palette, drop = FALSE) +
    labs(title = "Core target AUCell score by T-cell subtype", x = NULL, y = "Core target AUCell score") +
    box_theme
  save_plot_pdf_png(p_t_box, "15_Fig7_CoreTargets_AUC_by_TcellSubtype_InfectionStatus",
                    cfg$figure_dir, width = 10, height = 5.2)

  validation_features <- rownames(t_cells)[toupper(rownames(t_cells)) %in% c(toupper(th17_sig), toupper(treg_sig))]
  if (length(validation_features) > 0) {
    p_val_dot <- DotPlot(t_cells, features = validation_features, group.by = "T_subtype_paper",
                         dot.scale = 5, col.min = -2, col.max = 2) +
      RotatedAxis() +
      labs(title = "Th17/Treg marker validation", x = NULL, y = NULL) +
      theme_bw(base_size = 10, base_family = "sans") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

    t_sig_long <- t_meta %>%
      select(cell, T_subtype_paper, Infection_Status, Th17, Treg) %>%
      pivot_longer(cols = c(Th17, Treg), names_to = "Signature", values_to = "Score")
    p_val_score <- ggplot(t_sig_long, aes(T_subtype_paper, Score, fill = T_subtype_paper)) +
      geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.85, color = "#2F2F2F") +
      facet_wrap(~Signature, ncol = 1, scales = "free_y") +
      labs(title = "Signature score", x = NULL, y = "Average log-normalized expression") +
      theme_bw(base_size = 10, base_family = "sans") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "none")

    save_plot_pdf_png(p_val_dot + p_val_score + plot_layout(widths = c(1.4, 1)),
                      "16_Fig8_Tcell_Th17_Treg_marker_signature_validation",
                      cfg$figure_dir, width = 14, height = 6.5)
  }
}

report_lines <- c(
  "# GSE252663 AUCell Result Report",
  "",
  "## Dataset",
  "",
  "- GEO: GSE252663, single-cell RNA-seq on Klebsiella-infected mouse lung.",
  "- Samples: KO-KP1, KO-KP2, WT-KP1, WT-KP2, WT-Control, KO-Control.",
  "- Main analysis grouping: KP infection versus Control; genotype is retained as sample metadata.",
  "",
  "## Analysis Notes",
  "",
  "- This pipeline uses mouse lung marker/signature annotation instead of the human PBMC Monaco/STCAT workflow.",
  "- Core target AUCell uses the nine target genes after mouse gene-symbol matching by uppercase symbols.",
  "- Th17/Treg labels in mouse lung T cells are marker/signature-supported exploratory labels; sparse Il17a/Foxp3 expression should be interpreted cautiously.",
  "- Cell-level Wilcoxon tests are provided as descriptive screening statistics, not as fully independent biological replicates.",
  "",
  "## Key Output Files",
  "",
  "- 07_Fig1: global mouse lung cell-type UMAP.",
  "- 09_Fig3: core target AUCell FeaturePlot split by KP infection status.",
  "- 10_Fig4: core target AUCell boxplot by lung cell type and infection status.",
  "- 11_Fig5: single-gene DotPlot for the nine core targets.",
  "- 14_Fig6 to 16_Fig8: T-cell subtype and Th17/Treg marker-signature validation, if enough T cells are present."
)
writeLines(report_lines, file.path(cfg$figure_dir, "00_Result_report.md"))

message("[GSE252663] Completed. Outputs saved in: ", cfg$figure_dir)
quit(save = "no", status = 0, runLast = FALSE)
