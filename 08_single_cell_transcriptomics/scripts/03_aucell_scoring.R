suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

source("scripts/utils.R")

message("[Step3] Loading annotated objects...")
combined <- readRDS(file.path(cfg$step2_dir, "step2_annotated_combined.rds"))
t_cells  <- readRDS(file.path(cfg$step2_dir, "step2_tcells_subclustered.rds"))

message(paste0("[Step3] Reading ", cfg$target_name, " gene set..."))
hub_genes <- read_hub_genes(cfg)

message("[Step3] Joining RNA layers...")
combined <- JoinLayers(combined, assay = "RNA")

expr_mat <- get_matrix_layer(combined, assay = "RNA", layer = "data")
rownames(expr_mat) <- toupper(rownames(expr_mat))

genes_use <- intersect(hub_genes, rownames(expr_mat))
if (length(genes_use) < 5) {
  stop("Too few hub genes matched in expression matrix (<5). Please check gene symbols.")
}
write.csv(
  data.frame(signature = "CoreTargets",
             n_matched = length(genes_use),
             matched_genes = paste(genes_use, collapse = ";"),
             stringsAsFactors = FALSE),
  file.path(cfg$step3_dir, "AUCell_matched_core_targets.csv"),
  row.names = FALSE
)

message("[Step3] Running AUCell...")
rankings <- AUCell_buildRankings(expr_mat, plotStats = FALSE, verbose = FALSE)
gene_sets <- setNames(list(genes_use), "CoreTargets")
auc <- AUCell_calcAUC(gene_sets, rankings, aucMaxRank = ceiling(0.05 * nrow(rankings)), verbose = FALSE)

auc_score <- as.numeric(getAUC(auc)[1, ])
names(auc_score) <- colnames(expr_mat)

score_name <- cfg$aucell_score_name
combined[[score_name]] <- auc_score[colnames(combined)]

common_cells <- intersect(colnames(t_cells), names(auc_score))
t_cells[[score_name]] <- NA_real_
t_cells[[score_name]][match(common_cells, colnames(t_cells))] <- auc_score[common_cells]

p_auc_global <- FeaturePlot(
  combined,
  features = score_name,
  reduction = "umap",
  split.by = "Phase",
  keep.scale = "all",
  min.cutoff = "q5",
  max.cutoff = "q95"
) &
  scale_color_gradientn(colors = c("#D3D3D3", "#C6DBEF", "#6BAED6", "#2171B5", "#084594")) &
  theme_classic(base_size = 11, base_family = "sans")
save_plot_pdf_png(p_auc_global, "CoreTargets_AUCell_FeaturePlot_Global", cfg$step3_dir, width = 13, height = 6)

p_auc_t <- FeaturePlot(
  t_cells,
  features = score_name,
  reduction = "umap",
  split.by = "Phase",
  keep.scale = "all",
  min.cutoff = "q5",
  max.cutoff = "q95"
) &
  scale_color_gradientn(colors = c("#D3D3D3", "#C6DBEF", "#6BAED6", "#2171B5", "#084594")) &
  theme_classic(base_size = 11, base_family = "sans")
save_plot_pdf_png(p_auc_t, "CoreTargets_AUCell_FeaturePlot_Tcell", cfg$step3_dir, width = 13, height = 6)

auc_meta <- data.frame(
  cell = names(auc_score),
  stringsAsFactors = FALSE
)
auc_meta[[score_name]] <- auc_score
write.csv(auc_meta, file.path(cfg$step3_dir, "AUCell_scores_all_cells.csv"), row.names = FALSE)

saveRDS(combined, file.path(cfg$step3_dir, "step3_annotated_with_aucell.rds"))
saveRDS(t_cells,  file.path(cfg$step3_dir, "step3_tcells_with_aucell.rds"))

celltype_order <- c("T_cells", "Myeloid", "B_cells", "NK_cells",
                    "Plasma", "Dendritic", "Epithelial",
                    "Neutrophils", "Basophils", "Progenitors", "Unknown")
meta_vln <- combined@meta.data
meta_vln$CellType_f <- factor(meta_vln$CellType,
                               levels = intersect(celltype_order, unique(meta_vln$CellType)))

p_vln_celltype <- ggplot(
  meta_vln[!is.na(meta_vln[[score_name]]), ],
  aes(x = CellType_f, y = .data[[score_name]], fill = CellType_f)
) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.5) +
  labs(
    x = "Cell type",
    y = score_name
  ) +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
save_plot_pdf_png(p_vln_celltype, "AUCell_VlnPlot_by_CellType", cfg$step3_dir, width = 10, height = 6)

message("[Step3] Completed. AUCell score added to metadata.")
