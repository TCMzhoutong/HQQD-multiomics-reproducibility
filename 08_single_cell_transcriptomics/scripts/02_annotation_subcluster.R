suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SingleR)
  library(celldex)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

source("scripts/utils.R")

message("[Step2] Loading integrated object...")
combined <- readRDS(file.path(cfg$step1_dir, "step1_harmony_combined.rds"))

message("[Step2] Loading Monaco Immune reference...")
ref <- celldex::MonacoImmuneData()

message("[Step2] Running SingleR for major cell type annotation...")
sce_combined <- SingleCellExperiment(
  assays = list(
    counts    = GetAssayData(combined, assay = "RNA", layer = "counts"),
    logcounts = GetAssayData(combined, assay = "RNA", layer = "data")
  )
)
pred_main <- SingleR(
  test   = sce_combined,
  ref    = ref,
  labels = ref$label.main
)
combined$SingleR_main  <- pred_main$labels
combined$SingleR_pruned <- pred_main$pruned.labels

main_map <- c(
  "CD8+ T cells" = "T_cells",
  "CD4+ T cells" = "T_cells",
  "T cells"      = "T_cells",
  "B cells"      = "B_cells",
  "Monocytes"    = "Myeloid",
  "NK cells"     = "NK_cells",
  "Dendritic cells" = "Dendritic",
  "Neutrophils"  = "Neutrophils",
  "Basophils"    = "Basophils",
  "Progenitors"  = "Progenitors"
)
combined$CellType <- dplyr::recode(combined$SingleR_main, !!!main_map, .default = "Unknown")
combined$CellType <- factor(combined$CellType)

p_umap_clust <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  theme_classic(base_size = 11, base_family = "sans")
save_plot_pdf_png(p_umap_clust, "Global_UMAP_ClusterNumbers", cfg$step2_dir, width = 9, height = 7)

p_global <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "CellType",
  label = TRUE,
  repel = TRUE
) +
  theme_classic(base_size = 11, base_family = "sans")
save_plot_pdf_png(p_global, "Global_UMAP_CellType", cfg$step2_dir, width = 9, height = 7)

if ("tsne" %in% names(combined@reductions)) {
  p_tsne_global <- DimPlot(combined, reduction = "tsne", group.by = "CellType",
                           label = TRUE, repel = TRUE) +
    theme_classic(base_size = 11, base_family = "sans")
  save_plot_pdf_png(p_tsne_global, "Global_tSNE_CellType", cfg$step2_dir, width = 9, height = 7)
}

t_cells <- subset(combined, subset = CellType == "T_cells")

if (ncol(t_cells) == 0) {
  stop("No cells found in CellType == 'T_cells'. Please check Step2 SingleR main labels.")
}

message("[Step2] T cell sub-clustering (with Harmony batch correction)...")
t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 2000)
t_cells[["RNA"]]$scale.data <- NULL
t_cells <- ScaleData(t_cells)
t_cells <- RunPCA(t_cells, npcs = 30, verbose = FALSE)
# Re-run Harmony in the T-cell subset by patient rather than sample. Sample is close
# to disease state in GSE263380, so sample-level correction can remove biology of interest.
t_cells <- RunHarmony(t_cells, group.by.vars = "Patient_ID", reduction = "pca",
                      reduction.save = "harmony", verbose = FALSE)
t_cells <- RunUMAP(t_cells, reduction = "harmony", dims = 1:20)
t_cells <- FindNeighbors(t_cells, reduction = "harmony", dims = 1:20)
t_cells <- FindClusters(t_cells, resolution = cfg$seurat$tcell_cluster_resolution)
t_cells <- JoinLayers(t_cells)

n_tcell_clusters <- nlevels(t_cells$seurat_clusters)
p_tcell_clust <- DimPlot(
  t_cells,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label     = TRUE,
  repel     = TRUE
) +
  theme_classic(base_size = 11, base_family = "sans")
save_plot_pdf_png(p_tcell_clust, "Tcell_UMAP_ClusterNumbers", cfg$step2_dir, width = 9, height = 7)
message(sprintf("[Step2] T cell clusters: %d (resolution=%.1f)",
                n_tcell_clusters, cfg$seurat$tcell_cluster_resolution))


message("[Step2] Exporting T-cell matrices for STCAT (via Matrix::writeMM)...")
export_dir <- file.path(cfg$step2_dir, "stcat_export")
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

Matrix::writeMM(
  GetAssayData(t_cells, assay = "RNA", layer = "counts"),
  file.path(export_dir, "counts.mtx")
)
Matrix::writeMM(
  GetAssayData(t_cells, assay = "RNA", layer = "data"),
  file.path(export_dir, "data.mtx")
)
writeLines(colnames(t_cells), file.path(export_dir, "barcodes.tsv"))
writeLines(rownames(t_cells), file.path(export_dir, "features.tsv"))

obs_meta <- t_cells@meta.data
obs_meta$Barcode <- rownames(obs_meta)
write.csv(obs_meta, file.path(export_dir, "obs_metadata.csv"))

# UMAP embedding
umap_df <- as.data.frame(Embeddings(t_cells, "umap"))
umap_df$Barcode <- rownames(umap_df)
write.csv(umap_df, file.path(export_dir, "umap.csv"))

message("[Step2] Matrix export done: ", export_dir)

if (!file.exists(file.path(export_dir, "counts.mtx"))) {
  stop("Matrix export failed: counts.mtx not found in ", export_dir)
}

py_script <- file.path(cfg$root_dir, "scripts", "stcat_annotate.py")
if (!file.exists(py_script)) {
  stop("STCAT script not found: ", py_script)
}

pred_path <- file.path(cfg$step2_dir, "stcat_predictions.csv")
stcat_input <- export_dir
conda_bin <- Sys.getenv("ANALYSIS_CONDA_EXE", "")
if (!nzchar(conda_bin)) {
  conda_bin <- Sys.which("conda")
}
if (!nzchar(conda_bin) && .Platform$OS.type == "windows") {
  conda_bin <- Sys.which("conda.bat")
}
if (!nzchar(conda_bin)) {
  stop("Conda executable not found. Please set ANALYSIS_CONDA_EXE to conda/conda.bat path.")
}

message("[Step2] Running STCAT in conda env: stcat_env ...")
stcat_log <- file.path(cfg$step2_dir, "stcat_run.log")
stcat_cmd <- c(
  "run", "-n", "stcat_env",
  "python", py_script,
  "--input", stcat_input,
  "--output", pred_path
)
# Write both standard output and standard error to the same log.
stcat_code <- system2(conda_bin, args = stcat_cmd,
                      stdout = stcat_log, stderr = stcat_log)
if (is.null(stcat_code)) stcat_code <- 0
# Print the log to simplify troubleshooting.
if (file.exists(stcat_log)) {
  message("----- STCAT log (start) -----")
  message(paste(readLines(stcat_log), collapse = "\n"))
  message("----- STCAT log (end) -----")
}
if (stcat_code != 0 || !file.exists(pred_path)) {
  stop("STCAT failed (exit code ", stcat_code, "). See log above for Python traceback.")
}

message("[Step2] Reading STCAT predictions...")
stcat_pred <- read.csv(pred_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
if (!"Prediction" %in% colnames(stcat_pred)) {
  stop("Column 'Prediction' not found in STCAT output: ", pred_path)
}

pred_aligned <- stcat_pred[rownames(t_cells@meta.data), , drop = FALSE]
pred_label <- pred_aligned$Prediction

# Use the original STCAT labels without marker-based remapping.
pred_label[is.na(pred_label) | pred_label == ""] <- "Unassigned"
t_cells@meta.data$T_subtype <- factor(pred_label, levels = sort(unique(pred_label)))

cluster_anno_tbl <- as.data.frame.matrix(table(
  Cluster = t_cells$seurat_clusters,
  T_subtype = t_cells$T_subtype
))
cluster_anno_tbl$cluster <- rownames(cluster_anno_tbl)
cluster_anno_tbl <- cluster_anno_tbl[, c("cluster", setdiff(colnames(cluster_anno_tbl), "cluster"))]
write.csv(cluster_anno_tbl,
          file.path(cfg$step2_dir, "Tcell_cluster_subtype_annotation.csv"),
          row.names = FALSE)

if ("Uncertainty score" %in% colnames(pred_aligned)) {
  t_cells@meta.data$STCAT_Uncertainty <- as.numeric(pred_aligned$`Uncertainty score`)
}
message("[Step2] STCAT annotation completed (using raw STCAT subtype labels).")

combined@meta.data$T_subtype <- "Non_T"
common_cells <- intersect(rownames(combined@meta.data), rownames(t_cells@meta.data))
combined@meta.data[common_cells, "T_subtype"] <-
  as.character(t_cells@meta.data[common_cells, "T_subtype"])
message("[Step2] T_subtype distribution in combined (Non_T = non-T-cell populations):")
print(table(combined@meta.data$T_subtype))

p_t <- DimPlot(
  t_cells,
  reduction = "umap",
  group.by = "T_subtype",
  label = TRUE,
  repel = TRUE
) +
  theme_classic(base_size = 11, base_family = "sans")
save_plot_pdf_png(p_t, "Tcell_Subcluster_UMAP", cfg$step2_dir, width = 9, height = 7)

t_sub_counts <- as.data.frame(table(t_cells$T_subtype), stringsAsFactors = FALSE)
colnames(t_sub_counts) <- c("T_subtype", "n_cells")
t_sub_counts <- t_sub_counts[order(t_sub_counts$n_cells, decreasing = TRUE), ]
t_sub_counts$T_subtype <- factor(t_sub_counts$T_subtype, levels = t_sub_counts$T_subtype)

p_bar_t <- ggplot(t_sub_counts, aes(x = T_subtype, y = n_cells)) +
  geom_col(fill = "#4C78A8", width = 0.75) +
  labs(x = "STCAT subtype", y = "Cell count") +
  theme_bw(base_size = 11, base_family = "sans") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot_pdf_png(p_bar_t, "Tcell_STCAT_Subtype_Counts", cfg$step2_dir, width = 10, height = 5)

# ---- Core target genes vs major cell types DotPlot ----
hub_genes_cfg <- read_hub_genes(cfg)
genes_dot <- rownames(combined)[toupper(rownames(combined)) %in% hub_genes_cfg]
genes_dot <- head(genes_dot, 30)
if (length(genes_dot) > 0) {
  Idents(combined) <- "CellType"
  p_dot_core <- DotPlot(combined, features = genes_dot, group.by = "CellType",
                        dot.scale = 5, col.min = -2, col.max = 2) +
    RotatedAxis() +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 10, base_family = "sans") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  save_plot_pdf_png(p_dot_core, "CoreTarget_Gene_DotPlot_by_CellType",
                    cfg$step2_dir,
                    width  = min(4 + length(genes_dot) * 0.35, 18),
                    height = 6)
}

# ---- Core target genes vs T-cell subtypes DotPlot ----
genes_dot_t <- rownames(t_cells)[toupper(rownames(t_cells)) %in% hub_genes_cfg]
genes_dot_t <- head(genes_dot_t, 30)
if (length(genes_dot_t) > 0) {
  Idents(t_cells) <- "T_subtype"
  p_dot_core_t <- DotPlot(t_cells, features = genes_dot_t, group.by = "T_subtype",
                           dot.scale = 5, col.min = -2, col.max = 2) +
    RotatedAxis() +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 10, base_family = "sans") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  save_plot_pdf_png(p_dot_core_t, "CoreTarget_Gene_DotPlot_by_Tsubtype",
                    cfg$step2_dir,
                    width  = min(4 + length(genes_dot_t) * 0.35, 18),
                    height = 6)
}

saveRDS(combined, file.path(cfg$step2_dir, "step2_annotated_combined.rds"))
saveRDS(t_cells,  file.path(cfg$step2_dir, "step2_tcells_subclustered.rds"))

write.csv(
  as.data.frame(table(combined$CellType, combined$Phase)),
  file.path(cfg$step2_dir, "celltype_phase_counts.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(combined$CellType, combined$Disease_Group)),
  file.path(cfg$step2_dir, "celltype_disease_group_counts.csv"),
  row.names = FALSE
)

if ("Disease_State" %in% colnames(combined@meta.data)) {
  write.csv(
    as.data.frame(table(combined$CellType, combined$Disease_State)),
    file.path(cfg$step2_dir, "celltype_disease_state_counts.csv"),
    row.names = FALSE
  )
}

write.csv(
  as.data.frame(table(t_cells$T_subtype, t_cells$Phase)),
  file.path(cfg$step2_dir, "tcell_subtype_phase_counts.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(t_cells$T_subtype, t_cells$Disease_Group)),
  file.path(cfg$step2_dir, "tcell_subtype_disease_group_counts.csv"),
  row.names = FALSE
)

if ("Disease_State" %in% colnames(t_cells@meta.data)) {
  write.csv(
    as.data.frame(table(t_cells$T_subtype, t_cells$Disease_State)),
    file.path(cfg$step2_dir, "tcell_subtype_disease_state_counts.csv"),
    row.names = FALSE
  )
}

message("[Step2] Running FindAllMarkers (this may take a few minutes)...")
Idents(combined) <- "CellType"
all_markers <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

top50_markers <- all_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 50) %>%
  dplyr::ungroup()

write.csv(top50_markers, file.path(cfg$step2_dir, "Top50_markers_per_celltype.csv"), row.names = FALSE)

meta_export <- combined@meta.data
meta_export$Barcode <- rownames(meta_export)
write.csv(
  meta_export[, intersect(c("Barcode", "Sample", "Patient_ID", "Severity", "Phase",
                             "Disease_State", "Disease_Group", "GEO_Accession",
                             "seurat_clusters", "CellType", "T_subtype"),
                           colnames(meta_export))],
  file.path(cfg$step2_dir, "cell_metadata_annotated.csv"),
  row.names = FALSE
)

message("[Step2] Completed. Outputs saved in results/.")
