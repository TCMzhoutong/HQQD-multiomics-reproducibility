suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

source("scripts/utils.R")

message("[Step1] Loading 10X samples...")
seu_list <- lapply(seq_len(nrow(cfg$sample_info)), function(i) {
  smp <- cfg$sample_info[i, ]
  data_dir <- file.path(cfg$raw_data_dir, smp$sample_name)

  counts <- Read10X(data.dir = data_dir)
  seu <- CreateSeuratObject(counts = counts, project = smp$sample_name, min.cells = 3, min.features = 100)

  seu$Sample <- smp$sample_name
  seu$GEO_Accession <- smp$geo_accession
  seu$Patient_ID <- smp$patient_id
  seu$Severity <- smp$severity
  seu$Phase <- smp$phase
  seu$Disease_State <- smp$disease_state
  seu$Disease_Group <- smp$disease_group

  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu
})

message("[Step1] Merge samples...")
combined <- Reduce(function(x, y) merge(x, y), seu_list)

p_qc_pre <- VlnPlot(
  combined,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "Sample",
  pt.size = 0
) + plot_annotation(title = "QC metrics BEFORE filtering",
                    theme = theme(plot.title = element_text(size = 11)))
save_plot_pdf_png(p_qc_pre, "QC_violin_before_filter", cfg$step1_dir, width = 14, height = 5)

message("[Step1] QC filtering...")
qc <- cfg$qc
combined <- subset(
  combined,
  subset = nFeature_RNA >= qc$min_features &
    nFeature_RNA <= qc$max_features &
    percent.mt <= qc$max_percent_mt
)

message("[Step1] Normalization and feature selection...")
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- JoinLayers(combined, assay = "RNA")
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = cfg$seurat$npcs, verbose = FALSE)

top20 <- head(VariableFeatures(combined), 20)
p_hvg <- VariableFeaturePlot(combined)
p_hvg <- LabelPoints(plot = p_hvg, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)
save_plot_pdf_png(p_hvg, "HVG_VariableFeaturePlot", cfg$step1_dir, width = 9, height = 6)

# ---- ElbowPlot -------------------------------------------------------
p_elbow <- ElbowPlot(combined, ndims = cfg$seurat$npcs) +
  ggtitle("ElbowPlot: PC selection reference") +
  theme_classic(base_size = 11, base_family = "sans")
save_plot_pdf_png(p_elbow, "PCA_ElbowPlot", cfg$step1_dir, width = 7, height = 5)

p_pca <- DimPlot(combined, reduction = "pca", group.by = "Sample") +
  ggtitle("PCA (by Sample)") +
  theme_classic(base_size = 11, base_family = "sans")
save_plot_pdf_png(p_pca, "PCA_DimPlot_by_sample", cfg$step1_dir, width = 8, height = 6)

message("[Step1] Harmony batch correction by Patient_ID...")
combined <- RunHarmony(
  object = combined,
  group.by.vars = "Patient_ID",
  reduction = "pca",
  reduction.save = "harmony",
  max_iter = 20
)

dims_use <- cfg$seurat$dims_use
combined <- RunUMAP(combined, reduction = "harmony", dims = dims_use)
combined <- RunTSNE(combined, reduction = "harmony", dims = dims_use, check_duplicates = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = dims_use)
combined <- FindClusters(combined, resolution = cfg$seurat$cluster_resolution)

saveRDS(combined, file.path(cfg$step1_dir, "step1_harmony_combined.rds"))

p_qc_post <- VlnPlot(
  combined,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "Sample",
  pt.size = 0
) + plot_annotation(title = "QC metrics AFTER filtering",
                    theme = theme(plot.title = element_text(size = 11)))
save_plot_pdf_png(p_qc_post, "QC_violin_after_filter", cfg$step1_dir, width = 14, height = 5)

p_umap_clust <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Global UMAP – Cluster numbers (Harmony)") +
  theme_classic(base_size = 11, base_family = "sans")
save_plot_pdf_png(p_umap_clust, "Global_UMAP_preliminary_clusters", cfg$step1_dir, width = 8, height = 6)

p_tsne_clust <- DimPlot(combined, reduction = "tsne", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Global tSNE – Cluster numbers (Harmony)") +
  theme_classic(base_size = 11, base_family = "sans")
save_plot_pdf_png(p_tsne_clust, "Global_tSNE_preliminary_clusters", cfg$step1_dir, width = 8, height = 6)

message("[Step1] Completed. Output: results/rds/step1_harmony_combined.rds")
