suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

source("scripts/utils.R")

message("[Step4] Loading AUCell object...")
combined <- readRDS(file.path(cfg$step3_dir, "step3_annotated_with_aucell.rds"))

score_name <- cfg$aucell_score_name
group_var <- cfg$comparison$group_var
group_order <- cfg$comparison$group_order

if (!(score_name %in% colnames(combined@meta.data))) {
  stop("AUCell score not found in metadata: ", score_name)
}
if (!(group_var %in% colnames(combined@meta.data))) {
  stop("Group variable not found in metadata: ", group_var)
}

meta <- combined@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(!is.na(.data[[score_name]]), !is.na(.data[[group_var]]))

meta[[group_var]] <- factor(as.character(meta[[group_var]]), levels = group_order)

summarise_score <- function(df, by_cols) {
  df %>%
    group_by(across(all_of(by_cols))) %>%
    summarise(
      n_cells = n(),
      mean = mean(.data[[score_name]], na.rm = TRUE),
      median = median(.data[[score_name]], na.rm = TRUE),
      q25 = quantile(.data[[score_name]], 0.25, na.rm = TRUE),
      q75 = quantile(.data[[score_name]], 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

test_by_label <- function(df, label_col) {
  df %>%
    group_by(.data[[label_col]]) %>%
    summarise(
      n_groups = n_distinct(.data[[group_var]]),
      p_value = tryCatch({
        if (n_distinct(.data[[group_var]]) < 2) {
          NA_real_
        } else {
          kruskal.test(.data[[score_name]] ~ .data[[group_var]])$p.value
        }
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))
}

plot_grouped_box <- function(df, label_col, title, max_labels = 30) {
  counts <- df %>%
    count(.data[[label_col]], sort = TRUE)
  keep <- head(counts[[label_col]], max_labels)
  plot_df <- df %>%
    filter(.data[[label_col]] %in% keep)
  label_levels <- counts[[label_col]][counts[[label_col]] %in% keep]
  plot_df[[label_col]] <- factor(as.character(plot_df[[label_col]]), levels = label_levels)

  ggplot(plot_df, aes(x = .data[[label_col]], y = .data[[score_name]], fill = .data[[group_var]])) +
    geom_boxplot(position = position_dodge(width = 0.72), width = 0.55,
                 outlier.shape = NA, alpha = 0.82) +
    geom_point(aes(color = .data[[group_var]]),
               position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.72),
               size = 0.22, alpha = 0.22, show.legend = FALSE) +
    labs(title = title, x = NULL, y = score_name, fill = NULL) +
    theme_bw(base_size = 10, base_family = "sans") +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
      panel.grid.minor = element_blank()
    )
}

cell_df <- meta %>%
  filter(!is.na(CellType))

write.csv(
  summarise_score(cell_df, c("CellType", group_var)),
  file.path(cfg$step4_dir, "AUCell_summary_CellType_by_DiseaseGroup.csv"),
  row.names = FALSE
)
write.csv(
  test_by_label(cell_df, "CellType"),
  file.path(cfg$step4_dir, "AUCell_kruskal_CellType_by_DiseaseGroup.csv"),
  row.names = FALSE
)

p_cell <- plot_grouped_box(cell_df, "CellType", "Core target AUCell by major immune cell type")
save_plot_pdf_png(p_cell, "AUCell_Boxplot_CellType_DiseaseGroup", cfg$step4_dir, width = 10, height = 5)

exclude_t <- cfg$comparison$t_exclude_labels
t_df <- meta %>%
  filter(CellType == "T_cells",
         !is.na(T_subtype),
         !T_subtype %in% exclude_t)

if (nrow(t_df) > 0) {
  write.csv(
    summarise_score(t_df, c("T_subtype", group_var)),
    file.path(cfg$step4_dir, "AUCell_summary_Tsubtype_by_DiseaseGroup.csv"),
    row.names = FALSE
  )
  write.csv(
    test_by_label(t_df, "T_subtype"),
    file.path(cfg$step4_dir, "AUCell_kruskal_Tsubtype_by_DiseaseGroup.csv"),
    row.names = FALSE
  )

  p_t <- plot_grouped_box(t_df, "T_subtype", "Core target AUCell by T cell subtype")
  save_plot_pdf_png(p_t, "AUCell_Boxplot_Tsubtype_DiseaseGroup", cfg$step4_dir, width = 13, height = 5.5)
}

meta_full <- combined@meta.data
meta_full$Barcode <- rownames(meta_full)
cols_keep <- c("Barcode", "Sample", "GEO_Accession", "Patient_ID", "Severity", "Phase",
               "Disease_Group", "seurat_clusters", "CellType", "T_subtype", score_name)
cols_keep <- intersect(cols_keep, colnames(meta_full))
write.csv(
  meta_full[, cols_keep],
  file.path(cfg$step4_dir, "cell_metadata_full_with_AUC.csv"),
  row.names = FALSE
)

message("[Step4] Completed. Disease-group summaries saved.")
