suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

source("config/analysis_config.R")

cfg <- get_analysis_config()

dir.create(cfg$results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$step1_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$step2_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$step3_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$step4_dir,   showWarnings = FALSE, recursive = TRUE)

save_plot_pdf_png <- function(plot_obj, file_stem, out_dir, width = 8, height = 6, dpi = 320) {
  pdf_file <- file.path(out_dir, paste0(file_stem, ".pdf"))
  png_file <- file.path(out_dir, paste0(file_stem, ".png"))
  svg_file <- file.path(out_dir, paste0(file_stem, ".svg"))

  ggsave(pdf_file, plot = plot_obj, width = width, height = height, device = cairo_pdf)
  ggsave(png_file, plot = plot_obj, width = width, height = height, dpi = dpi)
  ggsave(svg_file, plot = plot_obj, width = width, height = height, device = grDevices::svg)

  invisible(list(pdf = pdf_file, png = png_file, svg = svg_file))
}

get_matrix_layer <- function(seu, assay = "RNA", layer = "counts") {
  GetAssayData(seu, assay = assay, layer = layer)
}

safe_features <- function(seu, genes) {
  genes[genes %in% rownames(seu)]
}

assign_by_max_score <- function(df, prefix, labels, unknown_label = "Unknown") {
  score_cols <- paste0(prefix, labels)
  score_cols <- score_cols[score_cols %in% colnames(df)]

  if (length(score_cols) == 0) {
    return(rep(unknown_label, nrow(df)))
  }

  mat <- as.matrix(df[, score_cols, drop = FALSE])
  idx <- max.col(mat, ties.method = "first")
  best <- colnames(mat)[idx]
  best <- gsub(paste0("^", prefix), "", best)
  best <- gsub("1$", "", best)

  max_score <- apply(mat, 1, max)
  best[max_score <= 0] <- unknown_label
  best
}

read_hub_genes <- function(cfg) {
  hub <- read.csv(cfg$hub_gene_file, stringsAsFactors = FALSE)

  if (!(cfg$hub_gene_column %in% colnames(hub))) {
    stop("Hub gene column not found: ", cfg$hub_gene_column)
  }

  if (!is.null(cfg$hub_filter_column) && !is.null(cfg$hub_filter_value)) {
    if (!(cfg$hub_filter_column %in% colnames(hub))) {
      stop("Hub filter column not found: ", cfg$hub_filter_column)
    }
    hub <- hub[as.character(hub[[cfg$hub_filter_column]]) == cfg$hub_filter_value, , drop = FALSE]
  }

  raw_vals <- as.character(na.omit(hub[[cfg$hub_gene_column]]))

  if (!is.null(cfg$hub_gene_split) && nzchar(cfg$hub_gene_split)) {
    parts <- strsplit(raw_vals, split = cfg$hub_gene_split, fixed = TRUE)
    raw_vals <- unlist(parts, use.names = FALSE)
  }

  genes <- unique(toupper(trimws(raw_vals)))
  genes[genes != "" & !is.na(genes)]
}
