# ============================================================================
# Step 06: Mantel Network Analysis
# Omics features vs disease indicator modules
#
# Mantel test and network heatmap use Hy4m/linkET source code:
# https://github.com/Hy4m/linkET
# The color palette is the reference palette for steps 05 and 06.
# ============================================================================

rm(list = ls())

required_packages <- c("dplyr", "ggplot2", "stringr", "tibble", "purrr", "magrittr", "vegan", "grid", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)
  library(vegan)
  library(grid)
  library(readr)
})

source_linket <- function() {
  if (requireNamespace("linkET", quietly = TRUE)) {
    suppressPackageStartupMessages(library(linkET))
    return("installed package")
  }
  stop("Missing required R package: linkET. Install Hy4m/linkET before running this script.")
}

save_plot_svg <- function(filename, plot, width, height) {
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = filename, width = width, height = height)
  } else {
    grDevices::svg(filename = filename, width = width, height = height, onefile = FALSE)
  }
  print(plot)
  grDevices::dev.off()
}

linket_source <- source_linket()

output_dir <- "06_mantel_network_analysis"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

read_csv_flexible <- function(path) {
  encodings <- c("UTF-8", "GB18030", "GBK")
  last_error <- NULL
  for (enc in encodings) {
    dat <- tryCatch(
      suppressWarnings(
        readr::read_csv(
          path,
          locale = readr::locale(encoding = enc),
          show_col_types = FALSE,
          name_repair = "minimal"
        ) %>%
          as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)
      ),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(dat) && ncol(dat) > 1) return(dat)
  }
  stop("Could not read CSV file: ", path, "; last error: ", conditionMessage(last_error))
}

write_csv_utf8 <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

short_label <- function(x, max_chars = 36) {
  x <- as.character(x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

make_sample_ids <- function(df) {
  target_groups <- c("Control", "Model", "HQQD-L")
  group_to_prefix <- c("Control" = "K", "Model" = "M", "HQQD-L" = "Z")
  df %>%
    filter(Group %in% target_groups) %>%
    group_by(Group) %>%
    mutate(.rep = row_number(), SampleID = paste0(group_to_prefix[as.character(Group)], .rep)) %>%
    ungroup() %>%
    select(-.rep)
}

clean_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  as.numeric(gsub("%", "", as.character(x)))
}

color_palette <- colorRampPalette(c(
  "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
  "#FFFFFF",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
))(100)

plot_mantel_network <- function(omics_mat, disease_df, disease_group_list, title_prefix, save_prefix) {
  if (ncol(omics_mat) == 0 || ncol(disease_df) == 0) return(NULL)

  complete_idx <- complete.cases(disease_df)
  omics_mat <- omics_mat[complete_idx, , drop = FALSE]
  disease_df <- disease_df[complete_idx, , drop = FALSE]

  keep_omics <- vapply(omics_mat, function(x) length(unique(x[!is.na(x)])) > 1, logical(1))
  keep_disease <- vapply(disease_df, function(x) length(unique(x[!is.na(x)])) > 1, logical(1))
  omics_mat <- omics_mat[, keep_omics, drop = FALSE]
  disease_df <- disease_df[, keep_disease, drop = FALSE]
  disease_group_list <- lapply(disease_group_list, function(cols) intersect(cols, colnames(disease_df)))
  disease_group_list <- disease_group_list[vapply(disease_group_list, length, integer(1)) > 0]

  if (nrow(disease_df) < 3 || ncol(omics_mat) == 0 || length(disease_group_list) == 0) return(NULL)
  is_metabolite_plot <- identical(save_prefix, "metabolites")

  set.seed(ifelse(save_prefix == "metabolites", 20260909, 20260908))
  mantel_res <- mantel_test(
    spec = disease_df,
    env = omics_mat,
    spec_select = disease_group_list,
    mantel_fun = "mantel",
    spec_dist = dist_func(.FUN = "dist", method = "euclidean"),
    env_dist = dist_func(.FUN = "dist", method = "euclidean"),
    permutations = 999
  )

  mantel_res <- mantel_res %>%
    mutate(
      rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
      pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")),
      link_type = ifelse(r >= 0, "Positive", "Negative")
    )

  write_csv_utf8(mantel_res, file.path(output_dir, paste0(save_prefix, "_mantel_results.csv")))

  # Correlation heatmap among omics features, with Mantel links to disease modules.
  p <- qcorrplot(correlate(omics_mat, method = "spearman"), type = "lower", diag = FALSE) +
    geom_square() +
    geom_couple(
      aes(colour = pd, size = rd),
      data = mantel_res,
      curvature = 0.15,
      nudge_x = 0.2,
      label.size = 3.6,
      label.family = "sans",
      label.fontface = 1
    ) +
    scale_fill_gradientn(colours = color_palette, limits = c(-1, 1), name = "Spearman rho", guide = guide_colorbar(display = "rectangles")) +
    scale_colour_manual(values = c("< 0.01" = "#B2182B", "0.01 - 0.05" = "#EF8A62", ">= 0.05" = "grey70"), name = "Mantel p") +
    scale_size_manual(values = c("< 0.2" = 0.3, "0.2 - 0.4" = 0.7, ">= 0.4" = 1.2), name = "Mantel r") +
    labs(title = title_prefix, x = NULL, y = NULL) +
    theme_minimal(base_family = "sans", base_size = 10) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = ifelse(is_metabolite_plot, 8, 10)),
      axis.text.y = element_text(color = "black", size = ifelse(is_metabolite_plot, 8, 10)),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10)
    )

  height <- max(6, min(16, 0.26 * ncol(omics_mat) + 3))
  ggsave(file.path(output_dir, paste0(save_prefix, "_mantel_network_heatmap.png")), p, width = 10, height = height, dpi = 300)
  ggsave(file.path(output_dir, paste0(save_prefix, "_mantel_network_heatmap.pdf")), p, width = 10, height = height)
  save_plot_svg(file.path(output_dir, paste0(save_prefix, "_mantel_network_heatmap.svg")), p, width = 10, height = height)

  mantel_res
}

cat("\n=== Step 06: Data Import ===\n")
cat("Using linkET from: ", linket_source, "\n", sep = "")

disease_file <- file.path("..", "01_animal_experiment", "data", "animal_outcomes.csv")
if (!file.exists(disease_file)) stop("Disease indicator input not found: ", disease_file)

disease_raw <- read_csv_flexible(disease_file)
disease_data <- make_sample_ids(disease_raw)

selected_cols <- c(
  "Serum TNF-α", "Serum IL-17A", "Serum IL-10", "Serum IL-6", "Serum IL-1β",
  "Ileum TNF-α", "Ileum IL-17A", "Ileum IL-10", "Ileum IL-6", "Ileum IL-1β",
  "Lung TNF-α", "Lung IL-17A", "Lung IL-10", "Lung IL-6", "Lung IL-1β",
  "Lung Bacterial Load (log10 CFU)", "Lung Pathology Score",
  "MLN_treg_CD25+FOXP3+(%)", "MLN_TH17_CD8-IL-17A+(%)",
  "Lung treg_CD25+FOXP3+(%)", "Lung TH17_CD8-IL-17A+(%)"
)

cols_to_keep <- intersect(selected_cols, colnames(disease_data))
disease_clean <- disease_data[, cols_to_keep, drop = FALSE] %>%
  mutate(across(everything(), clean_numeric)) %>%
  as.data.frame()
rownames(disease_clean) <- disease_data$SampleID

disease_group_list <- list(
  "Serum Cytokines" = c("Serum TNF-α", "Serum IL-17A", "Serum IL-10", "Serum IL-6", "Serum IL-1β"),
  "Ileum Cytokines" = c("Ileum TNF-α", "Ileum IL-17A", "Ileum IL-10", "Ileum IL-6", "Ileum IL-1β"),
  "Lung Cytokines" = c("Lung TNF-α", "Lung IL-17A", "Lung IL-10", "Lung IL-6", "Lung IL-1β"),
  "Lung Pathology" = c("Lung Bacterial Load (log10 CFU)", "Lung Pathology Score"),
  "Immune Cells" = c("MLN_treg_CD25+FOXP3+(%)", "MLN_TH17_CD8-IL-17A+(%)", "Lung treg_CD25+FOXP3+(%)", "Lung TH17_CD8-IL-17A+(%)")
)

disease_group_list <- lapply(disease_group_list, function(cols) intersect(cols, colnames(disease_clean)))
disease_group_list <- disease_group_list[vapply(disease_group_list, length, integer(1)) > 0]

taxa_raw <- read.csv(
  "04_differential_species_metabolites/differential_taxa_abundance_Model_vs_HQQD.csv",
  row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)
metab_raw <- read.csv(
  "04_differential_species_metabolites/differential_metabolites_abundance_Model_vs_HQQD.csv",
  check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)

taxa_mat <- taxa_raw %>% select(matches("^(M|Z|K)\\d+")) %>% as.matrix()
metab_name_col <- if ("Name" %in% colnames(metab_raw)) "Name" else if ("name" %in% colnames(metab_raw)) "name" else NA
metab_labels <- ifelse(is.na(metab_raw[[metab_name_col]]) | metab_raw[[metab_name_col]] == "", metab_raw$ID, metab_raw[[metab_name_col]])
metab_labels <- make.unique(short_label(metab_labels, 34))
metab_mat <- metab_raw %>% select(matches("^(M|Z|K)\\d+")) %>% as.matrix()
rownames(metab_mat) <- metab_labels

common_samples <- Reduce(intersect, list(rownames(disease_clean), colnames(taxa_mat), colnames(metab_mat)))
common_samples <- c(
  sort(common_samples[grepl("^K\\d+", common_samples)]),
  sort(common_samples[grepl("^M\\d+", common_samples)]),
  sort(common_samples[grepl("^Z\\d+", common_samples)])
)
if (length(common_samples) < 3) stop("Too few matched samples for Mantel analysis.")

disease_use <- disease_clean[common_samples, , drop = FALSE]
taxa_use <- as.data.frame(t(taxa_mat[, common_samples, drop = FALSE]))
metab_use <- as.data.frame(t(metab_mat[, common_samples, drop = FALSE]))

cat("Matched samples: ", length(common_samples), "\n", sep = "")
cat("Disease modules: ", paste(names(disease_group_list), collapse = ", "), "\n", sep = "")
cat("Taxa features: ", ncol(taxa_use), "; metabolite features: ", ncol(metab_use), "\n", sep = "")

cat("\n=== Mantel Tests and Network Heatmaps via linkET ===\n")

taxa_results <- plot_mantel_network(taxa_use, disease_use, disease_group_list, "Taxa vs Disease Indicator Modules", "taxa")
metab_results <- plot_mantel_network(metab_use, disease_use, disease_group_list, "Metabolites vs Disease Indicator Modules", "metabolites")
all_results <- bind_rows(
  if (!is.null(taxa_results)) mutate(taxa_results, Omics = "Taxa") else NULL,
  if (!is.null(metab_results)) mutate(metab_results, Omics = "Metabolites") else NULL
)

write_csv_utf8(all_results, file.path(output_dir, "combined_mantel_results.csv"))

cat("Taxa Mantel rows: ", ifelse(is.null(taxa_results), 0, nrow(taxa_results)), "\n", sep = "")
cat("Metabolite Mantel rows: ", ifelse(is.null(metab_results), 0, nrow(metab_results)), "\n", sep = "")
cat("Significant Mantel tests (p < 0.05): ", sum(all_results$p < 0.05, na.rm = TRUE), "\n", sep = "")
cat("Step 06 complete. Outputs saved to: ", output_dir, "\n", sep = "")
