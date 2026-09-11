################################################################################
# KEGG over-representation analysis for pairwise differential metabolites
# Input:
#   - 06_merge_diff_revised/06_merged_*_diff.csv
#   - meta_kegg_anno_neg.csv / meta_kegg_anno_pos.csv
# Output:
#   - 09_KEGG_enrichment_revised/
#
# Statistical definition:
#   Query set    = differential metabolites for each comparison.
#   Background   = all detected metabolites with KEGG pathway annotation.
#   Test         = hypergeometric over-representation analysis.
#   Correction   = Benjamini-Hochberg FDR.
################################################################################

source("00_revised_config.R", encoding = "UTF-8")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

out_dir <- "09_KEGG_enrichment_revised"
min_annotated_query <- 3

reset_output_dir <- function(path) {
  if (dir.exists(path)) {
    unlink(list.files(path, full.names = TRUE), recursive = TRUE, force = TRUE)
  }
  ensure_dir(path)
}

integer_count_breaks <- function(x) {
  sort(unique(as.integer(x[is.finite(x)])))
}

read_kegg_anno <- function(path, ion_mode) {
  x <- read_csv_mode(path)
  colnames(x) <- gsub("^#", "", colnames(x))
  if (!"ID" %in% colnames(x)) colnames(x)[1] <- "ID"
  if (!"name" %in% colnames(x)) x$name <- x$ID
  if (!"KEGG_pathway_annotation" %in% colnames(x)) {
    stop("Missing KEGG_pathway_annotation in ", path)
  }
  if (!"KEGG_annotation" %in% colnames(x)) x$KEGG_annotation <- NA_character_

  x %>%
    transmute(
      ID = as.character(ID),
      Name = as.character(name),
      Ion_Mode = ion_mode,
      KEGG_annotation = as.character(KEGG_annotation),
      KEGG_pathway_annotation = as.character(KEGG_pathway_annotation)
    )
}

prepare_kegg_background <- function(output_dir = out_dir) {
  background_raw <- bind_rows(
    read_kegg_anno("meta_kegg_anno_neg.csv", "NEG"),
    read_kegg_anno("meta_kegg_anno_pos.csv", "POS")
  ) %>%
    mutate(
      KEGG_annotation = str_replace_all(KEGG_annotation, "\\s+", ""),
      KEGG_annotation = na_if(KEGG_annotation, ""),
      KEGG_annotation = na_if(KEGG_annotation, "--"),
      KEGG_pathway_annotation = na_if(str_trim(KEGG_pathway_annotation), ""),
      KEGG_pathway_annotation = na_if(KEGG_pathway_annotation, "--")
    ) %>%
    filter(!is.na(KEGG_pathway_annotation))

  term2gene <- background_raw %>%
    separate_rows(KEGG_pathway_annotation, sep = ";;|;|\\|") %>%
    mutate(
      KEGG_pathway_annotation = str_trim(KEGG_pathway_annotation),
      Pathway_ID = str_extract(KEGG_pathway_annotation, "ko[0-9]+"),
      Description = str_trim(str_remove(KEGG_pathway_annotation, "\\(ko[0-9]+\\)")),
      Description = ifelse(is.na(Description) | Description == "", Pathway_ID, Description)
    ) %>%
    filter(!is.na(Pathway_ID), Pathway_ID != "") %>%
    transmute(Pathway_ID, Description, ID, Name, Ion_Mode, KEGG_annotation) %>%
    distinct()

  write_csv_utf8(term2gene, file.path(output_dir, "KEGG_ORA_background_TERM2GENE.csv"))
  term2gene
}

empty_ora_result <- function(comp_name, set_name) {
  data.frame(
    Comparison = character(),
    Set = character(),
    Pathway_ID = character(),
    Description = character(),
    Count = integer(),
    BgCount = integer(),
    GeneRatio = character(),
    BgRatio = character(),
    RichFactor = numeric(),
    FoldEnrichment = numeric(),
    pvalue = numeric(),
    p.adjust = numeric(),
    geneID = character(),
    metabolite_name = character(),
    stringsAsFactors = FALSE
  )
}

run_ora <- function(query_ids, comp_name, set_name, term2gene, min_query = min_annotated_query) {
  universe <- unique(term2gene$ID)
  query <- intersect(unique(query_ids[!is.na(query_ids) & query_ids != ""]), universe)
  n_universe <- length(universe)
  n_query <- length(query)

  if (n_query < min_query || n_universe == 0) {
    return(empty_ora_result(comp_name, set_name))
  }

  term2gene %>%
    group_by(Pathway_ID, Description) %>%
    summarise(
      BgCount = n_distinct(ID),
      QueryIDs = list(intersect(unique(ID), query)),
      .groups = "drop"
    ) %>%
    mutate(Count = lengths(QueryIDs)) %>%
    filter(Count > 0) %>%
    rowwise() %>%
    mutate(
      QueryNames = paste(unique(term2gene$Name[term2gene$ID %in% unlist(QueryIDs)]), collapse = "/")
    ) %>%
    ungroup() %>%
    mutate(
      Comparison = comp_name,
      Set = set_name,
      GeneRatio = paste0(Count, "/", n_query),
      BgRatio = paste0(BgCount, "/", n_universe),
      RichFactor = Count / BgCount,
      FoldEnrichment = (Count / n_query) / (BgCount / n_universe),
      pvalue = phyper(Count - 1, BgCount, n_universe - BgCount, n_query, lower.tail = FALSE),
      p.adjust = p.adjust(pvalue, method = "BH"),
      geneID = vapply(QueryIDs, paste, character(1), collapse = "/"),
      metabolite_name = QueryNames
    ) %>%
    select(Comparison, Set, Pathway_ID, Description, Count, BgCount, GeneRatio, BgRatio,
           RichFactor, FoldEnrichment, pvalue, p.adjust, geneID, metabolite_name) %>%
    arrange(p.adjust, pvalue, desc(Count))
}

plot_ora <- function(res, comp_name, set_name, path_prefix) {
  if (nrow(res) == 0) return(invisible(FALSE))
  subtitle_text <- str_wrap(
    "Query: differential metabolites; background: all KEGG-annotated detected metabolites",
    width = 78
  )

  plot_df <- res %>%
    filter(!is.na(p.adjust), p.adjust < 0.05) %>%
    arrange(p.adjust, pvalue) %>%
    slice_head(n = 20) %>%
    mutate(
      Pathway_Label = str_wrap(paste0(Description, " (", Pathway_ID, ")"), width = 46),
      Pathway_Label = factor(Pathway_Label, levels = rev(unique(Pathway_Label))),
      GeneRatioValue = Count / as.numeric(str_extract(GeneRatio, "(?<=/)\\d+")),
      p_adjust_plot = pmin(pmax(p.adjust, 1e-16), 1)
    )

  if (nrow(plot_df) == 0) return(invisible(FALSE))

  p <- ggplot(plot_df, aes(GeneRatioValue, Pathway_Label)) +
    geom_point(aes(size = Count, color = p_adjust_plot), alpha = 0.9) +
    scale_color_gradient(
      low = "#D73027",
      high = "#4575B4",
      limits = c(0, 0.05),
      breaks = c(0.01, 0.03, 0.05),
      name = "p.adjust"
    ) +
    scale_size_continuous(
      range = c(2.2, 7),
      breaks = integer_count_breaks(plot_df$Count),
      labels = integer_count_breaks(plot_df$Count),
      name = "Count"
    ) +
    labs(
      x = "Gene ratio (Count / KEGG-annotated query metabolites)",
      y = NULL,
      title = sprintf("KEGG ORA: %s (%s)", comp_name, set_name),
      subtitle = paste(subtitle_text, "Only BH-FDR < 0.05 pathways are shown.", sep = "\n")
    ) +
    theme_bw(base_size = 12, base_family = "sans") +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  ggsave(paste0(path_prefix, ".pdf"), p, width = 9.2, height = max(4.5, 0.34 * nrow(plot_df) + 2))
  ggsave(paste0(path_prefix, ".png"), p, width = 9.2, height = max(4.5, 0.34 * nrow(plot_df) + 2), dpi = 300)
  invisible(TRUE)
}

reset_output_dir(out_dir)
term2gene <- prepare_kegg_background(out_dir)
universe <- unique(term2gene$ID)

diff_files <- list.files("06_merge_diff_revised", pattern = "^06_merged_.*_diff\\.csv$", full.names = FALSE)
if (length(diff_files) == 0) stop("No merged differential files found in 06_merge_diff_revised/")

summary_rows <- list()
all_results <- list()

for (f in diff_files) {
  comp_name <- str_remove(f, "^06_merged_")
  comp_name <- str_remove(comp_name, "_diff\\.csv$")
  message(sprintf("[KEGG ORA] %s", comp_name))

  comp_dir <- file.path(out_dir, comp_name)
  ensure_dir(comp_dir)

  diff_data <- read_csv_mode(file.path("06_merge_diff_revised", f))
  if (!"Sig_Status" %in% colnames(diff_data)) {
    diff_data$Sig_Status <- ifelse(as.numeric(diff_data$Log2FC) > 0, "Up", "Down")
  }

  query_sets <- list(
    All = diff_data$ID,
    Up = diff_data$ID[diff_data$Sig_Status == "Up"],
    Down = diff_data$ID[diff_data$Sig_Status == "Down"]
  )

  for (set_name in names(query_sets)) {
    ids <- unique(query_sets[[set_name]])
    res <- run_ora(ids, comp_name, set_name, term2gene)
    out_prefix <- file.path(comp_dir, paste0(comp_name, "_", set_name, "_KEGG_ORA"))
    write_csv_utf8(res, paste0(out_prefix, ".csv"))
    plot_ora(res, comp_name, set_name, out_prefix)

    summary_rows[[paste(comp_name, set_name, sep = "_")]] <- data.frame(
      Comparison = comp_name,
      Set = set_name,
      QueryCount = length(unique(ids[!is.na(ids) & ids != ""])),
      KEGGAnnotatedQueryCount = length(intersect(unique(ids), universe)),
      AnalysisStatus = ifelse(
        length(intersect(unique(ids), universe)) < min_annotated_query,
        paste0("Not tested: fewer than ", min_annotated_query, " KEGG-annotated query metabolites"),
        "Tested"
      ),
      EnrichedPathways = nrow(res),
      RawP_lt_0.05 = sum(res$pvalue < 0.05, na.rm = TRUE),
      FDR_lt_0.05 = sum(res$p.adjust < 0.05, na.rm = TRUE),
      FDR_lt_0.10 = sum(res$p.adjust < 0.10, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if (nrow(res) > 0) all_results[[paste(comp_name, set_name, sep = "_")]] <- res
  }
}

summary_df <- bind_rows(summary_rows)
write_csv_utf8(summary_df, file.path(out_dir, "KEGG_ORA_pairwise_summary.csv"))

all_res <- bind_rows(all_results)
if (nrow(all_res) > 0) {
  write_csv_utf8(all_res, file.path(out_dir, "KEGG_ORA_pairwise_all_results.csv"))

  overview_df <- all_res %>%
    filter(Set == "All", !is.na(p.adjust), p.adjust < 0.05) %>%
    group_by(Comparison) %>%
    arrange(p.adjust, pvalue, .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    mutate(
      Pathway_Label = str_wrap(paste0(Description, " (", Pathway_ID, ")"), width = 40),
      Comparison = factor(Comparison, levels = unique(Comparison)),
      GeneRatioValue = Count / as.numeric(str_extract(GeneRatio, "(?<=/)\\d+")),
      p_adjust_plot = pmin(pmax(p.adjust, 1e-16), 1)
    )

  if (nrow(overview_df) > 0) {
    p_overview <- ggplot(overview_df, aes(Comparison, Pathway_Label)) +
      geom_point(aes(size = Count, color = p_adjust_plot)) +
      scale_color_gradient(
        low = "#D73027",
        high = "#4575B4",
        limits = c(0, 0.05),
        breaks = c(0.01, 0.03, 0.05),
        name = "p.adjust"
      ) +
      scale_size_continuous(
        range = c(2.2, 7),
        breaks = integer_count_breaks(overview_df$Count),
        labels = integer_count_breaks(overview_df$Count),
        name = "Count"
      ) +
      labs(
        x = NULL,
        y = NULL,
        title = "Top KEGG ORA pathways across pairwise comparisons",
        subtitle = "Only BH-FDR < 0.05 pathways are shown"
      ) +
      theme_bw(base_size = 12, base_family = "sans") +
      theme(
        panel.grid.major = element_line(color = "grey90"),
        axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      )

    ggsave(file.path(out_dir, "KEGG_ORA_pairwise_overview.pdf"), p_overview, width = 9, height = 7)
    ggsave(file.path(out_dir, "KEGG_ORA_pairwise_overview.png"), p_overview, width = 9, height = 7, dpi = 300)
  }
}

writeLines(c(
  "KEGG ORA method note",
  "Query set: differential metabolites for each comparison.",
  "Background: all detected metabolites with KEGG pathway annotation in NEG and POS modes.",
  "Test: hypergeometric over-representation analysis.",
  "Multiple testing correction: Benjamini-Hochberg FDR.",
  paste0("Sets with fewer than ", min_annotated_query, " KEGG-annotated query metabolites are not tested to avoid unstable ORA."),
  "VIP is not used in KEGG ORA itself; it only affects upstream differential-metabolite selection when OPLS-DA is validated."
), con = file.path(out_dir, "KEGG_ORA_method_note.txt"))

cat(sprintf("Pairwise KEGG ORA complete. Results saved to %s\n", out_dir))
