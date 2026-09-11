################################################################################
# KEGG over-representation analysis for HQQD-reversed metabolites
# Input:
#   - 04a_neg_merged_revised/*Control_vs_Model* and *Model_vs_HQQD*
#   - 04b_pos_merged_revised/*Control_vs_Model* and *Model_vs_HQQD*
#   - meta_kegg_anno_neg.csv / meta_kegg_anno_pos.csv
# Output:
#   - 09b_KEGG_enrichment_reverse_revised/
#
# Statistical definition:
#   Query set    = HQQD-reversed metabolites.
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

out_dir <- "09b_KEGG_enrichment_reverse_revised"
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

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES", "Y")
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
    )

  lookup <- background_raw %>% distinct(ID, .keep_all = TRUE)

  term2gene <- background_raw %>%
    filter(!is.na(KEGG_pathway_annotation)) %>%
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
  list(lookup = lookup, term2gene = term2gene)
}

collect_hqqd_reverse_candidates <- function() {
  collect_mode <- function(mode) {
    cfg <- mode_config[[mode]]
    cm_file <- file.path(cfg$merge_dir, paste0(cfg$file_prefix, "_Control_vs_Model_merged.csv"))
    mh_file <- file.path(cfg$merge_dir, paste0(cfg$file_prefix, "_Model_vs_HQQD_merged.csv"))
    if (!file.exists(cm_file) || !file.exists(mh_file)) {
      warning(sprintf("Missing Control_vs_Model or Model_vs_HQQD merged file for %s", mode))
      return(NULL)
    }

    cm <- read_csv_mode(cm_file)
    mh <- read_csv_mode(mh_file)

    normalise_types <- function(d) {
      numeric_cols <- intersect(
        c("QC_CV", "VIP", "Q2_cum", "pQ2", "FC", "Log2FC", "p_value", "p_adj", "t_stat",
          grep("^Mean_", colnames(d), value = TRUE)),
        colnames(d)
      )
      d[numeric_cols] <- lapply(d[numeric_cols], as.numeric)
      if ("VIP_available" %in% colnames(d)) d$VIP_available <- as_bool(d$VIP_available)
      if ("VIP_valid" %in% colnames(d)) d$VIP_valid <- as_bool(d$VIP_valid)
      d
    }
    cm <- normalise_types(cm)
    mh <- normalise_types(mh)

    cm_pref <- cm %>% rename_with(~ paste0("CM_", .x), .cols = -ID)
    mh_pref <- mh %>% rename_with(~ paste0("MH_", .x), .cols = -ID)
    joined <- inner_join(cm_pref, mh_pref, by = "ID")
    joined$Ion_Mode <- cfg$ion_mode
    joined
  }

  joined <- bind_rows(collect_mode("NEG"), collect_mode("POS"))
  if (nrow(joined) == 0) stop("No evaluable HQQD reversal rows.")

  joined <- annotate_metabolite_exclusions(joined)

  mean_control_col <- "CM_Mean_Control"
  mean_model_col <- "CM_Mean_Model"
  mean_hqqd_col <- "MH_Mean_HQQD"

  joined <- joined %>%
    mutate(
      across(
        any_of(c("CM_QC_CV", "MH_QC_CV", "CM_VIP", "MH_VIP", "CM_FC", "MH_FC",
                 "CM_Log2FC", "MH_Log2FC", "CM_p_value", "MH_p_value", "CM_p_adj", "MH_p_adj")),
        ~ as.numeric(.x)
      ),
      QC_CV = coalesce(CM_QC_CV, MH_QC_CV),
      Pass_QC = is.na(QC_CV) | QC_CV <= analysis_cutoffs$qc_cv,
      CM_VIP_required = CM_VIP_available & CM_VIP_valid,
      CM_Pass_VIP = ifelse(CM_VIP_required, CM_VIP > analysis_cutoffs$vip, TRUE),
      CM_Pass = Pass_QC &
        abs(CM_Log2FC) >= log2(analysis_cutoffs$model_fc) &
        CM_p_value < analysis_cutoffs$p &
        CM_Pass_VIP,
      MH_Pass = Pass_QC &
        abs(MH_Log2FC) >= log2(analysis_cutoffs$treatment_fc) &
        MH_p_value < analysis_cutoffs$p,
      Opposite_Direction = CM_Log2FC * MH_Log2FC < 0,
      Control_log2_mean = log2(pmax(.data[[mean_control_col]], 1e-12)),
      Model_log2_mean = log2(pmax(.data[[mean_model_col]], 1e-12)),
      HQQD_log2_mean = log2(pmax(.data[[mean_hqqd_col]], 1e-12)),
      RecoveryScore = 1 - abs(HQQD_log2_mean - Control_log2_mean) / abs(Model_log2_mean - Control_log2_mean),
      RecoveryScore = ifelse(is.finite(RecoveryScore), RecoveryScore, NA_real_),
      Recovered_Toward_Control = !is.na(RecoveryScore) & RecoveryScore > 0,
    Reverse_Type = case_when(
        CM_Log2FC > 0 & MH_Log2FC < 0 ~ "Model_Up_HQQD_Down",
        CM_Log2FC < 0 & MH_Log2FC > 0 ~ "Model_Down_HQQD_Up",
        TRUE ~ "Not_Opposite"
      ),
      ReverseCandidate_Raw = CM_Pass & MH_Pass & Opposite_Direction & Recovered_Toward_Control,
      ReverseCandidate = ReverseCandidate_Raw & !Excluded_Annotation,
      ReverseConfidence = -log10(pmax(CM_p_value, 1e-300)) +
        -log10(pmax(MH_p_value, 1e-300)) +
        pmax(RecoveryScore, 0)
    ) %>%
    arrange(desc(ReverseCandidate), desc(RecoveryScore), CM_p_value, MH_p_value)

  list(
    all = joined,
    candidates = joined %>% filter(ReverseCandidate),
    excluded = joined %>% filter(ReverseCandidate_Raw, Excluded_Annotation)
  )
}

empty_ora_result <- function(label) {
  data.frame(
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

run_ora <- function(query_ids, label, term2gene, min_query = min_annotated_query) {
  universe <- unique(term2gene$ID)
  query <- intersect(unique(query_ids[!is.na(query_ids) & query_ids != ""]), universe)
  n_universe <- length(universe)
  n_query <- length(query)

  if (n_query < min_query || n_universe == 0) {
    return(empty_ora_result(label))
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
      Set = label,
      GeneRatio = paste0(Count, "/", n_query),
      BgRatio = paste0(BgCount, "/", n_universe),
      RichFactor = Count / BgCount,
      FoldEnrichment = (Count / n_query) / (BgCount / n_universe),
      pvalue = phyper(Count - 1, BgCount, n_universe - BgCount, n_query, lower.tail = FALSE),
      p.adjust = p.adjust(pvalue, method = "BH"),
      geneID = vapply(QueryIDs, paste, character(1), collapse = "/"),
      metabolite_name = QueryNames
    ) %>%
    select(Set, Pathway_ID, Description, Count, BgCount, GeneRatio, BgRatio,
           RichFactor, FoldEnrichment, pvalue, p.adjust, geneID, metabolite_name) %>%
    arrange(p.adjust, pvalue, desc(Count))
}

plot_ora <- function(res, label) {
  if (is.null(res) || nrow(res) == 0) return(invisible(FALSE))
  clean_label <- str_replace_all(label, "_", " ")
  subtitle_text <- str_wrap(
    "Query: HQQD-reversed metabolites; background: all KEGG-annotated detected metabolites",
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
      title = paste0("KEGG ORA: ", clean_label),
      subtitle = paste(subtitle_text, "Only BH-FDR < 0.05 pathways are shown.", sep = "\n")
    ) +
    theme_bw(base_size = 12, base_family = "sans") +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  ggsave(file.path(out_dir, paste0(label, "_KEGG_ORA.pdf")), p, width = 9.2, height = max(4.5, 0.34 * nrow(plot_df) + 2))
  ggsave(file.path(out_dir, paste0(label, "_KEGG_ORA.png")), p, width = 9.2, height = max(4.5, 0.34 * nrow(plot_df) + 2), dpi = 300)
  invisible(TRUE)
}

reset_output_dir(out_dir)

bg <- prepare_kegg_background(out_dir)
term2gene <- bg$term2gene
universe <- unique(term2gene$ID)

reverse_data <- collect_hqqd_reverse_candidates()
reverse_all <- reverse_data$all
candidates <- reverse_data$candidates

write_csv_utf8(reverse_all, file.path(out_dir, "HQQD_reverse_all_evaluable_for_KEGG.csv"))
write_csv_utf8(candidates, file.path(out_dir, "HQQD_reverse_candidates_for_KEGG.csv"))
write_csv_utf8(reverse_data$excluded, file.path(out_dir, "HQQD_reverse_excluded_annotations_for_KEGG.csv"))

candidate_kegg <- candidates %>%
  left_join(bg$lookup %>% select(ID, KEGG_annotation, KEGG_pathway_annotation), by = "ID") %>%
  mutate(
    KEGG_annotation = na_if(KEGG_annotation, ""),
    KEGG_annotation = na_if(KEGG_annotation, "--"),
    KEGG_pathway_annotation = na_if(KEGG_pathway_annotation, ""),
    KEGG_pathway_annotation = na_if(KEGG_pathway_annotation, "--")
  )
write_csv_utf8(candidate_kegg, file.path(out_dir, "HQQD_reverse_candidates_with_KEGG.csv"))

summary_reverse <- reverse_all %>%
  summarise(
    Evaluable = n(),
    Pass_CM = sum(CM_Pass, na.rm = TRUE),
    Pass_MH = sum(MH_Pass, na.rm = TRUE),
    Opposite = sum(Opposite_Direction, na.rm = TRUE),
    Recovered = sum(Recovered_Toward_Control, na.rm = TRUE),
    ReverseCandidates = sum(ReverseCandidate, na.rm = TRUE),
    KEGGAnnotatedReverseCandidates = length(intersect(candidates$ID, universe))
  )
write_csv_utf8(summary_reverse, file.path(out_dir, "HQQD_reverse_candidate_summary.csv"))

query_sets <- list(
  All_reverse_candidates = candidates$ID,
  Model_Up_HQQD_Down = candidates$ID[candidates$Reverse_Type == "Model_Up_HQQD_Down"],
  Model_Down_HQQD_Up = candidates$ID[candidates$Reverse_Type == "Model_Down_HQQD_Up"]
)

results <- list()
summary_rows <- list()

for (label in names(query_sets)) {
  ids <- unique(query_sets[[label]])
  res <- run_ora(ids, label, term2gene)
  write_csv_utf8(res, file.path(out_dir, paste0(label, "_KEGG_ORA.csv")))
  plot_ora(res, label)
  results[[label]] <- res

  summary_rows[[label]] <- data.frame(
    Set = label,
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
}

ora_summary <- bind_rows(summary_rows)
write_csv_utf8(ora_summary, file.path(out_dir, "HQQD_reverse_KEGG_ORA_summary.csv"))

all_results <- bind_rows(results)
if (nrow(all_results) > 0) {
  write_csv_utf8(all_results, file.path(out_dir, "HQQD_reverse_KEGG_ORA_all_results.csv"))
}

writeLines(c(
  "HQQD reverse KEGG ORA method note",
  "Query set: HQQD-reversed metabolites defined by Control-vs-Model abnormality, Model-vs-HQQD opposite direction, p < 0.05, fold-change threshold, QC filter, and recovery toward Control.",
  "VIP handling: Control-vs-Model VIP is used only when the corresponding OPLS-DA model is validated; Model-vs-HQQD VIP is not used because its OPLS-DA validation failed.",
  "Background: all detected metabolites with KEGG pathway annotation in NEG and POS modes.",
  "Test: hypergeometric over-representation analysis.",
  "Multiple testing correction: Benjamini-Hochberg FDR.",
  paste0("Sets with fewer than ", min_annotated_query, " KEGG-annotated query metabolites are not tested to avoid unstable ORA."),
  "If no pathway remains significant after FDR correction, report the results as exploratory enrichment rather than definitive significant enrichment."
), con = file.path(out_dir, "HQQD_reverse_KEGG_ORA_method_note.txt"))

cat(sprintf(
  "HQQD reverse KEGG ORA complete: %d candidates, %d KEGG-annotated candidates. Results saved to %s\n",
  nrow(candidates),
  length(intersect(candidates$ID, universe)),
  out_dir
))
