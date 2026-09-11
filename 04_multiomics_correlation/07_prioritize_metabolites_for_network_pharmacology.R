# ============================================================================
# Step 07: Prioritize serum metabolites for network pharmacology
#
# Criteria:
#   1. Gut microbiota Spearman association: |rho| >= 0.6 and P < 0.05.
#   2. Disease-module Mantel association: |r| >= 0.4 and P < 0.05.
#   3. Candidate metabolites satisfy both criteria.
# ============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

output_dir <- "network_pharmacology_candidates"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

write_csv_utf8 <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

short_label <- function(x, max_chars = 34) {
  x <- as.character(x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

read_matrix_csv <- function(path) {
  x <- readr::read_csv(path, show_col_types = FALSE, name_repair = "minimal")
  row_ids <- as.character(x[[1]])
  x <- as.data.frame(x[-1], check.names = FALSE)
  rownames(x) <- row_ids
  as.matrix(x)
}

metabolites <- readr::read_csv(
  "04_differential_species_metabolites/differential_metabolites_abundance_Model_vs_HQQD.csv",
  show_col_types = FALSE,
  name_repair = "minimal"
) %>%
  transmute(
    ID,
    Name,
    Display_Name = Name,
    Short_Label = make.unique(short_label(Name)),
    Ion_Mode,
    Log2FC = as.numeric(Log2FC)
  )

rho <- read_matrix_csv("05_correlation_analysis/spearman_correlation_matrix.csv")
pval <- read_matrix_csv("05_correlation_analysis/spearman_pvalue_matrix.csv")
if (!identical(dim(rho), dim(pval)) || !identical(dimnames(rho), dimnames(pval))) {
  stop("Spearman correlation and P-value matrices are not aligned.")
}

spearman_edges <- as.data.frame(as.table(rho), stringsAsFactors = FALSE) %>%
  rename(Taxon = Var1, Short_Label = Var2, Spearman_rho = Freq) %>%
  mutate(
    Spearman_P = as.vector(pval),
    Abs_Spearman_rho = abs(Spearman_rho),
    Spearman_direction = ifelse(Spearman_rho >= 0, "Positive", "Negative")
  ) %>%
  filter(Abs_Spearman_rho >= 0.6, Spearman_P < 0.05) %>%
  left_join(metabolites, by = "Short_Label") %>%
  select(ID, Name, Display_Name, Short_Label, Ion_Mode, Log2FC,
         Taxon, Spearman_rho, Abs_Spearman_rho, Spearman_P, Spearman_direction)

mantel_edges <- readr::read_csv(
  "06_mantel_network_analysis/metabolites_mantel_results.csv",
  show_col_types = FALSE,
  name_repair = "minimal"
) %>%
  transmute(
    Short_Label = env,
    Disease_Module = spec,
    Mantel_r = as.numeric(r),
    Abs_Mantel_r = abs(as.numeric(r)),
    Mantel_P = as.numeric(p),
    Mantel_direction = ifelse(as.numeric(r) >= 0, "Positive", "Negative"),
    rd, pd, link_type
  ) %>%
  filter(Abs_Mantel_r >= 0.4, Mantel_P < 0.05) %>%
  left_join(metabolites, by = "Short_Label") %>%
  select(ID, Name, Display_Name, Short_Label, Ion_Mode, Log2FC,
         Disease_Module, Mantel_r, Abs_Mantel_r, Mantel_P,
         Mantel_direction, rd, pd, link_type)

spearman_metabolites <- spearman_edges %>% distinct(ID, .keep_all = TRUE) %>%
  select(ID, Name, Display_Name, Short_Label, Ion_Mode, Log2FC)
mantel_metabolites <- mantel_edges %>% distinct(ID, .keep_all = TRUE) %>%
  select(ID, Name, Display_Name, Short_Label, Ion_Mode, Log2FC)
candidate_ids <- intersect(spearman_metabolites$ID, mantel_metabolites$ID)

spearman_summary <- spearman_edges %>%
  filter(ID %in% candidate_ids) %>%
  group_by(ID) %>%
  summarise(
    Strong_Spearman_Taxa_Count = n_distinct(Taxon),
    Positive_Taxa_Count = sum(Spearman_rho > 0),
    Negative_Taxa_Count = sum(Spearman_rho < 0),
    Max_Abs_Spearman_rho = max(Abs_Spearman_rho),
    Min_Spearman_P = min(Spearman_P),
    Top_Spearman_Taxon = Taxon[which.max(Abs_Spearman_rho)][1],
    Top_Spearman_rho = Spearman_rho[which.max(Abs_Spearman_rho)][1],
    Top_Spearman_P = Spearman_P[which.max(Abs_Spearman_rho)][1],
    Strong_Spearman_Taxa = paste(unique(Taxon), collapse = "; "),
    .groups = "drop"
  )

mantel_summary <- mantel_edges %>%
  filter(ID %in% candidate_ids) %>%
  group_by(ID) %>%
  summarise(
    Strong_Mantel_Module_Count = n_distinct(Disease_Module),
    Positive_Module_Count = sum(Mantel_r > 0),
    Negative_Module_Count = sum(Mantel_r < 0),
    Max_Abs_Mantel_r = max(Abs_Mantel_r),
    Min_Mantel_P = min(Mantel_P),
    Top_Mantel_Module = Disease_Module[which.max(Abs_Mantel_r)][1],
    Top_Mantel_r = Mantel_r[which.max(Abs_Mantel_r)][1],
    Top_Mantel_P = Mantel_P[which.max(Abs_Mantel_r)][1],
    Strong_Mantel_Modules = paste(unique(Disease_Module), collapse = "; "),
    .groups = "drop"
  )

candidates <- metabolites %>%
  filter(ID %in% candidate_ids) %>%
  left_join(spearman_summary, by = "ID") %>%
  left_join(mantel_summary, by = "ID") %>%
  mutate(
    Candidate_Rank_Score = Strong_Spearman_Taxa_Count + Strong_Mantel_Module_Count +
      Max_Abs_Spearman_rho + Max_Abs_Mantel_r
  ) %>%
  arrange(desc(Candidate_Rank_Score), desc(Max_Abs_Spearman_rho), desc(Max_Abs_Mantel_r))

spearman_support <- spearman_edges %>%
  filter(ID %in% candidate_ids) %>%
  mutate(Edge_Type = "Gut microbiota Spearman")
mantel_support <- mantel_edges %>%
  filter(ID %in% candidate_ids) %>%
  mutate(Edge_Type = "Disease phenotype Mantel")
supporting_edges <- bind_rows(spearman_support, mantel_support)

write_csv_utf8(spearman_edges, file.path(output_dir, "gut_microbiota_spearman_strong_edges.csv"))
write_csv_utf8(spearman_metabolites, file.path(output_dir, "gut_microbiota_spearman_strong_metabolites.csv"))
write_csv_utf8(mantel_edges, file.path(output_dir, "disease_phenotype_mantel_strong_edges.csv"))
write_csv_utf8(mantel_metabolites, file.path(output_dir, "disease_phenotype_mantel_strong_metabolites.csv"))
write_csv_utf8(candidates, file.path(output_dir, "candidate_metabolites_intersection.csv"))
write_csv_utf8(supporting_edges, file.path(output_dir, "candidate_metabolites_supporting_edges.csv"))

report <- c(
  "# Candidate serum metabolites associated with gut microbiota and phenotypes",
  "",
  "## Selection criteria",
  "",
  "- Strong taxon-metabolite Spearman association: `|rho| >= 0.6` and `P < 0.05`.",
  "- Strong metabolite-phenotype Mantel association: `|r| >= 0.4` and `P < 0.05`.",
  "- Candidate metabolites satisfy both criteria.",
  "",
  "## Selection summary",
  "",
  sprintf("- Metabolites with strong Spearman support: %d; supported taxon-metabolite edges: %d.", nrow(spearman_metabolites), nrow(spearman_edges)),
  sprintf("- Metabolites with strong Mantel support: %d; supported metabolite-phenotype edges: %d.", nrow(mantel_metabolites), nrow(mantel_edges)),
  sprintf("- Candidate serum metabolites satisfying both criteria: %d.", nrow(candidates)),
  "",
  "## Candidate metabolites",
  "",
  "| Rank | ID | Name | Ion | Log2FC | Spearman taxa n | Max abs rho | Mantel modules n | Max abs r | Top taxa | Top phenotype |",
  "|---:|---|---|---|---:|---:|---:|---:|---:|---|---|"
)

if (nrow(candidates) > 0) {
  report_rows <- vapply(seq_len(nrow(candidates)), function(i) {
    x <- candidates[i, ]
    sprintf(
      "| %d | %s | %s | %s | %.3f | %d | %.3f | %d | %.3f | %s | %s |",
      i, x$ID, x$Name, x$Ion_Mode, x$Log2FC,
      x$Strong_Spearman_Taxa_Count, x$Max_Abs_Spearman_rho,
      x$Strong_Mantel_Module_Count, x$Max_Abs_Mantel_r,
      x$Top_Spearman_Taxon, x$Top_Mantel_Module
    )
  }, character(1))
  report <- c(report, report_rows)
}

report <- c(
  report,
  "",
  "## Notes",
  "",
  "- This list is generated from the outputs of steps 05 and 06 using the prespecified thresholds and is used as the metabolite input for network pharmacology.",
  "- Ranking is used for result organisation; candidate eligibility is defined by the intersection of the two criteria."
)
writeLines(report, file.path(output_dir, "candidate_metabolites_intersection_report.md"), useBytes = TRUE)

cat(sprintf(
  "Prioritization complete: %d Spearman metabolites, %d Mantel metabolites, %d intersection candidates.\n",
  nrow(spearman_metabolites), nrow(mantel_metabolites), nrow(candidates)
))
