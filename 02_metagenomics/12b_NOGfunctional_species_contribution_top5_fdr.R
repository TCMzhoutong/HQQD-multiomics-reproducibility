#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)

output_dir <- "12b_NOGfunctional_species_contribution_top5_fdr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

geneset_abundance <- read.table(
  "geneset.abundance.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  comment.char = "",
  check.names = FALSE
)

if (grepl("^#", colnames(geneset_abundance)[1])) {
  colnames(geneset_abundance)[1] <- gsub("^#", "", colnames(geneset_abundance)[1])
}

taxonomy_genes <- read.table(
  "nr.taxonomy.gene_list.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

eggnog_genes <- read.table(
  "eggNOG.NOG.gene_list.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

reverse_stats <- read.csv(
  file.path("08b_HQQD_reverse_eggNOG", "statistical_results_eggNOG_NOG_fdr.csv"),
  check.names = FALSE
) %>%
  arrange(desc(LDA_avg)) %>%
  slice_head(n = 5)

target_eggnog_ids <- reverse_stats$Taxon

target_metadata <- reverse_stats %>%
  select(
    eggNOG = Taxon,
    descrition,
    eggNOG_Class,
    Category,
    Trend_Model_vs_Control,
    Trend_HQQD_vs_Model,
    LDA_avg
  )

write.csv(
  target_metadata,
  file.path(output_dir, "target_eggNOG_NOG_top5.csv"),
  row.names = FALSE,
  quote = TRUE
)

taxonomy_genes$Genus <- vapply(strsplit(taxonomy_genes$Taxon, "; "), function(x) {
  if (length(x) >= 6) {
    genus <- gsub("^g__", "", x[6])
    trimws(genus)
  } else {
    "Unassigned"
  }
}, character(1))

target_gene_set <- eggnog_genes %>%
  filter(eggNOG %in% target_eggnog_ids) %>%
  pull(gene_list) %>%
  strsplit(";") %>%
  unlist(use.names = FALSE) %>%
  unique()
target_gene_set <- target_gene_set[target_gene_set %in% rownames(geneset_abundance)]

target_gene_genus_map <- bind_rows(lapply(seq_len(nrow(taxonomy_genes)), function(i) {
  genes <- unlist(strsplit(taxonomy_genes$gene_list[i], ";"), use.names = FALSE)
  genes <- intersect(genes, target_gene_set)
  if (length(genes) == 0) {
    return(NULL)
  }
  data.frame(Gene = genes, Genus = taxonomy_genes$Genus[i], stringsAsFactors = FALSE)
})) %>%
  distinct(Gene, .keep_all = TRUE)

message("Mapped target genes to genus: ", nrow(target_gene_genus_map))

calculate_contribution <- function(eggnog_id, eggnog_genes_df, gene_genus_map, abundance_df) {
  eggnog_row <- eggnog_genes_df[eggnog_genes_df$eggNOG == eggnog_id, , drop = FALSE]
  if (nrow(eggnog_row) == 0) {
    warning("No eggNOG gene list found for: ", eggnog_id)
    return(NULL)
  }

  function_genes <- unlist(strsplit(eggnog_row$gene_list[1], ";"), use.names = FALSE)
  function_genes <- function_genes[function_genes %in% rownames(abundance_df)]
  if (length(function_genes) == 0) {
    warning("No abundance-matched genes found for: ", eggnog_id)
    return(NULL)
  }

  function_abundance <- abundance_df[function_genes, , drop = FALSE]
  total_abundance_per_sample <- colSums(function_abundance)
  sample_names <- colnames(abundance_df)
  function_gene_map <- gene_genus_map %>%
    filter(Gene %in% function_genes)
  genus_names <- sort(unique(function_gene_map$Genus))

  contribution_matrix <- matrix(0, nrow = length(genus_names), ncol = length(sample_names))
  rownames(contribution_matrix) <- genus_names
  colnames(contribution_matrix) <- sample_names

  for (genus in genus_names) {
    overlap_genes <- function_gene_map$Gene[function_gene_map$Genus == genus]
    genus_abundance_per_sample <- colSums(function_abundance[overlap_genes, , drop = FALSE])
    contribution_matrix[genus, ] <- ifelse(
      total_abundance_per_sample > 0,
      genus_abundance_per_sample / total_abundance_per_sample,
      0
    )
  }

  contribution_df <- as.data.frame(contribution_matrix, stringsAsFactors = FALSE)
  contribution_df <- cbind(Genus = rownames(contribution_df), contribution_df, stringsAsFactors = FALSE)
  rownames(contribution_df) <- NULL
  mean_contribution <- rowMeans(contribution_df[, -1, drop = FALSE])
  contribution_df[order(mean_contribution, decreasing = TRUE), ]
}

coverage_rows <- list()
for (eggnog_id in target_eggnog_ids) {
  message("Calculating genus contribution for ", eggnog_id)
  eggnog_row <- eggnog_genes[eggnog_genes$eggNOG == eggnog_id, , drop = FALSE]
  total_genes <- length(unlist(strsplit(eggnog_row$gene_list[1], ";"), use.names = FALSE))
  matched_genes <- sum(unlist(strsplit(eggnog_row$gene_list[1], ";"), use.names = FALSE) %in% rownames(geneset_abundance))

  contribution_result <- calculate_contribution(
    eggnog_id = eggnog_id,
    eggnog_genes_df = eggnog_genes,
    gene_genus_map = target_gene_genus_map,
    abundance_df = geneset_abundance
  )

  if (!is.null(contribution_result)) {
    write.csv(
      contribution_result,
      file.path(output_dir, paste0(eggnog_id, ".csv")),
      row.names = FALSE,
      quote = FALSE
    )
  }

  coverage_rows[[eggnog_id]] <- data.frame(
    eggNOG = eggnog_id,
    total_genes = total_genes,
    matched_genes = matched_genes,
    matched_ratio = matched_genes / total_genes,
    stringsAsFactors = FALSE
  )
}

coverage_df <- bind_rows(coverage_rows)
write.csv(
  coverage_df,
  file.path(output_dir, "target_eggNOG_NOG_top5_gene_coverage.csv"),
  row.names = FALSE,
  quote = FALSE
)

message("Done. Results saved to ", output_dir)
