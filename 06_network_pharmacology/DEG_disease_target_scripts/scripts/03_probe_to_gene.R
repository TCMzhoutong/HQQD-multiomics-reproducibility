# ============================================================================
# ============================================================================

library(GEOquery)

script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0 || script_dir == "") {
  project_root <- getwd()
  if (basename(project_root) == "scripts") {
    project_root <- dirname(project_root)
  }
} else {
  project_root <- dirname(script_dir)
}

raw_data_dir <- file.path(project_root, "01_raw_data")
processed_data_dir <- file.path(project_root, "02_processed_data")

cat("=== Load data ===\n")

load(file.path(processed_data_dir, "GSE40012_selected.RData"))
load(file.path(processed_data_dir, "GSE20346_selected.RData"))

load(file.path(raw_data_dir, "GSE40012_data.RData"))
load(file.path(raw_data_dir, "GSE20346_data.RData"))

cat("\n=== Define probe-to-gene mapping function ===\n")

probe_to_gene <- function(expr_matrix, gene_symbols) {
  df <- data.frame(
    probe = rownames(expr_matrix),
    symbol = gene_symbols,
    stringsAsFactors = FALSE
  )
  
  valid_idx <- !is.na(df$symbol) & df$symbol != "" & df$symbol != "---"
  df <- df[valid_idx, ]
  expr_matrix <- expr_matrix[valid_idx, ]
  
  cat("Valid probes:", nrow(df), "\n")
  cat("Unique genes:", length(unique(df$symbol)), "\n")
  
  gene_expr <- aggregate(expr_matrix, 
                         by = list(Gene = df$symbol), 
                         FUN = mean)
  
  rownames(gene_expr) <- gene_expr$Gene
  gene_expr <- gene_expr[, -1]
  
  return(as.matrix(gene_expr))
}

cat("\n=== GSE40012 probe-to-gene mapping ===\n")

fdata_40012 <- fData(gse40012_data)
probe_ids_in_expr <- rownames(expr_40012_selected)

if (all(probe_ids_in_expr %in% rownames(fdata_40012))) {
  cat("All probes were matched in the GPL annotation.\n")
  gene_symbol_40012_selected <- fdata_40012[probe_ids_in_expr, "Symbol"]
} else {
  cat("Warning: some probes could not be matched\n")
  cat("Match rate:", sum(probe_ids_in_expr %in% rownames(fdata_40012)), "/", length(probe_ids_in_expr), "\n")
  gene_symbol_40012_selected <- rep(NA, length(probe_ids_in_expr))
  matched_idx <- probe_ids_in_expr %in% rownames(fdata_40012)
  gene_symbol_40012_selected[matched_idx] <- fdata_40012[probe_ids_in_expr[matched_idx], "Symbol"]
}

names(gene_symbol_40012_selected) <- probe_ids_in_expr

expr_40012_gene <- probe_to_gene(expr_40012_selected, gene_symbol_40012_selected)
cat("Mapped-matrix dimensions:", dim(expr_40012_gene), "\n")

cat("\n=== GSE20346 probe-to-gene mapping ===\n")

fdata_20346 <- fData(gse20346_data)
probe_ids_in_expr_20346 <- rownames(expr_20346_selected)

if (all(probe_ids_in_expr_20346 %in% rownames(fdata_20346))) {
  cat("All probes were matched in the GPL annotation.\n")
  gene_symbol_20346_selected <- fdata_20346[probe_ids_in_expr_20346, "Symbol"]
} else {
  cat("Warning: some probes could not be matched\n")
  cat("Match rate:", sum(probe_ids_in_expr_20346 %in% rownames(fdata_20346)), "/", length(probe_ids_in_expr_20346), "\n")
  gene_symbol_20346_selected <- rep(NA, length(probe_ids_in_expr_20346))
  matched_idx <- probe_ids_in_expr_20346 %in% rownames(fdata_20346)
  gene_symbol_20346_selected[matched_idx] <- fdata_20346[probe_ids_in_expr_20346[matched_idx], "Symbol"]
}

names(gene_symbol_20346_selected) <- probe_ids_in_expr_20346

expr_20346_gene <- probe_to_gene(expr_20346_selected, gene_symbol_20346_selected)
cat("Mapped-matrix dimensions:", dim(expr_20346_gene), "\n")

cat("\n=== Identify shared genes ===\n")

genes_40012 <- rownames(expr_40012_gene)
genes_20346 <- rownames(expr_20346_gene)

common_genes <- intersect(genes_40012, genes_20346)
cat("GSE40012Genes:", length(genes_40012), "\n")
cat("GSE20346Genes:", length(genes_20346), "\n")
cat("Shared genes:", length(common_genes), "\n")

expr_40012_common <- expr_40012_gene[common_genes, ]
expr_20346_common <- expr_20346_gene[common_genes, ]

cat("\nShared-gene expression-matrix dimensions:\n")
cat("GSE40012:", dim(expr_40012_common), "\n")
cat("GSE20346:", dim(expr_20346_common), "\n")

cat("\n=== Check data distributions ===\n")

cat("\nGSE40012expression-value range:\n")
cat("Min:", min(expr_40012_common), "\n")
cat("Max:", max(expr_40012_common), "\n")
cat("Mean:", mean(expr_40012_common), "\n")

cat("\nGSE20346expression-value range:\n")
cat("Min:", min(expr_20346_common), "\n")
cat("Max:", max(expr_20346_common), "\n")
cat("Mean:", mean(expr_20346_common), "\n")

if (max(expr_40012_common) > 50) {
  cat("\nGSE40012apply log2 transformation...\n")
  if (min(expr_40012_common) < 0) {
    cat("  Warning: negative values detected and replaced with 0.1\n")
    expr_40012_common[expr_40012_common < 0] <- 0.1
  }
  expr_40012_common <- log2(expr_40012_common + 1)
  cat("  Range after transformation: ", round(min(expr_40012_common), 2), " - ", 
      round(max(expr_40012_common), 2), "\n")
}

if (max(expr_20346_common) > 50) {
  cat("\nGSE20346apply log2 transformation...\n")
  if (min(expr_20346_common) < 0) {
    cat("  Warning: detected ", sum(expr_20346_common < 0), "negative values\n")
    cat("  Replacing negative values with 0.1\n")
    expr_20346_common[expr_20346_common < 0] <- 0.1
  }
  expr_20346_common <- log2(expr_20346_common + 1)
  cat("  Range after transformation: ", round(min(expr_20346_common), 2), " - ", 
      round(max(expr_20346_common), 2), "\n")
}

cat("\n=== Save data ===\n")

save(expr_40012_common, group_40012_selected,
     expr_20346_common, group_20346_selected,
     common_genes,
     file = file.path(processed_data_dir, "data_ready_for_merge.RData"))

cat("Data saved and ready for batch correction.\n")

cat("\n
==========================================================
                Data-preparation summary
==========================================================

GSE40012 (gene level):
  - Samples: ", ncol(expr_40012_common), "
  - Genes: ", nrow(expr_40012_common), "

GSE20346 (gene level):
  - Samples: ", ncol(expr_20346_common), "
  - Genes: ", nrow(expr_20346_common), "

Shared genes: ", length(common_genes), "

Next step: source('scripts/04_batch_correction.R')
==========================================================
", sep = "")
