#!/usr/bin/env Rscript
# ============================================================================
# ============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(readxl)

output_dir <- "11_KEGG3functional_species_contribution_fdr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# ============================================================================

cat("Reading data...\n")

geneset_abundance <- read.table("geneset.abundance.txt", 
                                header = TRUE, 
                                row.names = 1, 
                                sep = "\t", 
                                comment.char = "",
                                check.names = FALSE)

if (grepl("^#", colnames(geneset_abundance)[1])) {
  colnames(geneset_abundance)[1] <- gsub("^#", "", colnames(geneset_abundance)[1])
}

if (any(duplicated(colnames(geneset_abundance)))) {
  cat("Warning: duplicate sample-column names detected; resolving them...\n")
  dup_names <- colnames(geneset_abundance)[duplicated(colnames(geneset_abundance))]
  cat("Duplicate column names:", paste(unique(dup_names), collapse = ", "), "\n")
  colnames(geneset_abundance) <- make.unique(colnames(geneset_abundance), sep = "_")
  cat("Column names were made unique\n")
}

cat("Sample columns:", paste(colnames(geneset_abundance), collapse = ", "), "\n")
cat("Total: ", ncol(geneset_abundance), " samples\n")

taxonomy_genes <- read.table("nr.taxonomy.gene_list.txt", 
                             header = TRUE, 
                             sep = "\t", 
                             quote = "", 
                             comment.char = "",
                             stringsAsFactors = FALSE)

kegg_genes <- read.table("kegg.kegg_pathway_level3.gene_list.txt", 
                        header = TRUE, 
                        sep = "\t", 
                        quote = "", 
                        comment.char = "",
                        stringsAsFactors = FALSE)

cat("Data import complete.\n")

# ============================================================================
# ============================================================================

cat("Processing taxonomic annotations...\n")

taxonomy_genes$Genus <- sapply(strsplit(taxonomy_genes$Taxon, "; "), function(x) {
  if (length(x) >= 6) {
    genus <- x[6]
    genus <- gsub("^g__", "", genus)
    return(genus)
  } else {
    return("Unassigned")
  }
})

genus_to_genes <- list()
for (i in 1:nrow(taxonomy_genes)) {
  genus <- taxonomy_genes$Genus[i]
  genus <- trimws(genus)
  
  genes <- unlist(strsplit(taxonomy_genes$gene_list[i], ";"))
  
  if (genus %in% names(genus_to_genes)) {
    genus_to_genes[[genus]] <- c(genus_to_genes[[genus]], genes)
  } else {
    genus_to_genes[[genus]] <- genes
  }
}

cat("Taxonomic annotation complete. Identified ", length(genus_to_genes), " genera\n")

unique_genus_count <- length(unique(names(genus_to_genes)))
if (unique_genus_count != length(genus_to_genes)) {
  cat("Warning: duplicate genus names detected; merging entries...\n")
}
cat("Unique genera:", unique_genus_count, "\n")

# ============================================================================
# ============================================================================

calculate_contribution <- function(kegg_id, kegg_genes_df, genus_to_genes_list, abundance_df) {
  
  cat("\nCalculating ", kegg_id, " taxonomic contributions...\n")
  
  kegg_row <- kegg_genes_df[kegg_genes_df$kegg_pathway_id == kegg_id, ]
  
  if (nrow(kegg_row) == 0) {
    cat("Warning: KEGG ID not found: ", kegg_id, "\n")
    return(NULL)
  }
  
  kegg_function_name <- kegg_row$kegg_pathway[1]
  cat("Function:", kegg_function_name, "\n")
  
  function_genes <- unlist(strsplit(kegg_row$gene_list[1], ";"))
  cat("The function contains ", length(function_genes), " genes\n")
  
  function_genes <- function_genes[function_genes %in% rownames(abundance_df)]
  cat("Genes found in the abundance table: ", length(function_genes), " genes\n")
  
  if (length(function_genes) == 0) {
    cat("Warning: no genes for this function were found in the abundance table\n")
    return(NULL)
  }
  
  function_abundance <- abundance_df[function_genes, , drop = FALSE]
  
  total_abundance_per_sample <- colSums(function_abundance)
  
  sample_names <- colnames(abundance_df)
  genus_names <- unique(names(genus_to_genes_list))
  contribution_matrix <- matrix(0, nrow = length(genus_names), ncol = length(sample_names))
  rownames(contribution_matrix) <- genus_names
  colnames(contribution_matrix) <- sample_names
  
  cat("Calculating genus-level contributions...\n")
  for (genus in genus_names) {
    genus_genes <- genus_to_genes_list[[genus]]
    
    overlap_genes <- intersect(genus_genes, function_genes)
    
    if (length(overlap_genes) > 0) {
      overlap_abundance <- function_abundance[overlap_genes, , drop = FALSE]
      
      genus_abundance_per_sample <- colSums(as.matrix(overlap_abundance))
      
      contribution <- ifelse(total_abundance_per_sample > 0, 
                           genus_abundance_per_sample / total_abundance_per_sample, 
                           0)
      
      contribution_matrix[genus, ] <- contribution
    }
  }
  
  contribution_df <- as.data.frame(contribution_matrix, stringsAsFactors = FALSE)
  
  contribution_df <- cbind(
    Genus = rownames(contribution_df), 
    contribution_df, 
    stringsAsFactors = FALSE
  )
  rownames(contribution_df) <- NULL
  
  mean_contribution <- rowMeans(contribution_df[, -1, drop = FALSE])
  contribution_df <- contribution_df[order(mean_contribution, decreasing = TRUE), ]
  
  cat("Calculation complete: ", nrow(contribution_df), " genera\n")
  
  return(contribution_df)
}

# ============================================================================
# ============================================================================

target_kegg_ids <- c("ko02010", "ko02020", "ko02030", "ko02040", "ko00471")

for (kegg_id in target_kegg_ids) {
  cat("\n", strrep("=", 70), "\n")
  
  contribution_result <- calculate_contribution(
    kegg_id = kegg_id,
    kegg_genes_df = kegg_genes,
    genus_to_genes_list = genus_to_genes,
    abundance_df = geneset_abundance
  )
  
  if (!is.null(contribution_result)) {
    output_file <- file.path(output_dir, paste0(kegg_id, ".csv"))
    write.csv(contribution_result, output_file, row.names = FALSE, quote = FALSE)
    cat("Results saved to:", output_file, "\n")
    
    cat("\nTop 10 contributing genera:\n")
    print(head(contribution_result, 10))
  }
}

cat("\n", strrep("=", 70), "\n")
cat("All analyses complete.\n")
cat("Results directory:", output_dir, "\n")
