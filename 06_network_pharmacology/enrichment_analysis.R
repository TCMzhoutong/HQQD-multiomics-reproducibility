# ==============================================================================
# ==============================================================================

options(clusterProfiler.download.method = "auto")
options(timeout = 300)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(circlize)
  library(RColorBrewer)
})

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
if (script_dir != "") {
  setwd(script_dir)
}

if (!dir.exists("07.enrichment_analysis")) {
  dir.create("07.enrichment_analysis", recursive = TRUE)
}

# ==============================================================================
# ==============================================================================
cat("Reading hub-gene data...\n")
hub_genes <- read.csv("06.hub_genes\\03_hub_genes_intersection_all.csv", 
                      stringsAsFactors = FALSE)

gene_list <- hub_genes$name
cat("Identified ", length(gene_list), " hub genes\n")
cat("Gene list:\n")
print(gene_list)

# ==============================================================================
# ==============================================================================
cat("\nConverting gene identifiers...\n")
gene_df <- bitr(gene_list, 
                fromType = "SYMBOL",
                toType = c("ENTREZID", "ENSEMBL"),
                OrgDb = org.Hs.eg.db)

cat("Initial mapping result: ", nrow(gene_df), " rows (one-to-many mappings may be present)\n")

gene_df_unique <- gene_df %>%
  group_by(SYMBOL) %>%
  slice(1) %>%
  ungroup() %>%
  as.data.frame()

cat("After deduplication: ", nrow(gene_df_unique), " unique genes\n")
cat("Successfully mapped ", nrow(gene_df_unique), " gene identifiers\n")

write.csv(gene_df, 
          "07.enrichment_analysis/00_gene_id_conversion_full.csv", 
          row.names = FALSE)

write.csv(gene_df_unique, 
          "07.enrichment_analysis/00_gene_id_conversion.csv", 
          row.names = FALSE)

gene_df <- gene_df_unique

# ==============================================================================
# ==============================================================================
cat("\nRunning GO enrichment analysis...\n")

ego_BP <- enrichGO(gene         = gene_df$ENTREZID,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable     = TRUE)

ego_CC <- enrichGO(gene         = gene_df$ENTREZID,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable     = TRUE)

ego_MF <- enrichGO(gene         = gene_df$ENTREZID,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable     = TRUE)

ego_ALL <- enrichGO(gene         = gene_df$ENTREZID,
                    OrgDb        = org.Hs.eg.db,
                    ont          = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable     = TRUE)

if (!is.null(ego_BP)) {
  write.csv(as.data.frame(ego_BP), 
            "07.enrichment_analysis/01_GO_BP_enrichment.csv", 
            row.names = FALSE)
  cat("GO BP enrichment analysis complete; identified ", nrow(ego_BP), " significant terms\n")
}

if (!is.null(ego_CC)) {
  write.csv(as.data.frame(ego_CC), 
            "07.enrichment_analysis/02_GO_CC_enrichment.csv", 
            row.names = FALSE)
  cat("GO CC enrichment analysis complete; identified ", nrow(ego_CC), " significant terms\n")
}

if (!is.null(ego_MF)) {
  write.csv(as.data.frame(ego_MF), 
            "07.enrichment_analysis/03_GO_MF_enrichment.csv", 
            row.names = FALSE)
  cat("GO MF enrichment analysis complete; identified ", nrow(ego_MF), " significant terms\n")
}

if (!is.null(ego_ALL)) {
  write.csv(as.data.frame(ego_ALL), 
            "07.enrichment_analysis/04_GO_ALL_enrichment.csv", 
            row.names = FALSE)
  cat("GO ALL enrichment analysis complete; identified ", nrow(ego_ALL), " significant terms\n")
}

# ==============================================================================
# ==============================================================================
cat("\nRunning KEGG enrichment analysis...\n")

ekegg <- enrichKEGG(gene         = gene_df$ENTREZID,
                    organism     = 'hsa',
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)

if (!is.null(ekegg) && nrow(ekegg) > 0) {
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  write.csv(as.data.frame(ekegg), 
            "07.enrichment_analysis/05_KEGG_enrichment.csv", 
            row.names = FALSE)
  cat("KEGG enrichment analysis complete; identified ", nrow(ekegg), " significant pathways\n")
} else {
  cat("No significantly enriched KEGG pathways were identified\n")
}

# ==============================================================================
# ==============================================================================
cat("\nGenerating visualisations...\n")

if (!is.null(ego_ALL) && nrow(ego_ALL) > 0) {
  p_go <- dotplot(ego_ALL, showCategory = 5, split = "ONTOLOGY") + 
    facet_grid(ONTOLOGY~., scale = "free") +
    theme(
      plot.title = element_blank(),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(size = 10, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50))
  
  ggsave("07.enrichment_analysis/plot_01_GO_ALL_dotplot.pdf", 
         plot = p_go, width = 8, height = 10)
  ggsave("07.enrichment_analysis/plot_01_GO_ALL_dotplot.svg", 
         plot = p_go, width = 8, height = 10, device = "svg")
  cat("GO ALL dotplot saved\n")
}

if (!is.null(ekegg) && nrow(ekegg) > 0) {
  p_kegg <- dotplot(ekegg, showCategory = 20) +
    theme(
      plot.title = element_blank(),
      axis.text.y = element_text(size = 9),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45))
  
  ggsave("07.enrichment_analysis/plot_02_KEGG_dotplot.pdf", 
         plot = p_kegg, width = 8, height = 10)
  ggsave("07.enrichment_analysis/plot_02_KEGG_dotplot.svg", 
         plot = p_kegg, width = 8, height = 10, device = "svg")
  cat("KEGG dotplot saved\n")
}

if (!is.null(ekegg) && nrow(ekegg) >= 5) {
  tryCatch({
    library(circlize)
    library(RColorBrewer)
    library(grid)
    library(gridBase)
    
    n_category <- min(20, nrow(ekegg))
    
    kegg_df <- as.data.frame(ekegg)
    top_pathways <- head(kegg_df, n_category)
    
    links_list <- lapply(1:nrow(top_pathways), function(i) {
      genes <- strsplit(top_pathways$geneID[i], "/")[[1]]
      data.frame(
        gene = genes,
        pathway = top_pathways$Description[i],
        value = 1,
        stringsAsFactors = FALSE
      )
    })
    links_df <- do.call(rbind, links_list)
    
    unique_genes <- unique(links_df$gene)
    unique_pathways <- top_pathways$Description
    
    n_colors <- length(unique_pathways)
    base_colors <- c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  # Set1
      "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
      "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",  # Set2+Set3
      "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"   # Dark2
    )
    if (n_colors <= length(base_colors)) {
      pathway_colors <- base_colors[1:n_colors]
    } else {
      pathway_colors <- colorRampPalette(base_colors)(n_colors)
    }
    names(pathway_colors) <- unique_pathways
    
    gene_colors <- rep("#B0B0B0", length(unique_genes))
    names(gene_colors) <- unique_genes
    
    all_colors <- c(pathway_colors, gene_colors)
    
    sector_order <- c(unique_pathways, sort(unique_genes))
    
    gene_counts <- table(links_df$gene)
    pathway_counts <- table(links_df$pathway)
    
    n_genes <- length(unique_genes)
    n_pathways <- length(unique_pathways)
    total_sectors <- n_genes + n_pathways
    
    normal_gap <- 1
    gap_sizes <- c(rep(normal_gap, n_pathways - 1), normal_gap, 
                   rep(normal_gap, n_genes - 1), normal_gap)
    
    for (plot_device in c("pdf", "svg")) {
      if (plot_device == "pdf") {
        pdf("07.enrichment_analysis/plot_03_KEGG_chordplot.pdf",
            width = 18, height = 12)
      } else {
        svg("07.enrichment_analysis/plot_03_KEGG_chordplot.svg",
            width = 18, height = 12)
      }
    
    layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1.2))
    
    par(mar = c(1, 6, 1, 0))
    
    circos.clear()
    circos.par(
      gap.after = gap_sizes,
      start.degree = 90,
      clock.wise = TRUE,
      cell.padding = c(0.02, 0, 0.02, 0),
      canvas.xlim = c(-1.1, 1.1),
      canvas.ylim = c(-1.1, 1.1)
    )
    
    circos.initialize(
      factors = factor(sector_order, levels = sector_order),
      xlim = cbind(rep(0, length(sector_order)), 
                   c(as.numeric(pathway_counts[unique_pathways]),
                     as.numeric(gene_counts[sort(unique_genes)])))
    )
    
    circos.track(
      factors = factor(sector_order, levels = sector_order),
      ylim = c(0, 1),
      track.height = 0.05,
      bg.col = all_colors[sector_order],
      bg.border = NA,
      panel.fun = function(x, y) {
        sector_name <- get.cell.meta.data("sector.index")
        if (sector_name %in% unique_genes) {
          circos.text(
            CELL_META$xcenter, 
            CELL_META$cell.ylim[2] + 0.4,
            sector_name,
            facing = "clockwise",
            niceFacing = TRUE,
            adj = c(0, 0.5),
            cex = 1.1
          )
        }
      }
    )
    
    for (i in 1:nrow(links_df)) {
      gene <- links_df$gene[i]
      pathway <- links_df$pathway[i]
      
      gene_links <- links_df[links_df$gene == gene, ]
      gene_idx <- which(gene_links$pathway == pathway)
      
      pathway_links <- links_df[links_df$pathway == pathway, ]
      pathway_idx <- which(pathway_links$gene == gene)
      
      gene_start <- gene_idx - 1
      gene_end <- gene_idx
      
      pathway_start <- pathway_idx - 1
      pathway_end <- pathway_idx
      
      chord_color <- adjustcolor(pathway_colors[pathway], alpha.f = 0.6)
      
      circos.link(
        gene, c(gene_start, gene_end),
        pathway, c(pathway_start, pathway_end),
        col = chord_color,
        border = NA
      )
    }
    
    circos.clear()
    
    par(mar = c(1, 0, 1, 1))
    plot.new()
    
    wrapped_pathways <- sapply(unique_pathways, function(x) {
      if (nchar(x) > 35) {
        paste(strwrap(x, width = 35), collapse = "\n")
      } else {
        x
      }
    })
    
    legend(
      "left",
      legend = wrapped_pathways,
      fill = pathway_colors,
      border = NA,
      ncol = 1,
      cex = 1.1,
      bty = "n",
      x.intersp = 0.5,
      y.intersp = 1.5
    )
    
    dev.off()
    }
    cat("KEGG  chord diagram (top", n_category, ") PDF saved\n")
    
  }, error = function(e) {
    cat("Error generating the KEGG chord diagram:", e$message, "\n")
    print(e)
  })
}

# ==============================================================================
# ==============================================================================
cat("\nGenerating the analysis summary...\n")

summary_text <- paste0(
  "==============================================================================\n",
  "Functional enrichment summary for hub genes\n",
  "==============================================================================\n\n",
  "Analysis time: ", Sys.time(), "\n\n",
  "1. Input data\n",
  "   - Input genes: ", length(gene_list), "\n",
  "   - Successfully mapped genes: ", nrow(gene_df), "\n\n",
  "2. GO enrichment\n",
  "   - GO Biological Process terms: ", ifelse(!is.null(ego_BP), nrow(ego_BP), 0), "\n",
  "   - GO Cellular Component terms: ", ifelse(!is.null(ego_CC), nrow(ego_CC), 0), "\n",
  "   - GO Molecular Function terms: ", ifelse(!is.null(ego_MF), nrow(ego_MF), 0), "\n",
  "   - All GO terms: ", ifelse(!is.null(ego_ALL), nrow(ego_ALL), 0), "\n\n",
  "3. KEGG enrichment\n",
  "   - KEGG pathways: ", ifelse(!is.null(ekegg), nrow(ekegg), 0), "\n\n",
  "4. Output\n",
  "   - Results are written to: 07.enrichment_analysis/\n\n",
  "==============================================================================\n"
)

cat(summary_text)
writeLines(summary_text, "07.enrichment_analysis/00_analysis_summary.txt")

cat("\nAnalysis complete. Results were saved to 07.enrichment_analysis/.\n")

# ==============================================================================
# ==============================================================================
save.image("07.enrichment_analysis/enrichment_analysis_workspace.RData")
cat("R workspace saved\n")
