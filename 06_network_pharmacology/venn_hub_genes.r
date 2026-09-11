# ==============================================================================
# ==============================================================================

library(ggVennDiagram)
library(ggplot2)

if (.Platform$OS.type == "windows") {
  windowsFonts(sans = windowsFont("sans"))
}

cytohubba_df <- read.csv(file.path("06.hub_genes", "MCC.csv"), stringsAsFactors = FALSE)
cytonca_df <- read.csv(file.path("06.hub_genes", "cytoNCA.csv"), stringsAsFactors = FALSE)
mcode_df <- read.csv(file.path("06.hub_genes", "MCODE.csv"), stringsAsFactors = FALSE)

cytohubba_genes <- unique(cytohubba_df$name)
cytonca_genes <- unique(cytonca_df$Gene)
mcode_genes <- unique(mcode_df$name)

# Validate the standardized manual exports against the current STRING network.
string_edges <- read.delim(
  file.path("06.hub_genes", "string_edges_for_cytoscape.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
connected_genes <- unique(c(string_edges$source, string_edges$target))
if (length(connected_genes) != 28) stop("Expected 28 connected STRING nodes")
if (length(cytohubba_genes) != 10) stop("Expected 10 cytoHubba genes")
if (length(cytonca_genes) != 20) stop("Expected 20 CytoNCA genes after preprocessing")
if (length(mcode_genes) < 1) stop("MCODE export is empty")
if (any(is.na(c(cytohubba_genes, cytonca_genes, mcode_genes)))) stop("Hub-gene exports contain missing symbols")
if (length(setdiff(c(cytohubba_genes, cytonca_genes, mcode_genes), connected_genes)) > 0) {
  stop("At least one hub-gene export contains genes outside the current STRING network")
}

venn_list <- list(
  `CytoHubba` = cytohubba_genes,
  CytoNCA = cytonca_genes,
  MCODE = mcode_genes
)

cat("========== Dataset summary ==========\n")
cat("cytoHubba (MCC) genes: ", length(cytohubba_genes), "\n")
cat("CytoNCA genes: ", length(cytonca_genes), "\n")
cat("MCODE genes: ", length(mcode_genes), "\n\n")

intersect_CH_CN <- intersect(cytohubba_genes, cytonca_genes)
intersect_CH_M <- intersect(cytohubba_genes, mcode_genes)
intersect_CN_M <- intersect(cytonca_genes, mcode_genes)
intersect_ALL <- intersect(intersect_CH_CN, mcode_genes)

cat("========== Intersection summary ==========\n")
cat("cytoHubba (MCC) ∩ CytoNCA:", length(intersect_CH_CN), "\n")
cat("cytoHubba (MCC) ∩ MCODE:", length(intersect_CH_M), "\n")
cat("CytoNCA ∩ MCODE:", length(intersect_CN_M), "\n")
cat("cytoHubba (MCC) ∩ CytoNCA ∩ MCODE:", length(intersect_ALL), "\n\n")

if (!dir.exists("figures")) {
  dir.create("figures")
  cat("Created the figures directory\n")
}

# ==============================================================================
# ==============================================================================
cat("Generating the Venn diagram...\n")

venn <- Venn(venn_list)
data <- process_data(venn)


fill_colors <- c("1" = "#E64B3580",
                 "2" = "#4DBBD580",
                 "3" = "#00A08780",
                 "1/2" = "#F39B7F80",
                 "1/3" = "#8491B480",
                 "2/3" = "#91D1C280",
                 "1/2/3" = "#DC000080")

edge_colors <- c("1" = "#000000",
                 "2" = "#000000",
                 "3" = "#000000")

region_data <- venn_regionlabel(data)
total_unique <- sum(region_data$count)
region_data$percentage <- paste0(round(region_data$count / total_unique * 100, 1), "%")
region_data$label <- paste0(region_data$count, "\n(", region_data$percentage, ")")

p_custom <- ggplot() +
  geom_polygon(aes(X, Y, fill = id, group = id), 
               data = venn_regionedge(data),
               show.legend = FALSE) +
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(data), 
            linewidth = 0.8,
            show.legend = FALSE) +
  geom_text(aes(X, Y, label = name), 
            data = venn_setlabel(data),
            fontface = "bold",
            size = 6) +
  geom_text(aes(X, Y, label = label), 
            data = region_data,
            size = 5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = edge_colors) +
  coord_equal() +
  theme_void(base_family = "sans") +
  theme(text = element_text(family = "sans"))

print(p_custom)

ggsave("figures/venn_hub_genes.png", 
       plot = p_custom, 
       width = 10, 
       height = 8, 
       dpi = 300,
       bg = "white")

ggsave("figures/venn_hub_genes.pdf", 
       plot = p_custom, 
       width = 10, 
       height = 8,
       device = "pdf")

svg(
  filename = "figures/venn_hub_genes.svg",
  width = 10,
  height = 8,
  family = "sans",
  bg = "white"
)
print(p_custom)
dev.off()

# ==============================================================================
# ==============================================================================
cat("\nExporting intersection data...\n")

if (length(intersect_ALL) > 0) {
  write.csv(data.frame(name = sort(intersect_ALL)), 
            "06.hub_genes/03_hub_genes_intersection_all.csv", 
            row.names = FALSE)
  cat("Three-way intersection saved: 06.hub_genes/03_hub_genes_intersection_all.csv\n")
  cat("Intersection contains ", length(intersect_ALL), " genes\n")
} else {
  cat("Warning: the three-way intersection is empty\n")
}


# # cytoHubba (MCC) ∩ CytoNCA
# if (length(intersect_CH_CN) > 0) {
#   write.csv(data.frame(name = sort(intersect_CH_CN)), 
#             "06.hub_genes/hub_genes_cytoHubba_CytoNCA.csv", 
#             row.names = FALSE)
# }

# # cytoHubba (MCC) ∩ MCODE
# if (length(intersect_CH_M) > 0) {
#   write.csv(data.frame(name = sort(intersect_CH_M)), 
#             "06.hub_genes/hub_genes_cytoHubba_MCODE.csv", 
#             row.names = FALSE)
# }

# # CytoNCA ∩ MCODE
# if (length(intersect_CN_M) > 0) {
#   write.csv(data.frame(name = sort(intersect_CN_M)), 
#             "06.hub_genes/hub_genes_CytoNCA_MCODE.csv", 
#             row.names = FALSE)
# }

# ==============================================================================
# ==============================================================================
cat("\n========== Analysis complete ==========\n")
cat("Figure files:\n")
cat("  - figures/venn_hub_genes.png\n")
cat("  - figures/venn_hub_genes.pdf\n")
cat("\nIntersection data:\n")
cat("  - 06.hub_genes/03_hub_genes_intersection_all.csv (Three-way intersection:", length(intersect_ALL), " genes)\n")
