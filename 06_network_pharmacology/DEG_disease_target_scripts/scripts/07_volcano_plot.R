# ============================================================================
# ============================================================================

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("dplyr")) install.packages("dplyr")

library(ggplot2)
library(ggrepel)
library(dplyr)

script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0 || script_dir == "") {
  project_root <- getwd()
  if (basename(project_root) == "scripts") {
    project_root <- dirname(project_root)
  }
} else {
  project_root <- dirname(script_dir)
}

deg_results_dir <- file.path(project_root, "04_DEG_results")
figures_dir <- file.path(project_root, "05_figures")

cat("=== Load data ===\n")

load(file.path(deg_results_dir, "DEG_results.RData"))

cat("\n=== Prepare volcano-plot data ===\n")

volcano_data <- deg_results
volcano_data$neg_log10_pval <- -log10(volcano_data$adj.P.Val)

if (!"Significant" %in% colnames(volcano_data)) {
  volcano_data$Significant <- ifelse(
    abs(volcano_data$logFC) > logFC_threshold & volcano_data$adj.P.Val < pval_threshold,
    ifelse(volcano_data$logFC > 0, "Up", "Down"),
    "Not Sig"
  )
}

volcano_data$Significant <- factor(volcano_data$Significant, 
                                    levels = c("Up", "Down", "Not Sig"))

volcano_data$Gene <- rownames(volcano_data)

cat("DEG summary:\n")
print(table(volcano_data$Significant))

cat("\n=== Select labelled genes ===\n")

top_up <- volcano_data %>%
  filter(Significant == "Up") %>%
  arrange(desc(logFC)) %>%
  head(25)

top_down <- volcano_data %>%
  filter(Significant == "Down") %>%
  arrange(logFC) %>%
  head(25)

genes_to_label <- rbind(top_up, top_down)

cat("Labelled genes:", nrow(genes_to_label), "\n")
cat("- upregulated genes:", nrow(top_up), "\n")
cat("- downregulated genes:", nrow(top_down), "\n")

cat("\n=== Plot volcano plot ===\n")

volcano_colors <- c(
  "Up" = "#E64B35",
  "Down" = "#4DBBD5",
  "Not Sig" = "gray70"
)

y_max <- min(max(volcano_data$neg_log10_pval, na.rm = TRUE) * 1.1, 50)

p_volcano <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_pval)) +
  geom_point(aes(color = Significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = volcano_colors,
    labels = c(
      "Up" = paste0("Up-regulated (n=", sum(volcano_data$Significant == "Up"), ")"),
      "Down" = paste0("Down-regulated (n=", sum(volcano_data$Significant == "Down"), ")"),
      "Not Sig" = "Not significant"
    )
  ) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "gray40") +
  geom_text_repel(
    data = genes_to_label,
    aes(label = Gene),
    size = 2.5,
    max.overlaps = 30,
    segment.color = "gray50",
    segment.size = 0.3,
    box.padding = 0.3
  ) +
  scale_y_continuous(limits = c(0, y_max)) +
  labs(
    title = "Volcano Plot of DEGs",
    subtitle = paste0("Threshold: |log2FC| > ", logFC_threshold, " & adj.P < ", pval_threshold),
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~(adjusted~p-value)),
    color = "Regulation"
  ) +
  theme_bw(base_family = "sans") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_volcano)

ggsave(file.path(figures_dir, "Figure_C_Volcano_plot.pdf"), 
       p_volcano, width = 10, height = 8)
ggsave(file.path(figures_dir, "Figure_C_Volcano_plot.png"), 
       p_volcano, width = 10, height = 8, dpi = 300)
ggsave(file.path(figures_dir, "Figure_C_Volcano_plot.svg"), 
       p_volcano, width = 10, height = 8)

cat("\nVolcano plot saved to 05_figures.\n")

p_volcano_simple <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_pval)) +
  geom_point(aes(color = Significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = volcano_colors,
    labels = c(
      "Up" = paste0("Up (", sum(volcano_data$Significant == "Up"), ")"),
      "Down" = paste0("Down (", sum(volcano_data$Significant == "Down"), ")"),
      "Not Sig" = "NS"
    )
  ) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "gray40") +
  scale_y_continuous(limits = c(0, y_max)) +
  labs(
    title = "Volcano Plot",
    x = expression(log[2]~FC),
    y = expression(-log[10]~(adj.P)),
    color = ""
  ) +
  theme_bw(base_family = "sans") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(figures_dir, "Figure_C_Volcano_plot_simple.pdf"), 
       p_volcano_simple, width = 7, height = 7)
ggsave(file.path(figures_dir, "Figure_C_Volcano_plot_simple.png"), 
       p_volcano_simple, width = 7, height = 7, dpi = 300)
ggsave(file.path(figures_dir, "Figure_C_Volcano_plot_simple.svg"), 
       p_volcano_simple, width = 7, height = 7)

write.csv(genes_to_label, 
          file.path(figures_dir, "Top50_DEGs_for_volcano_labels.csv"), 
          row.names = FALSE)

cat("\n=== Volcano-plot visualisation complete ===\n")
cat("
Generated files (saved to 05_figures):
  - Figure_C_Volcano_plot.pdf/png: labelled volcano plot
  - Figure_C_Volcano_plot_simple.pdf/png: simplified volcano plot
  - Top50_DEGs_for_volcano_labels.csv: list of the top 50 labelled genes

Next step: source('scripts/08_heatmap.R')
")
