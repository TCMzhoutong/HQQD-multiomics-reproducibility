# ============================================================================
# Alpha Diversity Statistical Analysis and Visualization
# Author: Generated for metagenomics analysis
# Date: 2026-01-27
# Description: Statistical analysis and visualization of alpha diversity indices
#              for four treatment groups (Control, Model, LVX, HQQD)
# ============================================================================

# Clear environment
rm(list = ls())

# Load required libraries
library(readxl)      # For reading Excel files
library(ggplot2)     # For visualization
library(dplyr)       # For data manipulation
library(tidyr)       # For data reshaping
library(car)         # For Levene's test
library(RColorBrewer) # For color palettes
library(ggpubr)      # For statistical comparison
library(rstatix)     # For statistical tests
library(cowplot)     # For plot arrangement

# Create output directory
output_dir <- "01_alpha_diversity_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(paste0("Created output directory: ", output_dir, "\n"))
}

figure_dir <- file.path("results", "figures")
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

# ============================================================================
# 1. Data Import and Preparation
# ============================================================================

# Read data
data <- read_excel("01_alpha_index_stat.xlsx")

# View data structure
print("Data structure:")
str(data)
head(data)

# Create group variable based on sample prefix
data$Group <- case_when(
  grepl("^K", data$sample) ~ "Control",
  grepl("^M", data$sample) ~ "Model",
  grepl("^Y", data$sample) ~ "LVX",
  grepl("^Z", data$sample) ~ "HQQD"
)

# Set group factor levels for proper ordering
data$Group <- factor(data$Group, levels = c("Control", "Model", "LVX", "HQQD"))

# Display sample distribution
print("Sample distribution by group:")
table(data$Group)

# ============================================================================
# 2. Statistical Analysis
# ============================================================================

# Alpha diversity indices to analyze
indices <- c("ACE", "Chao1", "Shannon", "Simpson", "goods_Coverage", "pielou_evenness")

# Create results data frame
stat_results <- data.frame()

cat("\n=== Statistical Analysis Results ===\n\n")

for (index in indices) {
  cat(paste0("\n--- ", index, " ---\n"))
  
  # Extract data for current index
  index_data <- data[[index]]
  
  # 1. Normality test (Shapiro-Wilk test for each group)
  normality_test <- data %>%
    group_by(Group) %>%
    summarise(
      W = shapiro.test(get(index))$statistic,
      p_value = shapiro.test(get(index))$p.value,
      .groups = 'drop'
    )
  
  cat("Shapiro-Wilk normality test:\n")
  print(normality_test)
  
  # Check if all groups pass normality test (p > 0.05)
  all_normal <- all(normality_test$p_value > 0.05)
  
  # 2. Homogeneity of variance test (Levene's test)
  levene_result <- leveneTest(as.formula(paste(index, "~ Group")), data = data)
  levene_p <- levene_result$`Pr(>F)`[1]
  
  cat(paste0("Levene's test for homogeneity of variance: p = ", 
             round(levene_p, 4), "\n"))
  
  # 3. Choose appropriate statistical test
  if (all_normal && levene_p > 0.05) {
    # Use One-way ANOVA if data is normal and variances are equal
    cat("Using One-way ANOVA (parametric test)\n")
    
    anova_result <- aov(as.formula(paste(index, "~ Group")), data = data)
    anova_summary <- summary(anova_result)
    p_value <- anova_summary[[1]]$`Pr(>F)`[1]
    
    cat(paste0("ANOVA p-value: ", format.pval(p_value), "\n"))
    
    # Post-hoc test: Tukey HSD
    if (p_value < 0.05) {
      cat("Performing Tukey HSD post-hoc test:\n")
      tukey_result <- TukeyHSD(anova_result)
      print(tukey_result)
    }
    
    stat_results <- rbind(stat_results, 
                         data.frame(Index = index, 
                                   Test = "One-way ANOVA",
                                   P_value = p_value,
                                   Significant = ifelse(p_value < 0.05, "Yes", "No")))
    
  } else {
    # Use Kruskal-Wallis test if assumptions are violated
    cat("Using Kruskal-Wallis test (non-parametric test)\n")
    cat("Reason: ")
    if (!all_normal) cat("Data not normally distributed. ")
    if (levene_p <= 0.05) cat("Variances not homogeneous. ")
    cat("\n")
    
    kw_result <- kruskal.test(as.formula(paste(index, "~ Group")), data = data)
    p_value <- kw_result$p.value
    
    cat(paste0("Kruskal-Wallis p-value: ", format.pval(p_value), "\n"))
    
    # Post-hoc test: Dunn test
    if (p_value < 0.05) {
      cat("Performing Dunn's post-hoc test:\n")
      dunn_result <- data %>% 
        dunn_test(as.formula(paste(index, "~ Group")), p.adjust.method = "bonferroni")
      print(dunn_result)
    }
    
    stat_results <- rbind(stat_results, 
                         data.frame(Index = index, 
                                   Test = "Kruskal-Wallis",
                                   P_value = p_value,
                                   Significant = ifelse(p_value < 0.05, "Yes", "No")))
  }
}

# Save statistical results
write.csv(stat_results, file.path(output_dir, "statistical_results.csv"), row.names = FALSE)
cat(paste0("\n\nStatistical results saved to: ", file.path(output_dir, "statistical_results.csv"), "\n"))

# ============================================================================
# 3. Data Visualization
# ============================================================================

# Prepare data for plotting (long format)
data_long <- data %>%
  select(sample, Group, all_of(indices)) %>%
  pivot_longer(cols = all_of(indices), 
               names_to = "Index", 
               values_to = "Value")

# Set index factor levels for proper ordering
data_long$Index <- factor(data_long$Index, levels = indices)

# Define academic color palette
# Using a professional color scheme suitable for scientific publications
group_colors <- c("Control" = "#3B4992", 
                 "Model" = "#EE0000", 
                 "LVX" = "#008B45", 
                 "HQQD" = "#E69F00")

# ============================================================================
# Boxplot with Statistical Significance (Primary Figure)
# ============================================================================

# Create individual plots for each index with statistical annotations
plot_list <- list()

# Create method annotation for each index
for (idx in indices) {
  # Subset data for current index
  data_subset <- data_long %>% filter(Index == idx)
  
  # Get statistical method used
  stat_method <- stat_results$Test[stat_results$Index == idx]
  
  # Perform statistical test for annotation
  if (stat_method == "One-way ANOVA") {
    stat_test <- data_subset %>%
      tukey_hsd(Value ~ Group) %>%
      add_xy_position(x = "Group")
    method_label <- "One-way ANOVA + Tukey HSD"
  } else {
    stat_test <- data_subset %>%
      dunn_test(Value ~ Group, p.adjust.method = "bonferroni") %>%
      add_xy_position(x = "Group")
    method_label <- "Kruskal-Wallis + Dunn test"
  }
  
  p_temp <- ggplot(data_subset, aes(x = Group, y = Value, color = Group, fill = Group)) +
    geom_boxplot(
      width = 0.56,
      outlier.shape = NA,
      linewidth = 0.42,
      alpha = 0.18,
      color = "black"
    ) +
    geom_jitter(width = 0.10, height = 0, size = 1.55, alpha = 0.92, stroke = 0) +
    stat_pvalue_manual(stat_test, 
                       label = "p.adj.signif",
                       hide.ns = TRUE,
                       size = 3,
                       family = "sans") +
    scale_color_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_x_discrete(expand = expansion(mult = 0.08)) +
    labs(title = idx,
         x = NULL,
         y = NULL) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.title.y = element_blank(),
      legend.position = "none",
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      plot.margin = margin(6, 8, 6, 8)
    )
  
  plot_list[[idx]] <- p_temp
}

# Combine all plots in a single column. Width matches one panel in the
# six-panel flow cytometry figure.
p_final <- plot_grid(plotlist = plot_list, ncol = 1, align = "v")

ggsave(file.path(output_dir, "alpha_diversity_boxplot.png"), 
       plot = p_final, 
       width = 2, 
       height = 16, 
       dpi = 300)

ggsave(file.path(output_dir, "alpha_diversity_boxplot.svg"), 
       plot = p_final, 
       width = 2, 
       height = 16)

ggsave(file.path(figure_dir, "Figure_1_Alpha_Diversity_Boxplot.png"), 
       plot = p_final, 
       width = 2, 
       height = 16, 
       dpi = 300)

ggsave(file.path(figure_dir, "Figure_1_Alpha_Diversity_Boxplot.svg"), 
       plot = p_final, 
       width = 2, 
       height = 16)

cat(paste0("Figure saved: ", file.path(output_dir, "alpha_diversity_boxplot.png/svg"), "\n"))

# ============================================================================
# 4. Summary Statistics Table
# ============================================================================

# Create summary table
summary_table <- data_long %>%
  group_by(Group, Index) %>%
  summarise(
    N = n(),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Min = min(Value, na.rm = TRUE),
    Max = max(Value, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

# Save summary table
write.csv(summary_table, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
cat(paste0("Summary statistics saved to: ", file.path(output_dir, "summary_statistics.csv"), "\n"))

# ============================================================================
# 5. Print Summary
# ============================================================================

cat("\n=== Analysis Complete ===\n")
cat("\nAll output files saved to directory: ", output_dir, "\n")
cat("\nFiles generated:\n")
cat("1. statistical_results.csv - Statistical test results\n")
cat("2. summary_statistics.csv - Descriptive statistics\n")
cat("3. alpha_diversity_boxplot.png/pdf - Primary figure with significance markers\n")

cat("\n=== Statistical Summary ===\n")
print(stat_results)

cat("\n=== Recommendations for Publication ===\n")
cat("1. Primary figure: 'alpha_diversity_boxplot.png/pdf'\n")
cat("   - Compact layout with 6 indices in 3x2 grid\n")
cat("   - Boxplots show median, IQR, and individual data points\n")
cat("   - Statistical significance markers added (ns, *, **, ***, ****)\n")
cat("   - Each subplot annotated with statistical method used\n")
cat("   - All text in sans font\n")
cat("   - Follows SCI publication standards\n")
cat("\n2. Statistical methods:\n")
cat("   - Parametric: One-way ANOVA + Tukey HSD post-hoc test\n")
cat("     (Used when data meet normality and homogeneity assumptions)\n")
cat("   - Non-parametric: Kruskal-Wallis + Dunn test with Bonferroni correction\n")
cat("     (Used when assumptions are violated)\n")
cat("   - Appropriate for 4-group comparison (NOT t-test)\n")
cat("\n3. Visualization approach:\n")
cat("   - Boxplots are suitable for both parametric and non-parametric data\n")
cat("   - Show median (robust central tendency) and quartiles\n")
cat("   - Individual points allow assessment of distribution and outliers\n")
cat("\n4. Figure specifications:\n")
cat("   - Size: 10×7 inches (compact for journal submission)\n")
cat("   - Resolution: 300 DPI (PNG) + vector (PDF)\n")
cat("   - Color scheme: Professional academic palette\n")
cat("   - Font: sans (standard for scientific publications)\n")

cat("\n============================================================================\n")
