# ============================================================================
# LEfSe Differential eggNOG Class Boxplot with Statistical Annotation (HQQD vs Model Selection)
# Author: Generated for metagenomics analysis
# Date: 2026-02-24
# Description: Horizontal boxplot showing differential eggNOG Class functions
#              Selection: Features with significant callback trend (Model vs Control & HQQD vs Model)
#              Plotting: STAMP-like extended error bar plot with Wilcoxon tests
# ============================================================================

rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readxl, dplyr, tidyr, ggplot2, ggpubr, stringr,
  patchwork, microeco, magrittr, tibble
)

group_colors <- c(
  "Control" = "#3B4992",
  "Model" = "#EE0000",
  "LVX" = "#008B45",
  "HQQD" = "#FF8C00"
)

abundance_file <- "eggNOG.Class.relabundance.xlsx"
output_dir <- "08d_HQQD_reverse_eggNOG_Class"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

TOP_N <- 10
LDA_CUTOFF <- 2.0
PSEUDO_COUNT <- 0.0001
ALPHA <- 0.05
source("reverse_fdr_plot_helper.R")

class_levels <- list(
  list(name = "Class", column = "eggNOG_Class")
)

cat("\n=== Step 1: Loading data ===\n")
if (!file.exists(abundance_file)) stop(paste("Input file not found:", abundance_file))

abundance_data <- read_excel(abundance_file)
cat("eggNOG Class data loaded:", nrow(abundance_data), "features\n")

# Identify sample columns strictly by prefix to avoid non-numeric metadata columns
sample_cols <- grep("^[KMYZ][0-9]", colnames(abundance_data), value = TRUE)
cat("Number of samples:", length(sample_cols), "\n\n")

cat("=== Step 2: Creating sample metadata ===\n")
sample_metadata <- data.frame(
  Sample = sample_cols,
  Group = case_when(
    grepl("^K", sample_cols) ~ "Control",
    grepl("^M", sample_cols) ~ "Model",
    grepl("^Y", sample_cols) ~ "LVX",
    grepl("^Z", sample_cols) ~ "HQQD",
    TRUE ~ "Unknown"
  ),
  stringsAsFactors = FALSE
)
rownames(sample_metadata) <- sample_metadata$Sample
sample_metadata$Group <- factor(sample_metadata$Group, levels = c("HQQD", "LVX", "Model", "Control"))
cat("Sample distribution:\n")
print(table(sample_metadata$Group))
cat("\n")

get_clean_name <- function(feature_string, keep_prefix = FALSE) {
  levels <- strsplit(feature_string, "\\|")[[1]]
  last_level <- levels[length(levels)]
  if (keep_prefix) return(last_level)
  sub("^CLASS__", "", last_level)
}

match_feature_to_abundance <- function(feature_string, abundance_data) {
  target_feature <- sub("^CLASS__", "", get_clean_name(feature_string, keep_prefix = TRUE))
  matched_rows <- abundance_data[abundance_data[["eggNOG_Class"]] == target_feature, ]
  if (nrow(matched_rows) > 0) {
    matched_numeric <- matched_rows[, sample_cols, drop = FALSE]
    matched_numeric[] <- lapply(matched_numeric, function(x) as.numeric(x))
    return(colSums(matched_numeric, na.rm = TRUE))
  }
  rep(0, length(sample_cols))
}

get_lefse_lda_for_feature <- function(feature_string, abundance_data, group1_name, group2_name, class_level, sample_metadata, sample_cols) {
  cache_key <- paste(class_level, group1_name, group2_name, sep = "_")

  if (!exists("LEFSE_CACHE", envir = .GlobalEnv)) assign("LEFSE_CACHE", list(), envir = .GlobalEnv)
  LEFSE_CACHE <- get("LEFSE_CACHE", envir = .GlobalEnv)

  if (!(cache_key %in% names(LEFSE_CACHE))) {
    cat("    Running LEfSe for", cache_key, "...")

    comp_metadata <- sample_metadata %>%
      filter(Group %in% c(group1_name, group2_name)) %>%
      mutate(Group = factor(as.character(Group), levels = c(group1_name, group2_name))) %>%
      as.data.frame()
    rownames(comp_metadata) <- comp_metadata$Sample

    if (nrow(comp_metadata) < 2) {
      LEFSE_CACHE[[cache_key]] <- data.frame()
      assign("LEFSE_CACHE", LEFSE_CACHE, envir = .GlobalEnv)
      return(0)
    }

    comp_samples <- comp_metadata$Sample

    agg_abundance <- abundance_data %>%
      filter(!is.na(eggNOG_Class) & eggNOG_Class != "") %>%
      group_by(eggNOG_Class) %>%
      summarise(across(all_of(comp_samples), sum), .groups = "drop")

    if (nrow(agg_abundance) == 0) {
      cat(" no data\n")
      LEFSE_CACHE[[cache_key]] <- data.frame()
      assign("LEFSE_CACHE", LEFSE_CACHE, envir = .GlobalEnv)
      return(0)
    }

    colnames(agg_abundance)[1] <- "Feature_Name"
    agg_abundance$Feature_ID <- paste0("CLASS_", sprintf("%05d", seq_len(nrow(agg_abundance))))

    otu_table_comp <- agg_abundance %>%
      select(Feature_ID, all_of(comp_samples)) %>%
      tibble::column_to_rownames("Feature_ID") %>%
      as.data.frame()

    tax_table_comp <- data.frame(
      Kingdom = rep(class_level, nrow(agg_abundance)),
      Phylum = rep(class_level, nrow(agg_abundance)),
      Class = rep("Unassigned", nrow(agg_abundance)),
      Order = rep("Unassigned", nrow(agg_abundance)),
      Family = rep("Unassigned", nrow(agg_abundance)),
      Genus = rep("Unassigned", nrow(agg_abundance)),
      Species = agg_abundance$Feature_Name,
      row.names = agg_abundance$Feature_ID,
      stringsAsFactors = FALSE
    )

    feature_names_map <- setNames(agg_abundance$Feature_Name, agg_abundance$Feature_ID)

    tryCatch({
      comp_dataset <- microtable$new(sample_table = comp_metadata, otu_table = otu_table_comp, tax_table = tax_table_comp)
      comp_lefse <- trans_diff$new(
        dataset = comp_dataset,
        method = "lefse",
        group = "Group",
        alpha = ALPHA,
        lefse_subgroup = NULL,
        p_adjust_method = "fdr"
      )

      if (!is.null(comp_lefse$res_diff) && nrow(comp_lefse$res_diff) > 0) {
        comp_lefse$res_diff$Taxon_Name <- sapply(comp_lefse$res_diff$Taxa, function(x) {
          if (grepl("^CLASS_", x) && x %in% names(feature_names_map)) return(feature_names_map[[x]])
          parts <- strsplit(x, "\\|")[[1]]
          if (length(parts) > 0) return(parts[length(parts)])
          x
        })
      }

      LEFSE_CACHE[[cache_key]] <- comp_lefse$res_diff
      cat(" found", nrow(comp_lefse$res_diff), "features\n")
    }, error = function(e) {
      cat(" ERROR:", e$message, "\n")
      LEFSE_CACHE[[cache_key]] <- data.frame()
    })

    assign("LEFSE_CACHE", LEFSE_CACHE, envir = .GlobalEnv)
  }

  LEFSE_CACHE <- get("LEFSE_CACHE", envir = .GlobalEnv)
  lefse_res <- LEFSE_CACHE[[cache_key]]
  if (is.null(lefse_res) || nrow(lefse_res) == 0) return(0)

  target_feature <- sub("^CLASS__", "", get_clean_name(feature_string, keep_prefix = TRUE))
  matching_rows <- lefse_res[lefse_res$Taxon_Name == target_feature, ]
  if (nrow(matching_rows) > 0) return(matching_rows$LDA[1])
  0
}

check_callback_trend <- function(feature_string, abundance_data, class_level, sample_metadata, sample_cols) {
  abundances <- match_feature_to_abundance(feature_string, abundance_data)
  temp_df <- data.frame(Abundance = as.numeric(abundances), Group = sample_metadata$Group)
  temp_df$Log_Abundance <- log10(temp_df$Abundance + PSEUDO_COUNT)

  means <- temp_df %>% group_by(Group) %>% summarize(Mean = mean(Log_Abundance), .groups = "drop")
  mean_control <- means$Mean[means$Group == "Control"]
  mean_model <- means$Mean[means$Group == "Model"]
  mean_hqqd <- means$Mean[means$Group == "HQQD"]

  lda_model_control <- get_lefse_lda_for_feature(feature_string, abundance_data, "Model", "Control", class_level, sample_metadata, sample_cols)
  lda_hqqd_model <- get_lefse_lda_for_feature(feature_string, abundance_data, "HQQD", "Model", class_level, sample_metadata, sample_cols)

  control_vals <- temp_df$Log_Abundance[temp_df$Group == "Control"]
  model_vals   <- temp_df$Log_Abundance[temp_df$Group == "Model"]
  hqqd_vals    <- temp_df$Log_Abundance[temp_df$Group == "HQQD"]

  p_mc <- tryCatch(wilcox.test(model_vals, control_vals, exact = FALSE)$p.value, error = function(e) 1.0)
  p_hm <- tryCatch(wilcox.test(hqqd_vals,  model_vals,   exact = FALSE)$p.value, error = function(e) 1.0)

  trend_mc <- ifelse(!is.na(p_mc) && p_mc < ALPHA, ifelse(mean_model > mean_control, "Up", "Down"), "None")
  trend_hm <- ifelse(!is.na(p_hm) && p_hm < ALPHA, ifelse(mean_hqqd > mean_model,   "Up", "Down"), "None")

  if (lda_model_control < LDA_CUTOFF || lda_hqqd_model < LDA_CUTOFF) {
    return(list(
      is_callback = FALSE,
      p_model_control = p_mc,
      p_hqqd_model = p_hm,
      trend_model_vs_control = trend_mc,
      trend_hqqd_vs_model = trend_hm,
      lda_model_control = lda_model_control,
      lda_hqqd_model = lda_hqqd_model
    ))
  }

  is_callback <- (!is.na(p_mc) && p_mc < ALPHA) &&
    (!is.na(p_hm) && p_hm < ALPHA) &&
    (trend_mc != "None") && (trend_hm != "None") &&
    (trend_mc != trend_hm)

  list(
    is_callback = is_callback,
    p_model_control = p_mc,
    p_hqqd_model = p_hm,
    trend_model_vs_control = trend_mc,
    trend_hqqd_vs_model = trend_hm,
    lda_model_control = lda_model_control,
    lda_hqqd_model = lda_hqqd_model
  )
}

if (exists("LEFSE_CACHE", envir = .GlobalEnv)) rm(LEFSE_CACHE, envir = .GlobalEnv)

for (level_info in class_levels) {
  class_level <- level_info$name
  level_column <- level_info$column

  cat("\n", rep("=", 70), "\n", sep = "")
  cat("=== Processing eggNOG", class_level, "===\n")
  cat(rep("=", 70), "\n\n", sep = "")

  unique_features <- abundance_data %>%
    filter(!is.na(.data[[level_column]]) & .data[[level_column]] != "") %>%
    pull(.data[[level_column]]) %>%
    unique()

  cat("Total eggNOG features at", class_level, ":", length(unique_features), "\n")

  level_results <- data.frame(Taxa = character(), stringsAsFactors = FALSE)

  for (feature in unique_features) {
    feature_string <- paste0("CLASS__", feature)
    level_results <- rbind(level_results, data.frame(Taxa = feature_string, stringsAsFactors = FALSE))
  }

  if (nrow(level_results) == 0) {
    cat("No features found at this level. Skipping.\n")
    next
  }

  cat("Total candidates for filtering:", nrow(level_results), "\n")
  cat("Checking for significant callback trend (Model vs Control AND HQQD vs Model)...\n")
  cat("LDA threshold:", LDA_CUTOFF, "\n")

  callback_results <- vector("list", nrow(level_results))
  for (i in seq_len(nrow(level_results))) {
    callback_results[[i]] <- check_callback_trend(level_results$Taxa[i], abundance_data, class_level, sample_metadata, sample_cols)
    if (i %% 10 == 0) cat(".")
  }
  cat("\n")

  level_results$is_callback <- sapply(callback_results, function(x) x$is_callback)
  level_results$p_val_Model_Control <- sapply(callback_results, function(x) x$p_model_control)
  level_results$p_val_HQQD_Model <- sapply(callback_results, function(x) x$p_hqqd_model)
  level_results$Trend_Model_vs_Control <- sapply(callback_results, function(x) x$trend_model_vs_control)
  level_results$Trend_HQQD_vs_Model <- sapply(callback_results, function(x) x$trend_hqqd_vs_model)
  level_results$LDA_Model_Control <- sapply(callback_results, function(x) x$lda_model_control)
  level_results$LDA_HQQD_Model <- sapply(callback_results, function(x) x$lda_hqqd_model)
  level_results$q_val_Model_Control <- p.adjust(level_results$p_val_Model_Control, method = "BH")
  level_results$q_val_HQQD_Model <- p.adjust(level_results$p_val_HQQD_Model, method = "BH")
  level_results$is_callback_fdr <- level_results$is_callback &
    !is.na(level_results$q_val_Model_Control) &
    !is.na(level_results$q_val_HQQD_Model) &
    level_results$q_val_Model_Control < ALPHA &
    level_results$q_val_HQQD_Model < ALPHA

  sig_results <- level_results %>% filter(is_callback)
  cat("eggNOG features with significant callback trend:", nrow(sig_results), "\n")
  cat("eggNOG features retained after BH correction:", sum(level_results$is_callback_fdr, na.rm = TRUE), "\n")
  if (nrow(sig_results) == 0) next

  all_stat_results <- list()
  for (i in seq_len(nrow(sig_results))) {
    feature_string <- sig_results$Taxa[i]
    abundances <- match_feature_to_abundance(feature_string, abundance_data)

    temp_data <- data.frame(Abundance = as.numeric(abundances), Group = sample_metadata$Group)
    temp_data$Log_Abundance <- log10(temp_data$Abundance + PSEUDO_COUNT)

    mc_vals   <- temp_data$Log_Abundance[temp_data$Group == "Model"]
    ctrl_vals <- temp_data$Log_Abundance[temp_data$Group == "Control"]
    hqqd_vals_s <- temp_data$Log_Abundance[temp_data$Group == "HQQD"]
    lvx_vals  <- temp_data$Log_Abundance[temp_data$Group == "LVX"]

    p_mc <- tryCatch(wilcox.test(mc_vals,     ctrl_vals,   exact = FALSE)$p.value, error = function(e) NA)
    p_hm <- tryCatch(wilcox.test(hqqd_vals_s, mc_vals,     exact = FALSE)$p.value, error = function(e) NA)
    p_ml <- tryCatch(wilcox.test(mc_vals,     lvx_vals,    exact = FALSE)$p.value, error = function(e) NA)

    clean_name <- get_clean_name(feature_string)
    all_stat_results[[clean_name]] <- list(
      Taxa_Full          = feature_string,
      p_Model_vs_Control = p_mc,
      p_HQQD_vs_Model    = p_hm,
      p_Model_vs_LVX     = p_ml
    )
  }

  sig_results$LDA_avg <- (sig_results$LDA_Model_Control + sig_results$LDA_HQQD_Model) / 2
  top_features <- sig_results %>% arrange(desc(LDA_avg))

  plot_data_list <- vector("list", nrow(top_features))
  for (i in seq_len(nrow(top_features))) {
    feature_row <- top_features[i, ]
    feature_string <- feature_row$Taxa
    clean_name <- get_clean_name(feature_string)
    abundances <- match_feature_to_abundance(feature_string, abundance_data)

    feature_data <- data.frame(
      Sample = sample_cols,
      Abundance = as.numeric(abundances),
      Taxon = clean_name,
      Taxon_Full = feature_string,
      LDA_Score = feature_row$LDA_avg,
      stringsAsFactors = FALSE
    )

    feature_data <- feature_data %>% left_join(sample_metadata, by = "Sample")
    feature_data$Log_Abundance <- log10(feature_data$Abundance + PSEUDO_COUNT)
    plot_data_list[[i]] <- feature_data
  }

  plot_data_combined <- bind_rows(plot_data_list)
  feature_order <- top_features %>%
    arrange(desc(LDA_avg)) %>%
    mutate(Clean_Name = sapply(Taxa, function(x) get_clean_name(x)))
  plot_data_combined$Taxon <- factor(plot_data_combined$Taxon, levels = rev(feature_order$Clean_Name))

  n_features <- length(levels(plot_data_combined$Taxon))
  if (n_features >= 2) {
    alternating_bg <- data.frame(
      ymin = seq(1.5, n_features - 0.5, by = 1),
      ymax = seq(2.5, n_features + 0.5, by = 1)
    )
    alternating_bg <- alternating_bg[seq(1, nrow(alternating_bg), by = 2), , drop = FALSE]
  } else {
    alternating_bg <- data.frame(ymin = numeric(0), ymax = numeric(0))
  }

  p_left <- ggplot(plot_data_combined, aes(x = Log_Abundance, y = Taxon, color = Group, fill = Group)) +
    geom_rect(data = alternating_bg, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      fill = "gray85", alpha = 0.5, inherit.aes = FALSE) +
    geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA, width = 0.35, alpha = 0.3, linewidth = 0.5) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
    guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE)) +
    labs(
      x = expression(Log[10] ~ "Relative Abundance"),
      y = NULL,
      title = paste0(nrow(top_features), " Differential eggNOG Class Features at ", class_level),
      subtitle = "Filtered by HQQD-induced Reversal, Sorted by LDA"
    ) +
    theme_bw(base_family = "sans") +
    theme(
      axis.text.y = element_text(size = 10, hjust = 1, face = "bold", color = "black"),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12, face = "bold"),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(color = "black", linewidth = 0.5),
      axis.line.y.left = element_line(color = "black", linewidth = 0.5)
    ) +
    geom_segment(aes(x = -Inf, xend = Inf, y = Inf, yend = Inf), color = "black", linewidth = 0.5, inherit.aes = FALSE)

  p_middle <- ggplot(plot_data_combined, aes(x = 1, y = Taxon)) +
    geom_rect(data = alternating_bg, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      fill = "gray85", alpha = 0.5, inherit.aes = FALSE) +
    scale_y_discrete(limits = levels(plot_data_combined$Taxon)) +
    theme_void(base_family = "sans") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(10, 10, 10, 5), panel.border = element_blank()) +
    geom_segment(aes(x = Inf, xend = Inf, y = -Inf, yend = Inf), color = "black", linewidth = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = -Inf, xend = Inf, y = Inf, yend = Inf), color = "black", linewidth = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = -Inf, xend = Inf, y = -Inf, yend = -Inf), color = "black", linewidth = 0.5, inherit.aes = FALSE) +
    coord_cartesian(xlim = c(0, 10), clip = "off")

  # Perform three Wilcoxon rank-sum tests per feature (for plot annotations)
  # Comparisons: Model vs Control | HQQD vs Model | Model vs LVX
  plot_stat_results <- list()
  for (feature in levels(plot_data_combined$Taxon)) {
    fd <- plot_data_combined %>% filter(Taxon == feature)
    mc_v   <- fd$Log_Abundance[fd$Group == "Model"]
    ctrl_v <- fd$Log_Abundance[fd$Group == "Control"]
    hqqd_v <- fd$Log_Abundance[fd$Group == "HQQD"]
    lvx_v  <- fd$Log_Abundance[fd$Group == "LVX"]
    plot_stat_results[[feature]] <- list(
      p_Model_vs_Control = tryCatch(wilcox.test(mc_v, ctrl_v, exact = FALSE)$p.value, error = function(e) NA),
      p_HQQD_vs_Model    = tryCatch(wilcox.test(hqqd_v, mc_v, exact = FALSE)$p.value, error = function(e) NA),
      p_Model_vs_LVX     = tryCatch(wilcox.test(mc_v, lvx_v,  exact = FALSE)$p.value, error = function(e) NA)
    )
  }

  comparisons_fixed <- list(
    list(group1 = "Model", group2 = "Control", p_key = "p_Model_vs_Control"),
    list(group1 = "HQQD",  group2 = "Model",   p_key = "p_HQQD_vs_Model"),
    list(group1 = "Model", group2 = "LVX",     p_key = "p_Model_vs_LVX")
  )
  group_positions <- c("HQQD" = 0, "LVX" = 1, "Model" = 2, "Control" = 3)
  dodge_width <- 0.7
  n_groups    <- 4

  for (feature_name in levels(plot_data_combined$Taxon)) {
    feature_idx <- which(levels(plot_data_combined$Taxon) == feature_name)
    stat_res    <- plot_stat_results[[feature_name]]
    current_x   <- 2

    for (comp in comparisons_fixed) {
      p_val <- stat_res[[comp$p_key]]
      if (is.na(p_val) || p_val >= 0.05) {
        current_x <- current_x + 1.2
        next
      }

      sig_label <- dplyr::case_when(p_val < 0.001 ~ "***", p_val < 0.01 ~ "**", TRUE ~ "*")
      g1_idx <- group_positions[comp$group1]
      g2_idx <- group_positions[comp$group2]
      y1 <- feature_idx + (g1_idx - (n_groups - 1) / 2) * dodge_width / n_groups
      y2 <- feature_idx + (g2_idx - (n_groups - 1) / 2) * dodge_width / n_groups
      star_spacing <- dplyr::case_when(sig_label == "***" ~ 2.6, sig_label == "**" ~ 2.0, TRUE ~ 1.2)

      p_middle <- p_middle +
        annotate("segment", x = current_x, xend = current_x, y = y1, yend = y2, color = "gray30", linewidth = 0.6) +
        annotate("text", x = current_x + 0.3, y = (y1 + y2) / 2, label = sig_label,
                 size = 3.5, fontface = "bold", color = "gray20", hjust = 0)

      current_x <- current_x + star_spacing
    }
  }

  p_legend_temp <- ggplot(plot_data_combined, aes(x = Log_Abundance, y = Taxon, color = Group, fill = Group)) +
    geom_boxplot(alpha = 0.3) +
    scale_color_manual(values = group_colors, name = "Group") +
    scale_fill_manual(values = group_colors, name = "Group") +
    guides(color = guide_legend(reverse = TRUE, override.aes = list(fill = group_colors, alpha = 0.3)), fill = "none") +
    theme_minimal(base_family = "sans") +
    theme(legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 10), legend.key.size = unit(0.8, "cm"))

  p_legend <- as_ggplot(get_legend(p_legend_temp))
  p_combined <- p_left + p_middle + p_legend + plot_layout(widths = c(5, 1, 0.8))

  output_file_pdf <- file.path(output_dir, paste0("LEfSe_boxplot_eggNOG_Class_", class_level, ".pdf"))
  output_file_png <- file.path(output_dir, paste0("LEfSe_boxplot_eggNOG_Class_", class_level, ".png"))

  plot_height <- max(6, nrow(top_features) * 0.55 + 2)
  plot_width <- 8.5

  ggsave(filename = output_file_pdf, plot = p_combined, width = plot_width, height = plot_height, units = "in", dpi = 300)
  ggsave(filename = output_file_png, plot = p_combined, width = plot_width, height = plot_height, units = "in", dpi = 300)

  stat_summary <- data.frame(
    Taxon              = names(all_stat_results),
    Taxa_Full          = sapply(all_stat_results, function(x) x$Taxa_Full),
    p_Model_vs_Control = sapply(all_stat_results, function(x) x$p_Model_vs_Control),
    p_HQQD_vs_Model    = sapply(all_stat_results, function(x) x$p_HQQD_vs_Model),
    p_Model_vs_LVX     = sapply(all_stat_results, function(x) x$p_Model_vs_LVX),
    stringsAsFactors = FALSE
  )

  lda_info <- sig_results %>%
    mutate(Clean_Name = sapply(Taxa, function(x) get_clean_name(x))) %>%
    select(
      Clean_Name, LDA_Model_Control, LDA_HQQD_Model, LDA_avg,
      Trend_Model_vs_Control, Trend_HQQD_vs_Model,
      q_val_Model_Control, q_val_HQQD_Model, is_callback_fdr
    ) %>%
    rename(
      q_Model_vs_Control = q_val_Model_Control,
      q_HQQD_vs_Model = q_val_HQQD_Model,
      FDR_corrected_callback = is_callback_fdr
    )

  stat_summary <- stat_summary %>%
    left_join(lda_info, by = c("Taxon" = "Clean_Name")) %>%
    arrange(desc(LDA_avg))

  class_meta <- abundance_data %>%
    select(eggNOG_Class) %>%
    distinct(eggNOG_Class, .keep_all = TRUE)

  stat_summary <- stat_summary %>%
    left_join(class_meta, by = c("Taxon" = "eggNOG_Class"))

  final_cols <- c(
    "Taxon", "Taxa_Full",
    "p_Model_vs_Control", "p_HQQD_vs_Model", "p_Model_vs_LVX",
    "q_Model_vs_Control", "q_HQQD_vs_Model", "FDR_corrected_callback",
    "LDA_Model_Control", "LDA_HQQD_Model", "LDA_avg",
    "Trend_Model_vs_Control", "Trend_HQQD_vs_Model"
  )
  stat_summary <- stat_summary %>% select(all_of(final_cols[final_cols %in% colnames(stat_summary)]))

  stat_summary_fdr <- stat_summary %>% filter(FDR_corrected_callback == TRUE)

  write.csv(stat_summary,
    file.path(output_dir, paste0("statistical_results_eggNOG_Class_", class_level, ".csv")),
    row.names = FALSE
  )

  write.csv(stat_summary_fdr,
    file.path(output_dir, paste0("statistical_results_eggNOG_Class_", class_level, "_fdr.csv")),
    row.names = FALSE
  )

  if (nrow(stat_summary_fdr) > 0) {
    fdr_order <- stat_summary_fdr %>% arrange(desc(LDA_avg)) %>% pull(Taxon)
    plot_data_fdr <- plot_data_combined %>%
      filter(as.character(Taxon) %in% fdr_order)
    plot_data_fdr$Taxon <- factor(as.character(plot_data_fdr$Taxon), levels = rev(fdr_order))
    save_reverse_boxplot(
      plot_data_fdr,
      file.path(output_dir, paste0("LEfSe_boxplot_eggNOG_Class_", class_level, "_fdr.pdf")),
      file.path(output_dir, paste0("LEfSe_boxplot_eggNOG_Class_", class_level, "_fdr.png")),
      paste0(nrow(stat_summary_fdr), " FDR-corrected Differential eggNOG Class Features at ", class_level),
      plot_width = plot_width
    )
  }
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("All eggNOG Class boxplots generated successfully!\n")
cat("Output directory:", output_dir, "\n")
cat("\nGenerated files:\n")
cat("  - LEfSe_boxplot_eggNOG_Class_Class.pdf/png\n")
cat("  - LEfSe_boxplot_eggNOG_Class_Class_fdr.pdf/png\n")
cat("  - statistical_results_eggNOG_Class_Class.csv (Wilcoxon p-values: Model vs Control, HQQD vs Model, Model vs LVX)\n")
cat("  - statistical_results_eggNOG_Class_Class_fdr.csv (BH-corrected callback subset)\n")
cat(rep("=", 70), "\n\n", sep = "")
