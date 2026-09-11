# ============================================================================
# Statistical analysis of animal outcomes
# Suffix for all outputs: _corrected_6sample
#
# Analysis rules:
# 1. No automatic outlier removal.
# 2. Normal endpoints: one-way ANOVA, post hoc Bonferroni comparisons.
# 3. Non-normal endpoints: Kruskal-Wallis test, post hoc median test.
# 4. Manuscript figures are plotted as box plots overlaid with individual points.
# 5. Statistical summaries and pairwise tests are unchanged by plotting style.
# 6. Pairwise outputs focus on each group vs Model for manuscript readability.
# ============================================================================

rm(list = ls())

required_packages <- c("dplyr", "tidyr", "ggplot2", "cowplot", "stringr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(stringr)
})

suffix <- "_corrected_6sample"
alpha <- 0.05
group_levels <- c("Control", "Model", "LVX", "HQQD-L", "HQQD-M", "HQQD-H")
comparison_ref <- "Model"

group_colors <- c(
  "Control" = "#3B4992",
  "Model" = "#EE0000",
  "LVX" = "#008B45",
  "HQQD-L" = "#631879",
  "HQQD-M" = "#E69F00",
  "HQQD-H" = "#56B4E9"
)

figure_dir <- file.path("results", "figures")
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

panel_width <- 12 / 5
panel_height <- 8 / 3

save_plot_svg <- function(filename, plot, width, height) {
  grDevices::svg(filename = filename, width = width, height = height, onefile = FALSE)
  print(plot)
  grDevices::dev.off()
}

read_csv_flexible <- function(path) {
  encodings <- c("UTF-8-BOM", "UTF-8", "GB18030", "GBK")
  last_error <- NULL
  for (enc in encodings) {
    dat <- tryCatch(
      suppressWarnings(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = enc)),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(dat) && ncol(dat) > 1) {
      attr(dat, "source_encoding") <- enc
      return(dat)
    }
  }
  stop("Could not read CSV file: ", path, "; last error: ", conditionMessage(last_error))
}

write_csv_utf8 <- function(x, path, row.names = FALSE) {
  write.csv(x, path, row.names = row.names, fileEncoding = "UTF-8")
}

find_disease_data <- function() {
  six_sample_input <- file.path("data", "animal_outcomes.csv")
  if (file.exists(six_sample_input)) {
    dat <- read_csv_flexible(six_sample_input)
    if ("Group" %in% names(dat) && "Lung Bacterial Load (log10 CFU)" %in% names(dat)) {
      return(dat)
    }
  }

  stop("Required input not found: data/animal_outcomes.csv")
}

clean_percent <- function(x) {
  as.numeric(str_replace(as.character(x), "%", ""))
}

normalize_cytokine_marker <- function(marker) {
  marker <- str_replace_all(marker, fixed("?"), "")
  marker <- str_replace_all(marker, "¦Į|α", "alpha")
  marker <- str_replace_all(marker, "¦Ā|β", "beta")
  marker <- str_replace_all(marker, "[- ]", "")
  marker <- toupper(marker)
  case_when(
    marker %in% c("TNF", "TNFALPHA") ~ "TNF_ALPHA",
    marker %in% c("IL1", "IL1BETA") ~ "IL1_BETA",
    marker == "IL17A" ~ "IL17A",
    marker == "IL10" ~ "IL10",
    marker == "IL6" ~ "IL6",
    TRUE ~ marker
  )
}

format_cytokine_marker <- function(marker_key) {
  case_when(
    marker_key == "TNF_ALPHA" ~ "TNF-\u03b1",
    marker_key == "IL1_BETA" ~ "IL-1\u03b2",
    marker_key == "IL17A" ~ "IL-17A",
    marker_key == "IL10" ~ "IL-10",
    marker_key == "IL6" ~ "IL-6",
    TRUE ~ marker_key
  )
}

safe_shapiro_p <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3 || length(x) > 5000 || length(unique(x)) < 2) return(NA_real_)
  tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
}

safe_anova_p <- function(dat) {
  tryCatch(summary(aov(Value ~ Group, data = dat))[[1]]$`Pr(>F)`[1], error = function(e) NA_real_)
}

safe_kruskal_p <- function(dat) {
  tryCatch(kruskal.test(Value ~ Group, data = dat)$p.value, error = function(e) NA_real_)
}

format_p <- function(p) {
  ifelse(is.na(p), NA_character_,
    ifelse(p < 0.0001, "<0.0001", sprintf("%.4f", p))
  )
}

p_to_star <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.0001) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

median_ci <- function(x, conf = 0.95, nboot = 2000) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(c(Median = NA_real_, CI_low = NA_real_, CI_high = NA_real_))
  if (length(unique(x)) == 1) return(c(Median = median(x), CI_low = median(x), CI_high = median(x)))
  boot_medians <- replicate(nboot, median(sample(x, length(x), replace = TRUE)))
  c(
    Median = median(x),
    CI_low = as.numeric(quantile(boot_medians, (1 - conf) / 2, na.rm = TRUE)),
    CI_high = as.numeric(quantile(boot_medians, 1 - (1 - conf) / 2, na.rm = TRUE))
  )
}

pairwise_bonferroni_vs_model <- function(dat) {
  dat <- dat %>%
    filter(!is.na(Value), !is.na(Group)) %>%
    mutate(Group = factor(as.character(Group), levels = group_levels)) %>%
    filter(!is.na(Group)) %>%
    droplevels()

  groups <- setdiff(levels(dat$Group), comparison_ref)
  groups <- groups[groups %in% as.character(unique(dat$Group))]

  raw <- lapply(groups, function(g) {
    pair_dat <- dat %>% filter(Group %in% c(comparison_ref, g))
    if (length(unique(pair_dat$Group)) < 2) return(NULL)
    tt <- tryCatch(t.test(Value ~ Group, data = pair_dat, var.equal = TRUE), error = function(e) NULL)
    if (is.null(tt)) return(NULL)
    data.frame(
      group1 = g,
      group2 = comparison_ref,
      comparison = paste(g, comparison_ref, sep = " - "),
      estimate = diff(rev(tt$estimate)),
      statistic = unname(tt$statistic),
      p_raw = tt$p.value,
      stringsAsFactors = FALSE
    )
  })

  out <- bind_rows(raw)
  if (nrow(out) == 0) return(out)
  out %>%
    mutate(
      p_adjusted_within_endpoint = p.adjust(p_raw, method = "bonferroni"),
      p_adjust_method = "Bonferroni"
    )
}

median_test_pair <- function(x, g) {
  valid <- !is.na(x) & !is.na(g)
  x <- x[valid]
  g <- droplevels(factor(g[valid]))
  if (length(unique(g)) != 2 || length(unique(x)) < 2) return(c(statistic = NA_real_, p = NA_real_))

  pooled_median <- median(x)
  side <- ifelse(x <= pooled_median, "below_or_equal", "above")
  tab <- table(g, side)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(c(statistic = NA_real_, p = NA_real_))

  res <- tryCatch(suppressWarnings(chisq.test(tab, correct = FALSE)), error = function(e) NULL)
  if (is.null(res) || any(res$expected < 1)) {
    fisher <- tryCatch(fisher.test(tab), error = function(e) NULL)
    if (is.null(fisher)) return(c(statistic = NA_real_, p = NA_real_))
    return(c(statistic = NA_real_, p = fisher$p.value))
  }
  c(statistic = unname(res$statistic), p = res$p.value)
}

pairwise_median_vs_model <- function(dat) {
  dat <- dat %>%
    filter(!is.na(Value), !is.na(Group)) %>%
    mutate(Group = factor(as.character(Group), levels = group_levels)) %>%
    filter(!is.na(Group)) %>%
    droplevels()

  groups <- setdiff(levels(dat$Group), comparison_ref)
  groups <- groups[groups %in% as.character(unique(dat$Group))]

  raw <- lapply(groups, function(g) {
    pair_dat <- dat %>% filter(Group %in% c(comparison_ref, g))
    if (length(unique(pair_dat$Group)) < 2) return(NULL)
    mt <- median_test_pair(pair_dat$Value, pair_dat$Group)
    data.frame(
      group1 = g,
      group2 = comparison_ref,
      comparison = paste(g, comparison_ref, sep = " - "),
      estimate = median(pair_dat$Value[pair_dat$Group == g], na.rm = TRUE) -
        median(pair_dat$Value[pair_dat$Group == comparison_ref], na.rm = TRUE),
      statistic = unname(mt["statistic"]),
      p_raw = unname(mt["p"]),
      p_adjusted_within_endpoint = unname(mt["p"]),
      p_adjust_method = "Pairwise median test; unadjusted",
      stringsAsFactors = FALSE
    )
  })

  bind_rows(raw)
}

analyse_endpoint <- function(dat, variable_name, analysis_name, endpoint_type) {
  endpoint <- dat %>%
    filter(Variable == variable_name, !is.na(Value), !is.na(Group)) %>%
    mutate(Group = factor(as.character(Group), levels = group_levels)) %>%
    filter(!is.na(Group)) %>%
    droplevels()

  group_counts <- endpoint %>%
    count(Group, name = "N") %>%
    complete(Group = factor(group_levels, levels = group_levels), fill = list(N = 0)) %>%
    mutate(Group = as.character(Group))

  valid_group_count <- sum(group_counts$N > 0)
  if (valid_group_count < 2) {
    return(list(
      summary = data.frame(
        Analysis = analysis_name,
        Variable = variable_name,
        Endpoint_Type = endpoint_type,
        Recommended_Test = "Not tested",
        P_global = NA_real_,
        Normality_All_Groups = NA,
        Min_N = min(group_counts$N[group_counts$N > 0], na.rm = TRUE),
        Max_N = max(group_counts$N, na.rm = TRUE),
        Plot_Type = NA_character_,
        Notes = "Fewer than two non-empty groups",
        stringsAsFactors = FALSE
      ),
      pairwise = data.frame(),
      normality = data.frame()
    ))
  }

  normality <- endpoint %>%
    group_by(Group) %>%
    summarise(N = n(), Shapiro_P = safe_shapiro_p(Value), .groups = "drop") %>%
    mutate(Variable = variable_name, Analysis = analysis_name) %>%
    dplyr::select(Analysis, Variable, Group, N, Shapiro_P)

  all_shapiro_available <- all(!is.na(normality$Shapiro_P))
  all_normal <- all_shapiro_available && all(normality$Shapiro_P > alpha)

  if (all_normal) {
    recommended_test <- "One-way ANOVA + Bonferroni vs Model"
    global_p <- safe_anova_p(endpoint)
    pairwise <- pairwise_bonferroni_vs_model(endpoint)
    plot_type <- "scatter plot with mean +/- SD"
    notes <- "All group Shapiro p > 0.05; parametric workflow"
  } else {
    recommended_test <- "Kruskal-Wallis + pairwise median test vs Model"
    global_p <- safe_kruskal_p(endpoint)
    pairwise <- pairwise_median_vs_model(endpoint)
    plot_type <- "scatter plot with median and 95% CI"
    notes <- "Normality not established; nonparametric workflow"
  }

  if (nrow(pairwise) > 0) {
    pairwise <- pairwise %>%
      mutate(
        Analysis = analysis_name,
        Variable = variable_name,
        Endpoint_Type = endpoint_type,
        Recommended_Test = recommended_test,
        Significant_pairwise_within_endpoint = !is.na(p_adjusted_within_endpoint) &
          p_adjusted_within_endpoint < alpha,
        .before = 1
      )
  }

  summary <- data.frame(
    Analysis = analysis_name,
    Variable = variable_name,
    Endpoint_Type = endpoint_type,
    Recommended_Test = recommended_test,
    P_global = global_p,
    Normality_All_Groups = all_normal,
    Min_N = min(group_counts$N[group_counts$N > 0], na.rm = TRUE),
    Max_N = max(group_counts$N, na.rm = TRUE),
    Plot_Type = plot_type,
    Notes = notes,
    P_global_adjust_method = "none; one global test per endpoint",
    Significant_global = !is.na(global_p) & global_p < alpha,
    stringsAsFactors = FALSE
  )

  list(summary = summary, pairwise = pairwise, normality = normality)
}

descriptive_stats <- function(long_data) {
  set.seed(20260615)
  long_data %>%
    filter(!is.na(Value), !is.na(Group)) %>%
    mutate(Group = factor(as.character(Group), levels = group_levels)) %>%
    filter(!is.na(Group)) %>%
    group_by(Analysis, Variable, Group) %>%
    summarise(
      N = n(),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      SEM = SD / sqrt(N),
      Median = median(Value, na.rm = TRUE),
      Median_CI_low = median_ci(Value)["CI_low"],
      Median_CI_high = median_ci(Value)["CI_high"],
      Q1 = quantile(Value, 0.25, na.rm = TRUE),
      Q3 = quantile(Value, 0.75, na.rm = TRUE),
      Min = min(Value, na.rm = TRUE),
      Max = max(Value, na.rm = TRUE),
      .groups = "drop"
    )
}

make_plot <- function(endpoint_data, summary_row, pairwise_rows, y_label, panel_title = NULL) {
  endpoint_data <- endpoint_data %>%
    mutate(Group = factor(as.character(Group), levels = group_levels)) %>%
    filter(!is.na(Group), !is.na(Value))

  y_range <- range(endpoint_data$Value, na.rm = TRUE)
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span == 0) y_span <- max(abs(y_range), 1)

  sig_rows <- pairwise_rows %>%
    filter(Significant_pairwise_within_endpoint) %>%
    mutate(
      x1 = match(group2, group_levels),
      x2 = match(group1, group_levels),
      label = vapply(p_adjusted_within_endpoint, p_to_star, character(1)),
      y = y_range[2] + y_span * (0.12 + 0.10 * (row_number() - 1))
    ) %>%
    filter(!is.na(x1), !is.na(x2))

  p <- ggplot(endpoint_data, aes(x = Group, y = Value, color = Group, fill = Group)) +
    geom_boxplot(
      width = 0.56,
      outlier.shape = NA,
      linewidth = 0.42,
      alpha = 0.18,
      color = "black"
    ) +
    geom_jitter(width = 0.10, height = 0, size = 1.55, alpha = 0.92, stroke = 0)

  p <- p +
    scale_color_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE)

  p <- p +
    labs(x = NULL, y = y_label, title = panel_title) +
    theme_classic(base_family = "sans") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.title.y = element_text(size = 9, color = "black"),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      plot.title = element_text(size = 9, hjust = 0.5, face = "bold", color = "black"),
      plot.margin = margin(6, 8, 6, 8)
    )

  if (nrow(sig_rows) > 0) {
    p <- p +
      geom_segment(
        data = sig_rows,
        aes(x = x1, xend = x2, y = y, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      geom_segment(
        data = sig_rows,
        aes(x = x1, xend = x1, y = y - y_span * 0.025, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      geom_segment(
        data = sig_rows,
        aes(x = x2, xend = x2, y = y - y_span * 0.025, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      geom_text(
        data = sig_rows,
        aes(x = (x1 + x2) / 2, y = y + y_span * 0.025, label = label),
        inherit.aes = FALSE,
        size = 3,
        fontface = "bold"
      ) +
      expand_limits(y = max(sig_rows$y, na.rm = TRUE) + y_span * 0.10)
  }

  p
}

analyse_block <- function(long_data, endpoint_plan, analysis_name, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  descriptive <- descriptive_stats(long_data)

  results <- lapply(seq_len(nrow(endpoint_plan)), function(i) {
    analyse_endpoint(
      dat = long_data,
      variable_name = endpoint_plan$Variable[i],
      analysis_name = analysis_name,
      endpoint_type = endpoint_plan$Endpoint_Type[i]
    )
  })

  global_summary <- bind_rows(lapply(results, `[[`, "summary"))
  pairwise_summary <- bind_rows(lapply(results, `[[`, "pairwise"))
  normality_summary <- bind_rows(lapply(results, `[[`, "normality"))

  write_csv_utf8(descriptive, file.path(output_dir, paste0("descriptive_statistics", suffix, ".csv")))
  write_csv_utf8(global_summary, file.path(output_dir, paste0("statistical_summary", suffix, ".csv")))
  write_csv_utf8(pairwise_summary, file.path(output_dir, paste0("pairwise_comparisons_vs_model", suffix, ".csv")))
  write_csv_utf8(normality_summary, file.path(output_dir, paste0("normality_by_group", suffix, ".csv")))

  list(
    descriptive = descriptive,
    global_summary = global_summary,
    pairwise_summary = pairwise_summary,
    normality_summary = normality_summary,
    long_data = long_data
  )
}

disease_data <- find_disease_data()
disease_data <- disease_data %>%
  filter(!is.na(Group), Group != "") %>%
  mutate(Group = factor(Group, levels = group_levels))

# ---------------------------------------------------------------------------
# Step 01: Cytokines
# ---------------------------------------------------------------------------
cytokine_cols <- names(disease_data)[2:16]
cytokine_long <- disease_data %>%
  dplyr::select(Group, all_of(cytokine_cols)) %>%
  pivot_longer(cols = all_of(cytokine_cols), names_to = "Variable", values_to = "Value") %>%
  mutate(
    Value = as.numeric(Value),
    Analysis = "01_cytokine",
    Endpoint_Type = "continuous_cytokine",
    Tissue = word(Variable, 1),
    Marker = str_remove(Variable, paste0(Tissue, " ")),
    Tissue = factor(Tissue, levels = c("Serum", "Lung", "Ileum")),
    Marker_Key = normalize_cytokine_marker(Marker),
    Marker_Display = format_cytokine_marker(Marker_Key),
    Panel_Title = paste(as.character(Tissue), Marker_Display)
  )

cytokine_plan <- data.frame(
  Variable = cytokine_cols,
  Endpoint_Type = "continuous_cytokine",
  stringsAsFactors = FALSE
)

res_cytokine <- analyse_block(
  long_data = cytokine_long,
  endpoint_plan = cytokine_plan,
  analysis_name = "01_cytokine",
  output_dir = file.path("results", "01_cytokine_analysis")
)

cytokine_plot_order <- cytokine_long %>%
  distinct(Variable, Tissue, Marker, Marker_Key, Panel_Title) %>%
  mutate(
    Tissue_Order = match(as.character(Tissue), c("Serum", "Lung", "Ileum")),
    Marker_Order = match(Marker_Key, c("TNF_ALPHA", "IL17A", "IL10", "IL6", "IL1_BETA"))
  ) %>%
  arrange(Tissue_Order, Marker_Order)

cytokine_plots <- lapply(cytokine_plot_order$Variable, function(var) {
  plot_data <- cytokine_long %>% filter(Variable == var)
  panel_title <- cytokine_plot_order %>% filter(Variable == var) %>% pull(Panel_Title) %>% .[1]
  summary_row <- res_cytokine$global_summary %>% filter(Variable == var) %>% slice(1)
  pairwise_rows <- res_cytokine$pairwise_summary %>% filter(Variable == var)
  make_plot(plot_data, summary_row, pairwise_rows, "Cytokine concentration (pg/mL)", panel_title)
})

names(cytokine_plots) <- cytokine_plot_order$Variable
cytokine_panel <- plot_grid(plotlist = cytokine_plots, ncol = 5, align = "hv")
ggsave(file.path("01_cytokine_analysis_6sample", "Combined_Cytokine_Panel_corrected_6sample.png"), cytokine_panel, width = 12, height = 8, dpi = 300)
ggsave(file.path(figure_dir, "Figure_1_Inflammatory_cytokine_profiles.png"), cytokine_panel, width = 12, height = 8, dpi = 300)
save_plot_svg(file.path(figure_dir, "Figure_1_Inflammatory_cytokine_profiles.svg"), cytokine_panel, width = 12, height = 8)

# ---------------------------------------------------------------------------
# Step 02: Lung bacterial load and pathology score
# ---------------------------------------------------------------------------
pathology_cols <- c("Lung Bacterial Load (log10 CFU)", "Lung Pathology Score")
pathology_long <- disease_data %>%
  dplyr::select(Group, all_of(pathology_cols)) %>%
  pivot_longer(cols = all_of(pathology_cols), names_to = "Variable", values_to = "Value") %>%
  mutate(
    Value = as.numeric(Value),
    Analysis = "02_lung_pathology",
    Endpoint_Type = ifelse(Variable == "Lung Pathology Score", "ordinal_score", "continuous_log10_cfu")
  )

pathology_plan <- data.frame(
  Variable = pathology_cols,
  Endpoint_Type = c("continuous_log10_cfu", "ordinal_score"),
  stringsAsFactors = FALSE
)

res_pathology <- analyse_block(
  long_data = pathology_long,
  endpoint_plan = pathology_plan,
  analysis_name = "02_lung_pathology",
  output_dir = file.path("results", "02_lung_pathology_analysis")
)

pathology_plots <- lapply(pathology_cols, function(var) {
  plot_data <- pathology_long %>% filter(Variable == var)
  summary_row <- res_pathology$global_summary %>% filter(Variable == var) %>% slice(1)
  pairwise_rows <- res_pathology$pairwise_summary %>% filter(Variable == var)
  y_lab <- ifelse(
    var == "Lung Bacterial Load (log10 CFU)",
    "Bacterial load (log10 CFU)",
    "Histopathological score"
  )
  panel_title <- ifelse(
    var == "Lung Bacterial Load (log10 CFU)",
    "Lung Bacterial Load",
    "Lung Pathology Score"
  )
  p <- make_plot(plot_data, summary_row, pairwise_rows, y_lab, panel_title)
  safe_name <- str_replace_all(var, "[^[:alnum:]]", "_")
  ggsave(file.path("02_lung_pathology_analysis_6sample", paste0(safe_name, "_corrected_6sample.png")), p, width = 3.2, height = 4.2, dpi = 300)
  p
})
names(pathology_plots) <- pathology_cols
ggsave(file.path(figure_dir, "Figure_2A_Lung_bacterial_load.png"), pathology_plots[[1]], width = panel_width, height = panel_height, dpi = 300)
save_plot_svg(file.path(figure_dir, "Figure_2A_Lung_bacterial_load.svg"), pathology_plots[[1]], width = panel_width, height = panel_height)
ggsave(file.path(figure_dir, "Figure_2B_Lung_pathology_score.png"), pathology_plots[[2]], width = panel_width, height = panel_height, dpi = 300)
save_plot_svg(file.path(figure_dir, "Figure_2B_Lung_pathology_score.svg"), pathology_plots[[2]], width = panel_width, height = panel_height)
pathology_panel <- plot_grid(plotlist = pathology_plots, ncol = 2, align = "hv")
ggsave(file.path("02_lung_pathology_analysis_6sample", "Combined_Lung_Bacterial_Load_and_Pathology_corrected_6sample.png"), pathology_panel, width = panel_width * 2, height = panel_height, dpi = 300)
ggsave(file.path(figure_dir, "Figure_2_Lung_bacterial_load_and_pathology.png"), pathology_panel, width = panel_width * 2, height = panel_height, dpi = 300)
save_plot_svg(file.path(figure_dir, "Figure_2_Lung_bacterial_load_and_pathology.svg"), pathology_panel, width = panel_width * 2, height = panel_height)

# ---------------------------------------------------------------------------
# Step 03: Flow cytometry
# ---------------------------------------------------------------------------
flow_col_map <- c(
  "MLN Treg" = names(disease_data)[stringr::str_detect(names(disease_data), fixed("MLN_treg_CD25+FOXP3"))][1],
  "MLN Th17" = names(disease_data)[stringr::str_detect(names(disease_data), fixed("MLN_TH17_CD8-IL-17A"))][1],
  "Lung Treg" = names(disease_data)[stringr::str_detect(names(disease_data), fixed("Lung treg_CD25+FOXP3"))][1],
  "Lung Th17" = names(disease_data)[stringr::str_detect(names(disease_data), fixed("Lung TH17_CD8-IL-17A"))][1]
)

if (any(is.na(flow_col_map))) {
  stop("Could not identify all required flow cytometry columns.")
}

flow_wide <- disease_data %>%
  dplyr::select(Group, all_of(unname(flow_col_map))) %>%
  dplyr::rename(
    `MLN Treg` = all_of(flow_col_map[["MLN Treg"]]),
    `MLN Th17` = all_of(flow_col_map[["MLN Th17"]]),
    `Lung Treg` = all_of(flow_col_map[["Lung Treg"]]),
    `Lung Th17` = all_of(flow_col_map[["Lung Th17"]])
  ) %>%
  mutate(
    across(c(`MLN Treg`, `MLN Th17`, `Lung Treg`, `Lung Th17`), clean_percent),
    `MLN Ratio` = `MLN Th17` / `MLN Treg`,
    `Lung Ratio` = `Lung Th17` / `Lung Treg`
  )

flow_vars <- c("MLN Treg", "MLN Th17", "MLN Ratio", "Lung Treg", "Lung Th17", "Lung Ratio")
flow_long <- flow_wide %>%
  pivot_longer(cols = all_of(flow_vars), names_to = "Variable", values_to = "Value") %>%
  mutate(
    Analysis = "03_flow_cytometry",
    Endpoint_Type = ifelse(str_detect(Variable, "Ratio"), "ratio", "percentage")
  )

flow_plan <- data.frame(
  Variable = flow_vars,
  Endpoint_Type = ifelse(str_detect(flow_vars, "Ratio"), "ratio", "percentage"),
  stringsAsFactors = FALSE
)

res_flow <- analyse_block(
  long_data = flow_long,
  endpoint_plan = flow_plan,
  analysis_name = "03_flow_cytometry",
  output_dir = file.path("results", "03_flow_cytometry_analysis")
)

flow_plots <- lapply(flow_vars, function(var) {
  plot_data <- flow_long %>% filter(Variable == var)
  summary_row <- res_flow$global_summary %>% filter(Variable == var) %>% slice(1)
  pairwise_rows <- res_flow$pairwise_summary %>% filter(Variable == var)
  y_lab <- if (str_detect(var, "Ratio")) "Th17/Treg ratio" else "Cell proportion (%)"
  make_plot(plot_data, summary_row, pairwise_rows, y_lab, var)
})
names(flow_plots) <- flow_vars
flow_panel <- plot_grid(plotlist = flow_plots, ncol = 6, align = "hv")
flow_panel_width <- panel_width * length(flow_vars)
ggsave(file.path("03_flow_cytometry_analysis_6sample", "Combined_Flow_Panel_corrected_6sample.png"), flow_panel, width = flow_panel_width, height = panel_height, dpi = 300)
ggsave(file.path(figure_dir, "Figure_3_Flow_cytometry_profiles.png"), flow_panel, width = flow_panel_width, height = panel_height, dpi = 300)
save_plot_svg(file.path(figure_dir, "Figure_3_Flow_cytometry_profiles.svg"), flow_panel, width = flow_panel_width, height = panel_height)

combined_global <- bind_rows(
  res_cytokine$global_summary,
  res_pathology$global_summary,
  res_flow$global_summary
)

combined_pairwise <- bind_rows(
  res_cytokine$pairwise_summary,
  res_pathology$pairwise_summary,
  res_flow$pairwise_summary
)

combined_descriptive <- bind_rows(
  res_cytokine$descriptive,
  res_pathology$descriptive,
  res_flow$descriptive
)

combined_normality <- bind_rows(
  res_cytokine$normality_summary,
  res_pathology$normality_summary,
  res_flow$normality_summary
)

combined_raw_long <- bind_rows(
  res_cytokine$long_data,
  res_pathology$long_data,
  res_flow$long_data
) %>%
  filter(!is.na(Value), !is.na(Group)) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    Sample_ID = ave(as.character(Group), Analysis, Variable, Group, FUN = seq_along),
    Sample_Label = paste(as.character(Group), Sample_ID, sep = "_")
  ) %>%
  dplyr::select(Analysis, Variable, Endpoint_Type, Group, Sample_ID, Sample_Label, Value)

publication_descriptive <- combined_descriptive %>%
  left_join(
    combined_global %>% dplyr::select(Analysis, Variable, Endpoint_Type, Recommended_Test, Plot_Type, P_global, Significant_global),
    by = c("Analysis", "Variable")
  ) %>%
  mutate(
    Statistic_Display = ifelse(
      Plot_Type == "scatter plot with mean +/- SD",
      paste0("n=", N, "; ", sprintf("%.3f", Mean), " +/- ", sprintf("%.3f", SD)),
      paste0(
        "n=", N, "; ", sprintf("%.3f", Median), " (",
        sprintf("%.3f", Median_CI_low), "-", sprintf("%.3f", Median_CI_high), ")"
      )
    ),
    P_global_display = format_p(as.numeric(P_global)),
    Significant_global_display = ifelse(Significant_global, "Yes", "No")
  )

publication_summary_table <- publication_descriptive %>%
  dplyr::select(
    Analysis,
    Variable,
    Endpoint_Type,
    Recommended_Test,
    Plot_Type,
    Group,
    Statistic_Display,
    P_global_display,
    Significant_global_display
  ) %>%
  pivot_wider(names_from = Group, values_from = Statistic_Display) %>%
  arrange(Analysis, factor(Variable, levels = unique(combined_global$Variable)))

publication_pairwise_p_table <- combined_pairwise %>%
  mutate(
    p_display = format_p(as.numeric(p_adjusted_within_endpoint)),
    significance = vapply(as.numeric(p_adjusted_within_endpoint), p_to_star, character(1)),
    p_with_stars = paste0(p_display, " ", significance)
  ) %>%
  dplyr::select(
    Analysis,
    Variable,
    Endpoint_Type,
    Recommended_Test,
    comparison,
    p_raw,
    p_adjusted_within_endpoint,
    p_adjust_method,
    p_with_stars,
    Significant_pairwise_within_endpoint
  ) %>%
  arrange(Analysis, factor(Variable, levels = unique(combined_global$Variable)), comparison)

publication_pairwise_p_wide <- publication_pairwise_p_table %>%
  dplyr::select(Analysis, Variable, comparison, p_with_stars) %>%
  pivot_wider(names_from = comparison, values_from = p_with_stars) %>%
  left_join(
    combined_global %>%
      mutate(P_global_display = format_p(as.numeric(P_global))) %>%
      dplyr::select(Analysis, Variable, Recommended_Test, P_global_display, Significant_global),
    by = c("Analysis", "Variable")
  ) %>%
  arrange(Analysis, factor(Variable, levels = unique(combined_global$Variable)))

cat("\nAnimal-outcome analysis complete.\n")
cat("Outputs use suffix: ", suffix, "\n", sep = "")
cat("Step outputs are stored in their corresponding analysis folders only.\n")

q(save = "no", status = 0, runLast = FALSE)
