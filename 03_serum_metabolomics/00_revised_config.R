################################################################################
# Revised shared analysis configuration
# Purpose:
#   - Keep grouping, QC filtering, univariate tests, OPLS/PLS validation, and
#     HQQD reversal rules consistent across scripts.
#   - Preserve old result folders by writing revised outputs to *_revised /
#     *_validated folders.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggrepel)
})

save_plot_svg <- function(filename, plot, width, height) {
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = filename, width = width, height = height)
  } else {
    grDevices::svg(filename = filename, width = width, height = height, onefile = FALSE)
  }
  print(plot)
  grDevices::dev.off()
}

comparisons <- list(
  c("Control", "Model"),
  c("Control", "LVX"),
  c("Control", "HQQD"),
  c("Model", "LVX"),
  c("Model", "HQQD"),
  c("LVX", "HQQD")
)

group_levels <- c("Control", "Model", "LVX", "HQQD")
group_colors <- c(
  "Control" = "#3B4992",
  "Model" = "#EE0000",
  "LVX" = "#008B45",
  "HQQD" = "#FF8C00"
)

analysis_cutoffs <- list(
  qc_cv = 0.30,
  vip = 1.0,
  general_fc = 2.0,
  model_fc = 1.5,
  treatment_fc = 1.2,
  p = 0.05,
  fdr_exploratory = 0.10
)

mode_config <- list(
  NEG = list(
    label = "Neg",
    ion_mode = "NEG",
    data_file = "metabolites_exp_neg.csv",
    # The source CSV is UTF-8. Reading it as GB18030 silently corrupts Greek
    # characters in annotations (for example, alpha/beta in neg_5611).
    encoding = "UTF-8",
    opls_dir = "01a_OPLS-DA_neg_validated",
    fc_dir = "02a_FC_neg_revised",
    p_dir = "03a_p_neg_welch",
    merge_dir = "04a_neg_merged_revised",
    diff_dir = "05a_neg_diff_revised",
    file_prefix = "neg"
  ),
  POS = list(
    label = "Pos",
    ion_mode = "POS",
    data_file = "metabolites_exp_pos.csv",
    encoding = "UTF-8",
    opls_dir = "01b_OPLS-DA_pos_validated",
    fc_dir = "02b_FC_pos_revised",
    p_dir = "03b_p_pos_welch",
    merge_dir = "04b_pos_merged_revised",
    diff_dir = "05b_pos_diff_revised",
    file_prefix = "pos"
  )
)

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

sanitise_utf8 <- function(x) {
  iconv(as.character(x), from = "UTF-8", to = "UTF-8", sub = "\uFFFD")
}

write_csv_utf8 <- function(x, path) {
  ensure_dir(dirname(path))
  out <- x
  char_cols <- vapply(out, is.character, logical(1))
  out[char_cols] <- lapply(out[char_cols], sanitise_utf8)
  colnames(out) <- sanitise_utf8(colnames(out))
  write.csv(out, path, row.names = FALSE, fileEncoding = "UTF-8")
}

read_csv_mode <- function(path, mode = NULL, row.names = NULL) {
  encodings <- c()
  if (!is.null(mode)) encodings <- c(encodings, mode_config[[mode]]$encoding)
  encodings <- unique(c(encodings, "UTF-8", "GB18030", "latin1"))

  for (enc in encodings) {
    res <- tryCatch(
      suppressWarnings(
        readr::read_csv(
          path,
          locale = readr::locale(encoding = enc),
          show_col_types = FALSE,
          name_repair = "minimal"
        ) %>%
          as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)
      ),
      error = function(e) NULL
    )
    if (!is.null(res) && ncol(res) > 1) {
      char_cols <- vapply(res, is.character, logical(1))
      res[char_cols] <- lapply(res[char_cols], sanitise_utf8)
      colnames(res) <- sanitise_utf8(colnames(res))
      if (!is.null(row.names)) {
        row_col <- if (is.numeric(row.names)) colnames(res)[row.names] else row.names
        rownames(res) <- make.unique(as.character(res[[row_col]]))
        res[[row_col]] <- NULL
      }
      return(res)
    }
  }

  stop(sprintf("Failed to read CSV file: %s", path))
}

exclusion_manifest_file <- "excluded_metabolite_annotations.csv"

load_metabolite_exclusion_manifest <- function(path = exclusion_manifest_file) {
  if (!file.exists(path)) {
    stop("Metabolite exclusion manifest not found: ", path)
  }

  manifest <- read_csv_mode(path)
  required_cols <- c("ID", "Reported_Annotation", "Exclusion_Reason", "Decision", "Exclusion_Set_Version")
  missing_cols <- setdiff(required_cols, colnames(manifest))
  if (length(missing_cols) > 0) {
    stop("Metabolite exclusion manifest is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  manifest %>%
    filter(Decision == "exclude_from_biological_interpretation") %>%
    distinct(ID, .keep_all = TRUE)
}

annotate_metabolite_exclusions <- function(df) {
  manifest <- load_metabolite_exclusion_manifest()
  df %>%
    left_join(manifest, by = "ID") %>%
    mutate(Excluded_Annotation = !is.na(Decision))
}

sample_group <- function(sample_name) {
  if (grepl("^K", sample_name)) return("Control")
  if (grepl("^M", sample_name)) return("Model")
  if (grepl("^Y", sample_name)) return("LVX")
  if (grepl("^Z", sample_name)) return("HQQD")
  if (grepl("^QC", sample_name)) return("QC")
  return(NA_character_)
}

as_numeric_df <- function(df) {
  as.data.frame(lapply(df, function(x) as.numeric(as.character(x))), check.names = FALSE)
}

impute_zero_missing <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  min_val <- suppressWarnings(min(mat[mat > 0], na.rm = TRUE))
  if (!is.finite(min_val)) min_val <- 1
  mat[is.na(mat) | mat <= 0] <- min_val / 2
  mat
}

calc_qc_cv <- function(raw_df, id_col, qc_cols) {
  if (length(qc_cols) == 0) {
    return(data.frame(ID = raw_df[[id_col]], QC_CV = NA_real_))
  }
  qc_mat <- as.matrix(as_numeric_df(raw_df[, qc_cols, drop = FALSE]))
  qc_mean <- rowMeans(qc_mat, na.rm = TRUE)
  qc_sd <- apply(qc_mat, 1, sd, na.rm = TRUE)
  qc_cv <- qc_sd / qc_mean
  qc_cv[!is.finite(qc_cv)] <- NA_real_
  data.frame(ID = raw_df[[id_col]], QC_CV = qc_cv)
}

load_mode_data <- function(mode) {
  cfg <- mode_config[[mode]]
  raw <- read_csv_mode(cfg$data_file, mode = mode)

  id_col <- colnames(raw)[1]
  name_col <- if ("name" %in% colnames(raw)) "name" else NA_character_
  abundance_cols <- setdiff(colnames(raw), c(id_col, name_col))
  col_groups <- vapply(abundance_cols, sample_group, character(1))

  sample_cols <- abundance_cols[!is.na(col_groups) & col_groups %in% group_levels]
  qc_cols <- abundance_cols[!is.na(col_groups) & col_groups == "QC"]

  abundance_df <- as_numeric_df(raw[, sample_cols, drop = FALSE])
  abundance_mat <- as.matrix(abundance_df)
  rownames(abundance_mat) <- raw[[id_col]]

  feature_info <- data.frame(
    ID = raw[[id_col]],
    Name = if (!is.na(name_col)) raw[[name_col]] else raw[[id_col]],
    stringsAsFactors = FALSE
  )
  feature_info <- left_join(feature_info, calc_qc_cv(raw, id_col, qc_cols), by = "ID")

  sample_info <- data.frame(
    Sample = sample_cols,
    Group = factor(vapply(sample_cols, sample_group, character(1)), levels = group_levels),
    stringsAsFactors = FALSE
  )

  list(
    cfg = cfg,
    raw = raw,
    feature_info = feature_info,
    abundance_mat = abundance_mat,
    sample_info = sample_info,
    qc_cols = qc_cols,
    sample_cols = sample_cols
  )
}

comparison_name <- function(group1, group2) {
  paste0(group1, "_vs_", group2)
}

fc_direction <- function(fc, cutoff) {
  case_when(
    is.na(fc) ~ "NA",
    fc > cutoff ~ "Up",
    fc < 1 / cutoff ~ "Down",
    TRUE ~ "Stable"
  )
}

safe_metric <- function(df, candidates) {
  if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) return(NA_real_)
  for (cn in candidates) {
    if (cn %in% colnames(df)) return(as.numeric(df[[cn]][nrow(df)]))
  }
  NA_real_
}

valid_ropls_object <- function(model_obj) {
  if (is.null(model_obj)) return(FALSE)
  if (is.null(model_obj@scoreMN) || nrow(model_obj@scoreMN) == 0 || ncol(model_obj@scoreMN) == 0) return(FALSE)
  if (is.null(model_obj@loadingMN) || nrow(model_obj@loadingMN) == 0 || ncol(model_obj@loadingMN) == 0) return(FALSE)
  TRUE
}

extract_ropls_validation <- function(model_obj, model_type, fallback_used) {
  if (!valid_ropls_object(model_obj)) {
    return(data.frame(
      ModelType = model_type,
      FallbackUsed = fallback_used,
      ModelObjectValid = FALSE,
      R2X_cum = NA_real_,
      R2Y_cum = NA_real_,
      Q2_cum = NA_real_,
      RMSEE = NA_real_,
      pR2Y = NA_real_,
      pQ2 = NA_real_,
      Q2_intercept = NA_real_,
      R2Y_intercept = NA_real_,
      VIP_valid = FALSE,
      Reason = "empty score/loading matrix",
      stringsAsFactors = FALSE
    ))
  }

  sum_df <- tryCatch(getSummaryDF(model_obj), error = function(e) data.frame())
  r2x <- safe_metric(sum_df, c("R2X(cum)", "R2X"))
  r2y <- safe_metric(sum_df, c("R2Y(cum)", "R2Y"))
  q2 <- safe_metric(sum_df, c("Q2(cum)", "Q2"))
  rmsee <- safe_metric(sum_df, c("RMSEE"))
  pr2y <- safe_metric(sum_df, c("pR2Y"))
  pq2 <- safe_metric(sum_df, c("pQ2"))

  q2_intercept <- NA_real_
  r2_intercept <- NA_real_
  perm <- model_obj@suppLs[["permMN"]]
  if (!is.null(perm) && nrow(perm) > 2) {
    sim_col <- if ("sim" %in% colnames(perm)) "sim" else colnames(perm)[ncol(perm)]
    q2_col <- grep("^Q2", colnames(perm), value = TRUE)[1]
    r2_col <- grep("^R2Y", colnames(perm), value = TRUE)[1]
    if (!is.na(q2_col) && !is.na(r2_col)) {
      q2_intercept <- tryCatch(as.numeric(coef(lm(perm[, q2_col] ~ perm[, sim_col]))[1]), error = function(e) NA_real_)
      r2_intercept <- tryCatch(as.numeric(coef(lm(perm[, r2_col] ~ perm[, sim_col]))[1]), error = function(e) NA_real_)
    }
  }

  q2_pass <- !is.na(q2) && q2 > 0
  perm_pass <- !is.na(pq2) && pq2 <= analysis_cutoffs$p
  vip_valid <- q2_pass && perm_pass
  reason <- if (vip_valid) {
    "validated"
  } else {
    paste(c(
      if (!q2_pass) "Q2<=0_or_missing" else NULL,
      if (!perm_pass) "permutation_pQ2>=0.05_or_missing" else NULL
    ), collapse = "; ")
  }

  data.frame(
    ModelType = model_type,
    FallbackUsed = fallback_used,
    ModelObjectValid = TRUE,
    R2X_cum = r2x,
    R2Y_cum = r2y,
    Q2_cum = q2,
    RMSEE = rmsee,
    pR2Y = pr2y,
    pQ2 = pq2,
    Q2_intercept = q2_intercept,
    R2Y_intercept = r2_intercept,
    VIP_valid = vip_valid,
    Reason = reason,
    stringsAsFactors = FALSE
  )
}

save_model_summary <- function(model_obj, validation, path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines("=== Multivariate Model Summary ===", con)
  writeLines(sprintf("ModelType: %s", validation$ModelType[1]), con)
  writeLines(sprintf("FallbackUsed: %s", validation$FallbackUsed[1]), con)
  writeLines(sprintf("ModelObjectValid: %s", validation$ModelObjectValid[1]), con)
  writeLines(sprintf("VIP_valid: %s", validation$VIP_valid[1]), con)
  writeLines(sprintf("Reason: %s", validation$Reason[1]), con)
  writeLines("", con)
  writeLines(capture.output(print(validation)), con)
  writeLines("", con)
  if (valid_ropls_object(model_obj)) {
    writeLines(capture.output(print(model_obj)), con)
  } else {
    writeLines("No valid score/loading matrix was produced.", con)
  }
}

plot_ropls_score <- function(model_obj, validation, group_subset, output_prefix, title_text) {
  if (!valid_ropls_object(model_obj)) return(invisible(FALSE))

  score <- getScoreMN(model_obj)
  t1 <- score[, 1]
  ortho_score <- tryCatch(getScoreMN(model_obj, orthoL = TRUE), error = function(e) NULL)
  has_ortho <- !is.null(ortho_score) && ncol(ortho_score) >= 1 && nrow(ortho_score) == length(t1)
  y <- if (has_ortho) ortho_score[, 1] else rep(0, length(t1))

  model_df <- model_obj@modelDF
  x_r2 <- if (!is.null(model_df) && "R2X" %in% colnames(model_df) && nrow(model_df) >= 1) model_df$R2X[1] * 100 else NA_real_
  y_r2 <- if (has_ortho && !is.null(model_df) && "R2X" %in% colnames(model_df) && nrow(model_df) >= 2) model_df$R2X[2] * 100 else NA_real_
  y_label <- if (has_ortho) {
    sprintf("Orthogonal score to[1] (R2X %.1f%%)", y_r2)
  } else {
    "Second score unavailable"
  }

  plot_df <- data.frame(
    t1 = t1,
    to1 = y,
    Group = group_subset,
    Sample = rownames(score),
    stringsAsFactors = FALSE
  )

  p <- ggplot(plot_df, aes(t1, to1, color = Group, fill = Group)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "gray60") +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", color = "gray60") +
    stat_ellipse(geom = "polygon", alpha = 0.14, linewidth = 0.4, level = 0.95, show.legend = FALSE) +
    geom_point(size = 3.2, shape = 21, color = "black", stroke = 0.35) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    coord_equal() +
    labs(
      x = sprintf("Predictive score t[1] (R2X %.1f%%)", x_r2),
      y = y_label,
      title = title_text,
      subtitle = sprintf(
        "%s | R2Y=%.3f Q2=%.3f pQ2=%s VIP=%s",
        validation$ModelType[1],
        validation$R2Y_cum[1],
        validation$Q2_cum[1],
        ifelse(is.na(validation$pQ2[1]), "NA", sprintf("%.3f", validation$pQ2[1])),
        ifelse(validation$VIP_valid[1], "valid", "not used")
      )
    ) +
    theme_bw(base_size = 12, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right",
      aspect.ratio = 1
    )

  ggsave(paste0(output_prefix, "_scores_plot.pdf"), p, width = 6.5, height = 6)
  ggsave(paste0(output_prefix, "_scores_plot.png"), p, width = 6.5, height = 6, dpi = 300)
  invisible(TRUE)
}

plot_permutation <- function(model_obj, validation, output_prefix, title_text) {
  if (!valid_ropls_object(model_obj)) return(invisible(FALSE))
  perm <- model_obj@suppLs[["permMN"]]
  if (is.null(perm) || nrow(perm) < 3) return(invisible(FALSE))

  perm <- as.data.frame(perm)
  sim_col <- if ("sim" %in% colnames(perm)) "sim" else colnames(perm)[ncol(perm)]
  q2_col <- grep("^Q2", colnames(perm), value = TRUE)[1]
  r2_col <- grep("^R2Y", colnames(perm), value = TRUE)[1]
  if (is.na(q2_col) || is.na(r2_col)) return(invisible(FALSE))

  perm_export <- data.frame(
    Run = seq_len(nrow(perm)),
    RunType = ifelse(seq_len(nrow(perm)) == 1, "Original", "Permutation"),
    Similarity = perm[[sim_col]],
    R2Y = perm[[r2_col]],
    Q2 = perm[[q2_col]]
  )
  write_csv_utf8(perm_export, paste0(output_prefix, "_permutation_values.csv"))

  plot_df <- bind_rows(
    data.frame(
      Similarity = perm_export$Similarity,
      Value = perm_export$R2Y,
      Metric = "R2Y",
      RunType = perm_export$RunType
    ),
    data.frame(
      Similarity = perm_export$Similarity,
      Value = perm_export$Q2,
      Metric = "Q2",
      RunType = perm_export$RunType
    )
  )
  perm_points <- plot_df[plot_df$RunType == "Permutation", , drop = FALSE]
  original_points <- plot_df[plot_df$RunType == "Original", , drop = FALSE]

  lm_r2 <- lm(perm[[r2_col]] ~ perm[[sim_col]])
  lm_q2 <- lm(perm[[q2_col]] ~ perm[[sim_col]])

  p <- ggplot(plot_df, aes(Similarity, Value, color = Metric)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_abline(
      intercept = coef(lm_r2)[1],
      slope = coef(lm_r2)[2],
      color = "#3C5488",
      linewidth = 0.7,
      linetype = "dashed"
    ) +
    geom_abline(
      intercept = coef(lm_q2)[1],
      slope = coef(lm_q2)[2],
      color = "#E64B35",
      linewidth = 0.7,
      linetype = "dashed"
    ) +
    geom_point(data = perm_points, size = 1.8, alpha = 0.75) +
    geom_point(
      data = original_points,
      aes(fill = Metric, shape = Metric),
      color = "black",
      size = 4,
      stroke = 0.45
    ) +
    scale_color_manual(values = c("R2Y" = "#3C5488", "Q2" = "#E64B35")) +
    scale_fill_manual(values = c("R2Y" = "#3C5488", "Q2" = "#E64B35"), guide = "none") +
    scale_shape_manual(values = c("R2Y" = 21, "Q2" = 22), guide = "none") +
    coord_cartesian(xlim = c(0, 1.02)) +
    labs(
      x = "Y-permutation correlation",
      y = "R2Y(cum) / Q2(cum)",
      title = title_text,
      subtitle = sprintf(
        "Q2 intercept=%.3f; R2Y intercept=%.3f; pQ2=%s",
        validation$Q2_intercept[1],
        validation$R2Y_intercept[1],
        ifelse(is.na(validation$pQ2[1]), "NA", sprintf("%.3f", validation$pQ2[1]))
      ),
      caption = NULL
    ) +
    theme_bw(base_size = 12, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0, size = 8, color = "grey35")
    )

  ggsave(paste0(output_prefix, "_permutation_plot.pdf"), p, width = 6.5, height = 5)
  ggsave(paste0(output_prefix, "_permutation_plot.png"), p, width = 6.5, height = 5, dpi = 300)
  invisible(TRUE)
}

run_oplsda_pipeline <- function(mode, permI = 200) {
  if (!requireNamespace("ropls", quietly = TRUE)) stop("Package 'ropls' is required.")
  library(ropls)

  dat <- load_mode_data(mode)
  cfg <- dat$cfg
  ensure_dir(cfg$opls_dir)

  abundance <- impute_zero_missing(dat$abundance_mat)
  qc_keep <- is.na(dat$feature_info$QC_CV) | dat$feature_info$QC_CV <= analysis_cutoffs$qc_cv
  abundance <- abundance[qc_keep, , drop = FALSE]

  validation_rows <- list()

  for (comp_idx in seq_along(comparisons)) {
    comp <- comparisons[[comp_idx]]
    group1 <- comp[1]
    group2 <- comp[2]
    comp_name <- comparison_name(group1, group2)
    message(sprintf("\n[%s OPLS] %s", mode, comp_name))
    set.seed(20260203 + comp_idx + ifelse(mode == "POS", 1000, 0))

    selected <- dat$sample_info$Group %in% c(group1, group2)
    group_subset <- droplevels(dat$sample_info$Group[selected])
    data_subset <- t(abundance[, dat$sample_info$Sample[selected], drop = FALSE])

    nzv <- apply(data_subset, 2, var, na.rm = TRUE) < 1e-10
    if (any(nzv, na.rm = TRUE)) data_subset <- data_subset[, !nzv, drop = FALSE]
    data_log <- log2(impute_zero_missing(data_subset))

    fit_one <- function(orthoI, model_type) {
      tryCatch(
        opls(
          data_log,
          group_subset,
          predI = 1,
          orthoI = orthoI,
          permI = permI,
          scaleC = "pareto",
          crossvalI = min(7, nrow(data_log)),
          fig.pdfC = "none",
          info.txtC = "none"
        ),
        error = function(e) {
          message(sprintf("  %s failed: %s", model_type, e$message))
          NULL
        }
      )
    }

    model_obj <- fit_one(1, "OPLS-DA")
    model_type <- "OPLS-DA"
    fallback_used <- FALSE

    if (!valid_ropls_object(model_obj)) {
      message("  OPLS-DA did not produce valid scores/loadings; trying PLS-DA diagnostic fallback.")
      model_obj <- fit_one(0, "PLS-DA")
      model_type <- "PLS-DA"
      fallback_used <- TRUE
    }

    validation <- extract_ropls_validation(model_obj, model_type, fallback_used)
    validation$Comparison <- comp_name
    validation$Mode <- mode
    validation$FeaturesAfterQC <- ncol(data_log)
    validation$Samples <- nrow(data_log)
    validation <- validation[, c("Mode", "Comparison", "Samples", "FeaturesAfterQC", setdiff(colnames(validation), c("Mode", "Comparison", "Samples", "FeaturesAfterQC")))]
    validation_rows[[comp_name]] <- validation

    output_prefix <- file.path(cfg$opls_dir, comp_name)
    vip_path <- paste0(output_prefix, "_VIP_scores.csv")
    if (file.exists(vip_path)) file.remove(vip_path)
    save_model_summary(model_obj, validation, paste0(output_prefix, "_model_summary.txt"))

    if (valid_ropls_object(model_obj)) {
      plot_ropls_score(
        model_obj,
        validation,
        group_subset,
        output_prefix,
        sprintf("%s (%s mode): %s vs %s", validation$ModelType[1], cfg$label, group1, group2)
      )
      plot_permutation(
        model_obj,
        validation,
        output_prefix,
        sprintf("Permutation test (%s mode): %s vs %s", cfg$label, group1, group2)
      )
    }

    if (valid_ropls_object(model_obj) && isTRUE(validation$VIP_valid[1])) {
      vip <- getVipVn(model_obj, orthoL = FALSE)
      vip_df <- data.frame(
        ID = names(vip),
        VIP = as.numeric(vip),
        ModelType = validation$ModelType[1],
        VIP_valid = TRUE,
        Q2_cum = validation$Q2_cum[1],
        pQ2 = validation$pQ2[1],
        stringsAsFactors = FALSE
      ) %>% arrange(desc(VIP))
      write_csv_utf8(vip_df, vip_path)
      message(sprintf("  VIP exported: %d features", nrow(vip_df)))
    } else {
      message(sprintf("  VIP not exported: %s", validation$Reason[1]))
    }
  }

  validation_df <- bind_rows(validation_rows)
  write_csv_utf8(validation_df, file.path(cfg$opls_dir, "model_validation_summary.csv"))
  validation_df
}

run_fc_pipeline <- function(mode) {
  dat <- load_mode_data(mode)
  cfg <- dat$cfg
  ensure_dir(cfg$fc_dir)

  abundance <- impute_zero_missing(dat$abundance_mat)
  feature_info <- dat$feature_info

  summary_rows <- list()
  for (comp in comparisons) {
    group_ref <- comp[1]
    group_trg <- comp[2]
    comp_name <- comparison_name(group_ref, group_trg)
    message(sprintf("\n[%s FC] %s", mode, comp_name))

    ref_samples <- dat$sample_info$Sample[dat$sample_info$Group == group_ref]
    trg_samples <- dat$sample_info$Sample[dat$sample_info$Group == group_trg]

    mean_ref <- rowMeans(abundance[, ref_samples, drop = FALSE], na.rm = TRUE)
    mean_trg <- rowMeans(abundance[, trg_samples, drop = FALSE], na.rm = TRUE)
    fc <- mean_trg / mean_ref

    result <- feature_info %>%
      mutate(
        Comparison = comp_name,
        ReferenceGroup = group_ref,
        TargetGroup = group_trg,
        Mean_Ref = mean_ref,
        Mean_Trg = mean_trg,
        FC = fc,
        Log2FC = log2(fc),
        Regulation = fc_direction(fc, analysis_cutoffs$general_fc),
        Regulation_FC1.2 = fc_direction(fc, 1.2)
      )
    colnames(result)[colnames(result) == "Mean_Ref"] <- paste0("Mean_", group_ref)
    colnames(result)[colnames(result) == "Mean_Trg"] <- paste0("Mean_", group_trg)

    write_csv_utf8(result, file.path(cfg$fc_dir, paste0(comp_name, "_FC_results.csv")))
    summary_rows[[comp_name]] <- data.frame(
      Comparison = comp_name,
      Up_FC2 = sum(result$Regulation == "Up", na.rm = TRUE),
      Down_FC2 = sum(result$Regulation == "Down", na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- bind_rows(summary_rows)
  write_csv_utf8(summary_df, file.path(cfg$fc_dir, "Summary_FC_Counts.csv"))
  summary_df
}

run_hqqd_reverse_pipeline <- function() {
  out_dir <- "06_HQQD_reverse_revised"
  ensure_dir(out_dir)

  collect_mode <- function(mode) {
    cfg <- mode_config[[mode]]
    cm_file <- file.path(cfg$merge_dir, paste0(cfg$file_prefix, "_Control_vs_Model_merged.csv"))
    mh_file <- file.path(cfg$merge_dir, paste0(cfg$file_prefix, "_Model_vs_HQQD_merged.csv"))
    if (!file.exists(cm_file) || !file.exists(mh_file)) {
      warning(sprintf("Missing Control_vs_Model or Model_vs_HQQD merged file for %s", mode))
      return(NULL)
    }

    cm <- read_csv_mode(cm_file)
    mh <- read_csv_mode(mh_file)
    normalise_types <- function(d) {
      numeric_cols <- intersect(
        c("QC_CV", "VIP", "Q2_cum", "pQ2", "FC", "Log2FC", "p_value", "p_adj", "t_stat",
          grep("^Mean_", colnames(d), value = TRUE)),
        colnames(d)
      )
      d[numeric_cols] <- lapply(d[numeric_cols], as.numeric)
      bool_cols <- intersect(c("VIP_available", "VIP_valid"), colnames(d))
      d[bool_cols] <- lapply(d[bool_cols], function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y"))
      d
    }
    cm <- normalise_types(cm)
    mh <- normalise_types(mh)

    cm_pref <- cm %>% rename_with(~ paste0("CM_", .x), .cols = -ID)
    mh_pref <- mh %>% rename_with(~ paste0("MH_", .x), .cols = -ID)
    joined <- inner_join(cm_pref, mh_pref, by = "ID")
    joined$Ion_Mode <- cfg$ion_mode
    joined
  }

  joined <- bind_rows(collect_mode("NEG"), collect_mode("POS"))
  if (nrow(joined) == 0) stop("No evaluable HQQD reversal rows.")

  joined <- annotate_metabolite_exclusions(joined)

  mean_control_col <- "CM_Mean_Control"
  mean_model_col <- "CM_Mean_Model"
  mean_hqqd_col <- "MH_Mean_HQQD"

  joined <- joined %>%
    mutate(
      across(c(CM_QC_CV, MH_QC_CV, CM_VIP, MH_VIP, CM_FC, MH_FC, CM_Log2FC, MH_Log2FC, CM_p_value, MH_p_value, CM_p_adj, MH_p_adj), ~ as.numeric(.x)),
      QC_CV = coalesce(CM_QC_CV, MH_QC_CV),
      Pass_QC = is.na(QC_CV) | QC_CV <= analysis_cutoffs$qc_cv,
      CM_VIP_required = CM_VIP_available & CM_VIP_valid,
      CM_Pass_VIP = ifelse(CM_VIP_required, CM_VIP > analysis_cutoffs$vip, TRUE),
      CM_Pass = Pass_QC &
        abs(CM_Log2FC) >= log2(analysis_cutoffs$model_fc) &
        CM_p_value < analysis_cutoffs$p &
        CM_Pass_VIP,
      MH_Pass = Pass_QC &
        abs(MH_Log2FC) >= log2(analysis_cutoffs$treatment_fc) &
        MH_p_value < analysis_cutoffs$p,
      CM_Pass_padj = Pass_QC &
        abs(CM_Log2FC) >= log2(analysis_cutoffs$model_fc) &
        CM_p_adj < analysis_cutoffs$p &
        CM_Pass_VIP,
      MH_Pass_padj = Pass_QC &
        abs(MH_Log2FC) >= log2(analysis_cutoffs$treatment_fc) &
        MH_p_adj < analysis_cutoffs$p,
      Opposite_Direction = CM_Log2FC * MH_Log2FC < 0,
      Control_log2_mean = log2(pmax(.data[[mean_control_col]], 1e-12)),
      Model_log2_mean = log2(pmax(.data[[mean_model_col]], 1e-12)),
      HQQD_log2_mean = log2(pmax(.data[[mean_hqqd_col]], 1e-12)),
      RecoveryScore = 1 - abs(HQQD_log2_mean - Control_log2_mean) / abs(Model_log2_mean - Control_log2_mean),
      RecoveryScore = ifelse(is.finite(RecoveryScore), RecoveryScore, NA_real_),
      Recovered_Toward_Control = !is.na(RecoveryScore) & RecoveryScore > 0,
      Reverse_Type = case_when(
        CM_Log2FC > 0 & MH_Log2FC < 0 ~ "Model_Up_HQQD_Down",
        CM_Log2FC < 0 & MH_Log2FC > 0 ~ "Model_Down_HQQD_Up",
        TRUE ~ "Not_Opposite"
      ),
      ReverseCandidate_Raw = CM_Pass & MH_Pass & Opposite_Direction & Recovered_Toward_Control,
      ReverseCandidate_padj_Raw = CM_Pass_padj & MH_Pass_padj & Opposite_Direction & Recovered_Toward_Control,
      ReverseCandidate = ReverseCandidate_Raw & !Excluded_Annotation,
      ReverseCandidate_padj = ReverseCandidate_padj_Raw & !Excluded_Annotation,
      ReverseConfidence = -log10(pmax(CM_p_value, 1e-300)) +
        -log10(pmax(MH_p_value, 1e-300)) +
        pmax(RecoveryScore, 0),
      ReverseConfidence_padj = -log10(pmax(CM_p_adj, 1e-300)) +
        -log10(pmax(MH_p_adj, 1e-300)) +
        pmax(RecoveryScore, 0)
    ) %>%
    arrange(desc(ReverseCandidate), desc(ReverseCandidate_padj), desc(RecoveryScore), CM_p_value, MH_p_value)

  candidates <- joined %>% filter(ReverseCandidate)
  candidates_padj <- joined %>% filter(ReverseCandidate_padj)
  excluded_candidates <- joined %>% filter(ReverseCandidate_Raw, Excluded_Annotation)

  write_csv_utf8(joined, file.path(out_dir, "HQQD_reverse_all_evaluable.csv"))
  write_csv_utf8(candidates, file.path(out_dir, "HQQD_reverse_candidates.csv"))
  write_csv_utf8(candidates_padj, file.path(out_dir, "HQQD_reverse_candidates_padj.csv"))
  write_csv_utf8(excluded_candidates, file.path(out_dir, "HQQD_reverse_excluded_annotations.csv"))

  summary_df <- joined %>%
    summarise(
      Evaluable = n(),
      Pass_CM = sum(CM_Pass, na.rm = TRUE),
      Pass_MH = sum(MH_Pass, na.rm = TRUE),
      Opposite = sum(Opposite_Direction, na.rm = TRUE),
      Recovered = sum(Recovered_Toward_Control, na.rm = TRUE),
      ExcludedAnnotations = sum(Excluded_Annotation, na.rm = TRUE),
      ReverseCandidatesBeforeExclusion = sum(ReverseCandidate_Raw, na.rm = TRUE),
      ReverseCandidates = sum(ReverseCandidate, na.rm = TRUE),
      Pass_CM_padj = sum(CM_Pass_padj, na.rm = TRUE),
      Pass_MH_padj = sum(MH_Pass_padj, na.rm = TRUE),
      ReverseCandidates_padj = sum(ReverseCandidate_padj, na.rm = TRUE)
    )
  write_csv_utf8(summary_df, file.path(out_dir, "HQQD_reverse_summary.csv"))

  summary_df
}

run_welch_pipeline <- function(mode) {
  dat <- load_mode_data(mode)
  cfg <- dat$cfg
  ensure_dir(cfg$p_dir)

  abundance <- impute_zero_missing(dat$abundance_mat)
  log_data <- log2(abundance)
  feature_info <- dat$feature_info

  summary_rows <- list()
  for (comp in comparisons) {
    group1 <- comp[1]
    group2 <- comp[2]
    comp_name <- comparison_name(group1, group2)
    message(sprintf("\n[%s Welch] %s", mode, comp_name))

    samps1 <- dat$sample_info$Sample[dat$sample_info$Group == group1]
    samps2 <- dat$sample_info$Sample[dat$sample_info$Group == group2]

    p_values <- numeric(nrow(log_data))
    t_stats <- numeric(nrow(log_data))
    for (i in seq_len(nrow(log_data))) {
      vec1 <- as.numeric(log_data[i, samps1])
      vec2 <- as.numeric(log_data[i, samps2])
      if (isTRUE(var(vec1) == 0 && var(vec2) == 0)) {
        p_values[i] <- 1
        t_stats[i] <- 0
      } else {
        res <- tryCatch(t.test(vec1, vec2, var.equal = FALSE, paired = FALSE), error = function(e) NULL)
        if (is.null(res)) {
          p_values[i] <- 1
          t_stats[i] <- 0
        } else {
          p_values[i] <- res$p.value
          t_stats[i] <- as.numeric(res$statistic)
        }
      }
    }

    result <- feature_info %>%
      transmute(
        ID,
        Name,
        QC_CV,
        Comparison = comp_name,
        test_method = "Welch_t_test_on_log2",
        t_stat = t_stats,
        p_value = p_values,
        p_adj = p.adjust(p_values, method = "BH")
      ) %>%
      arrange(p_value)

    write_csv_utf8(result, file.path(cfg$p_dir, paste0(comp_name, "_ttest_results.csv")))
    summary_rows[[comp_name]] <- data.frame(
      Comparison = comp_name,
      p_lt_0.05 = sum(result$p_value < analysis_cutoffs$p, na.rm = TRUE),
      padj_lt_0.05 = sum(result$p_adj < analysis_cutoffs$p, na.rm = TRUE),
      padj_lt_0.10 = sum(result$p_adj < analysis_cutoffs$fdr_exploratory, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- bind_rows(summary_rows)
  write_csv_utf8(summary_df, file.path(cfg$p_dir, "Summary_Welch_Counts.csv"))
  summary_df
}

run_merge_pipeline <- function(mode) {
  cfg <- mode_config[[mode]]
  ensure_dir(cfg$merge_dir)

  fc_files <- list.files(cfg$fc_dir, pattern = "_FC_results\\.csv$", full.names = FALSE)
  if (length(fc_files) == 0) stop(sprintf("No FC files found in %s", cfg$fc_dir))

  summary_rows <- list()
  for (f_fc in fc_files) {
    comp_name <- sub("_FC_results\\.csv$", "", f_fc)
    message(sprintf("\n[%s merge] %s", mode, comp_name))

    path_fc <- file.path(cfg$fc_dir, f_fc)
    path_p <- file.path(cfg$p_dir, paste0(comp_name, "_ttest_results.csv"))
    path_vip <- file.path(cfg$opls_dir, paste0(comp_name, "_VIP_scores.csv"))

    if (!file.exists(path_p)) stop(sprintf("Missing Welch result: %s", path_p))

    df_fc <- read_csv_mode(path_fc)
    df_p <- read_csv_mode(path_p)
    df <- left_join(df_fc, df_p %>% select(ID, t_stat, p_value, p_adj, test_method), by = "ID")

    if (file.exists(path_vip)) {
      df_vip <- read_csv_mode(path_vip) %>% select(ID, VIP, ModelType, VIP_valid, Q2_cum, pQ2)
      df <- left_join(df, df_vip, by = "ID")
      df$VIP_available <- !is.na(df$VIP)
      message(sprintf("  Valid VIP file found: %s", basename(path_vip)))
    } else {
      df$VIP <- NA_real_
      df$ModelType <- NA_character_
      df$VIP_valid <- FALSE
      df$Q2_cum <- NA_real_
      df$pQ2 <- NA_real_
      df$VIP_available <- FALSE
      message("  No valid VIP file; keeping FC + Welch results with VIP_available=FALSE.")
    }

    priority <- c(
      "ID", "Name", "Comparison", "QC_CV", "VIP", "VIP_available", "VIP_valid",
      "ModelType", "Q2_cum", "pQ2", "FC", "Log2FC", "p_value", "p_adj",
      "Regulation", "Regulation_FC1.2", "test_method"
    )
    df <- df[, c(intersect(priority, colnames(df)), setdiff(colnames(df), priority))]

    out_file <- file.path(cfg$merge_dir, paste0(cfg$file_prefix, "_", comp_name, "_merged.csv"))
    write_csv_utf8(df, out_file)

    summary_rows[[comp_name]] <- data.frame(
      Comparison = comp_name,
      Rows = nrow(df),
      VIP_available_rows = sum(df$VIP_available, na.rm = TRUE),
      p_lt_0.05 = sum(df$p_value < analysis_cutoffs$p, na.rm = TRUE),
      padj_lt_0.05 = sum(df$p_adj < analysis_cutoffs$p, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- bind_rows(summary_rows)
  write_csv_utf8(summary_df, file.path(cfg$merge_dir, "Summary_Merge.csv"))
  summary_df
}

run_diff_pipeline <- function(mode) {
  cfg <- mode_config[[mode]]
  ensure_dir(cfg$diff_dir)

  files <- list.files(cfg$merge_dir, pattern = "_merged\\.csv$", full.names = FALSE)
  if (length(files) == 0) stop(sprintf("No merged files found in %s", cfg$merge_dir))

  summary_rows <- list()
  for (f in files) {
    comp_name <- str_remove(f, paste0("^", cfg$file_prefix, "_"))
    comp_name <- str_remove(comp_name, "_merged\\.csv$")
    message(sprintf("\n[%s diff] %s", mode, comp_name))

    df <- read_csv_mode(file.path(cfg$merge_dir, f)) %>%
      mutate(
        across(c(QC_CV, VIP, FC, Log2FC, p_value, p_adj), ~ as.numeric(.x)),
        Pass_QC = is.na(QC_CV) | QC_CV <= analysis_cutoffs$qc_cv,
        Pass_FC = FC > analysis_cutoffs$general_fc | FC < 1 / analysis_cutoffs$general_fc,
        Pass_P = p_value < analysis_cutoffs$p,
        Pass_FDR10 = p_adj < analysis_cutoffs$fdr_exploratory,
        Pass_VIP = ifelse(VIP_available & VIP_valid, VIP > analysis_cutoffs$vip, TRUE),
        VIP_rule = ifelse(VIP_available & VIP_valid, "VIP>1 required", "VIP not required: no validated multivariate model"),
        Sig_Status = case_when(
          Pass_QC & Pass_FC & Pass_P & Pass_VIP & FC > analysis_cutoffs$general_fc ~ "Up",
          Pass_QC & Pass_FC & Pass_P & Pass_VIP & FC < 1 / analysis_cutoffs$general_fc ~ "Down",
          TRUE ~ "NoSig"
        )
      )

    df_sig <- df %>% filter(Sig_Status != "NoSig")
    out_csv <- file.path(cfg$diff_dir, paste0(cfg$file_prefix, "_", comp_name, "_diff.csv"))
    write_csv_utf8(df_sig, out_csv)

    plot_df <- df %>%
      mutate(
        p_plot = pmax(p_value, min(p_value[p_value > 0], na.rm = TRUE) * 0.1),
        neg_log10_p = -log10(p_plot),
        VIP_plot = ifelse(is.na(VIP), 1, pmin(VIP, quantile(VIP, 0.98, na.rm = TRUE))),
        Label = NA_character_
      )
    vip_for_plot <- any(plot_df$VIP_available & plot_df$VIP_valid, na.rm = TRUE)
    vip_note <- ifelse(vip_for_plot, "VIP size shown", "VIP not used")
    volcano_caption <- ifelse(
      vip_for_plot,
      "Point size reflects validated VIP. Differential status requires QC, FC, p-value and VIP>1.",
      "Point size is fixed. VIP is not shown because this comparison has no validated OPLS-DA/VIP."
    )
    if (sum(plot_df$Sig_Status != "NoSig", na.rm = TRUE) > 0) {
      top_ids <- plot_df %>%
        filter(Sig_Status != "NoSig") %>%
        arrange(p_value) %>%
        head(10) %>%
        pull(ID)
      plot_df$Label[plot_df$ID %in% top_ids] <- plot_df$Name[plot_df$ID %in% top_ids]
    }

    p <- ggplot(plot_df, aes(Log2FC, neg_log10_p)) +
      geom_vline(xintercept = c(log2(analysis_cutoffs$general_fc), -log2(analysis_cutoffs$general_fc)), linetype = "dashed", color = "gray60") +
      geom_hline(yintercept = -log10(analysis_cutoffs$p), linetype = "dashed", color = "gray60")

    if (vip_for_plot) {
      p <- p +
        geom_point(aes(color = Sig_Status, size = VIP_plot), alpha = 0.75) +
        scale_size_continuous(range = c(1.2, 4), name = "VIP")
    } else {
      p <- p +
        geom_point(aes(color = Sig_Status), size = 1.8, alpha = 0.75)
    }

    p <- p +
      geom_text_repel(
        data = function(x) x[!is.na(x$Label), , drop = FALSE],
        aes(label = Label),
        size = 3,
        max.overlaps = 15,
        show.legend = FALSE
      ) +
      scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#3C5488", "NoSig" = "gray80")) +
      labs(
        x = expression(log[2]~"(fold change)"),
        y = expression(-log[10]~"(p-value)"),
        title = sprintf("Volcano plot (%s mode): %s", cfg$label, comp_name),
        subtitle = sprintf(
          "Up: %d | Down: %d | %s",
          sum(plot_df$Sig_Status == "Up"),
          sum(plot_df$Sig_Status == "Down"),
          vip_note
        ),
        color = "Status",
        caption = volcano_caption
      ) +
      theme_bw(base_size = 12, base_family = "sans") +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0, size = 8, color = "grey35")
      )

    ggsave(file.path(cfg$diff_dir, paste0(cfg$file_prefix, "_", comp_name, "_volcano.pdf")), p, width = 7, height = 6)
    ggsave(file.path(cfg$diff_dir, paste0(cfg$file_prefix, "_", comp_name, "_volcano.png")), p, width = 7, height = 6, dpi = 300)

    summary_rows[[comp_name]] <- data.frame(
      Comparison = comp_name,
      Rows = nrow(df_sig),
      Up = sum(df_sig$Sig_Status == "Up", na.rm = TRUE),
      Down = sum(df_sig$Sig_Status == "Down", na.rm = TRUE),
      VIP_required = any(df$VIP_available & df$VIP_valid, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- bind_rows(summary_rows)
  write_csv_utf8(summary_df, file.path(cfg$diff_dir, "Summary_Diff_Counts.csv"))
  summary_df
}

run_merge_diff_pipeline <- function() {
  out_dir <- "06_merge_diff_revised"
  ensure_dir(out_dir)
  excluded_label_ids <- load_metabolite_exclusion_manifest()$ID

  as_bool <- function(x) {
    if (is.logical(x)) return(x)
    tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
  }

  build_full_combined_volcano <- function(comp_name) {
    rows <- list()
    for (mode in c("NEG", "POS")) {
      cfg <- mode_config[[mode]]
      f <- file.path(cfg$merge_dir, paste0(cfg$file_prefix, "_", comp_name, "_merged.csv"))
      if (file.exists(f)) {
        d <- read_csv_mode(f)
        numeric_cols <- intersect(
          c("QC_CV", "VIP", "Q2_cum", "pQ2", "FC", "Log2FC", "p_value", "p_adj", "t_stat",
            grep("^Mean_", colnames(d), value = TRUE)),
          colnames(d)
        )
        d[numeric_cols] <- lapply(d[numeric_cols], as.numeric)
        bool_cols <- intersect(c("VIP_available", "VIP_valid"), colnames(d))
        d[bool_cols] <- lapply(d[bool_cols], as_bool)
        d$Ion_Mode <- cfg$ion_mode
        rows[[mode]] <- d
      }
    }
    if (length(rows) == 0) return(NULL)

    full <- bind_rows(rows)
    numeric_cols <- intersect(c("QC_CV", "VIP", "FC", "Log2FC", "p_value", "p_adj"), colnames(full))
    full <- full %>%
      mutate(
        across(all_of(numeric_cols), ~ as.numeric(.x)),
        VIP_available = if ("VIP_available" %in% colnames(.)) as_bool(VIP_available) else FALSE,
        VIP_valid = if ("VIP_valid" %in% colnames(.)) as_bool(VIP_valid) else FALSE,
        Pass_QC = is.na(QC_CV) | QC_CV <= analysis_cutoffs$qc_cv,
        Pass_FC = FC > analysis_cutoffs$general_fc | FC < 1 / analysis_cutoffs$general_fc,
        Pass_P = p_value < analysis_cutoffs$p,
        Pass_VIP = ifelse(VIP_available & VIP_valid, VIP > analysis_cutoffs$vip, TRUE),
        Sig_Status = case_when(
          Pass_QC & Pass_FC & Pass_P & Pass_VIP & FC > analysis_cutoffs$general_fc ~ "Up",
          Pass_QC & Pass_FC & Pass_P & Pass_VIP & FC < 1 / analysis_cutoffs$general_fc ~ "Down",
          TRUE ~ "NoSig"
        )
      )

    min_pos_p <- suppressWarnings(min(full$p_value[full$p_value > 0], na.rm = TRUE))
    if (!is.finite(min_pos_p)) min_pos_p <- .Machine$double.xmin
    valid_vip_vals <- full$VIP[full$VIP_available & full$VIP_valid & !is.na(full$VIP)]
    vip_cap <- if (length(valid_vip_vals) > 0) {
      as.numeric(quantile(valid_vip_vals, 0.98, na.rm = TRUE))
    } else {
      1
    }
    plot_df <- full %>%
      mutate(
        p_plot = pmax(p_value, min_pos_p * 0.1),
        neg_log10_p = -log10(p_plot),
        VIP_for_plot = ifelse(VIP_available & VIP_valid & !is.na(VIP), pmin(VIP, vip_cap), 1),
        Label = NA_character_
      )

    label_ids <- plot_df %>%
      filter(Sig_Status != "NoSig", !ID %in% excluded_label_ids) %>%
      arrange(p_value) %>%
      slice_head(n = 10) %>%
      pull(ID)
    plot_df$Label[plot_df$ID %in% label_ids] <- plot_df$Name[plot_df$ID %in% label_ids]

    vip_rules <- plot_df %>%
      group_by(Ion_Mode) %>%
      summarise(
        VIP_required = any(VIP_available & VIP_valid, na.rm = TRUE),
        VIP_rule = ifelse(VIP_required, "VIP used in filtering", "VIP not used"),
        .groups = "drop"
      )
    vip_caption <- paste(paste(vip_rules$Ion_Mode, vip_rules$VIP_rule, sep = ": "), collapse = "; ")
    vip_mode <- case_when(
      all(vip_rules$VIP_required) ~ "all",
      any(vip_rules$VIP_required) ~ "partial",
      TRUE ~ "none"
    )
    vip_subtitle <- case_when(
      vip_mode == "all" ~ "VIP size shown",
      vip_mode == "partial" ~ "VIP size shown only for ion mode(s) with validated OPLS-DA/VIP",
      TRUE ~ "VIP not used"
    )

    p <- ggplot(plot_df, aes(Log2FC, neg_log10_p)) +
      geom_vline(xintercept = c(log2(analysis_cutoffs$general_fc), -log2(analysis_cutoffs$general_fc)), linetype = "dashed", color = "gray60") +
      geom_hline(yintercept = -log10(analysis_cutoffs$p), linetype = "dashed", color = "gray60")

    if (vip_mode == "none") {
      p <- p +
        geom_point(aes(color = Sig_Status, shape = Ion_Mode), size = 2, alpha = 0.72)
    } else {
      p <- p +
        geom_point(aes(color = Sig_Status, size = VIP_for_plot, shape = Ion_Mode), alpha = 0.72) +
        scale_size_continuous(range = c(1, 4), name = ifelse(vip_mode == "all", "VIP", "Validated VIP"))
    }

    p <- p +
      geom_text_repel(
        data = function(x) x[!is.na(x$Label), , drop = FALSE],
        aes(label = Label),
        size = 3,
        box.padding = 0.5,
        point.padding = 0.3,
        max.overlaps = 15,
        show.legend = FALSE
      ) +
      scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#3C5488", "NoSig" = "gray80")) +
      scale_shape_manual(values = c("NEG" = 16, "POS" = 15)) +
      labs(
        x = expression(log[2]~"(Fold Change)"),
        y = expression(-log[10]~italic("(p-value)")),
        title = sprintf("Volcano Plot: %s", comp_name),
        subtitle = sprintf(
          "Up: %d | Down: %d | %s",
          sum(plot_df$Sig_Status == "Up"),
          sum(plot_df$Sig_Status == "Down"),
          vip_subtitle
        ),
        color = "Regulation",
        shape = "Ion Mode",
        caption = NULL
      ) +
      guides(
        color = guide_legend(override.aes = list(size = 4)),
        shape = guide_legend(override.aes = list(size = 4))
      ) +
      theme_bw(base_size = 12, base_family = "sans") +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0, size = 8, color = "grey35")
      )

    ggsave(file.path(out_dir, paste0("06_merged_", comp_name, "_volcano.pdf")), p, width = 8, height = 6)
    ggsave(file.path(out_dir, paste0("06_merged_", comp_name, "_volcano.png")), p, width = 8, height = 6, dpi = 300)
    save_plot_svg(file.path(out_dir, paste0("06_merged_", comp_name, "_volcano.svg")), p, width = 8, height = 6)

    list(
      full = full,
      vip_rules = vip_rules
    )
  }

  all_groups <- unique(c(
    str_match(list.files(mode_config$NEG$diff_dir, pattern = "_diff\\.csv$"), paste0(mode_config$NEG$file_prefix, "_(.*)_diff\\.csv"))[, 2],
    str_match(list.files(mode_config$POS$diff_dir, pattern = "_diff\\.csv$"), paste0(mode_config$POS$file_prefix, "_(.*)_diff\\.csv"))[, 2]
  ))
  all_groups <- all_groups[!is.na(all_groups)]

  summary_rows <- list()
  for (comp_name in all_groups) {
    message(sprintf("\n[merge diff] %s", comp_name))
    rows <- list()
    for (mode in c("NEG", "POS")) {
      cfg <- mode_config[[mode]]
      f <- file.path(cfg$diff_dir, paste0(cfg$file_prefix, "_", comp_name, "_diff.csv"))
      if (file.exists(f)) {
        d <- read_csv_mode(f)
        numeric_cols <- intersect(
          c("QC_CV", "VIP", "Q2_cum", "pQ2", "FC", "Log2FC", "p_value", "p_adj", "t_stat",
            grep("^Mean_", colnames(d), value = TRUE)),
          colnames(d)
        )
        d[numeric_cols] <- lapply(d[numeric_cols], as.numeric)
        bool_cols <- intersect(c("VIP_available", "VIP_valid", "Pass_QC", "Pass_FC", "Pass_P", "Pass_FDR10", "Pass_VIP"), colnames(d))
        d[bool_cols] <- lapply(d[bool_cols], as_bool)
        d$Ion_Mode <- cfg$ion_mode
        rows[[mode]] <- d
      }
    }
    if (length(rows) == 0) next
    merged <- bind_rows(rows) %>%
      select(ID, Name, Ion_Mode, QC_CV, everything())
    write_csv_utf8(merged, file.path(out_dir, paste0("06_merged_", comp_name, "_diff.csv")))

    plot_info <- build_full_combined_volcano(comp_name)
    if (!is.null(plot_info)) {
      vip_rule_wide <- plot_info$vip_rules %>%
        select(Ion_Mode, VIP_required)
      vip_neg <- vip_rule_wide$VIP_required[vip_rule_wide$Ion_Mode == "NEG"]
      vip_pos <- vip_rule_wide$VIP_required[vip_rule_wide$Ion_Mode == "POS"]
      if (length(vip_neg) == 0) vip_neg <- NA
      if (length(vip_pos) == 0) vip_pos <- NA
    } else {
      vip_neg <- NA
      vip_pos <- NA
    }

    summary_rows[[comp_name]] <- data.frame(
      Comparison = comp_name,
      Rows = nrow(merged),
      NEG_rows = sum(merged$Ion_Mode == "NEG", na.rm = TRUE),
      POS_rows = sum(merged$Ion_Mode == "POS", na.rm = TRUE),
      NEG_VIP_required = vip_neg[1],
      POS_VIP_required = vip_pos[1],
      stringsAsFactors = FALSE
    )
  }

  summary_df <- bind_rows(summary_rows)
  write_csv_utf8(summary_df, file.path(out_dir, "Summary_Merged_Diff.csv"))
  summary_df
}
