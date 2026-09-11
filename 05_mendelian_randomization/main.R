# ==============================================================================
# ==============================================================================
# ==============================================================================

rm(list = ls())
gc()

required_packages <- c(
  "data.table", "dplyr", "TwoSampleMR", "MRPRESSO",
  "MendelianRandomization", "ieugwasr", "LDlinkR",
  "pbapply", "parallel", "doParallel", "metafor"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

source("config.R")
source("data_sources_config.R")
log_info("Configuration loaded")
log_info(paste0("Active exposure sources: ", paste(ACTIVE_EXPOSURES, collapse = ", ")))
log_info(paste0("Active outcome: ", ACTIVE_OUTCOME))

source("functions/utils.R")
source("functions/clumping.R")
source("functions/proxy_snp.R")
source("functions/mr_analysis.R")
source("functions/plotting.R")

apply_exposure_filter <- function(ivs_data, source_name) {
  source_config <- get_exposure_source_config(source_name)
  filter_config <- EXPOSURE_FILTER_CONFIG[[source_name]]

  if (is.null(filter_config)) {
    return(ivs_data)
  }

  if (!is.null(filter_config$id_map)) {
    ivs_data <- rename_exposures(
      ivs_data,
      names(filter_config$id_map),
      filter_config$id_map
    )
  }

  if (identical(filter_config$filter_mode, "selected")) {
    selected_items <- filter_config$selected_items
    log_info(paste0("  Filter mode: selected exposures (", length(selected_items), ")"))
    ivs_data <- ivs_data %>%
      dplyr::filter(exposure %in% selected_items)
    log_info(paste0(
      "  Retained ", length(unique(ivs_data$exposure)),
      " exposures and ", nrow(ivs_data), " SNPs after selected filtering"
    ))
  } else {
    log_info("  Filter mode: all exposures")
  }

  ivs_data
}
log_info("Function modules loaded")

if (START_FROM_STEP <= 1) {
  log_info("Step1: Start data preprocessing and LD clumping")

all_ivs_list <- list()

for (source_name in ACTIVE_EXPOSURES) {
  source_config <- get_exposure_source_config(source_name)
  if (!is.null(source_config$preprocess_script)) {
    preprocess_script <- source_config$preprocess_script
    if (!file.exists(preprocess_script)) {
      stop(paste0("Preprocess script not found for ", source_name, ": ", preprocess_script))
    }
    log_info(paste0("  Running preprocess script: ", preprocess_script))
    source(preprocess_script, local = new.env(parent = globalenv()))
  }

  log_info(paste0("\nProcessing data source: ", source_name))
  
  source_config <- get_exposure_source_config(source_name)
  filter_config <- EXPOSURE_FILTER_CONFIG[[source_name]]
  
  ivs_data <- tryCatch({
    process_exposure_data(
      file_path = source_config$file_path,
      data_source = source_config$name,
      format_params = source_config$format_params,
      plink_path = PLINK_PATH,
      bfile_path = BFILE_PATH,
      maf_path = if(source_config$need_eaf) MAF_FILE_PATH else NULL,
      default_samplesize = if (!is.null(source_config$default_samplesize)) source_config$default_samplesize else NULL,
      pval_threshold = PVAL_THRESHOLD,
      f_threshold = F_THRESHOLD,
      clump_params = CLUMP_PARAMS,
      need_eaf = source_config$need_eaf,
      selected_exposures = if (!is.null(filter_config) && identical(filter_config$filter_mode, "selected")) filter_config$selected_items else NULL,
      id_map = if (!is.null(filter_config) && !is.null(filter_config$id_map)) filter_config$id_map else NULL,
      cis_filter_params = if (!is.null(source_config$cis_filter_params)) source_config$cis_filter_params else NULL
    )
  }, error = function(e) {
    log_error(paste0("Processing data source ", source_name, "  failed: ", e$message))
    return(NULL)
  })
  
  if (is.null(ivs_data) || nrow(ivs_data) == 0) {
    log_warning(paste0("Data source ", source_name, " produced no valid IVs; skipped"))
    next
  }
  
  # Apply the same source-level selection again as a harmless safeguard.
  ivs_data <- apply_exposure_filter(ivs_data, source_name)
  
  required_cols <- c(
    "id.exposure", "SNP", "effect_allele.exposure", "other_allele.exposure",
    "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure",
    "chr.exposure", "pos.exposure", "samplesize.exposure", "exposure",
    "mr_keep.exposure", "pval_origin.exposure", "R2", "F", "MAF"
  )
  missing_cols <- setdiff(required_cols, names(ivs_data))
  if (length(missing_cols) > 0) {
    for (col_name in missing_cols) {
      ivs_data[[col_name]] <- NA
    }
    log_warning(paste0("  Data source ", source_name, " is missing columns; filled with NA: ", paste(missing_cols, collapse = ", ")))
  }

  ivs_data <- ivs_data %>%
    dplyr::select(id.exposure, SNP, effect_allele.exposure, other_allele.exposure,
                  eaf.exposure, beta.exposure, se.exposure, pval.exposure,
                  chr.exposure, pos.exposure, samplesize.exposure, exposure,
                  mr_keep.exposure, pval_origin.exposure, R2, F, MAF)
  
  ivs_data <- ivs_data %>%
    dplyr::group_by(id.exposure, SNP) %>%
    dplyr::arrange(pval.exposure) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  output_file <- file.path(IVS_DIR, paste0("ivs_", source_name, ".csv"))
  safe_write_csv(ivs_data, output_file)
  
  all_ivs_list[[source_name]] <- ivs_data
}

} else {
  log_info("Skip step 1: load source-specific IVs from existing files")
}

if (START_FROM_STEP <= 2) {
  log_info("Step2: Combine all instrumental variables")
  
  if (START_FROM_STEP == 2) {
    all_ivs_list <- list()
    for (source_name in ACTIVE_EXPOSURES) {
      ivs_file <- file.path(IVS_DIR, paste0("ivs_", source_name, ".csv"))
      if (file.exists(ivs_file)) {
        ivs_data <- data.table::fread(ivs_file)
        ivs_data <- apply_exposure_filter(ivs_data, source_name)
        all_ivs_list[[source_name]] <- ivs_data
        log_info(paste0("  Loaded ", source_name, ": ", nrow(ivs_data), " SNPs"))
      } else {
        log_warning(paste0("File not found: ", ivs_file))
      }
    }
  }

  if (length(all_ivs_list) == 0) {
    stop("No valid IVs were obtained. Check the step-1 data paths and the p-value and F-statistic thresholds.")
  }

ivs_all <- data.table::rbindlist(all_ivs_list, fill = TRUE)

if (is.null(ivs_all) || nrow(ivs_all) == 0 || ncol(ivs_all) == 0) {
  stop("The merged IV table is empty. Check data import and the configured selection thresholds.")
}

safe_write_csv(ivs_all, OUTPUT_FILES$ivs_all)

log_info(paste0("Selected  ", nrow(ivs_all), "  instrumental variables"))
summarize_snp_counts(ivs_all)
summarize_f_statistics(ivs_all)

} else {
  log_info("Skip step 2: load merged IVs from the existing file")
  ivs_all <- data.table::fread(OUTPUT_FILES$ivs_all)
  log_info(paste0("  Loaded ivs_all.csv: ", nrow(ivs_all), " SNPs"))
}

if (START_FROM_STEP <= 3) {
  log_info("Step3: Prepare outcome data and identify proxy SNPs")

outcome_config <- get_outcome_source_config(ACTIVE_OUTCOME)
log_info(paste0("Outcome: ", outcome_config$name))

if (!is.null(outcome_config$file_path)) {
  outcome_file_path <- outcome_config$file_path
} else {
  outcome_file_path <- file.path(DATA_DIR, "diseases", outcome_config$file)
}
outcome_target_snps <- unique(ivs_all$SNP)
proxy_candidate_snps <- find_proxy_candidate_snps_plink(
  snps = outcome_target_snps,
  plink_path = PLINK_PATH,
  bfile_path = BFILE_PATH,
  r2_threshold = PROXY_R2_THRESHOLD,
  window_kb = 1000
)
outcome_lookup_snps <- unique(c(outcome_target_snps, proxy_candidate_snps))
oc_data <- read_outcome_subset_for_snps(
  file_path = outcome_file_path,
  target_snps = outcome_lookup_snps,
  snp_col = outcome_config$format_params$snp_col
)
if (nrow(oc_data) == 0) {
  stop("No outcome rows were extracted for IV/proxy SNPs")
}

if (!is.null(outcome_config$cases)) {
  oc_data[, cases := outcome_config$cases]
}
if (!is.null(outcome_config$controls)) {
  oc_data[, controls := outcome_config$controls]
}
if (!is.null(outcome_config$trait)) {
  oc_data[, trait := outcome_config$trait]
}
if (!is.null(outcome_config$id)) {
  oc_data[, id := outcome_config$id]
}

oc_data <- as.data.frame(oc_data)

oc_outcome <- do.call(TwoSampleMR::format_data, 
                      c(list(dat = oc_data, type = "outcome"), 
                        outcome_config$format_params))

proxy_result <- process_proxy_snps(
  exposure_data = ivs_all,
  outcome_data = oc_outcome,
  plink_path = PLINK_PATH,
  bfile_path = BFILE_PATH,
  r2_threshold = PROXY_R2_THRESHOLD,
  save_proxy_file = TRUE,
  output_file = OUTPUT_FILES$proxy_file
)

ivs_all_aftpxy <- proxy_result$exposure_with_proxy
safe_write_csv(ivs_all_aftpxy, OUTPUT_FILES$ivs_all_aftpxy)
safe_write_csv(proxy_result$outcome_filtered, OUTPUT_FILES$oc_outcome)

log_info(paste0("Proxy-SNP processing complete; final IVs: ", nrow(ivs_all_aftpxy)))

} else {
  log_info("Skip step 3: load proxy-processed data from existing files")
  ivs_all_aftpxy <- data.table::fread(OUTPUT_FILES$ivs_all_aftpxy)
  oc_outcome <- data.table::fread(OUTPUT_FILES$oc_outcome)
  log_info(paste0("  Loaded ivs_all_aftpxy.csv: ", nrow(ivs_all_aftpxy), " SNPs"))
  log_info(paste0("  Loaded oc_outcome.csv: ", nrow(oc_outcome), " SNPs"))
  proxy_result <- list(
    exposure_with_proxy = ivs_all_aftpxy,
    outcome_filtered = oc_outcome
  )
}

if (START_FROM_STEP <= 4) {
  log_info("Step4: Harmonise exposure and outcome data")

harmonise_data <- TwoSampleMR::harmonise_data(
  exposure_dat = ivs_all_aftpxy,
  outcome_dat = proxy_result$outcome_filtered,
  action = HARMONISE_STRICTNESS
)

harmonise_data_filtered <- harmonise_data %>%
  dplyr::filter(mr_keep == TRUE) %>%
  dplyr::group_by(id.exposure, id.outcome) %>%
  dplyr::filter(dplyr::n() >= 3) %>%
  dplyr::ungroup()

safe_write_csv(harmonise_data_filtered, OUTPUT_FILES$harmonise_main)
log_info(paste0("Harmonisation complete:  ", nrow(harmonise_data_filtered), "  SNPs"))

} else {
  log_info("Skip step 4: load harmonised data from the existing file")
  harmonise_data_filtered <- data.table::fread(OUTPUT_FILES$harmonise_main)
  log_info(paste0("  Loaded harmonised data: ", nrow(harmonise_data_filtered), " SNPs"))
}

if (START_FROM_STEP <= 5) {
  log_info("Step5: Run MR analysis")

mr_results <- parallel_mr_analysis(
  harmonise_data = harmonise_data_filtered,
  n_cores = N_CORES,
  methods = PRIMARY_MR_METHODS
)
safe_write_csv(mr_results, OUTPUT_FILES$results_main)

mr_significant <- select_mr_pairs(
  mr_results = mr_results,
  scope = "significant",
  significance_col = SIGNIFICANCE_COLUMN,
  significance_threshold = SIGNIFICANCE_THRESHOLD
)
log_info(paste0("Identified  ", nrow(mr_significant), "  results meeting  ", SIGNIFICANCE_COLUMN, " < ", SIGNIFICANCE_THRESHOLD, "  among the primary results"))

if (nrow(mr_results) > 0) {
  log_info("Step5.2: Run sensitivity analyses")
  
  sensitivity_pairs <- select_mr_pairs(
    mr_results = mr_results,
    scope = SENSITIVITY_SCOPE,
    significance_col = SIGNIFICANCE_COLUMN,
    significance_threshold = SIGNIFICANCE_THRESHOLD
  )
  log_info(paste0("Sensitivity-analysis pairs: ", nrow(sensitivity_pairs)))

  harmonise_sens <- harmonise_data_filtered %>%
    add_exposure_outcome_id() %>%
    dplyr::filter(exposure_outcome %in% paste(sensitivity_pairs$id.exposure, sensitivity_pairs$id.outcome, sep = "_"))
  
  weighted_median_results <- NULL
  conmix_results <- NULL
  mrpresso_results <- NULL
  sensitivity_tests <- NULL
  i2_stats <- NULL

  if (nrow(harmonise_sens) > 0) {
    # Weighted Median
    weighted_median_results <- TwoSampleMR::mr(harmonise_sens, method_list = c("mr_weighted_median"))
    safe_write_csv(weighted_median_results, OUTPUT_FILES$sens_weighted_median)
    
    # Contamination Mixture
    conmix_results <- run_conmix(harmonise_sens, sensitivity_pairs)
    if (!is.null(conmix_results)) {
      safe_write_csv(conmix_results, OUTPUT_FILES$sens_conmix)
    }
    
    # MR-PRESSO
    mrpresso_pairs <- select_mr_pairs(
      mr_results = mr_results,
      scope = MRPRESSO_SCOPE,
      significance_col = SIGNIFICANCE_COLUMN,
      significance_threshold = SIGNIFICANCE_THRESHOLD
    ) %>%
      dplyr::filter(nsnp >= MRPRESSO_MIN_NSNP)

    harmonise_presso <- harmonise_data_filtered %>%
      add_exposure_outcome_id() %>%
      dplyr::filter(exposure_outcome %in% paste(mrpresso_pairs$id.exposure, mrpresso_pairs$id.outcome, sep = "_"))

    if (nrow(mrpresso_pairs) > 0 && nrow(harmonise_presso) > 0) {
      mrpresso_results <- run_mrpresso(harmonise_presso, mrpresso_pairs, nb_dist = MRPRESSO_NB_DIST)
      if (!is.null(mrpresso_results) && nrow(mrpresso_results) > 0) {
        safe_write_csv(mrpresso_results, OUTPUT_FILES$sens_mrpresso)
      } else {
        log_warning("MR-PRESSO produced no valid result; output was not saved")
      }
    } else {
      log_info("No pair met the MR-PRESSO scope and SNP-count requirements; MR-PRESSO was skipped")
    }
    
    sensitivity_tests <- run_sensitivity_tests(harmonise_sens, outcome_prevalence = OUTCOME_PREVALENCE)
    safe_write_csv(sensitivity_tests$heterogeneity, OUTPUT_FILES$sens_heterogeneity)
    safe_write_csv(sensitivity_tests$pleiotropy, OUTPUT_FILES$sens_pleiotropy)
    safe_write_csv(sensitivity_tests$steiger, OUTPUT_FILES$sens_steiger)
    
    i2_stats <- calculate_i2_statistics(harmonise_sens)
    safe_write_csv(i2_stats, OUTPUT_FILES$sens_heterogeneity_isq)
  } else {
    log_warning("No harmonised data were available for sensitivity analysis")
  }
  
  final_results <- merge_mr_results(
    mr_results = mr_results,
    conmix_results = conmix_results,
    mrpresso_results = mrpresso_results,
    weighted_median_results = weighted_median_results,
    sensitivity_results = if (!is.null(sensitivity_tests)) list(
      pleiotropy = sensitivity_tests$pleiotropy,
      steiger = sensitivity_tests$steiger
    ) else NULL,
    i2_results = i2_stats,
    output_file = OUTPUT_FILES$results_sens
  )
  
  if ("b" %in% names(final_results)) {
    final_results_or <- TwoSampleMR::generate_odds_ratios(final_results)
    safe_write_csv(final_results_or, OUTPUT_FILES$results_or)
  }
}

} else {
  log_info("Skip step 5: MR analysis already complete")
}

if (START_FROM_STEP <= 6) {
  log_info("Step6: Filter results using sensitivity-analysis criteria")
  
  if (START_FROM_STEP == 6) {
    if (file.exists(OUTPUT_FILES$results_or)) {
      final_results_or <- data.table::fread(OUTPUT_FILES$results_or)
      log_info(paste0("  Loaded OR result file: ", nrow(final_results_or), " rows"))
    } else {
      log_error(paste0("OR result file not found: ", OUTPUT_FILES$results_or))
      stop("Run Step 5 before loading the OR result file")
    }
  }
  
  if (exists("final_results_or") && nrow(final_results_or) > 0) {
    final_results_filtered <- filter_by_sensitivity(
      mr_results = final_results_or,
      thresholds = SENSITIVITY_THRESHOLDS
    )
    
    if (nrow(final_results_filtered) > 0) {
      safe_write_csv(final_results_filtered, OUTPUT_FILES$results_or_filtered)
      log_info(paste0("Results passing sensitivity checks saved: ", OUTPUT_FILES$results_or_filtered))
      
      unique_pairs <- final_results_filtered %>%
        dplyr::select(exposure, outcome) %>%
        dplyr::distinct() %>%
        nrow()
      log_info(paste0("  A total of ", unique_pairs, " causal associations passed sensitivity checks"))
    } else {
      log_warning("No result passed the sensitivity checks")
    }
  } else {
    log_warning("OR results were not found; filtering was skipped")
  }
  
} else {
  log_info("Skip step 6: sensitivity filtering already complete")
}

if (START_FROM_STEP <= 7) {
  log_info("Step7: Generate forest plot")
  
  if (START_FROM_STEP == 7) {
    if (!file.exists(OUTPUT_FILES$results_or_filtered)) {
      log_error(paste0("Filtered OR result file not found: ", OUTPUT_FILES$results_or_filtered))
      log_warning("Attempting to plot the unfiltered OR results...")
      
      if (file.exists(OUTPUT_FILES$results_or)) {
        data_for_plot <- OUTPUT_FILES$results_or
        log_info("  Using unfiltered OR results")
      } else {
        log_error(paste0("OR result file not found: ", OUTPUT_FILES$results_or))
        stop("Run Step 5 or 6 before loading the OR result file")
      }
    } else {
      data_for_plot <- OUTPUT_FILES$results_or_filtered
      log_info("  Using filtered OR results")
    }
  } else {
    if (exists("final_results_filtered") && nrow(final_results_filtered) > 0) {
      data_for_plot <- OUTPUT_FILES$results_or_filtered
      log_info("  Using filtered OR results")
    } else if (file.exists(OUTPUT_FILES$results_or)) {
      data_for_plot <- OUTPUT_FILES$results_or
      log_info("  Using unfiltered OR results")
    } else {
      log_warning("OR result file not found; plotting skipped")
      data_for_plot <- NULL
    }
  }
  
  if (!is.null(data_for_plot) && file.exists(data_for_plot)) {
    tryCatch({
      plot_forestplot_smart(
        data_file = data_for_plot,
        output_file = OUTPUT_FILES$plot_forestplot,
        font_family = "sans"
      )
      log_info(paste0("Forest plot saved: ", OUTPUT_FILES$plot_forestplot))
    }, error = function(e) {
      log_error(paste0("Generate forest plot failed: ", e$message))
      log_warning("If a font error occurs, ensure that a sans-serif font is installed")
    })
  } else {
    log_warning("Data file not found; plotting skipped")
  }
  
} else {
  log_info("Skip step 7: plotting already complete")
}

log_info(strrep("=", 80))
log_info("Analysis complete.")
log_info("Principal result files:")
if (file.exists(OUTPUT_FILES$results_main)) {
  log_info(paste0("  - IVW results: ", OUTPUT_FILES$results_main))
}
if (file.exists(OUTPUT_FILES$results_sens)) {
  log_info(paste0("  - Sensitivity analyses: ", OUTPUT_FILES$results_sens))
}
if (file.exists(OUTPUT_FILES$results_or)) {
  log_info(paste0("  - OR values (all): ", OUTPUT_FILES$results_or))
}
if (file.exists(OUTPUT_FILES$results_or_filtered)) {
  log_info(paste0("  - OR values (filtered): ", OUTPUT_FILES$results_or_filtered))
}
if (file.exists(OUTPUT_FILES$plot_forestplot)) {
  log_info(paste0("  - Forest plot: ", OUTPUT_FILES$plot_forestplot))
}
log_info(strrep("=", 80))
