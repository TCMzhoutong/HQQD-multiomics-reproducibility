# ==============================================================================
# ==============================================================================
# ==============================================================================

#' 
#' @export
parallel_mr_analysis <- function(harmonise_data,
                                  n_cores = 4,
                                  methods = c("mr_ivw_mre")) {
  
  message(paste0("Starting parallel MR analysis with  ", n_cores, "  cores..."))
  
  dat_list <- split(harmonise_data, harmonise_data$id.exposure)
  message(paste0("A total of ", length(dat_list), " exposure-outcome pairs"))
  
  mr_speed <- function(dat) {
    TwoSampleMR::mr(dat, method_list = methods)
  }
  
  if (n_cores <= 1) {
    start_time <- Sys.time()
    results <- lapply(dat_list, mr_speed)
    end_time <- Sys.time()

    results <- do.call(rbind, results)
    results$fdr <- p.adjust(results$pval, method = "fdr")

    message(paste0("MR analysis completed in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds"))
    message(paste0("Results: ", nrow(results)))
    message(paste0("p < 0.05 results: ", sum(results$pval < 0.05)))

    return(results)
  }

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, library(TwoSampleMR))
  parallel::clusterExport(cl = cl, varlist = c("dat_list", "mr_speed", "methods"), envir = environment())
  
  start_time <- Sys.time()
  results <- parallel::parLapply(cl = cl, X = dat_list, fun = mr_speed)
  end_time <- Sys.time()
  
  parallel::stopCluster(cl)
  
  results <- do.call(rbind, results)
  results$fdr <- p.adjust(results$pval, method = "fdr")
  
  message(paste0("MR analysis complete; elapsed time: ", round(difftime(end_time, start_time, units = "secs"), 2), "  seconds"))
  message(paste0("A total of ", nrow(results), " results"))
  message(paste0("Results with p < 0.05: ", sum(results$pval < 0.05)))
  
  return(results)
}


select_mr_pairs <- function(mr_results,
                            scope = "significant",
                            significance_col = "fdr",
                            significance_threshold = 0.05,
                            method_pattern = "Inverse variance weighted") {
  if (is.null(mr_results) || nrow(mr_results) == 0) {
    return(mr_results[0, ])
  }

  if (scope == "all") {
    return(mr_results)
  }

  if (!significance_col %in% names(mr_results)) {
    warning(paste0("Significance column does not exist: ", significance_col, ", using pval instead"))
    significance_col <- "pval"
  }

  mr_results %>%
    dplyr::filter(
      grepl(method_pattern, method, ignore.case = TRUE),
      !is.na(.data[[significance_col]]),
      .data[[significance_col]] < significance_threshold
    )
}


#' 
#' @export
run_conmix <- function(harmonise_data, significant_pairs) {
  
  message("Starting contamination-mixture analysis...")
  
  harmonise_data <- harmonise_data %>%
    dplyr::mutate(exposure_outcome = paste(id.exposure, id.outcome, sep = "_")) %>%
    dplyr::filter(exposure_outcome %in% paste(significant_pairs$id.exposure, significant_pairs$id.outcome, sep = "_"))
  
  split_data <- split(harmonise_data, harmonise_data$exposure_outcome)
  results <- list()
  
  for (name in names(split_data)) {
    cat("Processing: ", name, "\n")
    current_data <- split_data[[name]]
    
    tryCatch({
      conmix_result <- MendelianRandomization::mr_conmix(
        MendelianRandomization::mr_input(
          bx = current_data$beta.exposure,
          bxse = current_data$se.exposure,
          by = current_data$beta.outcome,
          byse = current_data$se.outcome,
          snps = current_data$SNP
        )
      )
      
      se_value <- (conmix_result@CIUpper - conmix_result@CILower) / (2 * 1.96)
      
      results[[name]] <- data.frame(
        id.exposure = unique(current_data$id.exposure),
        id.outcome = unique(current_data$id.outcome),
        exposure = unique(current_data$exposure),
        outcome = unique(current_data$outcome),
        b = conmix_result@Estimate,
        se = se_value,
        pval = conmix_result@Pvalue,
        method = "Contamination mixture"
      )
    }, error = function(e) {
      cat("Processing", name, " failed:", e$message, "\n")
    })
  }
  
  if (length(results) > 0) {
    results_df <- do.call(rbind, results)
    message(paste0("ConMix analysis complete:  ", nrow(results_df), "  results"))
    return(results_df)
  } else {
    message("ConMix produced no valid result")
    return(NULL)
  }
}


#' 
#' @export
run_mrpresso <- function(harmonise_data, significant_pairs, nb_dist = 1000) {
  
  message(paste0("Starting MR-PRESSO analysis (simulations: ", nb_dist, ")..."))
  
  harmonise_data <- harmonise_data %>%
    dplyr::mutate(exposure_outcome = paste(id.exposure, id.outcome, sep = "_")) %>%
    dplyr::filter(exposure_outcome %in% paste(significant_pairs$id.exposure, significant_pairs$id.outcome, sep = "_"))
  
  harmonise_data <- as.data.frame(harmonise_data)
  split_data <- split(harmonise_data, harmonise_data$exposure_outcome)
  
  MRPRESSO <- data.frame()
  set.seed(1234)
  total <- length(split_data)
  
  for (i in seq_along(split_data)) {
    name <- names(split_data)[i]
    cat("Processing:", name, "(", i, "/", total, ")\n")
    current_data <- split_data[[name]]
    
    tryCatch({
      mr_presso <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                                       BetaExposure = "beta.exposure",
                                       SdOutcome = "se.outcome",
                                       SdExposure = "se.exposure",
                                       OUTLIERtest = TRUE,
                                       DISTORTIONtest = TRUE,
                                       data = current_data,
                                       NbDistribution = nb_dist,
                                       SignifThreshold = 0.05)
      
      if (i == 1) {
        cat("\n=== MR-PRESSO return-structure diagnostics ===\n")
        cat("mr_pressoClass:", class(mr_presso), "\n")
        cat("mr_pressoLength:", length(mr_presso), "\n")
        cat("\nFirst element:\n")
        cat("  Class:", class(mr_presso[[1]]), "\n")
        if (is.data.frame(mr_presso[[1]])) {
          cat("  Column names:", names(mr_presso[[1]]), "\n")
          cat("  rows:", nrow(mr_presso[[1]]), "\n")
          print(mr_presso[[1]])
        }
        cat("\nSecond element:\n")
        cat("  Class:", class(mr_presso[[2]]), "\n")
        if (is.list(mr_presso[[2]])) {
          cat("  Element names:", names(mr_presso[[2]]), "\n")
        }
        cat("===============================\n\n")
      }
      
      main_mr_results <- mr_presso[[1]]
      
      if (length(mr_presso) >= 2 && is.list(mr_presso[[2]])) {
        presso_pval <- mr_presso[[2]][["Global Test"]][["Pvalue"]]
        if (is.null(presso_pval)) {
          presso_pval <- NA
        }
      } else {
        presso_pval <- NA
      }
      
      if (is.null(main_mr_results) || !is.data.frame(main_mr_results) || nrow(main_mr_results) == 0) {
        stop("MR-PRESSO did not return valid Main MR results")
      }
      
      n_rows <- nrow(main_mr_results)
      has_corrected <- FALSE
      
      if (!is.na(n_rows) && n_rows >= 2) {
        causal_est_2 <- main_mr_results$`Causal Estimate`[2]
        if (length(causal_est_2) > 0 && !is.na(causal_est_2)) {
          has_corrected <- TRUE
        }
      }
      
      if (has_corrected) {
        b <- main_mr_results$`Causal Estimate`[2]
        se <- main_mr_results$Sd[2]
        pval <- main_mr_results$`P-value`[2]
      } else {
        b <- main_mr_results$`Causal Estimate`[1]
        se <- main_mr_results$Sd[1]
        pval <- main_mr_results$`P-value`[1]
      }
      
      MRPRESSO <- rbind(MRPRESSO, data.frame(
        id.exposure = unique(current_data$id.exposure),
        id.outcome = unique(current_data$id.outcome),
        exposure = unique(current_data$exposure),
        outcome = unique(current_data$outcome),
        b = b,
        se = se,
        pval = pval,
        PRESSO_pval = presso_pval,
        method = "MR-PRESSO (Outlier-corrected)"
      ))
    }, error = function(e) {
      cat("Processing", name, " failed:", e$message, "\n")
    })
  }
  
  message(paste0("MR-PRESSO analysis complete:  ", nrow(MRPRESSO), "  results"))
  
  if (nrow(MRPRESSO) == 0) {
    message("Warning: all MR-PRESSO analyses failed; returning NULL")
    return(NULL)
  }
  
  return(MRPRESSO)
}


#' 
#' @export
run_sensitivity_tests <- function(harmonise_data, outcome_prevalence = NULL) {
  
  message("Starting sensitivity analyses...")
  
  message("  - Calculate heterogeneity...")
  heterogeneity <- TwoSampleMR::mr_heterogeneity(harmonise_data)
  
  message("  - Calculate pleiotropy...")
  pleiotropy <- TwoSampleMR::mr_pleiotropy_test(harmonise_data)
  
  message("  - Run Steiger directionality test...")
  
  harmonise_data$r.exposure <- TwoSampleMR::get_r_from_bsen(
    b = harmonise_data$beta.exposure,
    se = harmonise_data$se.exposure,
    n = harmonise_data$samplesize.exposure
  )
  
  has_binary_outcome_info <- all(c("ncase.outcome", "ncontrol.outcome") %in% names(harmonise_data))
  
  if (has_binary_outcome_info) {
    if (!is.null(outcome_prevalence) && !is.na(outcome_prevalence)) {
      harmonise_data$prevalence <- outcome_prevalence
      message(paste0("  - Using the configured population prevalence to calculate r.outcome: ", outcome_prevalence))
    } else {
      harmonise_data$prevalence <- harmonise_data$ncase.outcome / harmonise_data$samplesize.outcome
      warning("OUTCOME_PREVALENCE is not configured; Steiger will use case/(case+control), which is not equivalent to population prevalence in a case-control GWAS")
    }
    harmonise_data$r.outcome <- mapply(
      TwoSampleMR::get_r_from_lor,
      lor = harmonise_data$beta.outcome,
      af = harmonise_data$eaf.outcome,
      ncase = harmonise_data$ncase.outcome,
      ncontrol = harmonise_data$ncontrol.outcome,
      prevalence = harmonise_data$prevalence
    )
  } else {
    message("  - Continuous outcome detected (no ncase/ncontrol); calculating r.outcome using the BSEN formula")
    harmonise_data$r.outcome <- TwoSampleMR::get_r_from_bsen(
      b = harmonise_data$beta.outcome,
      se = harmonise_data$se.outcome,
      n = harmonise_data$samplesize.outcome
    )
  }
  
  steiger <- tryCatch({
    TwoSampleMR::directionality_test(harmonise_data)
  }, error = function(e) {
    warning(paste0("Steiger directionality test failed: ", e$message))
    data.frame()
  })
  
  message("Sensitivity analyses completed")
  
  return(list(
    heterogeneity = heterogeneity,
    pleiotropy = pleiotropy,
    steiger = steiger
  ))
}


#' 
#' @export
calculate_i2_statistics <- function(harmonise_data) {
  
  message("Calculating I-squared statistics...")
  
  rs_single <- TwoSampleMR::mr_singlesnp(harmonise_data)
  rs_single2 <- rs_single[grep("^rs", rs_single$SNP), ]
  
  unique_combinations <- unique(rs_single2[c("id.exposure", "id.outcome", "exposure", "outcome")])
  rs_isq <- data.frame()
  
  for (i in 1:nrow(unique_combinations)) {
    current_combination <- unique_combinations[i, ]
    
    subset_rs_single2 <- subset(
      rs_single2,
      id.exposure == current_combination$id.exposure &
        id.outcome == current_combination$id.outcome &
        exposure == current_combination$exposure &
        outcome == current_combination$outcome
    )
    
    subset_harmonise_data <- subset(
      harmonise_data,
      id.exposure == current_combination$id.exposure &
        id.outcome == current_combination$id.outcome &
        exposure == current_combination$exposure &
        outcome == current_combination$outcome
    )
    
    rs_meta <- metafor::rma(
      yi = subset_rs_single2$b,
      sei = subset_rs_single2$se,
      weights = 1 / subset_harmonise_data$se.outcome^2,
      data = subset_rs_single2,
      method = "FE"
    )
    
    temp_df <- data.frame(
      id.exposure = current_combination$id.exposure,
      id.outcome = current_combination$id.outcome,
      exposure = current_combination$exposure,
      outcome = current_combination$outcome,
      I2 = as.numeric(rs_meta$I2),
      H2 = as.numeric(rs_meta$H2),
      QE = as.numeric(rs_meta$QE),
      QEp = as.numeric(rs_meta$QEp)
    )
    
    rs_isq <- rbind(rs_isq, temp_df)
  }
  
  message(paste0("I-squared calculation complete:  ", nrow(rs_isq), "  results"))
  return(rs_isq)
}


#' 
#' @export
merge_mr_results <- function(mr_results,
                              conmix_results = NULL,
                              mrpresso_results = NULL,
                              weighted_median_results = NULL,
                              sensitivity_results = NULL,
                              i2_results = NULL,
                              output_file = NULL) {
  
  message("Merging MR results...")
  
  all_methods <- list(mr_results)
  if (!is.null(conmix_results)) all_methods <- c(all_methods, list(conmix_results))
  if (!is.null(mrpresso_results)) all_methods <- c(all_methods, list(mrpresso_results))
  if (!is.null(weighted_median_results)) all_methods <- c(all_methods, list(weighted_median_results))
  
  merged_data <- data.table::rbindlist(all_methods, fill = TRUE)

  method_keys <- intersect(
    c("id.exposure", "id.outcome", "exposure", "outcome", "method"),
    names(merged_data)
  )
  if (length(method_keys) > 0) {
    duplicated_methods <- merged_data[, .N, by = method_keys][N > 1]
    if (nrow(duplicated_methods) > 0) {
      message(paste0(
        "  Removed ", sum(duplicated_methods$N - 1),
        " duplicate MR method rows before merging sensitivity metrics"
      ))
      merged_data <- unique(merged_data, by = method_keys)
    }
  }

  prepare_sensitivity_join <- function(dat, keep_cols, label) {
    if (is.null(dat) || nrow(dat) == 0) {
      return(NULL)
    }

    dt <- data.table::as.data.table(dat)
    keys <- intersect(c("id.exposure", "id.outcome", "exposure", "outcome"), names(dt))
    keys <- intersect(keys, names(merged_data))
    if (length(keys) == 0) {
      stop(paste0("No shared merge keys found for ", label))
    }

    cols <- unique(c(keys, keep_cols))
    cols <- intersect(cols, names(dt))
    dt <- dt[, ..cols]

    duplicated_keys <- dt[, .N, by = keys][N > 1]
    if (nrow(duplicated_keys) > 0) {
      message(paste0(
        "  Removed ", sum(duplicated_keys$N - 1),
        " duplicate ", label, " rows before merge"
      ))
      dt <- unique(dt, by = keys)
    }

    list(data = dt, keys = keys)
  }
  
  if (!is.null(sensitivity_results)) {
    if (!is.null(sensitivity_results$pleiotropy)) {
      pleiotropy_data <- sensitivity_results$pleiotropy
      pleiotropy_data$egger_pval <- pleiotropy_data$pval
      pleiotropy_data$pval <- NULL
      join_data <- prepare_sensitivity_join(
        pleiotropy_data,
        c("egger_intercept", "egger_pval"),
        "pleiotropy"
      )
      merged_data <- merge(merged_data, join_data$data, by = join_data$keys, all.x = TRUE)
    }
    
    if (!is.null(i2_results)) {
      join_data <- prepare_sensitivity_join(
        i2_results,
        c("QE", "QEp", "I2"),
        "I2"
      )
      merged_data <- merge(merged_data, join_data$data, by = join_data$keys, all.x = TRUE)
    }
    
    if (!is.null(sensitivity_results$steiger)) {
      join_data <- prepare_sensitivity_join(
        sensitivity_results$steiger,
        c("correct_causal_direction", "steiger_pval"),
        "Steiger"
      )
      merged_data <- merge(merged_data, join_data$data, by = join_data$keys, all.x = TRUE)
    }
    
    if (!is.null(mrpresso_results) && nrow(mrpresso_results) > 0 && "PRESSO_pval" %in% names(mrpresso_results)) {
      join_data <- prepare_sensitivity_join(
        mrpresso_results,
        "PRESSO_pval",
        "MR-PRESSO"
      )
      merged_data <- merge(merged_data, join_data$data, by = join_data$keys, all.x = TRUE)
    }
  }
  
  if (!is.null(output_file)) {
    safe_write_csv(merged_data, output_file, row.names = FALSE)
  }
  
  message(paste0("Result merging complete: ", nrow(merged_data), " rows"))
  return(merged_data)
}

#' 
#' @export
filter_by_sensitivity <- function(mr_results, thresholds = NULL) {
  
  if (is.null(thresholds)) {
    message("No filtering threshold supplied; using defaults")
    thresholds <- list(
      egger_pval_min = 0.05,
      heterogeneity_qep_min = 0.05,
      i2_max = 50,
      presso_pval_min = 0.05,
      require_correct_direction = TRUE
    )
  }
  
  message("Filtering results using sensitivity-analysis criteria...")
  message(paste0("  Before filtering: ", nrow(mr_results), " rows"))
  
  if ("PRESSO_pval.x" %in% names(mr_results) && !"PRESSO_pval" %in% names(mr_results)) {
    mr_results$PRESSO_pval <- mr_results$PRESSO_pval.x
  }
  if ("PRESSO_pval.y" %in% names(mr_results) && !"PRESSO_pval" %in% names(mr_results)) {
    mr_results$PRESSO_pval <- mr_results$PRESSO_pval.y
  }
  if ("PRESSO_pval.x" %in% names(mr_results) && "PRESSO_pval.y" %in% names(mr_results)) {
    mr_results$PRESSO_pval <- ifelse(!is.na(mr_results$PRESSO_pval.y), 
                                     mr_results$PRESSO_pval.y, 
                                     mr_results$PRESSO_pval.x)
  }
  
  original_count <- nrow(mr_results)
  require_primary <- isTRUE(thresholds$require_primary_significance)
  significance_col <- thresholds$significance_col
  if (is.null(significance_col) || !significance_col %in% names(mr_results)) {
    significance_col <- "pval"
  }
  significance_threshold <- thresholds$significance_threshold
  if (is.null(significance_threshold)) {
    significance_threshold <- 0.05
  }
  primary_method_pattern <- thresholds$primary_method_pattern
  if (is.null(primary_method_pattern)) {
    primary_method_pattern <- "Inverse variance weighted"
  }
  
  mr_results <- mr_results %>%
    dplyr::group_by(exposure, outcome) %>%
    dplyr::mutate(
      pass_egger = if("egger_pval" %in% names(.)) {
        is.na(egger_pval) | egger_pval > thresholds$egger_pval_min
      } else { TRUE },
      
      pass_heterogeneity = if("QEp" %in% names(.)) {
        is.na(QEp) | QEp > thresholds$heterogeneity_qep_min
      } else { TRUE },
      
      pass_i2 = if("I2" %in% names(.)) {
        is.na(I2) | I2 < thresholds$i2_max
      } else { TRUE },
      
      pass_presso = if("PRESSO_pval" %in% names(.)) {
        is.na(PRESSO_pval) | PRESSO_pval > thresholds$presso_pval_min
      } else { TRUE },
      
      pass_steiger = if(thresholds$require_correct_direction && "correct_causal_direction" %in% names(.)) {
        is.na(correct_causal_direction) | correct_causal_direction == TRUE
      } else { TRUE },

      pass_primary_significance = if(require_primary) {
        any(
          grepl(primary_method_pattern, method, ignore.case = TRUE) &
            !is.na(.data[[significance_col]]) &
            .data[[significance_col]] < significance_threshold
        )
      } else { TRUE },
      
      pass_all = pass_egger & pass_heterogeneity & pass_i2 & pass_presso & pass_steiger & pass_primary_significance
    ) %>%
    dplyr::filter(all(pass_all)) %>%
    dplyr::select(-pass_egger, -pass_heterogeneity, -pass_i2, -pass_presso, -pass_steiger, -pass_primary_significance, -pass_all) %>%
    dplyr::ungroup()
  
  filtered_count <- nrow(mr_results)
  removed_count <- original_count - filtered_count
  
  message(paste0("  After filtering: ", filtered_count, " rows"))
  message(paste0("  Removed: ", removed_count, " rows"))
  
  if (removed_count > 0) {
    message("  Filtering criteria:")
    message(paste0("    - Egger P value > ", thresholds$egger_pval_min))
    message(paste0("    - Heterogeneity QEp > ", thresholds$heterogeneity_qep_min))
    message(paste0("    - I² < ", thresholds$i2_max, "%"))
    message(paste0("    - PRESSO P value > ", thresholds$presso_pval_min))
    if (thresholds$require_correct_direction) {
      message("    - Steiger direction consistent")
    }
  }
  
  unique_pairs <- mr_results %>%
    dplyr::select(exposure, outcome) %>%
    dplyr::distinct() %>%
    nrow()
  
  message(paste0("  Causal pairs passing sensitivity checks: ", unique_pairs))
  
  return(mr_results)
}

# ==============================================================================
# ==============================================================================
