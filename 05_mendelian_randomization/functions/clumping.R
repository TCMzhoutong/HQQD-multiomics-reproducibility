# ==============================================================================
# ==============================================================================
# ==============================================================================

#' 
#' @export
perform_ld_clump <- function(dat, 
                              plink_path, 
                              bfile_path,
                              clump_kb = 10000,
                              clump_r2 = 0.01,
                              clump_p = 1) {
  
  message(paste0("Starting LD clumping; parameters: kb=", clump_kb, ", r²=", clump_r2, ", p=", clump_p))
  message(paste0("Input SNPs: ", nrow(dat)))
  
  unique_exposures <- unique(dat$id.exposure)
  message(paste0("Contains  ", length(unique_exposures), "  independent exposures"))
  
  all_clump_results <- list()
  success_count <- 0
  fail_count <- 0
  
  for (exp_id in unique_exposures) {
    exp_data <- dat[dat$id.exposure == exp_id, ]
    n_snps <- nrow(exp_data)
    
    message(paste0("\nProcessing ", exp_id, " (", n_snps, " SNPs)..."))
    
    tryCatch({
      clump_result <- ieugwasr::ld_clump(
        dplyr::tibble(
          rsid = exp_data$SNP,
          pval = exp_data$pval.exposure,
          id = exp_data$id.exposure
        ),
        clump_kb = clump_kb,
        clump_r2 = clump_r2,
        clump_p = clump_p,
        pop = "EUR",
        plink_bin = plink_path,
        bfile = bfile_path
      )
      
      if (!is.null(clump_result) && nrow(clump_result) > 0) {
        all_clump_results[[exp_id]] <- clump_result
        success_count <- success_count + 1
        message(paste0("  ✓ Retained  ", nrow(clump_result), "  independent SNPs"))
      } else {
        warning(paste0("  ✗ ", exp_id, "  clumping produced no result"))
        fail_count <- fail_count + 1
      }
      
    }, error = function(e) {
      warning(paste0("  ✗ ", exp_id, "  clumping failed: ", e$message))
      fail_count <<- fail_count + 1
    })
  }
  
  if (length(all_clump_results) == 0) {
    warning("Clumping failed for all exposures.")
    message(paste0("\nSummary: successful=0, failed=", fail_count))
    message(paste0("SNPs after clumping: 0"))
    message(paste0("Removed SNPs: ", nrow(dat)))
    return(dplyr::tibble(rsid = character(0), pval = numeric(0), id = character(0)))
  }
  
  final_result <- dplyr::bind_rows(all_clump_results)
  
  message(paste0("\nSummary: successful=", success_count, ", failed=", fail_count))
  message(paste0("Total SNPs after clumping: ", nrow(final_result)))
  message(paste0("Removed SNPs: ", nrow(dat) - nrow(final_result)))
  
  return(final_result)
}


filter_cis_variants <- function(dat, cis_filter_params) {
  if (is.null(cis_filter_params)) {
    return(dat)
  }

  window <- cis_filter_params$window
  if (is.null(window)) {
    window <- 1000000
  }

  required_cols <- c(
    cis_filter_params$chr_col,
    cis_filter_params$pos_col,
    cis_filter_params$gene_chr_col,
    cis_filter_params$gene_start_col,
    cis_filter_params$gene_end_col
  )
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop(paste0("cis filtering is missing required columns: ", paste(missing_cols, collapse = ", ")))
  }

  n_before <- nrow(dat)
  snp_chr <- as.character(dat[[cis_filter_params$chr_col]])
  gene_chr <- as.character(dat[[cis_filter_params$gene_chr_col]])
  snp_chr <- sub("^chr", "", snp_chr, ignore.case = TRUE)
  gene_chr <- sub("^chr", "", gene_chr, ignore.case = TRUE)

  snp_pos <- suppressWarnings(as.numeric(dat[[cis_filter_params$pos_col]]))
  gene_start <- suppressWarnings(as.numeric(dat[[cis_filter_params$gene_start_col]]))
  gene_end <- suppressWarnings(as.numeric(dat[[cis_filter_params$gene_end_col]]))

  keep <- !is.na(snp_pos) & !is.na(gene_start) & !is.na(gene_end) &
    snp_chr == gene_chr &
    snp_pos >= (gene_start - window) &
    snp_pos <= (gene_end + window)

  dat <- dat[keep, ]
  message(paste0("Cis filtering: ", n_before, " -> ", nrow(dat), " rows (window=", window, " bp)"))
  return(dat)
}


#' 
#' @param selected_exposures Optional exposure names to keep before LD clumping
#' @param id_map Optional named vector for exposure renaming before selection
#' @export
process_exposure_data <- function(file_path,
                                   data_source,
                                   format_params,
                                   plink_path,
                                   bfile_path,
                                   maf_path = NULL,
                                   default_samplesize = NULL,
                                   pval_threshold = 1e-5,
                                    f_threshold = 10,
                                    clump_params = list(kb = 10000, r2 = 0.01, p = 1),
                                    need_eaf = FALSE,
                                    selected_exposures = NULL,
                                    id_map = NULL,
                                    cis_filter_params = NULL) {
  
  message("\n", strrep("=", 80))
  message("Processing data source: ", data_source)
  message(strrep("=", 80))
  
  message("Step1: Read data file...")
  exp_data <- data.table::fread(file_path, header = TRUE)
  message(paste0("Raw-data rows: ", nrow(exp_data)))
  
  message(paste0("Step2: Filter SNPs at p < ", pval_threshold, "  SNPs..."))
  
  if (!is.null(format_params$log_pval) && format_params$log_pval == TRUE) {
    threshold_log10 <- -log10(pval_threshold)
    message(paste0("  (LOG10P format detected; threshold: LOG10P > ", threshold_log10, ")"))
    exp_data <- subset(exp_data, get(format_params$pval_col) > threshold_log10)
  } else {
    exp_data <- subset(exp_data, get(format_params$pval_col) < pval_threshold)
  }
  
  message(paste0("Rows after filtering: ", nrow(exp_data)))

  if (is.null(format_params$snp_col) || !format_params$snp_col %in% names(exp_data)) {
    stop(paste0("SNP column not found: ", format_params$snp_col))
  }
  snp_values <- exp_data[[format_params$snp_col]]
  snp_ids <- trimws(as.character(snp_values))
  invalid_snp <- is.na(snp_values) | snp_ids %in% c("", ".", "NA", "N/A", "NULL")
  if (any(invalid_snp)) {
    message(paste0("  Removed ", sum(invalid_snp), " rows with blank/placeholder SNP IDs before formatting"))
    exp_data <- exp_data[!invalid_snp, ]
  }

  exp_data <- filter_cis_variants(exp_data, cis_filter_params)

  if (!is.null(default_samplesize)) {
    if (is.null(format_params$samplesize_col) || !format_params$samplesize_col %in% names(exp_data)) {
      exp_data$samplesize <- default_samplesize
      format_params$samplesize_col <- "samplesize"
      message(paste0("  Sample-size column missing; using fixed sample size: ", default_samplesize))
    }
  }
  
  message("Step3: Format data...")
  
  if ("data.table" %in% class(exp_data)) {
    exp_data <- as.data.frame(exp_data)
  }
  
  exp_formatted <- do.call(TwoSampleMR::format_data, c(list(dat = exp_data, type = "exposure"), format_params))
  
  message(paste0("  Rows after formatting: ", nrow(exp_formatted)))
  message(paste0("  Unique SNPs: ", length(unique(exp_formatted$SNP))))
  message(paste0("  Unique exposures: ", length(unique(exp_formatted$id.exposure))))

  if (!is.null(id_map)) {
    for (i in seq_along(id_map)) {
      exp_formatted$exposure[exp_formatted$exposure == names(id_map)[i]] <- unname(id_map[i])
    }
  }

  if (!is.null(selected_exposures)) {
    n_before <- nrow(exp_formatted)
    n_exp_before <- length(unique(exp_formatted$exposure))
    exp_formatted <- exp_formatted[exp_formatted$exposure %in% selected_exposures, ]
    message(paste0(
      "  Selected exposure filter before clumping: ",
      n_exp_before, " exposures/", n_before, " rows -> ",
      length(unique(exp_formatted$exposure)), " exposures/", nrow(exp_formatted), " rows"
    ))
    if (nrow(exp_formatted) == 0) {
      warning(paste0("Data source ", data_source, ": no data remained after selected-exposure filtering"))
      return(NULL)
    }
  }
  
  dup_check <- duplicated(exp_formatted[, c("SNP", "id.exposure")])
  if (any(dup_check)) {
    n_dup <- sum(dup_check)
    message(paste0("  Warning: Identified  ", n_dup, "  duplicate (SNP, exposure) pairs."))
  }
  
  if (need_eaf && !is.null(maf_path)) {
    message("Step4: Supplement EAF from 1000 Genomes...")
    source("functions/utils.R")
    exp_formatted <- get_eaf_from_1000G(exp_formatted, maf_path, type = "exposure")
  } else {
    message("Step4: Skip EAF supplementation")
  }
  
  # 5. LD Clumping
  message("Step5: Run LD clumping...")
  clump_result <- perform_ld_clump(
    dat = exp_formatted,
    plink_path = plink_path,
    bfile_path = bfile_path,
    clump_kb = clump_params$kb,
    clump_r2 = clump_params$r2,
    clump_p = clump_params$p
  )
  
  if (is.null(clump_result) || nrow(clump_result) == 0) {
    warning(paste0("Data source ", data_source, ": no SNPs remained after clumping; subsequent processing skipped"))
    return(NULL)
  }
  
  message(paste0("  exp_formattedrows: ", nrow(exp_formatted)))
  message(paste0("  exp_formattedUnique SNPs: ", length(unique(exp_formatted$SNP))))
  message(paste0("  exp_formattedUnique exposures: ", length(unique(exp_formatted$id.exposure))))
  message(paste0("  clump_resultrows: ", nrow(clump_result)))
  message(paste0("  clump_resultunique rsIDs: ", length(unique(clump_result$rsid))))
  message(paste0("  clump_resultunique IDs: ", length(unique(clump_result$id))))
  
  ivs_data <- merge(exp_formatted, clump_result[, c("rsid", "id")], 
                    by.x = c("SNP", "id.exposure"), 
                    by.y = c("rsid", "id"),
                    all.x = FALSE, all.y = FALSE)
  
  message(paste0("  After merging: rows: ", nrow(ivs_data)))
  message(paste0("  After merging: Unique SNPs: ", length(unique(ivs_data$SNP))))
  message(paste0("  After merging: Unique exposures: ", length(unique(ivs_data$id.exposure))))
  
  if (nrow(ivs_data) == 0) {
    warning(paste0("Data source ", data_source, ": no data remained after matching SNPs"))
    return(NULL)
  }
  
  message(paste0("Step6: Calculate F statistics (threshold: F > ", f_threshold, ")..."))
  source("functions/utils.R")
  ivs_data <- get_f(ivs_data, F_value = f_threshold)
  
  message("Step7: Calculate MAF...")
  ivs_data <- add_MAF_from_EAF(ivs_data)
  
  message("\n--- Final summary ---")
  message(paste0("Final SNPs: ", nrow(ivs_data)))
  if ("F" %in% names(ivs_data)) {
    message(paste0("F-statistic range: ", round(min(ivs_data$F), 2), " - ", round(max(ivs_data$F), 2)))
    message(paste0("Median F statistic: ", round(median(ivs_data$F), 2)))
  }
  if ("R2" %in% names(ivs_data)) {
    message(paste0("R-squared range: ", round(min(ivs_data$R2), 4), " - ", round(max(ivs_data$R2), 4)))
  }
  
  message(strrep("=", 80), "\n")
  
  return(ivs_data)
}


#' 
#' @export
batch_process_exposures <- function(data_configs,
                                     plink_path,
                                     bfile_path,
                                     maf_path = NULL,
                                     output_dir = "02_IVs") {
  
  results <- list()
  
  for (i in seq_along(data_configs)) {
    config <- data_configs[[i]]
    
    tryCatch({
      ivs_data <- process_exposure_data(
        file_path = config$file_path,
        data_source = config$data_source,
        format_params = config$format_params,
        plink_path = plink_path,
        bfile_path = bfile_path,
        maf_path = maf_path,
        default_samplesize = config$default_samplesize %||% NULL,
        pval_threshold = config$pval_threshold %||% 1e-5,
        f_threshold = config$f_threshold %||% 10,
        clump_params = config$clump_params %||% list(kb = 10000, r2 = 0.01, p = 1),
        need_eaf = config$need_eaf %||% FALSE,
        selected_exposures = config$selected_exposures %||% NULL,
        id_map = config$id_map %||% NULL,
        cis_filter_params = config$cis_filter_params %||% NULL
      )
      
      output_file <- file.path(output_dir, paste0("ivs_", config$output_name, ".csv"))
      safe_write_csv(ivs_data, output_file, row.names = FALSE)
      
      results[[config$data_source]] <- ivs_data
      
    }, error = function(e) {
      message(paste0("Processing ", config$data_source, "  failed: ", e$message))
      results[[config$data_source]] <- NULL
    })
  }
  
  return(results)
}


#' @keywords internal
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

# ==============================================================================
# ==============================================================================
