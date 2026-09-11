# ==============================================================================
# ==============================================================================
# ==============================================================================

map_proxy_allele_from_phase <- function(phase, target_allele) {
  if (is.na(phase) || is.na(target_allele) || nchar(target_allele) != 1) {
    return(NA_character_)
  }

  haplotypes <- strsplit(phase, "/", fixed = TRUE)[[1]]
  for (haplotype in haplotypes) {
    if (nchar(haplotype) < 2) {
      next
    }
    target_phase_allele <- substr(haplotype, 1, 1)
    proxy_phase_allele <- substr(haplotype, 2, 2)
    if (target_phase_allele == target_allele) {
      return(proxy_phase_allele)
    }
  }

  return(NA_character_)
}

get_proxy_eaf_from_outcome <- function(proxy_allele,
                                       effect_allele_outcome,
                                       other_allele_outcome,
                                       eaf_outcome) {
  if (is.na(proxy_allele) || is.na(eaf_outcome)) {
    return(NA_real_)
  }
  if (!is.na(effect_allele_outcome) && proxy_allele == effect_allele_outcome) {
    return(as.numeric(eaf_outcome))
  }
  if (!is.na(other_allele_outcome) && proxy_allele == other_allele_outcome) {
    return(1 - as.numeric(eaf_outcome))
  }
  return(NA_real_)
}

#' 
#' @export
find_proxy_candidate_snps_plink <- function(snps,
                                            plink_path,
                                            bfile_path,
                                            r2_threshold = 0.8,
                                            window_kb = 1000) {
  snps <- unique(snps[!is.na(snps) & snps != ""])
  if (length(snps) == 0) {
    return(character(0))
  }

  message(paste0("Precomputing proxy candidates for ", length(snps), " SNPs"))

  temp_dir <- tempdir()
  snp_list_file <- file.path(temp_dir, "proxy_candidate_snps.txt")
  output_prefix <- file.path(temp_dir, "proxy_candidates")
  ld_file <- paste0(output_prefix, ".ld")

  writeLines(snps, snp_list_file)

  plink_cmd <- paste(
    shQuote(plink_path),
    "--bfile", shQuote(bfile_path),
    "--r2", "in-phase", "with-freqs",
    "--ld-snp-list", shQuote(snp_list_file),
    "--ld-window-kb", window_kb,
    "--ld-window-r2", r2_threshold,
    "--out", shQuote(output_prefix)
  )

  result <- system(plink_cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE)
  on.exit({
    unlink(c(snp_list_file, paste0(output_prefix, ".log"),
             paste0(output_prefix, ".nosex"), ld_file))
  }, add = TRUE)

  if (result != 0 || !file.exists(ld_file)) {
    message("  No proxy candidate file was produced by PLINK")
    return(character(0))
  }

  ld_data <- data.table::fread(ld_file, header = TRUE)
  if (nrow(ld_data) == 0 || !"SNP_B" %in% names(ld_data)) {
    message("  PLINK found no proxy candidates")
    return(character(0))
  }

  proxy_snps <- unique(setdiff(as.character(ld_data$SNP_B), snps))
  proxy_snps <- proxy_snps[!is.na(proxy_snps) & proxy_snps != ""]
  message(paste0("  Proxy candidate SNPs: ", length(proxy_snps)))
  return(proxy_snps)
}

find_proxy_snps_plink <- function(missing_snps,
                                   outcome_data,
                                   plink_path,
                                   bfile_path,
                                   r2_threshold = 0.8,
                                   window_kb = 1000) {
  
  if (length(missing_snps) == 0) {
    message("No missing SNPs; proxy lookup skipped")
    return(data.frame())
  }
  
  message("\n", strrep("=", 80))
  message("Local PLINK proxy-SNP search")
  message(strrep("=", 80))
  message(paste0("Missing SNPs: ", length(missing_snps)))
  message(paste0("R-squared threshold: ", r2_threshold))
  message(paste0("LD window: ", window_kb, " kb"))
  
  temp_dir <- tempdir()
  snp_list_file <- file.path(temp_dir, "missing_snps.txt")
  output_prefix <- file.path(temp_dir, "proxy_ld")
  ld_file <- paste0(output_prefix, ".ld")
  
  writeLines(missing_snps, snp_list_file)
  
  message("Calling PLINK to calculate LD...")
  plink_cmd <- paste(
    shQuote(plink_path),
    "--bfile", shQuote(bfile_path),
    "--r2", "in-phase", "with-freqs",
    "--ld-snp-list", shQuote(snp_list_file),
    "--ld-window-kb", window_kb,
    "--ld-window-r2", r2_threshold,
    "--out", shQuote(output_prefix)
  )
  
  result <- system(plink_cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE)
  
  if (result != 0) {
    warning("PLINK execution failed; returning an empty result")
    return(data.frame())
  }
  
  if (!file.exists(ld_file)) {
    message("No proxy SNP was found; the SNPs may be absent from the reference panel")
    return(data.frame())
  }
  
  message("Reading LD results...")
  ld_data <- data.table::fread(ld_file, header = TRUE)
  
  if (nrow(ld_data) == 0) {
    message("PLINK found no LD relationships")
    return(data.frame())
  }
  
  message(paste0("PLINK found  ", nrow(ld_data), "  LD relationships"))
  
  proxy_data <- ld_data %>%
    dplyr::rename(
      target_snp = SNP_A,
      proxy_snp = SNP_B,
      proxy_chr = CHR_B,
      proxy_pos = BP_B,
      phase = PHASE,
      target_maf = MAF_A,
      proxy_maf = MAF_B,
      r2 = R2
    ) %>%
    dplyr::filter(target_snp != proxy_snp) %>%
    dplyr::select(target_snp, proxy_snp, proxy_chr, proxy_pos, phase, target_maf, proxy_maf, r2)
  
  message(paste0("Remaining after excluding self-pairs:  ", nrow(proxy_data), "  valid LD relationships"))
  
  if (nrow(proxy_data) == 0) {
    message("All LD relationships were self-pairs; no true proxy SNP was found")
    return(data.frame())
  }
  
  message("Filtering proxies present in the outcome data...")
  proxy_in_outcome <- proxy_data %>%
    dplyr::filter(proxy_snp %in% outcome_data$SNP)
  
  message(paste0("  Initial proxies: ", nrow(proxy_data)))
  message(paste0("  Present in outcome: ", nrow(proxy_in_outcome)))
  
  if (nrow(proxy_in_outcome) == 0) {
    message("No proxy is present in the outcome data")
    return(data.frame())
  }
  
  proxy_best <- proxy_in_outcome %>%
    dplyr::group_by(target_snp) %>%
    dplyr::slice_max(r2, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  message(paste0("Best outcome-available proxies found for ", length(unique(proxy_best$target_snp)), " SNPs"))
  if (nrow(proxy_best) > 0) {
    message("Proxy details:")
    for (i in 1:min(5, nrow(proxy_best))) {
      message(sprintf("  %s -> %s (r² = %.3f)", 
                     proxy_best$target_snp[i], 
                     proxy_best$proxy_snp[i], 
                     proxy_best$r2[i]))
    }
  }
  
  message("Adding outcome information...")
  optional_cols <- c("effect_allele.outcome", "other_allele.outcome",
                     "beta.outcome", "se.outcome", "pval.outcome", "eaf.outcome")
  available_cols <- c("SNP", intersect(optional_cols, names(outcome_data)))
  
  proxy_with_outcome <- proxy_best %>%
    dplyr::left_join(
      outcome_data %>% dplyr::select(dplyr::all_of(available_cols)),
      by = c("proxy_snp" = "SNP")
    )
  
  message(paste0("Proxies successfully assigned to ", nrow(proxy_with_outcome), " target SNPs"))
  
  tryCatch({
    unlink(c(snp_list_file, paste0(output_prefix, ".log"), 
             paste0(output_prefix, ".nosex"), ld_file))
  }, error = function(e) {
  })
  
  return(as.data.frame(proxy_with_outcome))
}


#' 
#' @export
process_proxy_snps <- function(exposure_data,
                                outcome_data,
                                plink_path,
                                bfile_path,
                                r2_threshold = 0.8,
                                save_proxy_file = FALSE,
                                output_file = NULL) {
  
  message("\n", strrep("=", 80))
  message("Starting the local PLINK proxy-SNP workflow")
  message(strrep("=", 80))
  
  missing_snps <- setdiff(exposure_data$SNP, outcome_data$SNP)
  n_missing <- length(missing_snps)
  
  message(paste0("Exposure data contain ", n_missing, " SNPs missing from the outcome data"))
  
  if (n_missing == 0) {
    message("All SNPs are present in the outcome data; no proxy search is required")
    return(list(
      exposure_with_proxy = exposure_data,
      outcome_filtered = outcome_data,
      proxy_info = data.frame()
    ))
  }
  
  proxy_results <- find_proxy_snps_plink(
    missing_snps = missing_snps,
    outcome_data = outcome_data,
    plink_path = plink_path,
    bfile_path = bfile_path,
    r2_threshold = r2_threshold,
    window_kb = 1000
  )
  
  if (nrow(proxy_results) == 0) {
    message("No valid proxy SNP was found")
    
    exposure_available <- exposure_data %>%
      dplyr::filter(SNP %in% outcome_data$SNP)
    
    outcome_filtered <- outcome_data %>%
      dplyr::filter(SNP %in% exposure_data$SNP)
    
    message(paste0("Retained ", nrow(exposure_available), " available SNPs"))
    
    return(list(
      exposure_with_proxy = exposure_available,
      outcome_filtered = outcome_filtered,
      proxy_info = data.frame()
    ))
  }
  
  message("\nPreparing proxy-SNP data...")
  
  exposure_proxy <- exposure_data %>%
    dplyr::filter(SNP %in% proxy_results$target_snp) %>%
    dplyr::left_join(proxy_results, by = c("SNP" = "target_snp")) %>%
    dplyr::mutate(
      target_snp = SNP,
      target_a1 = effect_allele.exposure,
      target_a2 = other_allele.exposure,
      proxy_effect_allele = mapply(map_proxy_allele_from_phase, phase, effect_allele.exposure),
      proxy_other_allele = mapply(map_proxy_allele_from_phase, phase, other_allele.exposure),
      proxy_effect_eaf = mapply(
        get_proxy_eaf_from_outcome,
        proxy_effect_allele,
        effect_allele.outcome,
        other_allele.outcome,
        eaf.outcome
      )
    ) %>%
    dplyr::filter(
      !is.na(proxy_effect_allele),
      !is.na(proxy_other_allele),
      proxy_effect_allele != proxy_other_allele
    ) %>%
    dplyr::mutate(
      SNP = proxy_snp,
      effect_allele.exposure = proxy_effect_allele,
      other_allele.exposure = proxy_other_allele,
      eaf.exposure = ifelse(!is.na(proxy_effect_eaf), proxy_effect_eaf, eaf.exposure),
      proxy.outcome = TRUE
    ) %>%
    dplyr::select(-proxy_snp, -proxy_chr, -proxy_pos,
                  -proxy_effect_allele, -proxy_other_allele, -proxy_effect_eaf,
                  -effect_allele.outcome, -other_allele.outcome,
                  -beta.outcome, -se.outcome, -pval.outcome, -eaf.outcome)
  
  exposure_available <- exposure_data %>%
    dplyr::filter(SNP %in% outcome_data$SNP) %>%
    dplyr::mutate(
      proxy.outcome = FALSE,
      target_snp = SNP,
      target_a1 = effect_allele.exposure,
      target_a2 = other_allele.exposure,
      r2 = NA_real_
    )
  
  exposure_combined <- dplyr::bind_rows(exposure_available, exposure_proxy)
  
  outcome_filtered <- outcome_data %>%
    dplyr::filter(SNP %in% exposure_combined$SNP)
  
  message("\n", strrep("=", 80))
  message("Proxy-SNP processing complete")
  message(strrep("=", 80))
  message(paste0("Originally available SNPs: ", nrow(exposure_available)))
  message(paste0("Proxy SNPs: ", nrow(exposure_proxy)))
  message(paste0("Total: ", nrow(exposure_combined)))
  message(paste0("SNPs without proxies: ", n_missing - nrow(exposure_proxy)))
  message(paste0("Data-loss rate: ", round((n_missing - nrow(exposure_proxy)) / nrow(exposure_data) * 100, 2), "%"))
  
  if (save_proxy_file && !is.null(output_file) && nrow(proxy_results) > 0) {
    proxy_info <- proxy_results %>%
      dplyr::select(target_snp, proxy_snp, r2, phase, target_maf, proxy_maf, proxy_chr, proxy_pos)
    
    tryCatch({
      data.table::fwrite(proxy_info, output_file)
      message(paste0("Proxy-SNP information saved: ", output_file))
    }, error = function(e) {
      warning(paste0("Failed to save proxy-SNP file: ", e$message))
    })
  }
  
  return(list(
    exposure_with_proxy = exposure_combined,
    outcome_filtered = outcome_filtered,
    proxy_info = if(nrow(proxy_results) > 0) proxy_results else data.frame()
  ))
}

# ==============================================================================
# ==============================================================================
