# ==============================================================================
# ==============================================================================
# ==============================================================================

#' 
#' @export
get_eaf_from_1000G <- function(dat, path, type = "exposure") {
  
  corrected_eaf_expo <- function(data_MAF) {
    effect <- data_MAF$effect_allele.exposure
    other <- data_MAF$other_allele.exposure
    A1 <- data_MAF$A1
    A2 <- data_MAF$A2
    MAF_num <- data_MAF$MAF
    EAF_num <- 1 - MAF_num
    
    harna <- is.na(data_MAF$A1)
    harna <- data_MAF$SNP[which(harna == TRUE)]
    
    cor1 <- which(data_MAF$effect_allele.exposure != data_MAF$A1)
    data_MAF$eaf.exposure <- data_MAF$MAF
    data_MAF$type <- "raw"
    data_MAF$eaf.exposure[cor1] <- EAF_num[cor1]
    data_MAF$type[cor1] <- "corrected"
    
    cor2 <- which(data_MAF$other_allele.exposure == data_MAF$A1)
    cor21 <- setdiff(cor2, cor1)
    cor12 <- setdiff(cor1, cor2)
    error <- c(cor12, cor21)
    data_MAF$eaf.exposure[error] <- NA
    data_MAF$type[error] <- "error"
    
    data_MAF <- list(
      data_MAF = data_MAF,
      cor1 = cor1,
      harna = harna,
      error = error
    )
    return(data_MAF)
  }
  
  corrected_eaf_out <- function(data_MAF) {
    effect <- data_MAF$effect_allele.outcome
    other <- data_MAF$other_allele.outcome
    A1 <- data_MAF$A1
    A2 <- data_MAF$A2
    MAF_num <- data_MAF$MAF
    EAF_num <- 1 - MAF_num
    
    harna <- is.na(data_MAF$A1)
    harna <- data_MAF$SNP[which(harna == TRUE)]
    
    cor1 <- which(data_MAF$effect_allele.outcome != data_MAF$A1)
    data_MAF$eaf.outcome <- data_MAF$MAF
    data_MAF$type <- "raw"
    data_MAF$eaf.outcome[cor1] <- EAF_num[cor1]
    data_MAF$type[cor1] <- "corrected"
    
    cor2 <- which(data_MAF$other_allele.outcome == data_MAF$A1)
    cor21 <- setdiff(cor2, cor1)
    cor12 <- setdiff(cor1, cor2)
    error <- c(cor12, cor21)
    data_MAF$eaf.outcome[error] <- NA
    data_MAF$type[error] <- "error"
    
    data_MAF <- list(
      data_MAF = data_MAF,
      cor1 = cor1,
      harna = harna,
      error = error
    )
    return(data_MAF)
  }
  
  if (type == "exposure" && (!"eaf.exposure" %in% names(dat) || any(is.na(dat$eaf.exposure)))) {
    if (!"eaf.exposure" %in% names(dat)) {
      dat$eaf.exposure <- NA_real_
    }
    r <- nrow(dat)
    MAF <- data.table::fread(file.path(path, "fileFrequency.frq"), header = TRUE)
    dat <- merge(dat, MAF, by.x = "SNP", by.y = "SNP", all.x = TRUE)
    dat <- corrected_eaf_expo(dat)
    cor1 <- dat$cor1
    harna <- dat$harna
    error <- dat$error
    dat <- dat$data_MAF
    
    print(paste0("A total of ", (r - length(harna) - length(error)), " SNPs were successfully matched to EAF values (", 
                 round((r - length(harna) - length(error)) / r * 100, 2), "%"))
    print(paste0("A total of ", length(cor1), " SNPs used the major allele; EAF was calculated as 1-MAF (", 
                 round(length(cor1) / (r - length(harna) - length(error)) * 100, 2), "%"))
    print(paste0("A total of ", length(harna), " SNPs were not found in 1000 Genomes (", 
                 round(length(harna) / r * 100, 2), "%"))
    print(paste0("A total of ", length(error), " SNPs had inconsistent effect and reference alleles; EAF was set to missing (", 
                 round(length(error) / r * 100, 2), "%"))
    return(dat)
  }
  
  if (type == "outcome" && (!"eaf.outcome" %in% names(dat) || any(is.na(dat$eaf.outcome)))) {
    if (!"eaf.outcome" %in% names(dat)) {
      dat$eaf.outcome <- NA_real_
    }
    r <- nrow(dat)
    MAF <- data.table::fread(file.path(path, "fileFrequency.frq"), header = TRUE)
    dat <- merge(dat, MAF, by.x = "SNP", by.y = "SNP", all.x = TRUE)
    dat <- corrected_eaf_out(dat)
    cor1 <- dat$cor1
    harna <- dat$harna
    error <- dat$error
    dat <- dat$data_MAF
    
    print(paste0("A total of ", (r - length(harna) - length(error)), " SNPs were successfully matched to EAF values (", 
                 round((r - length(harna) - length(error)) / r * 100, 2), "%"))
    print(paste0("A total of ", length(cor1), " SNPs used the major allele; EAF was calculated as 1-MAF (", 
                 round(length(cor1) / (r - length(harna) - length(error)) * 100, 2), "%"))
    print(paste0("A total of ", length(harna), " SNPs were not found in 1000 Genomes (", 
                 round(length(harna) / r * 100, 2), "%"))
    print(paste0("A total of ", length(error), " SNPs had inconsistent effect and reference alleles; EAF was set to missing (", 
                 round(length(error) / r * 100, 2), "%"))
    return(dat)
  } else {
    return(dat)
  }
}


#' 
#' @export
get_f <- function(dat, F_value = 10) {

  if (!"beta.exposure" %in% names(dat) || all(is.na(dat$beta.exposure))) {
    print("The data do not contain beta values; the F statistic cannot be calculated")
    return(dat)
  }
  
  if (!"se.exposure" %in% names(dat) || all(is.na(dat$se.exposure))) {
    print("The data do not contain standard errors; the F statistic cannot be calculated")
    return(dat)
  }

  has_eaf <- "eaf.exposure" %in% names(dat) && !all(is.na(dat$eaf.exposure))
  has_samplesize <- "samplesize.exposure" %in% names(dat) && !all(is.na(dat$samplesize.exposure))

  if (has_eaf && has_samplesize) {
    R2 <- (2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$beta.exposure^2)) /
      ((2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$beta.exposure^2)) +
         (2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$se.exposure^2) * dat$samplesize.exposure))

    F <- (dat$samplesize.exposure - 2) * R2 / (1 - R2)
    dat$R2 <- R2
    dat$F <- F
  } else {
    warning("Valid EAF or sample size is missing; using approximate F = (beta/se)^2 and setting R2 to NA")
    dat$R2 <- NA_real_
    dat$F <- (dat$beta.exposure / dat$se.exposure)^2
  }

  original_n <- nrow(dat)
  dat <- subset(dat, !is.na(F) & F > F_value)

  print(paste0("SNPs before filtering: ", original_n))
  print(paste0("Number of SNPs with F > ", F_value, ": ", nrow(dat)))
  if (nrow(dat) > 0) {
    print(paste0("F-statistic range: ", round(min(dat$F, na.rm = TRUE), 2), " - ", round(max(dat$F, na.rm = TRUE), 2)))
  }

  return(dat)
}


#' 
#' @export
add_MAF_from_EAF <- function(data) {
  if (!"eaf.exposure" %in% names(data) || all(is.na(data$eaf.exposure))) {
    warning("eaf.exposure is missing; MAF will be set to NA")
    data$MAF <- NA_real_
    return(data)
  }
  data$MAF <- ifelse(data$eaf.exposure > 0.5, 1 - data$eaf.exposure, data$eaf.exposure)
  return(data)
}


#' 
#' @export
add_MAF_from_EAF_outcome <- function(data) {
  if (!"eaf.outcome" %in% names(data)) {
    stop("The eaf.outcome column is absent")
  }
  data$MAF.outcome <- ifelse(data$eaf.outcome > 0.5, 1 - data$eaf.outcome, data$eaf.outcome)
  return(data)
}


#' 
#' @export
summarize_snp_counts <- function(dat) {
  if (is.null(dat) || nrow(dat) == 0) {
    warning("Input data are empty; SNP counts cannot be summarised")
    return(invisible(NULL))
  }
  if (!"exposure" %in% names(dat)) {
    warning("The exposure column is missing; SNP counts cannot be summarised")
    return(invisible(NULL))
  }

  snp_counts <- dat %>%
    dplyr::group_by(exposure) %>%
    dplyr::summarise(SNP_count = dplyr::n(), .groups = "drop")
  
  median_count <- median(snp_counts$SNP_count)
  range_count <- range(snp_counts$SNP_count)
  
  cat("\n=== SNP-count summary ===\n")
  cat("Median:", median_count, "\n")
  cat("Range:", range_count[1], "to", range_count[2], "\n")
  print(snp_counts)
  
  return(invisible(snp_counts))
}


#' 
#' @export
summarize_f_statistics <- function(dat) {
  if (!"F" %in% names(dat)) {
    warning("The F column is absent")
    return(invisible(NULL))
  }
  
  cat("\n=== F-statistic summary ===\n")
  print(summary(dat$F))
  
  if ("R2" %in% names(dat)) {
    cat("\n=== R-squared summary ===\n")
    print(summary(dat$R2))
  }
}


#' 
#' @export
rename_exposures <- function(dat, old_names, new_names) {
  if (length(old_names) != length(new_names)) {
    stop("old_names and new_names must have the same length")
  }
  
  for (i in seq_along(old_names)) {
    dat$exposure[dat$exposure == old_names[i]] <- new_names[i]
  }
  
  return(dat)
}


#' 
#' @export
add_exposure_outcome_id <- function(dat) {
  if (!all(c("id.exposure", "id.outcome") %in% names(dat))) {
    stop("The data must contain id.exposure and id.outcome columns")
  }
  dat$exposure_outcome <- paste(dat$id.exposure, dat$id.outcome, sep = "_")
  return(dat)
}


#' 
#' @export
safe_read_csv <- function(file_path, ...) {
  if (!file.exists(file_path)) {
    warning(paste0("File does not exist: ", file_path))
    return(NULL)
  }
  
  tryCatch({
    data.table::fread(file_path, ...)
  }, error = function(e) {
    warning(paste0("Failed to read file: ", file_path, "\nError: ", e$message))
    return(NULL)
  })
}


#' 
#' @export
safe_write_csv <- function(dat, file_path, ...) {
  if (is.null(file_path) || !is.character(file_path) || length(file_path) == 0) {
    stop("file_path must be a non-empty character vector")
  }
  
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message(paste0("Creating directory: ", dir_path))
  }
  
  tryCatch({
    data.table::fwrite(dat, file_path, ...)
    message(paste0("File written successfully: ", file_path))
  }, error = function(e) {
    warning(paste0("Failed to write file: ", file_path, "\nError: ", e$message))
  })
}

read_outcome_subset_for_snps <- function(file_path,
                                         target_snps,
                                         snp_col,
                                         chunk_size = 200000) {
  target_snps <- unique(target_snps[!is.na(target_snps) & target_snps != ""])
  if (length(target_snps) == 0) {
    stop("No target SNPs were provided for outcome extraction")
  }
  if (!file.exists(file_path)) {
    stop(paste0("Outcome file does not exist: ", file_path))
  }

  message(paste0("Fast outcome extraction for ", length(target_snps), " target SNPs"))

  con <- if (grepl("\\.gz$", file_path, ignore.case = TRUE)) {
    gzfile(file_path, open = "rt")
  } else {
    file(file_path, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  header <- readLines(con, n = 1, warn = FALSE)
  if (length(header) == 0) {
    stop(paste0("Outcome file is empty: ", file_path))
  }

  cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
  if (!snp_col %in% cols) {
    stop(paste0("Outcome SNP column not found: ", snp_col))
  }

  escaped_snps <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", target_snps, perl = TRUE)
  snp_pattern <- paste0("(^|[\t,;])(", paste(escaped_snps, collapse = "|"), ")([\t,;]|$)")

  matched_chunks <- list()
  n_read <- 0
  n_matched <- 0

  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (length(lines) == 0) {
      break
    }

    keep <- grepl(snp_pattern, lines, perl = TRUE)
    if (any(keep)) {
      matched_lines <- lines[keep]
      matched_chunks[[length(matched_chunks) + 1]] <- matched_lines
      n_matched <- n_matched + length(matched_lines)
    }

    n_read <- n_read + length(lines)
    if (n_read %% 5000000 < chunk_size) {
      message(paste0("  Scanned ", n_read, " outcome rows; matched ", n_matched))
    }
  }

  empty_dt <- data.table::as.data.table(
    stats::setNames(rep(list(character()), length(cols)), cols)
  )

  if (length(matched_chunks) == 0) {
    message("  No outcome rows matched target SNPs")
    return(empty_dt)
  }

  dt <- data.table::fread(
    text = paste(c(header, unlist(matched_chunks, use.names = FALSE)), collapse = "\n"),
    header = TRUE
  )

  if (nrow(dt) == 0) {
    return(empty_dt)
  }

  target_lookup <- new.env(hash = TRUE, parent = emptyenv())
  for (snp in target_snps) {
    target_lookup[[snp]] <- TRUE
  }

  expanded_rows <- vector("list", nrow(dt))
  row_count <- 0
  snp_values <- as.character(dt[[snp_col]])

  for (i in seq_len(nrow(dt))) {
    tokens <- trimws(strsplit(snp_values[i], "[,;]")[[1]])
    hits <- tokens[vapply(tokens, exists, logical(1), envir = target_lookup, inherits = FALSE)]
    if (length(hits) == 0) {
      next
    }

    for (hit in unique(hits)) {
      row_count <- row_count + 1
      row <- dt[i]
      row[[snp_col]] <- hit
      expanded_rows[[row_count]] <- row
    }
  }

  if (row_count == 0) {
    message("  Regex matched lines, but no exact SNP IDs remained after parsing")
    return(empty_dt)
  }

  result <- data.table::rbindlist(expanded_rows[seq_len(row_count)], fill = TRUE)
  message(paste0("  Extracted ", nrow(result), " exact outcome rows"))
  return(result)
}

# ==============================================================================
# ==============================================================================
