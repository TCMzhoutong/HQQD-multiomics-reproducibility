# ============================================================================
# Step 04: Differential Species and Metabolites Analysis
#
# Current inputs:
# 1. Differential taxa list:
#    omics data/Phylum_Genus.csv, using columns Level and Taxon.
# 2. Taxa abundance tables:
#    omics data/taxonomy.genus.relabundance.csv
#    omics data/taxonomy.phylum.relabundance.csv
# 3. Differential metabolites list and expression matrices from module 03.
#
# Outputs are kept compatible with steps 05 and 06.
# ============================================================================

rm(list = ls())

required_packages <- c("dplyr", "stringr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
})

output_dir <- "04_differential_species_metabolites"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory: ", output_dir, "\n", sep = "")
}

read_csv_flexible <- function(path) {
  encodings <- c("UTF-8", "GB18030", "GBK")
  last_error <- NULL
  for (enc in encodings) {
    dat <- tryCatch(
      suppressWarnings(
        readr::read_csv(
          path,
          locale = readr::locale(encoding = enc),
          show_col_types = FALSE,
          name_repair = "minimal"
        ) %>%
          as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)
      ),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(dat) && ncol(dat) > 1) return(dat)
  }
  stop("Could not read CSV file: ", path, "; last error: ", conditionMessage(last_error))
}

write_csv_utf8 <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

pick_first_nonempty <- function(df, cols) {
  available <- intersect(cols, colnames(df))
  if (length(available) == 0) return(rep(NA_character_, nrow(df)))
  out <- as.character(df[[available[1]]])
  if (length(available) > 1) {
    for (col in available[-1]) {
      v <- as.character(df[[col]])
      idx <- is.na(out) | out == ""
      out[idx] <- v[idx]
    }
  }
  out
}

# ============================================================================
# Part 1: Differential Taxa
# ============================================================================

cat("\n=== Part 1: Differential Taxa Analysis ===\n\n")

taxa_list_file <- "omics data/Phylum_Genus.csv"
if (!file.exists(taxa_list_file)) {
  stop("Differential taxa file not found: ", taxa_list_file)
}

taxa_diff_source <- read_csv_flexible(taxa_list_file) %>%
  select(any_of(c("Level", "Taxon"))) %>%
  filter(!is.na(Level), !is.na(Taxon), Level != "", Taxon != "") %>%
  mutate(Level = str_trim(Level), Taxon = str_trim(Taxon)) %>%
  distinct(Level, Taxon)

cat("Differential taxa list loaded: ", nrow(taxa_diff_source), " rows\n", sep = "")
print(table(taxa_diff_source$Level))

taxa_levels <- list(
  list(level = "Genus", prefix = "g", abundance_file = "omics data/taxonomy.genus.relabundance.csv"),
  list(level = "Phylum", prefix = "p", abundance_file = "omics data/taxonomy.phylum.relabundance.csv")
)

extract_taxa_abundance <- function(abundance_df, taxa_names, taxa_col_name,
                                   model_cols, hqqd_cols, k_cols, prefix) {
  filtered <- abundance_df %>%
    filter(.data[[taxa_col_name]] %in% taxa_names)

  missing_taxa <- setdiff(taxa_names, filtered[[taxa_col_name]])
  if (length(missing_taxa) > 0) {
    cat("  Warning: missing taxa in abundance table: ", paste(missing_taxa, collapse = ", "), "\n", sep = "")
  }

  if (nrow(filtered) == 0) return(NULL)

  result <- filtered %>%
    select(all_of(c(taxa_col_name, model_cols, hqqd_cols, k_cols)))

  colnames(result)[1] <- "Taxon"
  result$Taxon <- paste0(prefix, "_", result$Taxon)
  result
}

all_diff_taxa_list <- list()
all_diff_taxa_summary <- list()

for (cfg in taxa_levels) {
  lvl <- cfg$level
  pfx <- cfg$prefix
  af <- cfg$abundance_file

  if (!file.exists(af)) {
    cat("Skipping ", lvl, " level: abundance file not found (", af, ")\n", sep = "")
    next
  }

  diff_df <- taxa_diff_source %>%
    filter(Level == lvl) %>%
    select(Level, Taxon)

  cat("Processing ", lvl, " level: ", nrow(diff_df), " taxa listed\n", sep = "")
  if (nrow(diff_df) == 0) next

  abund_df <- read_csv_flexible(af)
  taxa_col_name <- colnames(abund_df)[1]
  all_cols <- colnames(abund_df)
  model_cols <- all_cols[str_detect(all_cols, "^M\\d+")]
  hqqd_cols <- all_cols[str_detect(all_cols, "^Z\\d+")]
  k_cols <- all_cols[str_detect(all_cols, "^K\\d+")]

  cat("  Abundance table: ", nrow(abund_df), " rows; ",
      length(model_cols), " Model, ", length(hqqd_cols), " HQQD, ",
      length(k_cols), " Control samples\n", sep = "")

  abund_filtered <- extract_taxa_abundance(
    abund_df, diff_df$Taxon, taxa_col_name, model_cols, hqqd_cols, k_cols, pfx
  )

  if (!is.null(abund_filtered)) {
    all_diff_taxa_list[[lvl]] <- abund_filtered
  }

  all_diff_taxa_summary[[lvl]] <- diff_df %>%
    mutate(Taxon = paste0(pfx, "_", Taxon))
}

combined_taxa <- if (length(all_diff_taxa_list) > 0) bind_rows(all_diff_taxa_list) else data.frame()
all_diff_taxa <- if (length(all_diff_taxa_summary) > 0) bind_rows(all_diff_taxa_summary) else data.frame()

cat("Combined taxa extracted: ", nrow(combined_taxa), "\n", sep = "")

if (nrow(combined_taxa) > 0) {
  write_csv_utf8(combined_taxa, file.path(output_dir, "differential_taxa_abundance_Model_vs_HQQD.csv"))
  write_csv_utf8(all_diff_taxa, file.path(output_dir, "differential_taxa_list.csv"))
  cat("Saved differential taxa outputs.\n")
} else {
  cat("Warning: No taxa were extracted.\n")
}

# ============================================================================
# Part 2: Differential Metabolites
# ============================================================================

cat("\n=== Part 2: Differential Metabolites Analysis ===\n\n")

metabolite_diff_file <- file.path(
  "..", "03_serum_metabolomics", "results_reference",
  "HQQD_reverse_candidates_32.csv"
)
if (!file.exists(metabolite_diff_file)) {
  stop("Differential metabolite file not found: ", metabolite_diff_file)
}

metabolite_diff <- read_csv_flexible(metabolite_diff_file)
if (!"ID" %in% colnames(metabolite_diff)) {
  stop("Differential metabolite table must contain an ID column.")
}

metabolite_filtered <- data.frame(
  ID = metabolite_diff$ID,
  Name = pick_first_nonempty(metabolite_diff, c("Name", "name", "MH_Name", "CM_Name")),
  Log2FC = pick_first_nonempty(metabolite_diff, c("Log2FC", "MH_Log2FC", "CM_Log2FC")),
  Ion_Mode = pick_first_nonempty(metabolite_diff, c("Ion_Mode", "ion_mode")),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(ID), ID != "") %>%
  distinct(ID, .keep_all = TRUE)

cat("Differential metabolites listed: ", nrow(metabolite_filtered), "\n", sep = "")
cat("By ion mode:\n")
print(table(metabolite_filtered$Ion_Mode, useNA = "ifany"))

extract_metabolite_abundance <- function(path, ion_mode_label, metabolite_filtered) {
  if (!file.exists(path)) {
    cat("Warning: metabolite expression file not found: ", path, "\n", sep = "")
    return(NULL)
  }

  metabolite_exp <- read_csv_flexible(path)
  metab_cols <- colnames(metabolite_exp)
  id_col <- if ("ID" %in% metab_cols) "ID" else metab_cols[1]
  name_col <- if ("name" %in% metab_cols) "name" else if ("Name" %in% metab_cols) "Name" else metab_cols[2]

  model_cols <- metab_cols[str_detect(metab_cols, "^M\\d+")]
  hqqd_cols <- metab_cols[str_detect(metab_cols, "^Z\\d+")]
  k_cols <- metab_cols[str_detect(metab_cols, "^K\\d+")]

  cat("  ", ion_mode_label, " matrix: ", nrow(metabolite_exp), " rows; ",
      length(model_cols), " Model, ", length(hqqd_cols), " HQQD, ",
      length(k_cols), " Control samples\n", sep = "")

  extracted <- metabolite_exp %>%
    filter(.data[[id_col]] %in% metabolite_filtered$ID) %>%
    select(all_of(c(id_col, name_col, model_cols, hqqd_cols, k_cols))) %>%
    rename(ID = all_of(id_col), Matrix_Name = all_of(name_col)) %>%
    left_join(metabolite_filtered, by = "ID") %>%
    mutate(
      Name = ifelse(is.na(Name) | Name == "", Matrix_Name, Name),
      Ion_Mode = ifelse(is.na(Ion_Mode) | Ion_Mode == "", ion_mode_label, Ion_Mode)
    ) %>%
    select(ID, Name, Log2FC, Ion_Mode, all_of(model_cols), all_of(hqqd_cols), all_of(k_cols))

  missing_ids <- setdiff(
    metabolite_filtered$ID[str_starts(metabolite_filtered$ID, paste0(tolower(ion_mode_label), "_"))],
    extracted$ID
  )
  if (length(missing_ids) > 0) {
    cat("    Warning: ", length(missing_ids), " ", ion_mode_label, " metabolite IDs not found in expression matrix\n", sep = "")
  }

  extracted
}

metabolite_exp_files <- c(
  NEG = file.path("..", "03_serum_metabolomics", "metabolites_exp_neg.csv"),
  POS = file.path("..", "03_serum_metabolomics", "metabolites_exp_pos.csv")
)

metabolite_abundance_list <- lapply(names(metabolite_exp_files), function(mode) {
  extract_metabolite_abundance(metabolite_exp_files[[mode]], mode, metabolite_filtered)
})
metabolite_abundance_list <- metabolite_abundance_list[!vapply(metabolite_abundance_list, is.null, logical(1))]

metabolite_abundance <- if (length(metabolite_abundance_list) > 0) {
  bind_rows(metabolite_abundance_list) %>%
    distinct(ID, .keep_all = TRUE)
} else {
  data.frame()
}

cat("Extracted differential metabolites: ", nrow(metabolite_abundance), "\n", sep = "")

if (nrow(metabolite_abundance) > 0) {
  write_csv_utf8(metabolite_abundance, file.path(output_dir, "differential_metabolites_abundance_Model_vs_HQQD.csv"))
  write_csv_utf8(metabolite_filtered, file.path(output_dir, "differential_metabolites_list.csv"))
  cat("Saved differential metabolite outputs.\n")
} else {
  cat("Warning: No metabolites were extracted.\n")
}

# ============================================================================
# Summary
# ============================================================================

cat("\n=== Analysis Summary ===\n")
cat("Taxa extracted: ", nrow(combined_taxa), "\n", sep = "")
for (lvl in names(all_diff_taxa_summary)) {
  cat("  - ", lvl, " level listed: ", nrow(all_diff_taxa_summary[[lvl]]), "\n", sep = "")
}
cat("Metabolites extracted: ", nrow(metabolite_abundance), "\n", sep = "")
cat("All results saved to: ", output_dir, "\n", sep = "")
cat("\n=== Analysis Complete ===\n")
