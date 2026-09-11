# ============================================================================
# ============================================================================

library(GEOquery)

cat("=== Set data paths ===\n\n")

script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0 || script_dir == "") {
  project_root <- getwd()
  if (basename(project_root) == "scripts") {
    project_root <- dirname(project_root)
  }
} else {
  project_root <- dirname(script_dir)
}

raw_data_dir <- file.path(project_root, "01_raw_data")
processed_data_dir <- file.path(project_root, "02_processed_data")

cat("Project root:", project_root, "\n")
cat("Raw-data directory:", raw_data_dir, "\n")
cat("Processed-data directory:", processed_data_dir, "\n\n")

if (!dir.exists(processed_data_dir)) {
  dir.create(processed_data_dir, recursive = TRUE)
}

cat("=== Load GSE40012 data ===\n")

file_40012 <- list.files(raw_data_dir, 
                         pattern = "GSE40012.*\\.txt\\.gz$", 
                         full.names = TRUE)

if (length(file_40012) == 0) {
  stop("Error: GSE40012 data file not found.
Place GSE40012_series_matrix.txt.gz in 01_raw_data.")
}

cat("Loading file:", basename(file_40012[1]), "\n")

gse40012_data <- getGEO(filename = file_40012[1], getGPL = FALSE)

cat("GSE40012data dimensions:", dim(exprs(gse40012_data)), "\n")
cat("- Probes:", nrow(exprs(gse40012_data)), "\n")
cat("- Samples:", ncol(exprs(gse40012_data)), "\n\n")

cat("=== Load GSE20346 data ===\n")

file_20346 <- list.files(raw_data_dir, 
                         pattern = "GSE20346.*\\.txt\\.gz$", 
                         full.names = TRUE)

if (length(file_20346) == 0) {
  stop("Error: GSE20346 data file not found.
Place GSE20346_series_matrix.txt.gz in 01_raw_data.")
}

cat("Loading file:", basename(file_20346[1]), "\n")

gse20346_data <- getGEO(filename = file_20346[1], getGPL = FALSE)

cat("GSE20346data dimensions:", dim(exprs(gse20346_data)), "\n")
cat("- Probes:", nrow(exprs(gse20346_data)), "\n")
cat("- Samples:", ncol(exprs(gse20346_data)), "\n\n")

cat("=== Load GPL6947 platform annotation ===\n")

gpl_file_txt <- list.files(raw_data_dir, 
                           pattern = "GPL6947.*\\.txt$", 
                           full.names = TRUE)

if (length(gpl_file_txt) == 0) {
  stop("Error: GPL6947 annotation file not found.
Download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6947
Download the GPL6947 annotation using 'Download full table...'
Place the file in 01_raw_data.")
}

cat("Loading GPL annotation:", basename(gpl_file_txt[1]), "\n")

gpl_table <- read.delim(
  gpl_file_txt[1],
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE,
  quote = ""
)

cat("GPL6947 annotation-table dimensions:", dim(gpl_table), "\n")
cat("Column names:\n")
print(colnames(gpl_table))

if (!"Symbol" %in% colnames(gpl_table)) {
  stop("Error: the Symbol column was not found in the GPL annotation.")
}

cat("\nGene-symbol examples (first 10 non-empty values):\n")
non_empty <- gpl_table$Symbol[gpl_table$Symbol != "" & !is.na(gpl_table$Symbol)]
print(head(non_empty, 10))

cat("\n=== Merge GPL annotation with expression data ===\n")

rownames(gpl_table) <- gpl_table$ID

probe_ids_40012 <- featureNames(gse40012_data)
probe_ids_20346 <- featureNames(gse20346_data)

matched_40012 <- probe_ids_40012 %in% rownames(gpl_table)
matched_20346 <- probe_ids_20346 %in% rownames(gpl_table)

cat("GSE40012 matched probes:", sum(matched_40012), "/", length(probe_ids_40012), "\n")
cat("GSE20346 matched probes:", sum(matched_20346), "/", length(probe_ids_20346), "\n")

fdata_40012 <- gpl_table[probe_ids_40012, ]
fdata_20346 <- gpl_table[probe_ids_20346, ]

fData(gse40012_data) <- fdata_40012
fData(gse20346_data) <- fdata_20346

cat("\nMerge complete.\n")

cat("\n=== Save data ===\n")

save(gse40012_data, file = file.path(raw_data_dir, "GSE40012_data.RData"))
save(gse20346_data, file = file.path(raw_data_dir, "GSE20346_data.RData"))

cat("Data saved to 01_raw_data.\n")

cat("\n
==========================================================
                    Data loading complete
==========================================================

GSE40012:
  - Probes: ", nrow(exprs(gse40012_data)), "
  - Samples: ", ncol(exprs(gse40012_data)), "
  - Annotation columns: ", ncol(fData(gse40012_data)), "

GSE20346:
  - Probes: ", nrow(exprs(gse20346_data)), "
  - Samples: ", ncol(exprs(gse20346_data)), "
  - Annotation columns: ", ncol(fData(gse20346_data)), "

Next step: source('scripts/02_sample_selection.R')
==========================================================
", sep = "")
