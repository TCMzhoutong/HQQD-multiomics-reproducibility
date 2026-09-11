# ============================================================================
# ============================================================================

library(GEOquery)

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

cat("=== Load data ===\n")

load(file.path(raw_data_dir, "GSE40012_data.RData"))
load(file.path(raw_data_dir, "GSE20346_data.RData"))

cat("\n=== GSE40012 sample-group details ===\n\n")

pdata_40012 <- pData(gse40012_data)

cat("Sample source:\n")
print(table(pdata_40012$source_name_ch1))

char_cols <- grep("characteristics", colnames(pdata_40012), value = TRUE)
cat("\nCharacteristics columns:\n")
for (col in char_cols) {
  cat(paste0("\n", col, ":\n"))
  print(table(pdata_40012[[col]]))
}

cat("\n\n=== GSE20346 sample-group details ===\n\n")

pdata_20346 <- pData(gse20346_data)

cat("Sample source:\n")
print(table(pdata_20346$source_name_ch1))

char_cols <- grep("characteristics", colnames(pdata_20346), value = TRUE)
cat("\nCharacteristics columns:\n")
for (col in char_cols) {
  cat(paste0("\n", col, ":\n"))
  print(table(pdata_20346[[col]]))
}

cat("\n=== Select GSE40012 samples ===\n")

day1_idx_40012 <- which(pdata_40012$characteristics_ch1.2 == "day: 1")
cat("Day 1 Samples:", length(day1_idx_40012), "\n")

sample_type_40012 <- pdata_40012$characteristics_ch1.1[day1_idx_40012]

bacterial_idx_40012 <- day1_idx_40012[sample_type_40012 == "sample type: bacterial pneumonia"]
cat("Bacterial pneumonia (day 1):", length(bacterial_idx_40012), "\n")

healthy_idx_40012 <- day1_idx_40012[sample_type_40012 == "sample type: healthy control"]
cat("Healthy controls (day 1):", length(healthy_idx_40012), "\n")

selected_idx_40012 <- c(bacterial_idx_40012, healthy_idx_40012)
cat("\nGSE40012 selected samples:", length(selected_idx_40012), "\n")

cat("\n=== Select GSE20346 samples ===\n")

source_20346 <- pdata_20346$source_name_ch1

bacterial_idx_20346 <- which(
  grepl("bacterial pneumonia", source_20346, ignore.case = TRUE) &
  grepl("day 1", source_20346, ignore.case = TRUE)
)
cat("Bacterial pneumonia (day 1):", length(bacterial_idx_20346), "\n")

healthy_idx_20346 <- which(grepl("Vaccine_baseline", source_20346, ignore.case = TRUE))
cat("Healthy controls (Vaccine_baseline):", length(healthy_idx_20346), "\n")

selected_idx_20346 <- c(bacterial_idx_20346, healthy_idx_20346)
cat("\nGSE20346 selected samples:", length(selected_idx_20346), "\n")

cat("\n=== Extract expression matrices ===\n")

# GSE40012
if (length(selected_idx_40012) > 0) {
  expr_40012_selected <- exprs(gse40012_data)[, selected_idx_40012]
  pdata_40012_selected <- pdata_40012[selected_idx_40012, ]
  
  group_40012_selected <- ifelse(
    pdata_40012_selected$characteristics_ch1.1 == "sample type: healthy control",
    "Control",
    "Bacterial_Pneumonia"
  )
  
  cat("GSE40012 selected-expression-matrix dimensions:", dim(expr_40012_selected), "\n")
  cat("  - Controls:", sum(group_40012_selected == "Control"), "\n")
  cat("  - Bacterial-pneumonia samples:", sum(group_40012_selected == "Bacterial_Pneumonia"), "\n")
} else {
  stop("Error: no GSE40012 samples were selected.")
}

# GSE20346
if (length(selected_idx_20346) > 0) {
  expr_20346_selected <- exprs(gse20346_data)[, selected_idx_20346]
  pdata_20346_selected <- pdata_20346[selected_idx_20346, ]
  
  group_20346_selected <- ifelse(
    grepl("Vaccine_baseline", pdata_20346_selected$source_name_ch1, ignore.case = TRUE),
    "Control",
    "Bacterial_Pneumonia"
  )
  
  cat("\nGSE20346 selected-expression-matrix dimensions:", dim(expr_20346_selected), "\n")
  cat("  - Controls:", sum(group_20346_selected == "Control"), "\n")
  cat("  - Bacterial-pneumonia samples:", sum(group_20346_selected == "Bacterial_Pneumonia"), "\n")
} else {
  stop("Error: no GSE20346 samples were selected.")
}

cat("\n=== Save results ===\n")

save(expr_40012_selected, pdata_40012_selected, group_40012_selected,
     file = file.path(processed_data_dir, "GSE40012_selected.RData"))

save(expr_20346_selected, pdata_20346_selected, group_20346_selected,
     file = file.path(processed_data_dir, "GSE20346_selected.RData"))

cat("Saved to 02_processed_data.\n")

cat("\n
==========================================================
                  Sample selection complete
==========================================================

GSE40012:
  - Selected samples: ", length(selected_idx_40012), "
  - Controls: ", sum(group_40012_selected == "Control"), "
  - Bacterial-pneumonia samples: ", sum(group_40012_selected == "Bacterial_Pneumonia"), "

GSE20346:
  - Selected samples: ", length(selected_idx_20346), "
  - Controls: ", sum(group_20346_selected == "Control"), "
  - Bacterial-pneumonia samples: ", sum(group_20346_selected == "Bacterial_Pneumonia"), "

Next step: source('scripts/03_probe_to_gene.R')
==========================================================
", sep = "")
