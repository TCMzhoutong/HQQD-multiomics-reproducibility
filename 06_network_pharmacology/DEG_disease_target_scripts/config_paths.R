# ============================================================================
# ============================================================================

current_dir <- getwd()

if (basename(current_dir) %in% c("01_data_download", "02_batch_correction", 
                                   "03_DEG_analysis", "04_visualization")) {
  PROJECT_ROOT <- dirname(current_dir)
  cat("Subdirectory detected; project root:", PROJECT_ROOT, "\n")
} else {
  PROJECT_ROOT <- current_dir
  cat("Project root:", PROJECT_ROOT, "\n")
}

DATA_DOWNLOAD_DIR <- file.path(PROJECT_ROOT, "01_data_download")
RAW_DATA_DIR <- file.path(DATA_DOWNLOAD_DIR, "raw_data")
PROCESSED_DATA_DIR <- file.path(DATA_DOWNLOAD_DIR, "processed_data")

BATCH_CORRECTION_DIR <- file.path(PROJECT_ROOT, "02_batch_correction")
DEG_ANALYSIS_DIR <- file.path(PROJECT_ROOT, "03_DEG_analysis")
VISUALIZATION_DIR <- file.path(PROJECT_ROOT, "04_visualization")

create_directories <- function() {
  dirs <- c(RAW_DATA_DIR, PROCESSED_DATA_DIR, 
            BATCH_CORRECTION_DIR, DEG_ANALYSIS_DIR, VISUALIZATION_DIR)
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Creating directory:", dir, "\n")
    }
  }
}

show_config <- function() {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════╗\n")
  cat("║                    Project path configuration                                   ║\n")
  cat("╚══════════════════════════════════════════════════════════════════╝\n")
  cat("Project root:      ", PROJECT_ROOT, "\n")
  cat("Data-download directory:    ", DATA_DOWNLOAD_DIR, "\n")
  cat("  - Raw data:    ", RAW_DATA_DIR, "\n")
  cat("  - Processed data:    ", PROCESSED_DATA_DIR, "\n")
  cat("Batch-correction directory:    ", BATCH_CORRECTION_DIR, "\n")
  cat("Differential-expression directory:    ", DEG_ANALYSIS_DIR, "\n")
  cat("Visualisation directory:      ", VISUALIZATION_DIR, "\n")
  cat("\n")
}

cat("Path configuration loaded\n")
