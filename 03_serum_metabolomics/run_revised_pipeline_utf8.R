# Reproducible Windows entry point for the revised metabolomics pipeline.

locale_ok <- try(
  Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"),
  silent = TRUE
)

if (inherits(locale_ok, "try-error") || is.na(locale_ok)) {
  stop("A UTF-8 LC_CTYPE locale is required to run the revised pipeline")
}

source("run_revised_pipeline.R", encoding = "UTF-8")
