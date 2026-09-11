# Run the Mantel analysis under an explicit UTF-8 locale on Windows.
# This wrapper prevents Greek characters in disease-indicator column names
# from being misread when R starts in the C locale.

locale_ok <- try(
  Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"),
  silent = TRUE
)

if (inherits(locale_ok, "try-error") || is.na(locale_ok)) {
  stop("A UTF-8 LC_CTYPE locale is required to run 06_mantel_analysis.R")
}

source("06_mantel_analysis.R", encoding = "UTF-8")
