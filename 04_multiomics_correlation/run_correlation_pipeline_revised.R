# Reproducible Windows entry point for steps 04-07.

locale_ok <- try(
  Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"),
  silent = TRUE
)

if (inherits(locale_ok, "try-error") || is.na(locale_ok)) {
  stop("A UTF-8 LC_CTYPE locale is required to run the correlation pipeline")
}

source("04_differential_species_metabolites.R", encoding = "UTF-8")
source("05_correlation_analysis.R", encoding = "UTF-8")
source("06_mantel_analysis.R", encoding = "UTF-8")
source("07_prioritize_metabolites_for_network_pharmacology.R", encoding = "UTF-8")
