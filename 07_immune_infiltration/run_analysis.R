# Standalone CIBERSORT workflow for bacterial pneumonia immune infiltration.

cat("\n============================================================\n")
cat("Running standalone CIBERSORT immune infiltration analysis\n")
cat("============================================================\n\n")

source("scripts/09_immune_infiltration_CIBERSORT.R")

cat("\nStandalone CIBERSORT workflow finished.\n")

quit(save = "no", status = 0)
