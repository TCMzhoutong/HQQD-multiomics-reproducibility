source("00_revised_config.R", encoding = "UTF-8")

cat("=== Revised differential filtering and volcano plots: NEG mode ===\n")
res <- run_diff_pipeline("NEG")
print(res)
cat("\nResults saved to: 05a_neg_diff_revised\n")

