source("00_revised_config.R", encoding = "UTF-8")

cat("=== Revised NEG/POS differential merge ===\n")
res <- run_merge_diff_pipeline()
print(res)
cat("\nResults saved to: 06_merge_diff_revised\n")

