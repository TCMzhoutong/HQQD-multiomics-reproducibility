source("00_revised_config.R", encoding = "UTF-8")

cat("=== Revised merge: NEG mode ===\n")
res <- run_merge_pipeline("NEG")
print(res)
cat("\nResults saved to: 04a_neg_merged_revised\n")

