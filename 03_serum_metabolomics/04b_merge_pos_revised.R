source("00_revised_config.R", encoding = "UTF-8")

cat("=== Revised merge: POS mode ===\n")
res <- run_merge_pipeline("POS")
print(res)
cat("\nResults saved to: 04b_pos_merged_revised\n")

