source("00_revised_config.R", encoding = "UTF-8")

cat("=== Revised differential filtering and volcano plots: POS mode ===\n")
res <- run_diff_pipeline("POS")
print(res)
cat("\nResults saved to: 05b_pos_diff_revised\n")

