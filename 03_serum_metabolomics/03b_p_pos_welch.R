source("00_revised_config.R", encoding = "UTF-8")

cat("=== Revised Welch t-test analysis: POS mode ===\n")
res <- run_welch_pipeline("POS")
print(res)
cat("\nResults saved to: 03b_p_pos_welch\n")

