source("00_revised_config.R", encoding = "UTF-8")

cat("=== Revised Welch t-test analysis: NEG mode ===\n")
res <- run_welch_pipeline("NEG")
print(res)
cat("\nResults saved to: 03a_p_neg_welch\n")

