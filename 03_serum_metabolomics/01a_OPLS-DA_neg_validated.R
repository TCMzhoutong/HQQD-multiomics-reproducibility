source("00_revised_config.R", encoding = "UTF-8")

cat("=== Revised validated OPLS/PLS analysis: NEG mode ===\n")
res <- run_oplsda_pipeline("NEG", permI = 200)
print(res)
cat("\nResults saved to: 01a_OPLS-DA_neg_validated\n")

