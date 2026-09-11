source("00_revised_config.R", encoding = "UTF-8")

cat("=== Revised validated OPLS/PLS analysis: POS mode ===\n")
res <- run_oplsda_pipeline("POS", permI = 200)
print(res)
cat("\nResults saved to: 01b_OPLS-DA_pos_validated\n")

