source("00_revised_config.R", encoding = "UTF-8")

cat("=== Revised fold-change analysis: NEG mode ===\n")
res <- run_fc_pipeline("NEG")
print(res)
cat("\nResults saved to: 02a_FC_neg_revised\n")

