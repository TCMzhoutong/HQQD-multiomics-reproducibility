source("00_revised_config.R", encoding = "UTF-8")

cat("=== Revised fold-change analysis: POS mode ===\n")
res <- run_fc_pipeline("POS")
print(res)
cat("\nResults saved to: 02b_FC_pos_revised\n")

