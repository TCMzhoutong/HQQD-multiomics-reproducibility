source("00_revised_config.R", encoding = "UTF-8")

cat("=== Revised HQQD reversal candidate analysis ===\n")
res <- run_hqqd_reverse_pipeline()
print(res)
cat("\nResults saved to: 06_HQQD_reverse_revised\n")

