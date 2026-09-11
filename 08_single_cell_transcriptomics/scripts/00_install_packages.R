required_pkgs <- c(
  "Seurat",
  "harmony",
  "AUCell",
  "ggplot2",
  "patchwork",
  "ggpubr",
  "dplyr",
  "tibble",
  "tidyr",
  "Matrix",
  "cowplot"
)

installed <- rownames(installed.packages())
need <- setdiff(required_pkgs, installed)

if (length(need) > 0) {
  install.packages(need, dependencies = TRUE)
}

if (!"BiocManager" %in% rownames(installed.packages())) {
  install.packages("BiocManager")
}

bioc_pkgs <- c("AUCell", "SingleR", "celldex")
for (pkg in bioc_pkgs) {
  if (!pkg %in% rownames(installed.packages())) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

message("Package check completed.")
