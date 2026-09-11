get_gse252663_config <- function() {
  root_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)

  list(
    root_dir = root_dir,
    raw_dir = file.path(root_dir, "data", "GSE252663", "raw"),
    tar_file = file.path(root_dir, "data", "GSE252663", "GSE252663_RAW.tar"),
    results_dir = file.path(root_dir, "results", "GSE252663"),
    intermediate_dir = file.path(root_dir, "results", "GSE252663", "intermediate"),
    figure_dir = file.path(root_dir, "results", "GSE252663", "source_tables"),
    hub_gene_file = file.path(root_dir, "01_raw_data", "03_hub_genes_intersection_all.csv"),
    aucell_score_name = "CoreTargets_AUC",

    sample_info = data.frame(
      sample_name = c("KO-KP1", "KO-KP2", "WT-KP1", "WT-KP2", "WT-Control", "KO-Control"),
      geo_accession = c("GSM8004969", "GSM8004970", "GSM8004971", "GSM8004972", "GSM8004973", "GSM8004974"),
      file_prefix = c(
        "GSM8004969_LEE1723A1",
        "GSM8004970_LEE1723A2",
        "GSM8004971_LEE1723A3",
        "GSM8004972_LEE1723A4",
        "GSM8004973_LEE1723A5",
        "GSM8004974_LEE1723A6"
      ),
      genotype = c("KO", "KO", "WT", "WT", "WT", "KO"),
      infection_status = c("KP", "KP", "KP", "KP", "Control", "Control"),
      condition = c("KO_KP", "KO_KP", "WT_KP", "WT_KP", "WT_Control", "KO_Control"),
      replicate = c("1", "2", "1", "2", "1", "1"),
      stringsAsFactors = FALSE
    ),

    qc = list(
      min_features = 200,
      max_features = 7500,
      max_percent_mt = 20
    ),

    seurat = list(
      nfeatures = 2500,
      npcs = 30,
      dims_use = 1:30,
      cluster_resolution = 0.6,
      tcell_cluster_resolution = 0.8
    ),

    group_order = c("Control", "KP"),
    condition_order = c("WT_Control", "KO_Control", "WT_KP", "KO_KP")
  )
}
