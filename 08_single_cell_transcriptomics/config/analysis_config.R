get_analysis_config <- function() {
  root_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
  output_suffix <- Sys.getenv("ANALYSIS_OUTPUT_SUFFIX", "")

  hub_gene_file_env       <- Sys.getenv("ANALYSIS_HUB_GENE_FILE", "")
  hub_gene_column_env     <- Sys.getenv("ANALYSIS_HUB_GENE_COLUMN", "")
  hub_filter_column_env   <- Sys.getenv("ANALYSIS_HUB_FILTER_COLUMN", "")
  hub_filter_value_env    <- Sys.getenv("ANALYSIS_HUB_FILTER_VALUE", "")
  hub_split_delim_env     <- Sys.getenv("ANALYSIS_HUB_GENE_SPLIT", "")
  target_name_env         <- Sys.getenv("ANALYSIS_TARGET_NAME", "")
  hub_name_prefix_env     <- Sys.getenv("ANALYSIS_HUB_NAME_PREFIX", "")

  hub_gene_file <- if (nzchar(hub_gene_file_env)) {
    file.path(root_dir, sub("^[./\\\\]+", "", hub_gene_file_env))
  } else {
    file.path(root_dir, "01_raw_data", "03_hub_genes_intersection_all.csv")
  }

  hub_gene_column <- if (nzchar(hub_gene_column_env)) hub_gene_column_env else "name"
  hub_filter_column <- if (nzchar(hub_filter_column_env)) hub_filter_column_env else NULL
  hub_filter_value  <- if (nzchar(hub_filter_value_env)) hub_filter_value_env else NULL
  hub_gene_split    <- if (nzchar(hub_split_delim_env)) hub_split_delim_env else NULL
  target_name <- if (nzchar(target_name_env)) target_name_env else "CoreTargets"
  hub_name_prefix <- if (nzchar(hub_name_prefix_env)) hub_name_prefix_env else target_name

  list(
    root_dir = root_dir,
    raw_data_dir = file.path(root_dir, "01_raw_data"),
    results_dir = file.path(root_dir, "results"),
    step1_dir = file.path(root_dir, "results", "step1"),
    step2_dir = file.path(root_dir, "results", paste0("step2", output_suffix)),
    step3_dir = file.path(root_dir, "results", paste0("step3", output_suffix)),
    step4_dir = file.path(root_dir, "results", paste0("step4", output_suffix)),

    sample_info = data.frame(
      sample_name = c(
        "Moderate_1_CRKP_Pneumonia",
        "Moderate_2_CRKP_Pneumonia",
        "Moderate_2_Recovery",
        "Severe_CRKP_Pneumonia",
        "Severe_Recovery"
      ),
      geo_accession = c(
        "GSM8190573",
        "GSM8190576",
        "GSM8190579",
        "GSM8190585",
        "GSM8190582"
      ),
      patient_id = c("Patient_1", "Patient_2", "Patient_2", "Patient_3", "Patient_3"),
      severity = c("Moderate", "Moderate", "Moderate", "Severe", "Severe"),
      phase = c("Acute", "Acute", "Convalescent", "Acute", "Convalescent"),
      disease_state = c("Moderate", "Moderate", "Convalescent", "Severe", "Convalescent"),
      disease_group = c(
        "Moderate",
        "Moderate",
        "Moderate_Convalescent",
        "Severe",
        "Severe_Convalescent"
      ),
      stringsAsFactors = FALSE
    ),

    qc = list(
      min_features = 200,
      max_features = 6500,
      max_percent_mt = 20
    ),

    seurat = list(
      nfeatures = 2000,
      npcs = 30,
      cluster_resolution = 0.6,
      tcell_cluster_resolution = 1.5,
      dims_use = 1:30
    ),

    markers = list(
      tcell_subtype = list()
    ),

    hub_gene_file = hub_gene_file,
    hub_gene_column = hub_gene_column,
    hub_filter_column = hub_filter_column,
    hub_filter_value = hub_filter_value,
    hub_gene_split = hub_gene_split,
    target_name = target_name,
    hub_name_prefix = hub_name_prefix,
    aucell_score_name = "CoreTargets_AUC",

    comparison = list(
      focus_celltypes = NULL,
      group_var = "Disease_State",
      group_order = c("Moderate", "Severe", "Convalescent"),
      t_exclude_labels = c("Non_T", "None T", "Unassigned", "other CD8 T", "CD8 quiescence")
    )
  )
}
