# ==============================================================================
# ==============================================================================
# ==============================================================================

# ==============================================================================
# ==============================================================================

get_exposure_source_config <- function(source_name) {
  
  configs <- list(
    
    chen = list(
      name = "Chen Plasma Metabolites",
      file_path = DATA_FILES$chen,
      format_params = list(
        snp_col = "variant_id",
        beta_col = "beta",
        se_col = "standard_error",
        eaf_col = "effect_allele_frequency",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        pval_col = "p_value",
        samplesize_col = "SampleSize",
        phenotype_col = "reportedTrait",
        id_col = "GCST",
        chr_col = "chromosome",
        pos_col = "base_pair_location",
        header = TRUE,
        min_pval = 1e-200,
        log_pval = FALSE
      ),
      need_eaf = FALSE
    ),

    # ---- METSIM metabolites ----
    metsim = list(
      name = "METSIM Metabolites",
      file_path = DATA_FILES$metsim,
      preprocess_script = file.path(DATA_DIR, "metabolites", "metsim", "build_metsim_exp_p1e5.R"),
      format_params = list(
        snp_col = "variant_id",
        beta_col = "beta",
        se_col = "standard_error",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        pval_col = "p_value",
        samplesize_col = "samplesize",
        phenotype_col = "reportedTrait",
        id_col = "phenocode",
        chr_col = "chromosome",
        pos_col = "base_pair_location",
        header = TRUE,
        min_pval = 1e-200,
        log_pval = FALSE
      ),
      need_eaf = TRUE
    ),

    # ---- Rhee metabolites ----
    rhee = list(
      name = "Rhee Plasma Metabolites",
      file_path = DATA_FILES$rhee,
      default_samplesize = 2076,
      format_params = list(
        snp_col = "rsID",
        beta_col = "beta",
        se_col = "se",
        eaf_col = "MAF",
        effect_allele_col = "min_all",
        other_allele_col = "maj_all",
        pval_col = "pval",
        phenotype_col = "trait",
        id_col = "trait",
        chr_col = "Chr",
        pos_col = "PhysPos",
        header = TRUE,
        min_pval = 1e-200,
        log_pval = FALSE
      ),
      need_eaf = FALSE
    ),
    
    mibiogen = list(
      name = "MiBioGen Microbiome",
      file_path = DATA_FILES$mibiogen,
      format_params = list(
        snp_col = "rsID",
        beta_col = "beta",
        se_col = "SE",
        effect_allele_col = "eff.allele",
        other_allele_col = "ref.allele",
        pval_col = "P.weightedSumZ",
        samplesize_col = "N",
        phenotype_col = "bac",
        z_col = "Z.weightedSumZ",
        chr_col = "chr",
        pos_col = "bp",
        header = TRUE,
        min_pval = 1e-200,
        log_pval = FALSE
      ),
      need_eaf = TRUE
    ),

    ruhlemann = list(
      name = "Ruhlemann Microbiome",
      file_path = DATA_FILES$ruhlemann,
      default_samplesize = 8956,
      format_params = list(
        snp_col = "variant_id",
        beta_col = "beta",
        se_col = "standard_error",
        eaf_col = "effect_allele_frequency",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        pval_col = "p_value",
        phenotype_col = "reportedTrait",
        id_col = "GCST",
        chr_col = "chromosome",
        pos_col = "base_pair_location",
        header = TRUE,
        min_pval = 1e-200,
        log_pval = FALSE
      ),
      need_eaf = TRUE
    ),

    qin = list(
      name = "Qin Microbiome",
      file_path = DATA_FILES$qin,
      default_samplesize = 5959,
      format_params = list(
        snp_col = "variant_id",
        beta_col = "beta",
        se_col = "standard_error",
        eaf_col = "effect_allele_frequency",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        pval_col = "p_value",
        phenotype_col = "reportedTrait",
        id_col = "GCST",
        chr_col = "chromosome",
        pos_col = "base_pair_location",
        header = TRUE,
        min_pval = 1e-200,
        log_pval = FALSE
      ),
      need_eaf = FALSE
    ),

    plasma_pqtl = list(
      name = "UKB-PPP Plasma Proteins",
      file_path = DATA_FILES$plasma_pqtl,
      format_params = list(
        snp_col = "SNP",
        beta_col = "BETA",
        se_col = "SE",
        eaf_col = "Eaf",
        effect_allele_col = "Effect_allele",
        other_allele_col = "Other_allele",
        pval_col = "LOG10P",
        samplesize_col = "N",
        phenotype_col = "Symbol",
        id_col = "UniProt",
        chr_col = "Chr",
        pos_col = "Pos",
        header = TRUE,
        min_pval = 1e-200,
        log_pval = TRUE
      ),
      cis_filter_params = list(
        chr_col = "Chr",
        pos_col = "Pos",
        gene_chr_col = "gene_chr",
        gene_start_col = "gene_start",
        gene_end_col = "gene_end",
        window = 1000000
      ),
      need_eaf = FALSE
    )
  )
  
  if (!source_name %in% names(configs)) {
    stop(paste0("Unknown exposure data source: ", source_name,
                "\nAvailable data sources: ", paste(names(configs), collapse = ", ")))
  }
  
  return(configs[[source_name]])
}

# ==============================================================================
# ==============================================================================

get_outcome_source_config <- function(outcome_name) {
  
  configs <- list(
    
    bp = list(
      name = "Bacterial pneumonia",
      short_name = "bp",
      file = "finngen_R12_J10_PNEUMOBACT.gz",
      cases = 21582,
      controls = 415538,
      trait = "Bacterial pneumoniae",
      id = "J10_PNEUMOBACT",
      format_params = list(
        snp_col = "rsids",
        beta_col = "beta",
        se_col = "sebeta",
        eaf_col = "af_alt",
        effect_allele_col = "alt",
        other_allele_col = "ref",
        pval_col = "pval",
        ncase_col = "cases",
        ncontrol_col = "controls",
        phenotype_col = "trait",
        id_col = "id",
        gene_col = "nearest_genes",
        chr_col = "#chrom",
        pos_col = "pos",
        header = TRUE,
        min_pval = 1e-200,
        log_pval = FALSE
      )
    ),

    bpother_r7 = list(
      name = "Bacterial pneumonia, not elsewhere classified FinnGen R7",
      short_name = "bpother_r7",
      file = "finngen_R7_J10_PNEUMOBACTEROTH.gz",
      cases = 11881,
      controls = 261719,
      trait = "Bacterial pneumonia, not elsewhere classified",
      id = "J10_PNEUMOBACTEROTH",
      format_params = list(
        snp_col = "rsids",
        beta_col = "beta",
        se_col = "sebeta",
        eaf_col = "af_alt",
        effect_allele_col = "alt",
        other_allele_col = "ref",
        pval_col = "pval",
        ncase_col = "cases",
        ncontrol_col = "controls",
        phenotype_col = "trait",
        id_col = "id",
        gene_col = "nearest_genes",
        chr_col = "#chrom",
        pos_col = "pos",
        header = TRUE,
        min_pval = 1e-200,
        log_pval = FALSE
      )
    )
  )
  
  if (!outcome_name %in% names(configs)) {
    stop(paste0("Unknown outcome data source: ", outcome_name,
                "\nAvailable outcomes: ", paste(names(configs), collapse = ", ")))
  }
  
  return(configs[[outcome_name]])
}

# ==============================================================================
# ==============================================================================
