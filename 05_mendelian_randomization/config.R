# ==============================================================================
# ==============================================================================
# ==============================================================================

PROJECT_ROOT <- getwd()

DATA_DIR <- file.path(PROJECT_ROOT, "01_raw_data")
IVS_DIR <- file.path(PROJECT_ROOT, "02_IVs")
HARMONISE_DIR <- file.path(PROJECT_ROOT, "03_harmonise")
RESULTS_DIR <- file.path(PROJECT_ROOT, "04_results")
SENSITIVITY_DIR <- file.path(PROJECT_ROOT, "05_sensitivity")
PLOTS_DIR <- file.path(PROJECT_ROOT, "06_plots")

PLINK_PATH <- Sys.getenv("PLINK_EXECUTABLE", unset = "plink")
BFILE_PATH <- Sys.getenv("PLINK_REFERENCE_PREFIX", unset = file.path(PROJECT_ROOT, "external_data", "EUR", "EUR"))
MAF_FILE_PATH <- Sys.getenv("LD_REFERENCE_DIRECTORY", unset = file.path(PROJECT_ROOT, "external_data"))

CLUMP_PARAMS <- list(
  kb = 1000,
  r2 = 0.01,
  p = 1
)

PVAL_THRESHOLD <- 1e-5
F_THRESHOLD <- 10
PROXY_R2_THRESHOLD <- 0.8

START_FROM_STEP <- 1

N_CORES <- 1
MRPRESSO_NB_DIST <- 1000
MRPRESSO_MIN_NSNP <- 4

PRIMARY_MR_METHODS <- c("mr_ivw_mre")
SIGNIFICANCE_COLUMN <- "pval"
SIGNIFICANCE_THRESHOLD <- 0.05

SENSITIVITY_SCOPE <- "all"
MRPRESSO_SCOPE <- "significant"

HARMONISE_STRICTNESS <- 2

OUTCOME_PREVALENCE <- NULL

COLOC_PP_H4_THRESHOLD <- 0.95

ALPHA <- 0.05

SENSITIVITY_THRESHOLDS <- list(
  egger_pval_min = 0.05,
  heterogeneity_qep_min = 0.05,
  i2_max = 50,
  presso_pval_min = 0.05,
  require_correct_direction = TRUE,
  require_primary_significance = TRUE,
  significance_col = SIGNIFICANCE_COLUMN,
  significance_threshold = SIGNIFICANCE_THRESHOLD,
  primary_method_pattern = "Inverse variance weighted"
)

# ==============================================================================
# ==============================================================================

#   shin, chen, metsim, rhee
#   mibiogen, ruhlemann, qin
#   tcell, plasma_pqtl
ACTIVE_EXPOSURES <- c("mibiogen", "ruhlemann", "qin")

#   "pna", "vp", "copd", "asthma", "tcell"
ACTIVE_OUTCOME <- "bp"

# filter_mode:
EXPOSURE_FILTER_CONFIG <- list(
  shin = list(
    filter_mode = "selected",
    selected_items = c(
      "Serotonin 5HT levels",
      "Isoleucine levels"
    )
  ),

  chen = list(
    filter_mode = "selected",
    selected_items = c(
      "Serotonin levels",
      "Isoleucine levels",
      "Thyroxine levels"
    )
  ),

  metsim = list(
    filter_mode = "selected",
    selected_items = c(
      "METSIM_C504_serotonin",
      "METSIM_C100003101_alpha_CEHC_glucuronide",
      "METSIM_C2029_azelate_C9_DC",
      "METSIM_C376_isoleucine",
      "METSIM_C1094_thyroxine"
    )
  ),

  rhee = list(
    filter_mode = "selected",
    selected_items = c(
      "serotonin",
      "isoleucine",
      "thyroxine"
    )
  ),

  mibiogen = list(
    filter_mode = "selected",
    selected_items = c(
      "p_Firmicutes",
      "p_Actinobacteria",
      "g_Clostridiuminnocuumgroup",
      "g_Clostridiumsensustricto1",
      "g_Escherichia_Shigella",
      "g_Roseburia",
      "g_Lachnoclostridium",
      "g_Dorea",
      "g_Butyrivibrio",
      "g_Subdoligranulum",
      "g_Coprococcus1",
      "g_Coprococcus2",
      "g_Coprococcus3",
      "g_Eisenbergiella",
      "g_Paraprevotella",
      "g_Anaerostipes"
    ),
    id_map = c(
      "phylum.Firmicutes.id.1672" = "p_Firmicutes",
      "phylum.Actinobacteria.id.400" = "p_Actinobacteria",
      "genus..Clostridiuminnocuumgroup.id.14397" = "g_Clostridiuminnocuumgroup",
      "genus.Clostridiumsensustricto1.id.1873" = "g_Clostridiumsensustricto1",
      "genus.Escherichia.Shigella.id.3504" = "g_Escherichia_Shigella",
      "genus.Roseburia.id.2012" = "g_Roseburia",
      "genus.Lachnoclostridium.id.11308" = "g_Lachnoclostridium",
      "genus.Dorea.id.1997" = "g_Dorea",
      "genus.Butyrivibrio.id.1993" = "g_Butyrivibrio",
      "genus.Subdoligranulum.id.2070" = "g_Subdoligranulum",
      "genus.Coprococcus1.id.11301" = "g_Coprococcus1",
      "genus.Coprococcus2.id.11302" = "g_Coprococcus2",
      "genus.Coprococcus3.id.11303" = "g_Coprococcus3",
      "genus.Eisenbergiella.id.11304" = "g_Eisenbergiella",
      "genus.Paraprevotella.id.962" = "g_Paraprevotella",
      "genus.Anaerostipes.id.1991" = "g_Anaerostipes"
    )
  ),

  ruhlemann = list(
    filter_mode = "selected",
    selected_items = c(
      "P_Firmicutes abundance",
      "P_Actinobacteria abundance",
      "G_Clostridium_XlVa abundance",
      "G_Roseburia abundance",
      "G_Subdoligranulum abundance",
      "G_Coprococcus abundance",
      "G_Paraprevotella abundance"
    )
  ),

  qin = list(
    filter_mode = "selected",
    selected_items = c(
      "Firmicutes A abundance in stool",
      "Firmicutes E abundance in stool",
      "Firmicutes I abundance in stool",
      "Clostridium I abundance in stool",
      "Clostridium P abundance in stool",
      "Dorea abundance in stool"
    )
  ),

  tcell = list(
    filter_mode = "all",
    selected_items = c()
  ),

  plasma_pqtl = list(
    filter_mode = "all",
    selected_items = c()
  )
)

DATA_FILES <- list(
  shin = file.path(DATA_DIR, "metabolites", "exp_Shin_p1e5.csv"),
  chen = file.path(DATA_DIR, "metabolites", "exp_Chen_p1e5.csv"),
  metsim = file.path(DATA_DIR, "metabolites", "exp_METSIM_p1e5.csv"),
  rhee = file.path(DATA_DIR, "metabolites", "exp_Rhee_modified.csv"),

  mibiogen = file.path(DATA_DIR, "Taxa_abundance", "MBG.allHits.p1e4.txt"),
  ruhlemann = file.path(DATA_DIR, "Taxa_abundance", "exp_Ruhlemann_p1e5.csv"),
  qin = file.path(DATA_DIR, "Taxa_abundance", "exp_Qin_p1e5.csv"),

  tcell = file.path(DATA_DIR, "Tcell", "exp_Tcell_p1e5.csv"),
  tcell_outcome = file.path(DATA_DIR, "Tcell", "oc_Tcell_all.csv"),

  plasma_pqtl = file.path(DATA_DIR, "plasma_pQTL", "exp_UKB-PPP_p1e5.csv"),

  lsi_confounder = file.path(IVS_DIR, "confounder_LSI.txt.gz"),
  bmi_confounder = file.path(IVS_DIR, "confounder_BMI.txt.gz")
)

# ==============================================================================
# ==============================================================================

get_exposure_tag <- function(active_sources = ACTIVE_EXPOSURES) {
  paste(active_sources, collapse = "&")
}

generate_output_files <- function(active_sources = ACTIVE_EXPOSURES, outcome_id = ACTIVE_OUTCOME) {
  exp_tag <- get_exposure_tag(active_sources)

  combined_files <- list(
    ivs_all = file.path(IVS_DIR, paste0("ivs_", exp_tag, ".csv")),
    ivs_all_aftpxy = file.path(IVS_DIR, paste0("ivs_", exp_tag, "_aftpxy.csv")),

    oc_outcome = file.path(IVS_DIR, paste0("oc_", outcome_id, ".csv")),
    proxy_file = file.path(IVS_DIR, paste0("ivs_proxysnps_for_", outcome_id, ".csv")),

    harmonise_main = file.path(HARMONISE_DIR, paste0("harmonise_", exp_tag, "_", outcome_id, ".csv")),

    results_main = file.path(RESULTS_DIR, paste0("Rs_", exp_tag, "_", outcome_id, "_onlyivw.csv")),
    results_or = file.path(RESULTS_DIR, paste0("Rs_", exp_tag, "_", outcome_id, "_or.csv")),
    results_or_filtered = file.path(RESULTS_DIR, paste0("Rs_", exp_tag, "_", outcome_id, "_or_filtered.csv")),
    results_sens = file.path(RESULTS_DIR, paste0("Rs_", exp_tag, "_", outcome_id, "_sens.csv")),

    sens_conmix = file.path(SENSITIVITY_DIR, paste0(exp_tag, "_", outcome_id, "_conmix.csv")),
    sens_mrpresso = file.path(SENSITIVITY_DIR, paste0(exp_tag, "_", outcome_id, "_MR-PRESSO.csv")),
    sens_weighted_median = file.path(SENSITIVITY_DIR, paste0(exp_tag, "_", outcome_id, "_weighted_median.csv")),
    sens_heterogeneity = file.path(SENSITIVITY_DIR, paste0(exp_tag, "_", outcome_id, "_heterogeneity.csv")),
    sens_heterogeneity_isq = file.path(SENSITIVITY_DIR, paste0(exp_tag, "_", outcome_id, "_heterogeneity_isq.csv")),
    sens_pleiotropy = file.path(SENSITIVITY_DIR, paste0(exp_tag, "_", outcome_id, "_pleiotropy.csv")),
    sens_steiger = file.path(SENSITIVITY_DIR, paste0(exp_tag, "_", outcome_id, "_steiger.csv")),

    plot_forestplot = file.path(PLOTS_DIR, paste0(exp_tag, "_", outcome_id, "_forestplot.pdf"))
  )

  return(combined_files)
}

OUTPUT_FILES <- generate_output_files(ACTIVE_EXPOSURES, ACTIVE_OUTCOME)

VERBOSE <- TRUE

log_info <- function(msg) {
  if (VERBOSE) {
    cat(paste0("[", Sys.time(), "] INFO: ", msg, "\n"))
  }
}

log_warning <- function(msg) {
  cat(paste0("[", Sys.time(), "] WARNING: ", msg, "\n"))
}

log_error <- function(msg) {
  cat(paste0("[", Sys.time(), "] ERROR: ", msg, "\n"))
}

check_directories <- function() {
  dirs <- c(DATA_DIR, IVS_DIR, HARMONISE_DIR, RESULTS_DIR, SENSITIVITY_DIR, PLOTS_DIR)
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      warning(paste0("Directory does not exist: ", dir))
    }
  }
}

check_plink <- function() {
  if (!file.exists(PLINK_PATH)) {
    stop("PLINK executable does not exist: ", PLINK_PATH)
  }
  if (!file.exists(paste0(BFILE_PATH, ".bed"))) {
    stop("Reference-genome files do not exist: ", BFILE_PATH)
  }
}

if (interactive()) {
  log_info("Configuration loaded")
  log_info(paste0("Project root: ", PROJECT_ROOT))
  check_directories()
}

# ==============================================================================
# ==============================================================================
