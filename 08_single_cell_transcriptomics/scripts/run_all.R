message("Pipeline start...")

source("scripts/01_preprocess_integrate.R")
source("scripts/02_annotation_subcluster.R")
source("scripts/03_aucell_scoring.R")
source("scripts/04_phase_comparison.R")
source("scripts/05_core_target_scoring.R")
