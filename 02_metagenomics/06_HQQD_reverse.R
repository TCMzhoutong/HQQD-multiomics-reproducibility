# ============================================================================
# LEfSe Differential Taxa Boxplot with Statistical Annotation (HQQD vs Model Selection)
# Author: Generated for metagenomics analysis
# Date: 2026-01-30
# Description: Horizontal boxplot showing top 10 differential taxa
#              Selection: Taxa with significant callback trend (Model vs Control & HQQD vs Model)
#              Plotting: STAMP-like extended error bar plot with Dunn's test
# ============================================================================

# Clear environment
rm(list = ls())

# Load required libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readxl,       # For reading Excel files
  dplyr,        # For data manipulation
  tidyr,        # For data reshaping
  ggplot2,      # For plotting
  ggpubr,       # For statistical annotation
  stringr,      # For string manipulation
  patchwork,    # For combining plots
  microeco,     # For LEfSe analysis
  magrittr,     # For pipe operations
  tibble        # For column_to_rownames
)

# ============================================================================
# Configuration
# ============================================================================

# Define color palette (same as previous analyses)
group_colors <- c(
  "Control" = "#3B4992", 
  "Model" = "#EE0000", 
  "LVX" = "#008B45", 
  "HQQD" = "#FF8C00"
)

# Input files
abundance_file <- "taxonomy.all_level.relabundance.xlsx"

# Output directory
output_dir <- "06_HQQD_reverse"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Parameters
TOP_N <- 10                 # Number of top taxa to display per taxonomic level
LDA_CUTOFF <- 2.0           # Minimum LDA score threshold
PSEUDO_COUNT <- 0.0001      # Small constant for log transformation
ALPHA <- 0.05               # Significance level
source("reverse_fdr_plot_helper.R")

# ============================================================================
# Step 1: Load Data
# ============================================================================

cat("\n=== Step 1: Loading data ===\n")

# Load relative abundance data
abundance_data <- read_excel(abundance_file)
cat("Abundance data loaded:", nrow(abundance_data), "taxa\n")

# Check column names
sample_cols <- grep("^[KMYZ][0-9]", colnames(abundance_data), value = TRUE)
cat("Number of samples:", length(sample_cols), "\n\n")

# ============================================================================
# Step 2: Create Sample Metadata
# ============================================================================

cat("=== Step 2: Creating sample metadata ===\n")

sample_metadata <- data.frame(
  Sample = sample_cols,
  Group = case_when(
    grepl("^K", sample_cols) ~ "Control",
    grepl("^M", sample_cols) ~ "Model",
    grepl("^Y", sample_cols) ~ "LVX",
    grepl("^Z", sample_cols) ~ "HQQD",
    TRUE ~ "Unknown"
  ),
  stringsAsFactors = FALSE
)

# Set rownames for microeco compatibility
rownames(sample_metadata) <- sample_metadata$Sample

# Convert to factor with proper ordering (reversed for top-to-bottom display)
sample_metadata$Group <- factor(sample_metadata$Group, 
                                levels = c("HQQD", "LVX", "Model", "Control"))

cat("Sample distribution:\n")
print(table(sample_metadata$Group))
cat("\n")

# ============================================================================
# Step 3: Helper Functions
# ============================================================================

# Function to extract taxonomic level from Taxa string
get_taxonomic_level <- function(taxa_string) {
  # Split by |
  levels <- strsplit(taxa_string, "\\|")[[1]]
  last_level <- levels[length(levels)]
  
  # Extract prefix
  if (grepl("^p__", last_level)) return("Phylum")
  if (grepl("^g__", last_level)) return("Genus")
  if (grepl("^s__", last_level)) return("Species")
  return("Other")
}

# Function to get clean taxon name
get_clean_name <- function(taxa_string, keep_prefix = FALSE) {
  levels <- strsplit(taxa_string, "\\|")[[1]]
  last_level <- levels[length(levels)]
  
  if (keep_prefix) {
    return(last_level)
  } else {
    # Remove prefix like p__, g__, s__
    return(sub("^[a-z]__", "", last_level))
  }
}

# Function to match LEfSe taxa to abundance data
match_taxa_to_abundance <- function(taxa_string, abundance_data) {
  # Parse taxonomic path
  levels <- strsplit(taxa_string, "\\|")[[1]]
  
  # Build filter conditions
  conditions <- rep(TRUE, nrow(abundance_data))
  
  for (level in levels) {
    if (grepl("^k__", level)) {
      conditions <- conditions & (abundance_data$kingdom == level)
    } else if (grepl("^p__", level)) {
      conditions <- conditions & (abundance_data$phylum == level)
    } else if (grepl("^c__", level)) {
      conditions <- conditions & (abundance_data$class == level)
    } else if (grepl("^o__", level)) {
      conditions <- conditions & (abundance_data$order == level)
    } else if (grepl("^f__", level)) {
      conditions <- conditions & (abundance_data$family == level)
    } else if (grepl("^g__", level)) {
      conditions <- conditions & (abundance_data$genus == level)
    } else if (grepl("^s__", level)) {
      conditions <- conditions & (abundance_data$species == level)
    }
  }
  
  matched_rows <- abundance_data[conditions, ]
  
  # Sum abundances if multiple rows match
  if (nrow(matched_rows) > 0) {
    abundance_values <- colSums(matched_rows[, sample_cols, drop = FALSE])
    return(abundance_values)
  } else {
    return(rep(0, length(sample_cols)))
  }
}

# Function to run LEfSe for two groups at specific taxonomic level and get LDA score for a taxon
get_lefse_lda_for_taxon <- function(taxa_string, abundance_data, group1_name, group2_name, tax_level, sample_metadata, sample_cols) {
  # Create cache key
  cache_key <- paste(tax_level, group1_name, group2_name, sep="_")
  
  # Use a global cache (we'll create this outside)
  if (!exists("LEFSE_CACHE", envir = .GlobalEnv)) {
    assign("LEFSE_CACHE", list(), envir = .GlobalEnv)
  }
  
  LEFSE_CACHE <- get("LEFSE_CACHE", envir = .GlobalEnv)
  
  # Check if we already computed LEfSe for this level and comparison
  if (!(cache_key %in% names(LEFSE_CACHE))) {
    # Need to run LEfSe for this level
    cat("    Running LEfSe for", cache_key, "...")
    
    # Filter sample metadata for these two groups
    comp_metadata <- sample_metadata %>%
      filter(Group %in% c(group1_name, group2_name)) %>%
      mutate(Group = factor(as.character(Group), levels = c(group1_name, group2_name))) %>%
      as.data.frame()
    
    # Set rownames to Sample column (critical for microeco)
    rownames(comp_metadata) <- comp_metadata$Sample
    
    if (nrow(comp_metadata) < 2) {
      LEFSE_CACHE[[cache_key]] <- data.frame()
      assign("LEFSE_CACHE", LEFSE_CACHE, envir = .GlobalEnv)
      return(0)
    }
    
    comp_samples <- comp_metadata$Sample
    
    # Aggregate abundance at this taxonomic level
    level_col <- tolower(tax_level)
    
    agg_abundance <- abundance_data %>%
      filter(!is.na(!!sym(level_col)) & !!sym(level_col) != "") %>%
      group_by(!!sym(level_col)) %>%
      summarise(across(all_of(comp_samples), sum), .groups = "drop")
    
    if (nrow(agg_abundance) == 0) {
      cat(" no data\n")
      LEFSE_CACHE[[cache_key]] <- data.frame()
      assign("LEFSE_CACHE", LEFSE_CACHE, envir = .GlobalEnv)
      return(0)
    }
    
    # Create OTU and tax table
    agg_abundance$OTU_ID <- paste0("OTU_", sprintf("%05d", 1:nrow(agg_abundance)))
    
    otu_table_comp <- agg_abundance %>%
      select(OTU_ID, all_of(comp_samples)) %>%
      tibble::column_to_rownames("OTU_ID") %>%
      as.data.frame()
    
    # Create tax table matching microeco format (with uppercase column names)
    tax_table_comp <- data.frame(
      Kingdom = rep(tax_level, nrow(agg_abundance)),
      Phylum = rep(tax_level, nrow(agg_abundance)),
      Class = rep(tax_level, nrow(agg_abundance)),
      Order = rep(tax_level, nrow(agg_abundance)),
      Family = rep(tax_level, nrow(agg_abundance)),
      Genus = rep(tax_level, nrow(agg_abundance)),
      Species = agg_abundance[[level_col]],
      row.names = agg_abundance$OTU_ID,
      stringsAsFactors = FALSE
    )
    
    # Store the actual taxon names for later matching
    taxon_names_map <- setNames(agg_abundance[[level_col]], agg_abundance$OTU_ID)
    # Store the actual taxon names for later matching
    taxon_names_map <- setNames(agg_abundance[[level_col]], agg_abundance$OTU_ID)
    
    tryCatch({
      comp_dataset <- microtable$new(
        sample_table = comp_metadata,
        otu_table = otu_table_comp,
        tax_table = tax_table_comp
      )
      
      comp_lefse <- trans_diff$new(
        dataset = comp_dataset,
        method = "lefse",
        group = "Group",
        alpha = ALPHA,
        lefse_subgroup = NULL,
        p_adjust_method = "fdr"
      )
      
      # Store results with taxa names
      if (!is.null(comp_lefse$res_diff) && nrow(comp_lefse$res_diff) > 0) {
        # Map OTU IDs back to taxon names (keep with prefix for matching)
        comp_lefse$res_diff$Taxon_Name <- sapply(comp_lefse$res_diff$Taxa, function(x) {
          # Extract the OTU ID if present
          if (grepl("^OTU_", x)) {
            # It's an OTU ID
            if (x %in% names(taxon_names_map)) {
              return(taxon_names_map[[x]])  # This returns the full name with prefix
            }
          }
          # Try to extract from the full taxonomic path
          parts <- strsplit(x, "\\|")[[1]]
          if (length(parts) > 0) {
            last_part <- parts[length(parts)]
            # Keep the prefix for matching
            return(last_part)
          }
          return(x)
        })
      }
      
      LEFSE_CACHE[[cache_key]] <- comp_lefse$res_diff
      cat(" found", nrow(comp_lefse$res_diff), "features\n")
    }, error = function(e) {
      cat(" ERROR:", e$message, "\n")
      LEFSE_CACHE[[cache_key]] <- data.frame()
    })
    
    assign("LEFSE_CACHE", LEFSE_CACHE, envir = .GlobalEnv)
  }
  
  # Get cached results
  LEFSE_CACHE <- get("LEFSE_CACHE", envir = .GlobalEnv)
  lefse_res <- LEFSE_CACHE[[cache_key]]
  
  if (is.null(lefse_res) || nrow(lefse_res) == 0) {
    return(0)
  }
  
  # Extract the target taxon name from taxa_string
  target_taxon <- get_clean_name(taxa_string, keep_prefix = TRUE)
  
  # Find matching rows
  matching_rows <- lefse_res[lefse_res$Taxon_Name == target_taxon, ]
  
  if (nrow(matching_rows) > 0) {
    return(matching_rows$LDA[1])
  } else {
    return(0)
  }
}

# Function to check callback trend (Model vs Control and HQQD vs Model)
check_callback_trend <- function(taxa_string, abundance_data, tax_level, sample_metadata, sample_cols) {
  # Get abundances
  abundances <- match_taxa_to_abundance(taxa_string, abundance_data)
  
  # Create temp dataframe
  temp_df <- data.frame(
    Abundance = as.numeric(abundances),
    Group = sample_metadata$Group
  )
  temp_df$Log_Abundance <- log10(temp_df$Abundance + PSEUDO_COUNT)
  
  # Extract log abundances for each group
  control_vals <- temp_df$Log_Abundance[temp_df$Group == "Control"]
  model_vals <- temp_df$Log_Abundance[temp_df$Group == "Model"]
  hqqd_vals <- temp_df$Log_Abundance[temp_df$Group == "HQQD"]
  
  # Calculate means for trend direction
  means <- temp_df %>% 
    group_by(Group) %>% 
    summarize(Mean = mean(Log_Abundance), .groups = "drop")
  
  mean_control <- means$Mean[means$Group == "Control"]
  mean_model <- means$Mean[means$Group == "Model"]
  mean_hqqd <- means$Mean[means$Group == "HQQD"]
  
   # Get LDA scores from LEfSe results
  lda_model_control <- get_lefse_lda_for_taxon(taxa_string, abundance_data, "Model", "Control", tax_level, sample_metadata, sample_cols)
  lda_hqqd_model <- get_lefse_lda_for_taxon(taxa_string, abundance_data, "HQQD", "Model", tax_level, sample_metadata, sample_cols)
  
  # Wilcoxon rank-sum tests for the two key comparisons. These are calculated
  # for all candidates so the BH-adjusted output uses real p-values across the
  # tested feature space; LDA still controls the original reversal selection.
  p_mc <- tryCatch(
    wilcox.test(model_vals, control_vals, exact = FALSE)$p.value,
    error = function(e) 1.0
  )
  p_hm <- tryCatch(
    wilcox.test(hqqd_vals, model_vals, exact = FALSE)$p.value,
    error = function(e) 1.0
  )
  
  trend_mc <- "None"
  if (!is.na(p_mc) && p_mc < ALPHA) {
    trend_mc <- ifelse(mean_model > mean_control, "Up", "Down")
  }
  
  trend_hm <- "None"
  if (!is.na(p_hm) && p_hm < ALPHA) {
    trend_hm <- ifelse(mean_hqqd > mean_model, "Up", "Down")
  }
  
  # Keep the original selection rule: both pairwise LEfSe LDA scores must pass.
  if (lda_model_control < LDA_CUTOFF || lda_hqqd_model < LDA_CUTOFF) {
    return(list(
      is_callback = FALSE,
      p_model_control = p_mc,
      p_hqqd_model = p_hm,
      trend_model_vs_control = trend_mc,
      trend_hqqd_vs_model = trend_hm,
      lda_model_control = lda_model_control,
      lda_hqqd_model = lda_hqqd_model
    ))
  }
  
  # Callback condition: both significant and trends are opposite
  is_callback <- (!is.na(p_mc) && p_mc < ALPHA) &&
                 (!is.na(p_hm) && p_hm < ALPHA) &&
                 (trend_mc != "None") && (trend_hm != "None") &&
                 (trend_mc != trend_hm)
  
  return(list(
    is_callback = is_callback,
    p_model_control = p_mc,
    p_hqqd_model = p_hm,
    trend_model_vs_control = trend_mc,
    trend_hqqd_vs_model = trend_hm,
    lda_model_control = lda_model_control,
    lda_hqqd_model = lda_hqqd_model
  ))
}

# ============================================================================
# Step 4: Process Each Taxonomic Level
# ============================================================================

# Initialize global LEfSe cache
if (exists("LEFSE_CACHE", envir = .GlobalEnv)) {
  rm(LEFSE_CACHE, envir = .GlobalEnv)
}

taxonomic_levels <- c("Phylum", "Genus", "Species")

for (tax_level in taxonomic_levels) {
  
  cat("\n", rep("=", 70), "\n", sep="")
  cat("=== Processing", tax_level, "Level ===\n")
  cat(rep("=", 70), "\n\n", sep="")
  
  # Get all taxa at this taxonomic level
  level_column <- tolower(tax_level)
  
  # Get unique taxa at this level (non-empty)
  unique_taxa <- abundance_data %>%
    filter(!is.na(!!sym(level_column)) & !!sym(level_column) != "") %>%
    pull(!!sym(level_column)) %>%
    unique()
  
  cat("Total taxa at", tax_level, "level:", length(unique_taxa), "\n")
  
  # Build taxa strings for each unique taxon
  level_results <- data.frame(
    Taxa = character(),
    stringsAsFactors = FALSE
  )
  
  for (taxon in unique_taxa) {
    # Find a representative row for this taxon
    repr_row <- abundance_data %>%
      filter(!!sym(level_column) == taxon) %>%
      slice(1)
    
    # Build taxonomic path
    taxa_path <- c()
    if (!is.na(repr_row$kingdom) && repr_row$kingdom != "") taxa_path <- c(taxa_path, repr_row$kingdom)
    if (!is.na(repr_row$phylum) && repr_row$phylum != "") taxa_path <- c(taxa_path, repr_row$phylum)
    if (tax_level == "Phylum") {
      # For Phylum, stop here
    } else {
      if (!is.na(repr_row$class) && repr_row$class != "") taxa_path <- c(taxa_path, repr_row$class)
      if (!is.na(repr_row$order) && repr_row$order != "") taxa_path <- c(taxa_path, repr_row$order)
      if (!is.na(repr_row$family) && repr_row$family != "") taxa_path <- c(taxa_path, repr_row$family)
      if (tax_level == "Genus") {
        if (!is.na(repr_row$genus) && repr_row$genus != "") taxa_path <- c(taxa_path, repr_row$genus)
      } else if (tax_level == "Species") {
        if (!is.na(repr_row$genus) && repr_row$genus != "") taxa_path <- c(taxa_path, repr_row$genus)
        if (!is.na(repr_row$species) && repr_row$species != "") taxa_path <- c(taxa_path, repr_row$species)
      }
    }
    
    taxa_string <- paste(taxa_path, collapse = "|")
    level_results <- rbind(level_results, data.frame(Taxa = taxa_string, stringsAsFactors = FALSE))
  }
  
  if (nrow(level_results) == 0) {
    cat("No taxa found at this level. Skipping.\n")
    next
  }

  cat("Total candidates for filtering:", nrow(level_results), "\n")
  cat("Checking for significant callback trend (Model vs Control AND HQQD vs Model)...\n")
  cat("LDA threshold:", LDA_CUTOFF, "\n")
  
  # Calculate callback trends for all candidates
  cat("Analyzing callback trends for all taxa at", tax_level, "level...\n") 
  callback_results <- list()
  for(i in 1:nrow(level_results)) {
     callback_results[[i]] <- check_callback_trend(level_results$Taxa[i], abundance_data, tax_level, sample_metadata, sample_cols)
     if(i %% 10 == 0) cat(".")
  }
  cat("\n")
  
  level_results$is_callback <- sapply(callback_results, function(x) x$is_callback)
  level_results$p_val_Model_Control <- sapply(callback_results, function(x) x$p_model_control)
  level_results$p_val_HQQD_Model <- sapply(callback_results, function(x) x$p_hqqd_model)
  level_results$Trend_Model_vs_Control <- sapply(callback_results, function(x) x$trend_model_vs_control)
  level_results$Trend_HQQD_vs_Model <- sapply(callback_results, function(x) x$trend_hqqd_vs_model)
  level_results$LDA_Model_Control <- sapply(callback_results, function(x) x$lda_model_control)
  level_results$LDA_HQQD_Model <- sapply(callback_results, function(x) x$lda_hqqd_model)
  level_results$q_val_Model_Control <- p.adjust(level_results$p_val_Model_Control, method = "BH")
  level_results$q_val_HQQD_Model <- p.adjust(level_results$p_val_HQQD_Model, method = "BH")
  level_results$is_callback_fdr <- level_results$is_callback &
    !is.na(level_results$q_val_Model_Control) &
    !is.na(level_results$q_val_HQQD_Model) &
    level_results$q_val_Model_Control < ALPHA &
    level_results$q_val_HQQD_Model < ALPHA
  
  # Filter significant ones
  sig_results <- level_results %>% 
    filter(is_callback == TRUE)
  
  cat("Taxa with significant callback trend:", nrow(sig_results), "\n")
  cat("Taxa retained after BH correction:", sum(level_results$is_callback_fdr, na.rm = TRUE), "\n")
  
  if (nrow(sig_results) == 0) {
    cat("No taxa found with significant callback trend. Skipping.\n")
    next
  }
  
  # ============================================================================
  # Perform statistical analysis on ALL significant results for CSV output
  # ============================================================================
  
  cat("\nPerforming Wilcoxon rank-sum tests on all significant taxa...\n")
  
  all_stat_results <- list()
  
  for (i in 1:nrow(sig_results)) {
    taxa_string <- sig_results$Taxa[i]
    
    # Get abundances
    abundances <- match_taxa_to_abundance(taxa_string, abundance_data)
    
    # Create data frame (include all 4 groups)
    temp_data <- data.frame(
      Abundance = as.numeric(abundances),
      Group = sample_metadata$Group
    )
    temp_data$Log_Abundance <- log10(temp_data$Abundance + PSEUDO_COUNT)
    
    # Extract per-group log-abundance vectors
    mc_vals   <- temp_data$Log_Abundance[temp_data$Group == "Model"]
    ctrl_vals <- temp_data$Log_Abundance[temp_data$Group == "Control"]
    hqqd_vals_s <- temp_data$Log_Abundance[temp_data$Group == "HQQD"]
    lvx_vals  <- temp_data$Log_Abundance[temp_data$Group == "LVX"]
    
    # Three Wilcoxon rank-sum tests
    p_mc <- tryCatch(wilcox.test(mc_vals, ctrl_vals,    exact = FALSE)$p.value, error = function(e) NA)
    p_hm <- tryCatch(wilcox.test(hqqd_vals_s, mc_vals, exact = FALSE)$p.value, error = function(e) NA)
    p_ml <- tryCatch(wilcox.test(mc_vals, lvx_vals,    exact = FALSE)$p.value, error = function(e) NA)
    
    clean_name <- get_clean_name(taxa_string, keep_prefix = FALSE)
    all_stat_results[[clean_name]] <- list(
      Taxa_Full          = taxa_string,
      p_Model_vs_Control = p_mc,
      p_HQQD_vs_Model    = p_hm,
      p_Model_vs_LVX     = p_ml
    )
    
    if(i %% 10 == 0) cat(".")
  }
  cat("\n")
  
  # Select all significant taxa, sorted by average LDA (descending)
  sig_results$LDA_avg <- (sig_results$LDA_Model_Control + sig_results$LDA_HQQD_Model) / 2
  
  top_taxa <- sig_results %>%
    arrange(desc(LDA_avg))
  
  cat("Displaying all", nrow(top_taxa), "significant taxa\n\n")
  
  # Prepare data for plotting
  plot_data_list <- list()
  
  for (i in 1:nrow(top_taxa)) {
    taxa_row <- top_taxa[i, ]
    taxa_string <- taxa_row$Taxa
    taxa_lda_avg <- taxa_row$LDA_avg
    
    # Get clean name
    clean_name <- get_clean_name(taxa_string, keep_prefix = FALSE)
    
    cat("Processing:", clean_name, "(LDA_avg:", round(taxa_lda_avg, 2), ")\n")
    
    # Match to abundance data
    abundances <- match_taxa_to_abundance(taxa_string, abundance_data)
    
    # Create data frame
    taxa_data <- data.frame(
      Sample = sample_cols,
      Abundance = as.numeric(abundances),
      Taxon = clean_name,
      Taxon_Full = taxa_string,
      LDA_Score = taxa_lda_avg,
      stringsAsFactors = FALSE
    )
    
    # Add group information
    taxa_data <- taxa_data %>%
      left_join(sample_metadata, by = "Sample")
    
    # Log transformation
    taxa_data$Log_Abundance <- log10(taxa_data$Abundance + PSEUDO_COUNT)
    
    plot_data_list[[i]] <- taxa_data
  }
  
  # Combine all taxa data
  plot_data_combined <- bind_rows(plot_data_list)
  
  # Order taxa by average LDA score
  taxa_order <- top_taxa %>%
    arrange(desc(LDA_avg)) %>%
    mutate(Clean_Name = sapply(Taxa, function(x) get_clean_name(x, keep_prefix = FALSE)))
  
  plot_data_combined$Taxon <- factor(plot_data_combined$Taxon, 
                                     levels = rev(taxa_order$Clean_Name))
  
  # ============================================================================
  # Step 5: Create Left Plot (Boxplot)
  # ============================================================================
  
  cat("\nGenerating boxplot...\n")
  
  # Prepare alternating background rectangles
  n_taxa <- length(levels(plot_data_combined$Taxon))
  alternating_bg <- data.frame(
    ymin = seq(1.5, n_taxa - 0.5, by = 1),
    ymax = seq(2.5, n_taxa + 0.5, by = 1)
  )
  alternating_bg <- alternating_bg[seq(1, nrow(alternating_bg), by = 2), ]
  
  # Create left plot: boxplot without statistical annotations
  p_left <- ggplot(plot_data_combined, aes(x = Log_Abundance, y = Taxon, color = Group, fill = Group)) +
    # Add alternating background
    geom_rect(data = alternating_bg, 
              aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              fill = "gray85", alpha = 0.5, inherit.aes = FALSE) +
    geom_boxplot(
                 position = position_dodge(width = 0.7),
                 outlier.shape = NA,
                 width = 0.35,
                 alpha = 0.3,
                 linewidth = 0.5) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE)) +
    labs(
      x = expression(Log[10]~"Relative Abundance"),
      y = NULL,
      title = paste0(
        nrow(top_taxa), " Differential Taxa at ", tax_level, " Level"
      ),
      subtitle = "Filtered by HQQD-induced Reversal, Sorted by LDA"
    ) +
    theme_bw(base_family = "sans") +
    theme(
      axis.text.y = element_text(size = 10, hjust = 1, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12, face = "bold"),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(color = "black", linewidth = 0.5),
      axis.line.y.left = element_line(color = "black", linewidth = 0.5)
    )
  
  # Set y-axis labels to black (unified color)
  p_left <- p_left + theme(axis.text.y = element_text(color = "black", face = "bold"))
  
  # Add top border to left plot
  p_left <- p_left + 
    geom_segment(aes(x = -Inf, xend = Inf, y = Inf, yend = Inf), 
                 color = "black", linewidth = 0.5, inherit.aes = FALSE)
  
  # ============================================================================
  # Step 6: Create Right Plot (Statistical Annotations Panel)
  # ============================================================================
  
  cat("Creating statistical annotations panel...\n")
  
  # Create a blank plot with same y-axis as left plot
  p_middle <- ggplot(plot_data_combined, aes(x = 1, y = Taxon)) +
    # Add alternating background (same as left plot)
    geom_rect(data = alternating_bg, 
              aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              fill = "gray85", alpha = 0.5, inherit.aes = FALSE) +
    scale_y_discrete(limits = levels(plot_data_combined$Taxon)) +
    theme_void(base_family = "sans") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(10, 10, 10, 5),
      panel.border = element_blank()
    ) +
    # Add borders: top, right, bottom
    geom_segment(aes(x = Inf, xend = Inf, y = -Inf, yend = Inf), 
                 color = "black", linewidth = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = -Inf, xend = Inf, y = Inf, yend = Inf), 
                 color = "black", linewidth = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = -Inf, xend = Inf, y = -Inf, yend = -Inf), 
                 color = "black", linewidth = 0.5, inherit.aes = FALSE) +
    coord_cartesian(xlim = c(0, 10), clip = "off")
  
  # Perform three Wilcoxon rank-sum tests for each taxon (for plot annotations)
  # Comparisons: Model vs Control | HQQD vs Model | Model vs LVX
  plot_stat_results <- list()
  
  for (taxon in levels(plot_data_combined$Taxon)) {
    taxon_data <- plot_data_combined %>% filter(Taxon == taxon)
    
    mc_v   <- taxon_data$Log_Abundance[taxon_data$Group == "Model"]
    ctrl_v <- taxon_data$Log_Abundance[taxon_data$Group == "Control"]
    hqqd_v <- taxon_data$Log_Abundance[taxon_data$Group == "HQQD"]
    lvx_v  <- taxon_data$Log_Abundance[taxon_data$Group == "LVX"]
    
    p_mc <- tryCatch(wilcox.test(mc_v, ctrl_v, exact = FALSE)$p.value, error = function(e) NA)
    p_hm <- tryCatch(wilcox.test(hqqd_v, mc_v,  exact = FALSE)$p.value, error = function(e) NA)
    p_ml <- tryCatch(wilcox.test(mc_v, lvx_v,  exact = FALSE)$p.value, error = function(e) NA)
    
    plot_stat_results[[as.character(taxon)]] <- list(
      p_Model_vs_Control = p_mc,
      p_HQQD_vs_Model    = p_hm,
      p_Model_vs_LVX     = p_ml
    )
  }
  
  # Fixed set of 3 comparisons to annotate (drawn left to right per taxon)
  comparisons_fixed <- list(
    list(group1 = "Model", group2 = "Control", p_key = "p_Model_vs_Control"),
    list(group1 = "HQQD",  group2 = "Model",   p_key = "p_HQQD_vs_Model"),
    list(group1 = "Model", group2 = "LVX",     p_key = "p_Model_vs_LVX")
  )
  
  group_positions <- c("HQQD" = 0, "LVX" = 1, "Model" = 2, "Control" = 3)
  dodge_width <- 0.7
  n_groups    <- 4
  
  cat("  Adding Wilcoxon significance annotations (Model vs Control, HQQD vs Model, Model vs LVX)...\n")
  
  for (taxon_name in levels(plot_data_combined$Taxon)) {
    taxon_idx <- which(levels(plot_data_combined$Taxon) == taxon_name)
    stat_res  <- plot_stat_results[[taxon_name]]
    current_x <- 2
    
    for (comp in comparisons_fixed) {
      p_val <- stat_res[[comp$p_key]]
      if (is.na(p_val) || p_val >= 0.05) {
        # Still advance x so columns stay consistent across taxa
        current_x <- current_x + 1.2
        next
      }
      
      sig_label <- dplyr::case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        TRUE          ~ "*"
      )
      
      g1_idx <- group_positions[comp$group1]
      g2_idx <- group_positions[comp$group2]
      
      y1 <- taxon_idx + (g1_idx - (n_groups - 1) / 2) * dodge_width / n_groups
      y2 <- taxon_idx + (g2_idx - (n_groups - 1) / 2) * dodge_width / n_groups
      
      star_spacing <- dplyr::case_when(
        sig_label == "***" ~ 2.6,
        sig_label == "**"  ~ 2.0,
        TRUE               ~ 1.2
      )
      
      p_middle <- p_middle +
        annotate("segment",
                 x = current_x, xend = current_x,
                 y = y1, yend = y2,
                 color = "gray30", linewidth = 0.6) +
        annotate("text",
                 x = current_x + 0.3,
                 y = (y1 + y2) / 2,
                 label = sig_label,
                 size = 3.5, fontface = "bold",
                 color = "gray20", hjust = 0)
      
      current_x <- current_x + star_spacing
    }
  }
  
  # ============================================================================
  # Step 7: Create Legend Plot (Right Panel)
  # ============================================================================
  
  cat("Creating legend panel...\n")
  
  # Extract legend from a temporary plot
  p_legend_temp <- ggplot(plot_data_combined, aes(x = Log_Abundance, y = Taxon, color = Group, fill = Group)) +
    geom_boxplot(alpha = 0.3) +
    scale_color_manual(values = group_colors, name = "Group") +
    scale_fill_manual(values = group_colors, name = "Group") +
    guides(
      color = guide_legend(reverse = TRUE, override.aes = list(fill = group_colors, alpha = 0.3)),
      fill = "none"
    ) +
    theme_minimal(base_family = "sans") +
    theme(
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.8, "cm")
    )
  
  # Extract just the legend
  library(ggpubr)
  p_legend <- as_ggplot(get_legend(p_legend_temp))
  
  # ============================================================================
  # Step 8: Combine Plots with Patchwork
  # ============================================================================
  
  cat("\nCombining plots...\n")
  
  # Combine: left (boxplot) + middle (annotation) + right (legend)
  # Widths: 5 : 1 : 0.8 (boxplot : annotation : legend)
  p_combined <- p_left + p_middle + p_legend + 
    plot_layout(widths = c(5, 1, 0.8))
  
  # ============================================================================
  # Step 9: Save Combined Plot
  # ============================================================================
  
  output_file_pdf <- file.path(output_dir, paste0("LEfSe_boxplot_", tax_level, ".pdf"))
  output_file_png <- file.path(output_dir, paste0("LEfSe_boxplot_", tax_level, ".png"))
  
  # Calculate dynamic height (compact spacing per taxon)
  plot_height <- max(6, nrow(top_taxa) * 0.55 + 2)
  # Adjust width based on taxonomic level (Species has longer names)
  plot_width <- ifelse(tax_level == "Species", 9, 8)
  
  ggsave(
    filename = output_file_pdf,
    plot = p_combined,
    width = plot_width,
    height = plot_height,
    units = "in",
    dpi = 300
  )
  cat("PDF saved to:", output_file_pdf, "\n")
  
  ggsave(
    filename = output_file_png,
    plot = p_combined,
    width = plot_width,
    height = plot_height,
    units = "in",
    dpi = 300
  )
  cat("PNG saved to:", output_file_png, "\n")
  
  # ============================================================================
  # Save ALL statistical results (not just top_taxa)
  # ============================================================================
  
  cat("\nSaving statistical results...\n")
  
  # Build summary with three Wilcoxon p-values
  stat_summary <- data.frame(
    Taxon              = names(all_stat_results),
    Taxa_Full          = sapply(all_stat_results, function(x) x$Taxa_Full),
    p_Model_vs_Control = sapply(all_stat_results, function(x) x$p_Model_vs_Control),
    p_HQQD_vs_Model    = sapply(all_stat_results, function(x) x$p_HQQD_vs_Model),
    p_Model_vs_LVX     = sapply(all_stat_results, function(x) x$p_Model_vs_LVX),
    stringsAsFactors = FALSE
  )
  
  # Merge with LDA scores and Trends from sig_results
  lda_info <- sig_results %>%
    mutate(Clean_Name = sapply(Taxa, function(x) get_clean_name(x, keep_prefix = FALSE))) %>%
    select(
      Clean_Name,
      LDA_Model_Control,
      LDA_HQQD_Model,
      LDA_avg,
      Trend_Model_vs_Control,
      Trend_HQQD_vs_Model,
      q_val_Model_Control,
      q_val_HQQD_Model,
      is_callback_fdr
    ) %>%
    rename(
      q_Model_vs_Control = q_val_Model_Control,
      q_HQQD_vs_Model = q_val_HQQD_Model,
      FDR_corrected_callback = is_callback_fdr
    )
  
  stat_summary <- stat_summary %>%
    left_join(lda_info, by = c("Taxon" = "Clean_Name")) %>%
    arrange(desc(LDA_avg))
  
  write.csv(stat_summary, 
            file.path(output_dir, paste0("statistical_results_", tax_level, ".csv")),
            row.names = FALSE)
  
  cat("Statistical results saved:", nrow(stat_summary), "taxa\n")
  
  stat_summary_fdr <- stat_summary %>%
    filter(FDR_corrected_callback == TRUE)
  
  write.csv(stat_summary_fdr,
            file.path(output_dir, paste0("statistical_results_", tax_level, "_FDR_corrected.csv")),
            row.names = FALSE)

  write.csv(stat_summary_fdr,
            file.path(output_dir, paste0("statistical_results_", tax_level, "_fdr.csv")),
            row.names = FALSE)

  if (nrow(stat_summary_fdr) > 0) {
    fdr_order <- stat_summary_fdr %>% arrange(desc(LDA_avg)) %>% pull(Taxon)
    plot_data_fdr <- plot_data_combined %>%
      filter(as.character(Taxon) %in% fdr_order)
    plot_data_fdr$Taxon <- factor(as.character(plot_data_fdr$Taxon), levels = rev(fdr_order))
    save_reverse_boxplot(
      plot_data_fdr,
      file.path(output_dir, paste0("LEfSe_boxplot_", tax_level, "_fdr.pdf")),
      file.path(output_dir, paste0("LEfSe_boxplot_", tax_level, "_fdr.png")),
      paste0(tax_level, "-level taxa"),
      plot_width = ifelse(tax_level == "Genus", 9, plot_width)
    )
  }
  
  cat("FDR-corrected statistical results saved:", nrow(stat_summary_fdr), "taxa\n")
  
  cat("\nCompleted", tax_level, "level boxplot\n")
}

# ============================================================================
# Analysis Complete
# ============================================================================

cat("\n", rep("=", 70), "\n", sep="")
cat("All boxplots generated successfully!\n")
cat("Output directory:", output_dir, "\n")
cat("\nGenerated files:\n")
cat("  - LEfSe_boxplot_Phylum.pdf/png\n")
cat("  - LEfSe_boxplot_Genus.pdf/png\n")
cat("  - LEfSe_boxplot_Species.pdf/png\n")
cat("  - LEfSe_boxplot_*_fdr.pdf/png\n")
cat("  - statistical_results_*.csv (Wilcoxon p-values: Model vs Control, HQQD vs Model, Model vs LVX)\n")
cat("  - statistical_results_*_fdr.csv (BH-corrected callback subset)\n")
cat("\nInterpretation notes:\n")
cat("  - Selection: Taxa with significant reversal trend (Model vs Control & HQQD vs Model, p<0.05), sorted by LDA\n")
cat("  - Y-axis labels are colored by enriched group\n")
cat("  - X-axis shows log10-transformed relative abundance\n")
cat("  - Right panel shows pairwise comparisons with significance levels\n")
cat("  - Boxplots show median and quartiles (no outliers displayed)\n")
cat("  - Each taxon shows 4 boxplots (one per group) side by side\n")
cat("  - Left: boxplot panel, Middle: statistical annotation panel, Right: legend\n")
cat(rep("=", 70), "\n\n", sep="")
