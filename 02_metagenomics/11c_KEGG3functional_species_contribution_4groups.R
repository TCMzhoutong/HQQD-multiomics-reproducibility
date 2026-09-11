library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(gridExtra)
library(grid)

# ==============================================================================
# ==============================================================================
work_dir <- getwd()
input_dir <- file.path(work_dir, "11_KEGG3functional_species_contribution_fdr")
output_dir <- file.path(work_dir, "11b_KEGG3functional_species_contribution_fdr")
map_file <- file.path(work_dir, "kegg.kegg_pathway_level3.gene_list.txt")
TOP_N_GENUS <- 10

if(!dir.exists(output_dir)) dir.create(output_dir)

# ==============================================================================
# ==============================================================================
mapping <- read_tsv(map_file, show_col_types = FALSE) %>%
  select(kegg_pathway_id, kegg_pathway) %>%
  distinct()

# ==============================================================================
# ==============================================================================
csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) {
  stop("No contribution CSV files found in: ", input_dir)
}

process_one_file <- function(fpath) {
  func_id <- str_remove(basename(fpath), "\\.csv$")
  df <- read_csv(fpath, show_col_types = FALSE)
  
  cols_keep <- c("Genus", grep("^[KMYZ]\\d+$", names(df), value = TRUE))
  df <- df[, cols_keep]
  
  df_long <- pivot_longer(df, cols = -Genus, names_to = "Sample", values_to = "Abundance")
  
  df_long <- df_long %>%
    mutate(GroupLabel = case_when(
      str_starts(Sample, "K") ~ "C", # Control
      str_starts(Sample, "M") ~ "M", # Model
      str_starts(Sample, "Y") ~ "L", # LVX
      str_starts(Sample, "Z") ~ "H"  # HQQD
    )) %>%
    filter(!is.na(GroupLabel))
  
  df_long <- df_long %>%
    group_by(Sample) %>%
    mutate(RelAb = Abundance / sum(Abundance) * 100) %>%
    ungroup()
  
  df_grouped <- df_long %>%
    group_by(GroupLabel, Genus) %>%
    summarise(MeanRelAb = mean(RelAb), .groups = "drop")
  
  df_grouped <- df_grouped %>%
    group_by(GroupLabel) %>%
    mutate(FinalRelAb = MeanRelAb / sum(MeanRelAb) * 100) %>%
    ungroup()

  df_plot <- df_grouped %>%
    transmute(GroupLabel, Genus, RelAb = FinalRelAb, FuncID = func_id)
  
  return(df_plot)
}

plot_data_list <- lapply(csv_files, process_one_file)
all_data <- bind_rows(plot_data_list)

special_taxa <- c("Unclassified", "Unassigned")
global_rank_df <- all_data %>%
  filter(!Genus %in% c(special_taxa, "Others")) %>%
  group_by(Genus) %>%
  summarise(Total = sum(RelAb, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total))

top_genera <- global_rank_df$Genus[1:min(TOP_N_GENUS, nrow(global_rank_df))]
cat("Displayed Genus:", paste(top_genera, collapse = ", "), "\n")

all_data <- all_data %>%
  mutate(Genus_Plot = case_when(
    Genus %in% special_taxa ~ Genus,
    Genus %in% top_genera ~ Genus,
    TRUE ~ "Others"
  )) %>%
  group_by(GroupLabel, FuncID, Genus_Plot) %>%
  summarise(RelAb = sum(RelAb), .groups = "drop")

# ==============================================================================
# ==============================================================================

global_rank_df <- all_data %>%
  filter(!Genus_Plot %in% c("Unclassified", "Unassigned", "Others")) %>%
  group_by(Genus_Plot) %>%
  summarise(Total = sum(RelAb)) %>%
  arrange(desc(Total))

sorted_taxa <- global_rank_df$Genus_Plot
special_exist <- intersect(c("Others", "Unassigned", "Unclassified"), unique(all_data$Genus_Plot))
genus_levels <- rev(c(sorted_taxa, "Others", "Unassigned", "Unclassified"))
genus_levels <- genus_levels[genus_levels %in% unique(all_data$Genus_Plot)]

all_data$Genus_Plot <- factor(all_data$Genus_Plot, levels = genus_levels)

all_data$GroupLabel <- factor(all_data$GroupLabel, levels = c("C", "M", "L", "H"))

pos_map <- c("C" = 1.0, "M" = 1.25, "L" = 1.5, "H" = 1.75)
all_data$x_pos <- pos_map[as.character(all_data$GroupLabel)]

other_special <- c("Others", "Unassigned", "Unclassified")
main_taxa <- setdiff(genus_levels, other_special)

nb_main <- length(main_taxa)
if (nb_main > 0) {
  main_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(nb_main)
  names(main_colors) <- main_taxa
} else {
  main_colors <- c()
}

special_colors <- c(
  "Others" = "#D3D3D3", 
  "Unassigned" = "#A9A9A9", 
  "Unclassified" = "#696969"
)
special_colors <- special_colors[names(special_colors) %in% genus_levels]
genus_colors <- c(main_colors, special_colors)

# ==============================================================================
# ==============================================================================

p_main <- ggplot(all_data, aes(x = x_pos, y = RelAb, fill = Genus_Plot)) +
  geom_bar(stat = "identity", width = 0.2, show.legend = TRUE) +
  facet_wrap(~FuncID, scales = "fixed", strip.position = "top", nrow = 1) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 101), expand = c(0,0)) +
  # Bar Width 0.2 -> Half 0.1. 
  # Limits = [1.0-0.1, 1.75+0.1] = [0.9, 1.85]
  scale_x_continuous(breaks = c(1.0, 1.25, 1.5, 1.75), 
                     labels = c("C", "M", "L", "H"),
                     limits = c(0.9, 1.85),
                     expand = c(0, 0)) + 
  scale_fill_manual(values = genus_colors) +
  theme_classic(base_family = "sans") +
  theme(
    strip.background = element_rect(fill = "#E0E0E0", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(color="black", size=10, face="bold"),
    axis.ticks.x = element_line(),
    axis.line.x = element_line(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(y = "Relative contribution (%)")

# 1. Genus Legend
tm_legend <- theme_void(base_family = "sans") +
  theme(
    legend.position = "left",
    legend.justification = "left",
    legend.title = element_text(face="bold", size=10, hjust=0),
    legend.text = element_text(size=9),
    legend.key.size = unit(0.5, "cm"),
    plot.margin = margin(0,0,0,0)
  )

p_legend_genus <- ggplot(all_data, aes(x=GroupLabel, y=RelAb, fill=Genus_Plot)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = genus_colors) +
  tm_legend +
  guides(fill = guide_legend(title = "Genus", ncol = 1, order = 1))
legend_genus <- cowplot::get_legend(p_legend_genus)

# 2. KEGG Pathway Text
func_list <- mapping %>%
  filter(kegg_pathway_id %in% unique(all_data$FuncID)) %>%
  mutate(text = stringr::str_wrap(paste0(kegg_pathway_id, ": ", kegg_pathway), width = 34))

desc_text <- paste(func_list$text, collapse = "\n")
n_lines <- max(sum(stringr::str_count(func_list$text, "\n") + 1), 1)

p_text <- ggplot() + 
  annotate("text", x = 0, y = n_lines, label = desc_text, hjust = 0, vjust = 1, size = 3, lineheight = 1.2) +
  labs(title = "KEGG pathway") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(face="bold", size=10, hjust=0, margin=margin(b=5)),
    plot.margin = margin(t=0, l=0, r=0, b=0)
  ) +
  xlim(0, 10) + 
  ylim(0, n_lines) +
  coord_cartesian(clip = "off")

if(requireNamespace("cowplot", quietly = TRUE)) {
  library(cowplot)
  
  h_genus <- 1 + length(genus_levels) * 1.1
  h_text <- 1 + n_lines * 1.1
  
  right_column <- plot_grid(
    legend_genus,
    p_text,
    NULL, 
    ncol = 1,
    rel_heights = c(h_genus, h_text, 1),
    align = "v", axis = "l"
  )
  
  final_plot <- plot_grid(
    p_main + theme(legend.position = "none"), 
    right_column,
    ncol = 2,
    rel_widths = c(0.8, 0.2)
  )
  
  ggsave(file.path(output_dir, "KEGG3_functional_species_contribution_4groups.pdf"), 
         final_plot, width = 12, height = 6, bg = "white")
  ggsave(file.path(output_dir, "KEGG3_functional_species_contribution_4groups.png"), 
         final_plot, width = 12, height = 6, dpi = 300, bg = "white")
  
} else {
  ggsave(file.path(output_dir, "KEGG3_functional_species_contribution_4groups_simple.pdf"), 
         p_main, width = 12, height = 6, bg = "white")
}

print("Four-group plotting complete. Files were saved to 11b_KEGG3functional_species_contribution_fdr.")
