library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(grid)

work_dir <- getwd()
input_dir <- file.path(work_dir, "12b_NOGfunctional_species_contribution_top5_fdr")
output_dir <- file.path(work_dir, "12c_NOGfunctional_species_contribution_top5_fdr")
map_file <- file.path(input_dir, "target_eggNOG_NOG_top5.csv")
TOP_N_GENUS <- 10

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

save_plot_svg <- function(filename, plot, width, height) {
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = filename, width = width, height = height)
  } else {
    grDevices::svg(filename = filename, width = width, height = height, onefile = FALSE)
  }
  print(plot)
  grDevices::dev.off()
}

mapping <- read.csv(map_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  select(eggNOG, descrition, eggNOG_Class, Category) %>%
  distinct()

# Keep the publication figure independent of stale upstream description text.
annotation_overrides <- c(
  "COG3385" = "IS4 transposase InsG",
  "COG0534" = "Na+-driven multidrug efflux pump, DinF/NorM/MATE family",
  "COG2972" = "Sensor histidine kinase YesM",
  "COG3279" = "DNA-binding response regulator, LytR/AlgR family",
  "COG0745" = "DNA-binding response regulator, OmpR family"
)
mapping <- mapping %>%
  mutate(
    descrition = if_else(
      eggNOG %in% names(annotation_overrides),
      unname(annotation_overrides[match(eggNOG, names(annotation_overrides))]),
      descrition
    )
  )

target_order <- mapping$eggNOG
csv_files <- file.path(input_dir, paste0(target_order, ".csv"))
csv_files <- csv_files[file.exists(csv_files)]
if (length(csv_files) == 0) {
  stop("No contribution CSV files found in: ", input_dir)
}

process_one_file <- function(fpath) {
  func_id <- str_remove(basename(fpath), "\\.csv$")
  df <- read.csv(fpath, check.names = FALSE, stringsAsFactors = FALSE)
  cols_keep <- c("Genus", grep("^[MZ][0-9]+$", names(df), value = TRUE))
  df <- df[, cols_keep]
  df[, setdiff(names(df), "Genus")] <- lapply(df[, setdiff(names(df), "Genus"), drop = FALSE], as.numeric)

  df_long <- pivot_longer(df, cols = -Genus, names_to = "Sample", values_to = "Abundance") %>%
    mutate(Group = case_when(
      str_starts(Sample, "M") ~ "Model",
      str_starts(Sample, "Z") ~ "HQQD",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Group))

  df_long %>%
    group_by(Sample) %>%
    mutate(
      sample_total = sum(Abundance, na.rm = TRUE),
      RelAb = if_else(sample_total > 0, Abundance / sample_total * 100, 0)
    ) %>%
    ungroup() %>%
    select(-sample_total) %>%
    mutate(FuncID = func_id)
}

all_data <- bind_rows(lapply(csv_files, process_one_file)) %>%
  mutate(FuncID = factor(FuncID, levels = target_order))

special_taxa <- c("Unclassified", "Unassigned")
global_rank_df <- all_data %>%
  filter(!Genus %in% c(special_taxa, "Others")) %>%
  group_by(Genus) %>%
  summarise(Total = sum(RelAb, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total))

top_genera <- global_rank_df$Genus[1:min(TOP_N_GENUS, nrow(global_rank_df))]

all_data <- all_data %>%
  mutate(Genus_Plot = case_when(
    Genus %in% special_taxa ~ Genus,
    Genus %in% top_genera ~ Genus,
    TRUE ~ "Others"
  )) %>%
  group_by(Sample, Group, FuncID, Genus_Plot) %>%
  summarise(RelAb = sum(RelAb), .groups = "drop")

global_rank_df <- all_data %>%
  filter(!Genus_Plot %in% c("Unclassified", "Unassigned", "Others")) %>%
  group_by(Genus_Plot) %>%
  summarise(Total = sum(RelAb), .groups = "drop") %>%
  arrange(desc(Total))

sorted_taxa <- global_rank_df$Genus_Plot
genus_levels <- rev(c(sorted_taxa, "Others", "Unassigned", "Unclassified"))
genus_levels <- genus_levels[genus_levels %in% unique(all_data$Genus_Plot)]
all_data$Genus_Plot <- factor(all_data$Genus_Plot, levels = genus_levels)
all_data$Group <- factor(all_data$Group, levels = c("Model", "HQQD"))

sample_info <- data.frame(Sample = unique(all_data$Sample), stringsAsFactors = FALSE) %>%
  mutate(
    Group = ifelse(str_starts(Sample, "M"), "Model", "HQQD"),
    Num = as.numeric(str_extract(Sample, "[0-9]+"))
  ) %>%
  arrange(factor(Group, levels = c("Model", "HQQD")), Num)

step_size <- 0.33
group_gap <- 0.07
current_x <- 1
last_group <- sample_info$Group[1]
sample_info$x_pos <- NA_real_

for (i in seq_len(nrow(sample_info))) {
  if (sample_info$Group[i] != last_group) {
    current_x <- current_x + group_gap
    last_group <- sample_info$Group[i]
  }
  sample_info$x_pos[i] <- current_x
  current_x <- current_x + step_size
}

all_data <- left_join(all_data, sample_info[, c("Sample", "x_pos")], by = "Sample")

other_special <- c("Others", "Unassigned", "Unclassified")
main_taxa <- setdiff(genus_levels, other_special)
if (length(main_taxa) > 0) {
  main_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(main_taxa))
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

group_colors <- c("Model" = "#EE0000", "HQQD" = "#FF8C00")
combined_colors <- c(genus_colors, group_colors)

strip_data <- sample_info %>%
  mutate(
    ymin = -5,
    ymax = -1,
    xmin = x_pos - 0.165,
    xmax = x_pos + 0.165
  )

p_main <- ggplot() +
  geom_rect(
    data = strip_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Group),
    show.legend = FALSE
  ) +
  geom_bar(
    data = all_data,
    aes(x = x_pos, y = RelAb, fill = Genus_Plot),
    stat = "identity",
    width = 0.3,
    show.legend = TRUE
  ) +
  facet_wrap(~FuncID, scales = "fixed", strip.position = "top", nrow = 1) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(-6, 101), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = combined_colors) +
  theme_classic(base_family = "sans") +
  theme(
    strip.background = element_rect(fill = "#E0E0E0", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    panel.spacing = unit(0.2, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(y = "Relative contribution (%)")

tm_legend <- theme_void(base_family = "sans") +
  theme(
    legend.position = "left",
    legend.justification = "left",
    legend.title = element_text(face = "bold", size = 10, hjust = 0),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.5, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )

p_legend_group <- ggplot(data.frame(g = factor(c("Model", "HQQD"), levels = c("Model", "HQQD"))), aes(x = 1, fill = g)) +
  geom_bar() +
  scale_fill_manual(values = group_colors) +
  tm_legend +
  guides(fill = guide_legend(title = "Group", order = 1, ncol = 1))
legend_group <- cowplot::get_legend(p_legend_group)

p_legend_genus <- ggplot(all_data, aes(x = x_pos, y = RelAb, fill = Genus_Plot)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genus_colors) +
  tm_legend +
  guides(fill = guide_legend(title = "Genus", ncol = 1, order = 2))
legend_genus <- cowplot::get_legend(p_legend_genus)

func_list <- mapping %>%
  filter(eggNOG %in% unique(as.character(all_data$FuncID))) %>%
  mutate(text = stringr::str_wrap(paste0(eggNOG, ": ", descrition), width = 38))

desc_text <- paste(func_list$text, collapse = "\n")
n_lines <- max(sum(stringr::str_count(func_list$text, "\n") + 1), 1)

p_text <- ggplot() +
  annotate("text", x = 0, y = n_lines, label = desc_text, hjust = 0, vjust = 1, size = 3, lineheight = 1.2) +
  labs(title = "eggNOG NOG") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 5)),
    plot.margin = margin(t = 0, l = 0, r = 0, b = 0)
  ) +
  xlim(0, 10) +
  ylim(0, n_lines) +
  coord_cartesian(clip = "off")

library(cowplot)
h_group <- 6
h_genus <- 1 + length(genus_levels) * 1.3
h_text <- 1 + n_lines * 1.1

right_column <- plot_grid(
  legend_group,
  ggplot() + theme_void(),
  legend_genus,
  ggplot() + theme_void(),
  p_text,
  ncol = 1,
  rel_heights = c(h_group, 0.6, h_genus, 0.8, h_text),
  align = "v",
  axis = "l"
)

final_plot <- plot_grid(
  p_main + theme(legend.position = "none"),
  right_column,
  ncol = 2,
  rel_widths = c(0.76, 0.24)
)

base_name <- "NOG_functional_species_contribution_Model_HQQD"
output_width <- 11
output_height <- 6
ggsave(file.path(output_dir, paste0(base_name, ".pdf")), final_plot, width = output_width, height = output_height, bg = "white")
ggsave(file.path(output_dir, paste0(base_name, ".png")), final_plot, width = output_width, height = output_height, dpi = 300, bg = "white")
save_plot_svg(file.path(output_dir, paste0(base_name, ".svg")), final_plot, width = output_width, height = output_height)

summary_table <- all_data %>%
  group_by(FuncID, Group, Genus_Plot) %>%
  summarise(mean_contribution = mean(RelAb, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Group, values_from = mean_contribution, values_fill = 0) %>%
  mutate(delta_HQQD_minus_Model = HQQD - Model) %>%
  arrange(FuncID, desc(abs(delta_HQQD_minus_Model)))

write.csv(
  summary_table,
  file.path(output_dir, "NOG_functional_species_contribution_Model_HQQD_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

message("Displayed Genus: ", paste(rev(setdiff(genus_levels, other_special)), collapse = ", "))
message("Plot and summary saved to ", output_dir)
