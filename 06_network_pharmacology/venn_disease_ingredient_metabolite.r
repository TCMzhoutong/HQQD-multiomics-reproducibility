# ==============================================================================
# ==============================================================================

library(ggVennDiagram)
library(ggplot2)
library(readr)

if (.Platform$OS.type == "windows") {
  windowsFonts(sans = windowsFont("sans"))
}

disease_df <- readr::read_csv("05.disease_targets/01_disease_targets.csv", show_col_types = FALSE)
ingredient_df <- readr::read_csv("05.disease_targets/02_ingredient_targets.csv", show_col_types = FALSE)
metabolite_df <- readr::read_csv("05.disease_targets/03_metabolite_targets.csv", show_col_types = FALSE)

disease_targets <- unique(disease_df$targets)
ingredient_targets <- unique(ingredient_df$targets)
metabolite_targets <- unique(metabolite_df$targets)

venn_list <- list(
  Disease = disease_targets,
  Ingredient = ingredient_targets,
  Metabolite = metabolite_targets
)

cat("========== Dataset summary ==========\n")
cat("Disease targets: ", length(disease_targets), "\n")
cat("Ingredient targets: ", length(ingredient_targets), "\n")
cat("Metabolite targets: ", length(metabolite_targets), "\n\n")

intersect_DI <- intersect(disease_targets, ingredient_targets)
intersect_DM <- intersect(disease_targets, metabolite_targets)
intersect_IM <- intersect(ingredient_targets, metabolite_targets)
intersect_DIM <- intersect(intersect_DI, metabolite_targets)

cat("========== Intersection summary ==========\n")
cat("Disease ∩ Ingredient:", length(intersect_DI), "\n")
cat("Disease ∩ Metabolite:", length(intersect_DM), "\n")
cat("Ingredient ∩ Metabolite:", length(intersect_IM), "\n")
cat("Disease ∩ Ingredient ∩ Metabolite:", length(intersect_DIM), "\n\n")

if (!dir.exists("figures")) {
  dir.create("figures")
  cat("Created the figures directory\n")
}

# ==============================================================================
# ==============================================================================
cat("Generating the Venn diagram...\n")

venn <- Venn(venn_list)
data <- process_data(venn)


fill_colors <- c("1" = "#E64B3580",
                 "2" = "#4DBBD580",
                 "3" = "#00A08780",
                 "1/2" = "#F39B7F80",
                 "1/3" = "#8491B480",
                 "2/3" = "#91D1C280",
                 "1/2/3" = "#DC000080")

edge_colors <- c("1" = "#000000",
                 "2" = "#000000",
                 "3" = "#000000")

region_data <- venn_regionlabel(data)
total_unique <- sum(region_data$count)
region_data$percentage <- paste0(round(region_data$count / total_unique * 100, 1), "%")
region_data$label <- paste0(region_data$count, "\n(", region_data$percentage, ")")

p_custom <- ggplot() +
  geom_polygon(aes(X, Y, fill = id, group = id), 
               data = venn_regionedge(data),
               show.legend = FALSE) +
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(data), 
            linewidth = 0.8,
            show.legend = FALSE) +
  geom_text(aes(X, Y, label = name), 
            data = venn_setlabel(data),
            fontface = "bold",
            size = 6) +
  geom_text(aes(X, Y, label = label), 
            data = region_data,
            size = 5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = edge_colors) +
  coord_equal() +
  theme_void(base_family = "sans") +
  theme(text = element_text(family = "sans"))

print(p_custom)

ggsave("figures/venn_disease_ingredient_metabolite.png", 
       plot = p_custom, 
       width = 10, 
       height = 8, 
       dpi = 300,
       bg = "white")

ggsave("figures/venn_disease_ingredient_metabolite.pdf", 
       plot = p_custom, 
       width = 10, 
       height = 8,
       device = "pdf")

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave("figures/venn_disease_ingredient_metabolite.svg",
         plot = p_custom,
         width = 10,
         height = 8,
         device = svglite::svglite)
} else {
  grDevices::svg("figures/venn_disease_ingredient_metabolite.svg",
                 width = 10,
                 height = 8,
                 onefile = FALSE)
  print(p_custom)
  grDevices::dev.off()
}

# ==============================================================================
# ==============================================================================
cat("\nExporting intersection data...\n")

if (length(intersect_DIM) > 0) {
  write.csv(data.frame(targets = sort(intersect_DIM)), 
            "05.disease_targets/04_disease_ingredient_metabolite_target.csv", 
            row.names = FALSE)
  cat("Three-way intersection saved: 05.disease_targets/04_disease_ingredient_metabolite_target.csv\n")
  cat("Intersection contains ", length(intersect_DIM), " targets\n")
} else {
  cat("Warning: the three-way intersection is empty\n")
}

# ==============================================================================
# ==============================================================================
cat("\n========== Analysis complete ==========\n")
cat("Figure files:\n")
cat("  - figures/venn_disease_ingredient_metabolite.png\n")
cat("  - figures/venn_disease_ingredient_metabolite.pdf\n")
cat("\nIntersection data:\n")
cat("  - 05.disease_targets/04_disease_ingredient_metabolite_target.csv\n")
cat("  - Contains ", length(intersect_DIM), " shared targets\n")
