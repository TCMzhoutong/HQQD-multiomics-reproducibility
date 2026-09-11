save_reverse_boxplot <- function(plot_data_combined, output_file_pdf, output_file_png, plot_title,
                                 plot_width = 8.5,
                                 plot_subtitle = NULL,
                                 output_file_svg = NULL,
                                 label_size = 3.5,
                                 title_size = 14,
                                 axis_text_size = 10,
                                 axis_title_size = 12,
                                 legend_title_size = 11,
                                 legend_text_size = 10,
                                 height_multiplier = 1,
                                 plot_height_override = NULL) {
  if (is.null(plot_data_combined) || nrow(plot_data_combined) == 0 || length(levels(plot_data_combined$Taxon)) == 0) {
    return(FALSE)
  }

  n_features <- length(levels(plot_data_combined$Taxon))
  if (n_features >= 2) {
    alternating_bg <- data.frame(
      ymin = seq(1.5, n_features - 0.5, by = 1),
      ymax = seq(2.5, n_features + 0.5, by = 1)
    )
    alternating_bg <- alternating_bg[seq(1, nrow(alternating_bg), by = 2), , drop = FALSE]
  } else {
    alternating_bg <- data.frame(ymin = numeric(0), ymax = numeric(0))
  }

  label_data <- data.frame(
    Taxon = factor(levels(plot_data_combined$Taxon), levels = levels(plot_data_combined$Taxon)),
    Label = stringr::str_wrap(levels(plot_data_combined$Taxon), width = 20)
  )

  p_labels <- ggplot(label_data, aes(x = 1, y = Taxon, label = Label)) +
    geom_text(hjust = 1, size = label_size, fontface = "plain", color = "black", lineheight = 0.9) +
    scale_y_discrete(limits = levels(plot_data_combined$Taxon)) +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    theme_void(base_family = "sans") +
    theme(
      plot.margin = margin(10, 4, 10, 5),
      panel.border = element_blank()
    )

  p_left <- ggplot(plot_data_combined, aes(x = Log_Abundance, y = Taxon, color = Group, fill = Group)) +
    geom_rect(data = alternating_bg, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      fill = "gray85", alpha = 0.5, inherit.aes = FALSE) +
    geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA, width = 0.35, alpha = 0.3, linewidth = 0.25) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    scale_y_discrete(limits = levels(plot_data_combined$Taxon)) +
    guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE)) +
    labs(
      x = expression(Log[10] ~ "Relative Abundance"),
      y = NULL,
      title = plot_title,
      subtitle = plot_subtitle
    ) +
    theme_bw(base_family = "sans") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = axis_text_size),
      axis.title.x = element_text(size = axis_title_size, face = "bold"),
      legend.position = "none",
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(color = "black", linewidth = 0.25),
      axis.line.y.left = element_line(color = "black", linewidth = 0.25)
    ) +
    geom_segment(aes(x = -Inf, xend = Inf, y = Inf, yend = Inf), color = "black", linewidth = 0.25, inherit.aes = FALSE)

  p_middle <- ggplot(plot_data_combined, aes(x = 1, y = Taxon)) +
    geom_rect(data = alternating_bg, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      fill = "gray85", alpha = 0.5, inherit.aes = FALSE) +
    scale_y_discrete(limits = levels(plot_data_combined$Taxon)) +
    theme_void(base_family = "sans") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(10, 10, 10, 5), panel.border = element_blank()) +
    geom_segment(aes(x = Inf, xend = Inf, y = -Inf, yend = Inf), color = "black", linewidth = 0.25, inherit.aes = FALSE) +
    geom_segment(aes(x = -Inf, xend = Inf, y = Inf, yend = Inf), color = "black", linewidth = 0.25, inherit.aes = FALSE) +
    geom_segment(aes(x = -Inf, xend = Inf, y = -Inf, yend = -Inf), color = "black", linewidth = 0.25, inherit.aes = FALSE) +
    coord_cartesian(xlim = c(0, 10), clip = "off")

  plot_stat_results <- list()
  for (feature in levels(plot_data_combined$Taxon)) {
    fd <- plot_data_combined %>% filter(Taxon == feature)
    mc_v   <- fd$Log_Abundance[fd$Group == "Model"]
    ctrl_v <- fd$Log_Abundance[fd$Group == "Control"]
    hqqd_v <- fd$Log_Abundance[fd$Group == "HQQD"]
    lvx_v  <- fd$Log_Abundance[fd$Group == "LVX"]
    plot_stat_results[[feature]] <- list(
      p_Model_vs_Control = tryCatch(wilcox.test(mc_v, ctrl_v, exact = FALSE)$p.value, error = function(e) NA),
      p_HQQD_vs_Model    = tryCatch(wilcox.test(hqqd_v, mc_v, exact = FALSE)$p.value, error = function(e) NA),
      p_Model_vs_LVX     = tryCatch(wilcox.test(mc_v, lvx_v, exact = FALSE)$p.value, error = function(e) NA)
    )
  }

  comparisons_fixed <- list(
    list(group1 = "Model", group2 = "Control", p_key = "p_Model_vs_Control"),
    list(group1 = "HQQD", group2 = "Model", p_key = "p_HQQD_vs_Model"),
    list(group1 = "Model", group2 = "LVX", p_key = "p_Model_vs_LVX")
  )
  group_positions <- c("HQQD" = 0, "LVX" = 1, "Model" = 2, "Control" = 3)
  dodge_width <- 0.7
  n_groups <- 4

  for (feature_name in levels(plot_data_combined$Taxon)) {
    feature_idx <- which(levels(plot_data_combined$Taxon) == feature_name)
    stat_res <- plot_stat_results[[feature_name]]
    current_x <- 2

    for (comp in comparisons_fixed) {
      p_val <- stat_res[[comp$p_key]]
      if (is.na(p_val) || p_val >= 0.05) {
        current_x <- current_x + 1.2
        next
      }

      sig_label <- dplyr::case_when(p_val < 0.001 ~ "***", p_val < 0.01 ~ "**", TRUE ~ "*")
      g1_idx <- group_positions[comp$group1]
      g2_idx <- group_positions[comp$group2]
      y1 <- feature_idx + (g1_idx - (n_groups - 1) / 2) * dodge_width / n_groups
      y2 <- feature_idx + (g2_idx - (n_groups - 1) / 2) * dodge_width / n_groups
      star_spacing <- dplyr::case_when(sig_label == "***" ~ 2.6, sig_label == "**" ~ 2.0, TRUE ~ 1.2)

      p_middle <- p_middle +
        annotate("segment", x = current_x, xend = current_x, y = y1, yend = y2, color = "gray30", linewidth = 0.25) +
        annotate("text", x = current_x + 0.3, y = (y1 + y2) / 2, label = sig_label,
                 size = 3.5, fontface = "bold", color = "gray20", hjust = 0)

      current_x <- current_x + star_spacing
    }
  }

  p_legend_temp <- ggplot(plot_data_combined, aes(x = Log_Abundance, y = Taxon, color = Group, fill = Group)) +
    geom_boxplot(alpha = 0.3) +
    scale_color_manual(values = group_colors, name = "Group") +
    scale_fill_manual(values = group_colors, name = "Group") +
    guides(color = guide_legend(reverse = TRUE, override.aes = list(fill = group_colors, alpha = 0.3)), fill = "none") +
    theme_minimal(base_family = "sans") +
    theme(legend.title = element_text(size = legend_title_size, face = "bold"), legend.text = element_text(size = legend_text_size), legend.key.size = unit(0.8, "cm"))

  p_legend <- as_ggplot(get_legend(p_legend_temp))
  p_combined <- p_labels + p_left + p_middle + p_legend + plot_layout(widths = c(1.65, 5, 1, 0.8))
  plot_height <- if (is.null(plot_height_override)) {
    (n_features * 0.55 + 2) * height_multiplier
  } else {
    plot_height_override
  }

  ggsave(filename = output_file_pdf, plot = p_combined, width = plot_width, height = plot_height, units = "in", dpi = 300)
  ggsave(filename = output_file_png, plot = p_combined, width = plot_width, height = plot_height, units = "in", dpi = 300)
  if (!is.null(output_file_svg)) {
    if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(file = output_file_svg, width = plot_width, height = plot_height)
    } else {
      grDevices::svg(filename = output_file_svg, width = plot_width, height = plot_height, onefile = FALSE)
    }
    print(p_combined)
    grDevices::dev.off()
  }
  TRUE
}
