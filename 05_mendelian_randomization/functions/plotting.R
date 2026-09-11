# ==============================================================================
# ==============================================================================
# ==============================================================================

make_ci_with_caps <- function(est, lower, upper, sizes, xlim, pch, gp, t_height, nudge_y) {
  
  grob_list <- list()
  
  if (is.null(sizes) || length(sizes) == 0) {
    sizes <- rep(0.4, length(est))
  } else if (length(sizes) == 1) {
    sizes <- rep(sizes, length(est))
  }
  
  for (i in seq_along(est)) {
    if (is.na(est[i])) {
      grob_list[[i]] <- nullGrob()
      next
    }
    
    lower_out <- !is.na(lower[i]) && lower[i] < xlim[1]
    upper_out <- !is.na(upper[i]) && upper[i] > xlim[2]
    
    draw_lower <- if (lower_out) xlim[1] else lower[i]
    draw_upper <- if (upper_out) xlim[2] else upper[i]
    
    x_range <- xlim[2] - xlim[1]
    est_x <- (est[i] - xlim[1]) / x_range
    lower_x <- (draw_lower - xlim[1]) / x_range
    upper_x <- (draw_upper - xlim[1]) / x_range
    
    est_x <- max(0, min(1, est_x))
    
    cap_height <- 0.15
    
    grobs <- gList()
    
    grobs <- gList(grobs, linesGrob(
      x = c(lower_x, upper_x), y = c(0.5, 0.5),
      gp = gpar(col = "black", lwd = 1.0)
    ))
    
    if (lower_out) {
      grobs <- gList(grobs, linesGrob(
        x = c(lower_x + 0.04, lower_x, lower_x + 0.04),
        y = c(0.5 + cap_height, 0.5, 0.5 - cap_height),
        gp = gpar(col = "black", lwd = 1.0)
      ))
    } else {
      grobs <- gList(grobs, linesGrob(
        x = c(lower_x, lower_x),
        y = c(0.5 - cap_height, 0.5 + cap_height),
        gp = gpar(col = "black", lwd = 1.0)
      ))
    }
    
    if (upper_out) {
      grobs <- gList(grobs, linesGrob(
        x = c(upper_x - 0.04, upper_x, upper_x - 0.04),
        y = c(0.5 + cap_height, 0.5, 0.5 - cap_height),
        gp = gpar(col = "black", lwd = 1.0)
      ))
    } else {
      grobs <- gList(grobs, linesGrob(
        x = c(upper_x, upper_x),
        y = c(0.5 - cap_height, 0.5 + cap_height),
        gp = gpar(col = "black", lwd = 1.0)
      ))
    }
    
    point_size <- if(length(sizes) >= i && !is.na(sizes[i])) sizes[i] else 0.4
    grobs <- gList(grobs, pointsGrob(
      x = est_x, y = 0.5,
      pch = 19,
      size = unit(point_size * 1.2, "char"),
      gp = gpar(col = "#377EB8", fill = "#377EB8")
    ))
    
    grob_list[[i]] <- gTree(children = grobs, vp = viewport(xscale = c(0, 1)))
  }
  
  return(grob_list)
}

#' 
#' @export
plot_forestplot_smart <- function(data_file, output_file, threshold = 100, font_family = "sans") {
  
  require(data.table)
  require(dplyr)
  require(forestploter)
  require(grid)
  
  message("Processing file: ", data_file)
  
  raw_data <- fread(data_file)
  
  data <- raw_data %>%
    filter(grepl("Inverse variance weighted", method, ignore.case = TRUE) | method == "IVW") %>%
    arrange(desc(or))
  
  if(nrow(data) == 0) stop("No IVW data found")
  
  N <- nrow(data)
  mid_point <- ceiling(N / 2)
  
  data_left <- data[1:mid_point, ]
  data_right_raw <- if (mid_point < N) data[(mid_point + 1):N, ] else data.frame()
  
  rows_needed <- mid_point - nrow(data_right_raw)
  if (rows_needed > 0 && nrow(data_right_raw) > 0) {
    empty_rows <- data_right_raw[1:rows_needed, ]
    empty_rows[] <- NA
    data_right <- bind_rows(data_right_raw, empty_rows)
  } else if (nrow(data_right_raw) == 0) {
    data_right <- data_left
    data_right[] <- NA
  } else {
    data_right <- data_right_raw
  }
  
  
  combined_df <- data.frame(
    Exposure_L = ifelse(nchar(as.character(data_left$exposure)) > 28, 
                        paste0(substr(data_left$exposure, 1, 26), "..."), 
                        as.character(data_left$exposure)),
    nsnp_L = as.character(data_left$nsnp),
    Method_L = "IVW",
    CI_L = paste(rep(" ", 18), collapse = " "),
    pval_L = ifelse(data_left$pval < 0.001, "<0.001", sprintf("%.3f", data_left$pval)),
    OR_L = sprintf("%.2f (%.2f to %.2f)", data_left$or, data_left$or_lci95, data_left$or_uci95),
    
    Sep = "|",
    
    Exposure_R = ifelse(is.na(data_right$exposure), "", 
                        ifelse(nchar(as.character(data_right$exposure)) > 28,
                               paste0(substr(data_right$exposure, 1, 26), "..."),
                               as.character(data_right$exposure))),
    nsnp_R = ifelse(is.na(data_right$nsnp), "", as.character(data_right$nsnp)),
    Method_R = ifelse(is.na(data_right$exposure), "", "IVW"),
    CI_R = paste(rep(" ", 18), collapse = " "),
    pval_R = ifelse(is.na(data_right$pval), "", 
                    ifelse(data_right$pval < 0.001, "<0.001", sprintf("%.3f", data_right$pval))),
    OR_R = ifelse(is.na(data_right$or), "",
                  sprintf("%.2f (%.2f to %.2f)", data_right$or, data_right$or_lci95, data_right$or_uci95)),
    
    stringsAsFactors = FALSE
  )
  
  colnames(combined_df) <- c(
    "exposure", "nsnp", "method", " ", "pval", "OR (95% CI)",
    "|",
    "exposure ", "nsnp ", "method ", "  ", "pval ", "OR (95% CI) "
  )
  
  est_l <- data_left$or
  low_l <- data_left$or_lci95
  upp_l <- data_left$or_uci95
  
  est_r <- ifelse(is.na(data_right$or), NA, data_right$or)
  low_r <- ifelse(is.na(data_right$or_lci95), NA, data_right$or_lci95)
  upp_r <- ifelse(is.na(data_right$or_uci95), NA, data_right$or_uci95)
  
  theme <- forest_theme(
    base_size = 8,
    base_family = font_family,
    colhead = list(
      fg_params = list(hjust = 0.5, x = 0.5, fontface = "bold")
    ),
    core = list(
      bg_params = list(fill = c("white", "#F5F5F5")),
      fg_params = list(hjust = 0.5, x = 0.5)
    ),
    refline_gp = gpar(col = "gray50", lty = 2, lwd = 1)
  )
  
  xlim_val <- c(0.5, 1.5)
  ticks_val <- c(0.5, 1.0, 1.5)
  sizes_l <- rep(0.4, length(est_l))
  sizes_r <- rep(0.4, length(est_r))
  
  p <- forest(
    data = combined_df,
    est = list(est_l, est_r),
    lower = list(low_l, low_r),
    upper = list(upp_l, upp_r),
    sizes = list(sizes_l, sizes_r),
    ci_column = c(4, 11),
    ref_line = 1,
    xlim = xlim_val,
    ticks_at = ticks_val,
    xlab = "",
    fn_ci = make_ci_with_caps,
    theme = theme
  )
  
  p <- add_border(p, part = "header", row = 1, where = "bottom", gp = gpar(lwd = 2))
  
  p_height <- max(6, nrow(combined_df) * 0.25 + 1.5)
  p_width <- 20
  
  cairo_pdf(output_file, width = p_width, height = p_height, family = font_family)
  print(p)
  dev.off()
  
  message("Plotting complete: ", output_file)
}
