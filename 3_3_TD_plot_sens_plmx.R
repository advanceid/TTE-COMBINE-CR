
###20260922
### Primary Outcome: 21-day ACM
### Secondary Outcomes: 7-day & 21-day ACM


# Clear
rm(list = ls())

library(ggplot2)
library(patchwork)

# ---- 1. Data --------------------------------------------------------------
df_bl <- data.frame(
  outcome  = c("14-day All-cause Mortality", "7-day All-cause Mortality", "21-day All-cause Mortality"),
  estimate = c(-1.7,  -2.6,  -3.9),
  lower    = c(-11.4, -11.6, -14),
  upper    = c(8.1,   6.3,   6.2)
)

df_lm <- data.frame(
  outcome  = c("14-day All-cause Mortality", "7-day All-cause Mortality", "21-day All-cause Mortality"),
  estimate = c(0.9,  -1.5,  -5),
  lower    = c(-11.4, -12, -18.1),
  upper    = c(13.3,   9,  8.1)
)

# ---- 2. Shared x-axis range (so both panels are directly comparable) ------
x_min    <- -25
x_max    <- 30
x_breaks <- seq(x_min, x_max, 10)
label_x  <- x_max + 6     # x-position of the CI text column

# ---- 3. Palette -------------------------------------------------------------
col_point   <- "#1B4965"   # deep navy for point estimates / CI lines
col_line    <- "#1B4965"
col_zero    <- "#8C8C8C"
col_band_l  <- "#EEF3F6"   # left half shading (favors combo)
col_band_r  <- "#FFFFFF"   # right half shading (favors mono)
col_text    <- "#1A1A1A"
col_axis    <- "#4D4D4D"

# ---- 4. Reusable forest-plot builder ---------------------------------------
# NOTE: the row-position logic now scales to however many outcomes are in
# `df`, instead of a hard-coded 2-row lookup table that silently dropped the
# 3rd row (the original bug that left "28-day" with no y-coordinate).
build_forest_plot <- function(df, header_label, show_axis = TRUE, show_arrows = TRUE) {

  n <- nrow(df)
  df$label <- sprintf("%.1f%% (%.1f to %.1f)", df$estimate, df$lower, df$upper)

  row_spacing <- 1.6
  df$y        <- rev(seq_len(n)) * row_spacing        # first row = top
  y_header    <- max(df$y) + 2.3
  y_top       <- max(df$y) + 1.35
  y_bg_bottom <- min(df$y) - 1.15
  y_axis      <- y_bg_bottom
  y_axis_text <- y_axis - 0.55
  y_arrow     <- y_axis_text - 0.95
  y_favor_text<- y_arrow - 0.65

  p <- ggplot() +

    # shaded background halves, split at 0
    annotate("rect", xmin = x_min, xmax = 0, ymin = y_bg_bottom, ymax = y_top,
             fill = col_band_l, colour = NA) +
    annotate("rect", xmin = 0, xmax = x_max, ymin = y_bg_bottom, ymax = y_top,
             fill = col_band_r, colour = NA) +

    # dashed reference line at 0
    geom_segment(aes(x = 0, xend = 0, y = y_bg_bottom, yend = y_top),
                 colour = col_zero, linetype = "dashed", linewidth = 0.5) +

    # header row, centered over the shaded band
    annotate("text", x = (x_min + x_max) / 2, y = y_header, label = header_label,
             fontface = "bold", size = 4.6, colour = col_text) +
    annotate("text", x = label_x + 7, y = y_header,
             label = "Risk Difference\n(95% CI)",
             hjust = 0.5, fontface = "bold", size = 4.1, lineheight = 0.95,
             colour = col_text) +

    # CI whiskers + point estimates
    geom_segment(data = df, aes(x = lower, xend = upper, y = y, yend = y),
                 colour = col_line, linewidth = 0.8) +
    geom_segment(data = df, aes(x = lower, xend = lower, y = y - 0.14, yend = y + 0.14),
                 colour = col_line, linewidth = 0.8) +
    geom_segment(data = df, aes(x = upper, xend = upper, y = y - 0.14, yend = y + 0.14),
                 colour = col_line, linewidth = 0.8) +
    geom_point(data = df, aes(x = estimate, y = y),
               size = 4.2, shape = 21, colour = col_point,
               fill = col_point, stroke = 0.8) +

    # row labels on the left, right-aligned into their own column so they
    # never bleed into the plot area (previously left-aligned text could
    # run rightward across the CI whiskers)
    annotate("text", x = x_min - 3, y = df$y, label = df$outcome,
             hjust = 1, fontface = "bold", size = 4.3, colour = col_text) +

    # treatment-difference values on the right
    geom_text(data = df, aes(x = label_x + 7, y = y, label = label),
              hjust = 0.5, size = 4, colour = col_text) +

    coord_cartesian(clip = "off") +
    scale_x_continuous(limits = c(x_min - 34, label_x + 16)) +
    theme_void(base_size = 13) +
    theme(
      text             = element_text(family = "sans"),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(10, 20, 5, 15)
    )

  y_limit_lo <- y_bg_bottom - 0.3   # sensible default if no axis/arrows

  if (show_axis) {
    p <- p +
      annotate("segment", x = x_min, xend = x_max, y = y_axis, yend = y_axis,
               colour = col_axis, linewidth = 0.5) +
      annotate("text", x = x_breaks, y = y_axis_text,
               label = paste0(x_breaks, "%"), size = 3.7, colour = col_axis)
    y_limit_lo <- y_axis_text - 0.5
  }

  if (show_arrows) {
    p <- p +
      annotate("segment", x = -1.5, xend = x_min + 1, y = y_arrow, yend = y_arrow,
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               colour = col_text, linewidth = 0.6) +
      annotate("segment", x = 1.5, xend = x_max - 1, y = y_arrow, yend = y_arrow,
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               colour = col_text, linewidth = 0.6) +
      annotate("text", x = (x_min + 1) / 2, y = y_favor_text,
               label = "Favors Polymyxin-based\nCombination therapy", size = 3.7,
               fontface = "italic", colour = col_text) +
      annotate("text", x = (x_max - 1) / 2, y = y_favor_text,
               label = "Favors Polymyxin Monotherapy", size = 3.7,
               fontface = "italic", colour = col_text)
    y_limit_lo <- y_favor_text - 0.6
  }

  p + scale_y_continuous(limits = c(y_limit_lo, y_header + 0.6), expand = c(0, 0))
}

# ---- 5. Build the two panels ------------------------------------------------
# Each panel now carries its own x-axis scale and its own favors-arrow strip,
# so the Baseline and Landmark cohorts can each be read on their own.
p_bl <- build_forest_plot(df_bl, "Intention-to-treat population",
                           show_axis = TRUE, show_arrows = TRUE)

p_lm <- build_forest_plot(df_lm, "Per-protocol population",
                           show_axis = TRUE, show_arrows = TRUE)

# ---- 6. Combine into a single figure ----------------------------------------
combined <- (p_bl / p_lm) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", colour = NA))
  )

print(combined)

# ---- 7. Save ------------------------------------------------------------------
ggsave("sens_forest_plot_plmx_21.png", combined, width = 12, height = 7, dpi = 300, bg = "white")

###-----------------------------------------------------------------------------
### ^.^ ###
