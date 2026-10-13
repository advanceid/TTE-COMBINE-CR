# ==================================================================
# Forest plot: Treatment Difference (%) in All-cause Mortality
# Combination therapy vs Monotherapy
#
# Layout: two side-by-side columns (Intention-to-treat population, Per-protocol population),
# each sharing the same three rows (14-day, 7-day, 21-day ACM).
#
# Theme: shaded quadrant background, colored "Favors" arrows that
# start exactly at x = 0, percent-formatted x-axis, sample sizes in
# the column headers, centered value-column header.
# ==================================================================
###20260904
### Primary Outcome: 14-day ACM
### Secondary Outcome: 7-day & 21-day ACM


# Clear
rm(list = ls())

#
library(ggplot2)
library(patchwork)   # install.packages("patchwork") if needed

# ---- 1. Data ---------------------------------------------------------
#n_baseline <- 632
#n_landmark <- 370

### Update RD(95% CI) if needed
df <- data.frame(
  cohort   = rep(c("Intention-to-treat population", "Per-protocol population"), each = 3),
  time     = rep(c("14-day All-cause Mortality", "7-day All-cause Mortality", "21-day All-cause Mortality"), times = 2),
  estimate = c(4.6,  1.4,  4.7,     # Baseline: 14d, 7d, 21d
               8.1,  2.2,  4),     # Landmark: 14d, 7d, 21d
  lower    = c(-2.9, -5.3, -3.2,
               -1.3, -5.7, -6.1),
  upper    = c(12.1,  8.1,  12.6,
               17.6,  10.1,  14.1)
)

df$time  <- factor(df$time, levels = rev(c("14-day All-cause Mortality", "7-day All-cause Mortality", "21-day All-cause Mortality")))
df$label <- sprintf("%.1f (%.1f to %.1f)", df$estimate, df$lower, df$upper)

x_limits <- c(-10, 20)

# ---- 2. Palette --------------------------------------------------------
col_combo   <- "#d95f02"   # orange  = favors combination therapy
col_mono    <- "#2c7fb8"   # blue    = favors monotherapy
col_ci      <- "#1F2937"   # dark navy CI bars
bg_combo    <- "#FCEFE1"   # pale orange background (x < 0)
bg_mono     <- "#E9F1F8"   # pale blue background   (x > 0)

# ---- 3. Shared theme ----------------------------------------------------
base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    axis.ticks          = element_blank()
  )

pct_labels <- function(x) paste0(x, "%")

# ---- 4. Header cell ------------------------------------------------------
make_header <- function(label_text) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = label_text,
             fontface = "bold", size = 3.6, hjust = 0.5, lineheight = 0.95,
             color = "grey15") +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

# ---- 5. Forest-plot cell (shaded background, no y-axis text) ------------
make_forest <- function(cohort_name) {
  d <- subset(df, cohort == cohort_name)
  n_rows <- length(unique(d$time))
  y_bottom <- 0.4
  y_top    <- n_rows + 0.6

  ggplot(d, aes(x = estimate, y = time)) +
    annotate("rect", xmin = x_limits[1], xmax = 0, ymin = y_bottom, ymax = y_top,
             fill = bg_combo, alpha = 0.6) +
    annotate("rect", xmin = 0, xmax = x_limits[2], ymin = y_bottom, ymax = y_top,
             fill = bg_mono, alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                    height = 0.18, linewidth = 0.7, color = col_ci) +
    geom_point(size = 3.2, color = col_mono, fill = col_mono, shape = 21, stroke = 0.4) +
    scale_x_continuous(limits = x_limits, breaks = seq(-10, 20, 5), labels = pct_labels) +
    coord_cartesian(ylim = c(y_bottom, y_top), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    base_theme +
    theme(
      axis.text.y = element_blank(),   # row labels shown once, in their own column
      axis.text.x = element_text(size = 9.5, color = "grey30"),
      plot.margin = margin(4, 8, 4, 2)
    )
}

# ---- 6. Value-label cell (centered to match the header above it) --------
# Add % to Estimate
df$label <- sprintf(
  "%.1f%% (%.1f to %.1f)",
  df$estimate,
  df$lower,
  df$upper
)

#
make_values <- function(cohort_name) {
  d <- subset(df, cohort == cohort_name)

  ggplot(d, aes(x = 0.5, y = time, label = label)) +
    geom_text(hjust = 0.5, size = 3.5, color = "grey15") +
    scale_x_continuous(limits = c(0, 1)) +
    theme_void() +
    theme(plot.margin = margin(4, 4, 4, 4))
}

# ---- 7. Row-label cell (14-day / 7-day/ 21-day, shown once, bold) ------
make_row_labels <- function() {
  d <- subset(df, cohort == "Intention-to-treat population")
  d$label <- sub(" ", "\n", d$time)

  ggplot(d, aes(x = 1, y = time, label = time)) +
    geom_text(hjust = 1, size = 3.6, fontface = "bold", color = "grey15") +
    scale_x_continuous(limits = c(0, 1)) +
    theme_void() +
    theme(plot.margin = margin(4, 2, 4, 4))
}

# ---- 8. Favors-arrow caption cell -----------------------------------------
# Uses the SAME x-scale as the forest panel above it, so 0 lines up exactly
# and each arrow starts precisely at x = 0.
make_caption <- function() {
  ggplot() +
    annotate("segment", x = 0, xend = x_limits[1] * 0.95, y = 0, yend = 0,
             color = col_combo, linewidth = 0.9,
             arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
    annotate("text", x = x_limits[1] * 0.45, y = -0.10, label = "Favors Combination therapy",
             color = col_combo, fontface = "italic", size = 3.0, hjust = 0.5) +
    annotate("segment", x = 0, xend = x_limits[2] * 0.95, y = 0, yend = 0,
             color = col_mono, linewidth = 0.9,
             arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
    annotate("text", x = x_limits[2] * 0.45, y = -0.10, label = "Favors Monotherapy",
             color = col_mono, fontface = "italic", size = 3.0, hjust = 0.5) +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.35, 0.25), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(2, 8, 4, 2))
}

blank <- plot_spacer()

# ---- 9. Assemble the grid -------------------------------------------------
col_widths <- c(0.14, 0.32, 0.16, 0.32, 0.16)

header_row  <- (blank | make_header(paste0("Intention-to-treat population")) |
                 make_header("Treatment Difference\n(95% CI)") |
                 make_header(paste0("Per-protocol population")) |
                 make_header("Treatment Difference\n(95% CI)")) +
                plot_layout(widths = col_widths)

content_row <- (make_row_labels() | make_forest("Intention-to-treat population") | make_values("Intention-to-treat population") |
                 make_forest("Per-protocol population") | make_values("Per-protocol population")) +
                plot_layout(widths = col_widths)

caption_row <- (blank | make_caption() | blank | make_caption() | blank) +
                plot_layout(widths = col_widths)

combined <- (header_row / content_row / caption_row) +
  plot_layout(heights = c(0.11, 0.78, 0.11)) +
  plot_annotation(
    title    = "Treatment Difference (%) in All-cause Mortality",
    subtitle = "Combination therapy vs Monotherapy",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12, hjust = 0.5,
                                    margin = margin(b = 2)),
      plot.subtitle = element_text(size = 12, color = "grey35", hjust = 0.5,
                                    margin = margin(b = 8))
    )
  )

print(combined)

# ---- 10. Save ---------------------------------------------------------------
ggsave("fig2_forest_plot_combo_vs_mono_21.png", combined, width = 15.2, height = 7, dpi = 320, bg = "white")
