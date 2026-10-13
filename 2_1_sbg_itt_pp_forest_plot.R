# =========================================================
# Subgroup analysis with interaction tests
# Target trial emulation: combo vs mono therapy, 14-day mortality
# Subgroups: infection syndrome (VAP/BSI), baseline pathogen (CRA/CRE/CRP)
#
# Expected columns in df (edit names to match your data):
#   id            - patient id
#   treatment     - factor: "Monotherapy" / "Combination therapy"
#   syndrome      - factor: "VAP" / "BSI"
#   pathogen      - factor: "CRA" / "CRE" / "CRP"
#   mortality_14d - 0/1 outcome
#   iptw          - stabilized IPTW weight (skip step 1 if you
#                   already have this)
#   [covariates]  - variables used in your PS model
# =========================================================
# Clear
rm(list = ls())

#
library(tidyverse)
library(dplyr)
library(survey)      # for weighted GLM with correct variance
library(broom)        # tidy() model output
library(ggplot2)
library(purrr)
library(scales)
library(ggtext)   # for styled title/caption (install.packages("ggtext") if needed)
library(patchwork)
library(scales)

# ---------------------------------------------------------
### Load data
df <- readRDS("data/clean data/clean/clean/final/tables_bl.RData")
df_pp <- readRDS("data/clean data/clean/clean/final/tables_lm.RData")

### Clean data
### recode types and org
df <- df %>%
  mutate(types = dplyr::recode(types, `0` = "VAP", `1` = "BSI"))

df <- df %>%
  mutate(org = dplyr::recode(org, `0` = "CRA", `1` = "CRE", `2` = "CRP"))

### factor
df$treat <- relevel(factor(df$treat), ref = "mono") # Combination vs Monotherapy
df$types <- relevel(factor(df$types), ref = "VAP")
df$org <- relevel(factor(df$org), ref = "CRA")

# ---------------------------------------------------------
### recode types and org
df_pp <- df_pp %>%
  mutate(types = dplyr::recode(types, `0` = "VAP", `1` = "BSI"))

df_pp <- df_pp %>%
  mutate(org = dplyr::recode(org, `0` = "CRA", `1` = "CRE", `2` = "CRP"))

### factor
df_pp$treat <- relevel(factor(df_pp$treat), ref = "mono") # Combination vs Monotherapy
df_pp$types <- relevel(factor(df_pp$types), ref = "VAP")
df_pp$org <- relevel(factor(df_pp$org), ref = "CRA")


# ===========================================================
# Assumes you've already run your existing cleaning block:
#   df$treat <- relevel(factor(df$treat), ref = "mono")
#   df$types <- relevel(factor(df$types), ref = "VAP")
#   df$org   <- relevel(factor(df$org),   ref = "CRA")
# and that run_subgroup_analysis() (from your original script)
# is already defined and unchanged.
# ===========================================================

# -----------------------------------------------------------
# STEP 0: build the PP dataset / design
# -----------------------------------------------------------
# OPTION A: you already have a PP-specific weight column on df
#           (e.g. re-derived IPTW restricted to the PP sample,
#           with weight = 0 / NA outside the PP population).
# df_pp <- df   # <-- ADJUST if PP is a subset, see Option B
# design_pp <- svydesign(ids = ~id, weights = ~sw_treat1_99_pp, data = df_pp)  # <-- ADJUST

# OPTION B: PP is a row-subset of df (e.g. df$pp == 1 flag for
#           patients who received >=X days of assigned therapy,
#           no early crossover, etc.), reusing the ITT weight col.
                                   # <-- ADJUST
design_pp <- svydesign(ids = ~id, weights = ~sw_treat2_99, data = df_pp)  # <-- ADJUST

design_itt <- svydesign(ids = ~id, weights = ~sw_treat1_99, data = df)

# -----------------------------------------------------------
# STEP 1: wrap steps 4-5 from your original script into one
# reusable function so ITT and PP go through identical logic
# -----------------------------------------------------------
###
run_subgroup_analysis <- function(design, subgroup_var, outcome = "death_14") {
  
  # --- interaction model (overall test of effect modification) ---
  f_int <- as.formula(paste0(outcome, " ~ treat * ", subgroup_var))
  m_int <- svyglm(f_int, design = design, family = quasibinomial())
  
  int_terms <- tidy(m_int, conf.int = TRUE, exponentiate = TRUE) %>%  #This is exp(β3)
    filter(grepl(":", term))   # interaction term(s)
  
  # joint test if subgroup has >2 levels (e.g., pathogen with 3 levels)
  # regtermtest gives a proper Wald joint test across all interaction df
  joint_test <- regTermTest(m_int, paste0("treat:", subgroup_var))
  
  # --- stratum-specific estimates (treatment effect within each level) ---
  levels_sub <- levels(design$variables[[subgroup_var]])
  
  ###
  strata_results <- map_dfr(levels_sub, function(lvl) {
    sub_design <- subset(design, design$variables[[subgroup_var]] == lvl)
    f_strat <- as.formula(paste0(outcome, " ~ treat"))
    m_strat <- svyglm(f_strat, design = sub_design, family = quasibinomial())
    
  ###  
    tidy(m_strat, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term != "(Intercept)") %>%
      mutate(subgroup_variable = subgroup_var, subgroup_level = lvl, n = nrow(sub_design))
  })
  
  list(
    interaction_terms = int_terms,
    joint_interaction_pvalue = joint_test$p,
    stratum_estimates = strata_results
  )
}


# ---------------------------------------------------------
###
run_full_forest <- function(design, outcome = "death_14") {

  res_syndrome <- run_subgroup_analysis(design, "types", outcome)
  res_pathogen <- run_subgroup_analysis(design, "org",   outcome)

  m_overall <- svyglm(as.formula(paste0(outcome, " ~ treat")),
                       design = design, family = quasibinomial())

  overall_est <- tidy(m_overall, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(subgroup_variable = "Overall", subgroup_level = "All patients",
           n = nrow(design$variables))

  forest_data <- bind_rows(
    overall_est,
    res_syndrome$stratum_estimates,
    res_pathogen$stratum_estimates
  ) %>%
    mutate(
      label = paste0(subgroup_level, " (n=", n, ")"),
      label = factor(label, levels = rev(label)),
      or_ci = sprintf("%.1f (%.1f–%.1f)", estimate, conf.low, conf.high)
    )

  list(
    forest_data = forest_data,
    p_syndrome  = res_syndrome$joint_interaction_pvalue,
    p_pathogen  = res_pathogen$joint_interaction_pvalue
  )
}

itt <- run_full_forest(design_itt)
pp  <- run_full_forest(design_pp)

# -----------------------------------------------------------
###




# -----------------------------------------------------------
# STEP 2: shared x-axis range/breaks across BOTH panels so the
# ITT and PP forests are visually comparable (same scale)
# -----------------------------------------------------------
all_ci <- bind_rows(itt$forest_data, pp$forest_data)
x_min <- min(all_ci$conf.low,  na.rm = TRUE) * 0.85
x_max <- max(all_ci$conf.high, na.rm = TRUE) * 1.15
x_breaks <- scales::breaks_log(n = 6)(c(x_min, x_max))
x_breaks <- sort(unique(c(x_breaks, 1)))

# -----------------------------------------------------------
# STEP 3: panel builder — one forest + one OR/CI text column,
# with a title on top (mirrors "Baseline cohort" / "Landmark
# cohort" headers in your reference image)
# -----------------------------------------------------------
make_forest_panel <- function(forest_data, title) {
  ggplot(forest_data, aes(x = estimate, y = label)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_log10(limits = c(x_min, x_max), breaks = x_breaks,
                  labels = scales::label_number(accuracy = 0.1)) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title  = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.text.y = element_text(face = "bold"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_blank(),
      plot.margin = margin(10, 4, -25, 10)
    )
}

make_text_panel <- function(forest_data) {
  ggplot(forest_data, aes(x = 0, y = label)) +
    geom_text(aes(label = or_ci), hjust = 0, size = 3.5) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = "OR (95% CI)") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title  = element_text(face = "bold", hjust = 0, size = 11),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks  = element_blank(),
      axis.title  = element_blank(),
      panel.grid  = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_blank(),
      plot.margin = margin(10, 10, -25, 4)
    )
}

make_arrow_panel <- function(caption_text, caption_gap = 6, block_gap = 6) {
  # caption_gap: space between the "Favors..." labels and the caption
  #              text right below them (the "blue" gap)
  # block_gap:   space between this caption and whatever comes next
  #              (the "red" gap, e.g. before the following title)
  ggplot() +
    annotate("segment",
             x = 1 * 0.97, xend = x_min * 1.02,
             y = 1.15, yend = 1.15,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             linewidth = 0.5) +
    annotate("segment",
             x = 1 * 1.03, xend = x_max * 0.98,
             y = 1.15, yend = 1.15,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             linewidth = 0.5) +
    annotate("text", x = sqrt(1 * x_min), y = 0.97,
             label = "Favors Combination\ntherapy",
             size = 3.0, fontface = "bold.italic", hjust = 0.5) +
    annotate("text", x = sqrt(1 * x_max), y = 0.97,
             label = "Favors Monotherapy",
             size = 3.0, fontface = "bold.italic", hjust = 0.5) +
    scale_x_log10(limits = c(x_min, x_max)) +
    scale_y_continuous(limits = c(0.75, 1.15)) +
    labs(caption = caption_text) +
    theme_void(base_size = 12) +
    theme(
      plot.caption = element_text(size = 8.5, colour = "grey30", hjust = 0.5,
                                  margin = margin(t = caption_gap)),
      plot.margin = margin(-25, 4, block_gap, 10)
    )
}

# -----------------------------------------------------------
# STEP 4: assemble panels
# -----------------------------------------------------------
p_forest_itt <- make_forest_panel(itt$forest_data, "Intention-to-treat population")
p_text_itt   <- make_text_panel(itt$forest_data)
p_arrow_itt  <- make_arrow_panel(
  paste0("Interaction p — syndrome: ", round(itt$p_syndrome, 3),
         " | pathogen: ", round(itt$p_pathogen, 3)),
  caption_gap = 2,   # blue: tighter gap under "Favors..." labels
  block_gap   = 30   # red: more room before the PP title below
)

p_forest_pp <- make_forest_panel(pp$forest_data, "Per-protocol population")
p_text_pp   <- make_text_panel(pp$forest_data)
p_arrow_pp  <- make_arrow_panel(
  paste0("Interaction p — syndrome: ", round(pp$p_syndrome, 3),
         " | pathogen: ", round(pp$p_pathogen, 3)),
  caption_gap = 2,   # keep the two blocks visually consistent
  block_gap   = 6    # last row — no following block, keep it small
)

p_spacer <- plot_spacer()

# Each population gets its own forest+text row, immediately
# followed by its own arrow+caption row, then the next
# population starts below — arrow sits right under its forest
# instead of both arrows being pushed to the very bottom.
final_plot <-
  (p_forest_itt + p_text_itt) /
  (p_arrow_itt  + p_spacer)   /
  (p_forest_pp  + p_text_pp)  /
  (p_arrow_pp   + p_spacer)   +
  plot_layout(heights = c(6, 1.1, 6, 1.1), widths = c(4, 1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

###
final_plot

###
ggsave("subg_itt_pp_forest_plot.png", final_plot, width = 9, height = 12.4, dpi = 300)

#-------------------------------------#----------------------------------------#





###-----------------------------------------------------------------------------
### ^.^ ###