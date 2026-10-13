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
### OR interaction

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
library(tibble)

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
# Helper: extract ratio-of-ORs between two levels of a subgroup
# from the JOINT interaction model (correct SE via covariance)
# -----------------------------------------------------------
get_ror <- function(model, subgroup_var, level1, level2, treat_term = "treat") {
  b  <- coef(model)
  V  <- vcov(model)
  nm <- names(b)
  
  find_name <- function(lvl) {
    hits <- nm[grepl(paste0("^", treat_term, "[A-Za-z0-9_]*:", subgroup_var, lvl, "$"), nm)]
    if (length(hits) == 0) return(NA_character_)
    hits[1]
  }
  
  n1 <- find_name(level1)   # NA if level1 is the reference level
  n2 <- find_name(level2)   # NA if level2 is the reference level
  
  if (is.na(n1) && is.na(n2)) {
    warning(paste0(level1, " and ", level2, " both look like reference levels — check subgroup_var/levels"))
    return(tibble(comparison = paste0(level1, " vs ", level2),
                  ror = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  
  # contrast vector: reference level contributes 0 (no coefficient, no variance)
  cvec <- setNames(rep(0, length(b)), nm)
  if (!is.na(n1)) cvec[n1] <- cvec[n1] + 1
  if (!is.na(n2)) cvec[n2] <- cvec[n2] - 1
  
  est <- sum(cvec * b)
  se  <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
  
  tibble(
    comparison = paste0(level1, " vs ", level2),
    ror     = exp(est),
    ci_low  = exp(est - 1.96 * se),
    ci_high = exp(est + 1.96 * se)
  )
}
  

# -----------------------------------------------------------
# Wrapper: all pairwise ROR comparisons across levels of a subgroup
# -----------------------------------------------------------
get_all_pairwise_ror <- function(model, subgroup_var, levels_sub, treat_term = "treat") {
  combn(levels_sub, 2, simplify = FALSE) %>%
    map_dfr(~ get_ror(model, subgroup_var, .x[1], .x[2], treat_term)) %>%
    mutate(subgroup_variable = subgroup_var, .before = 1)
}

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
  
  strata_results <- map_dfr(levels_sub, function(lvl) {
    sub_design <- subset(design, design$variables[[subgroup_var]] == lvl)
    f_strat <- as.formula(paste0(outcome, " ~ treat"))
    m_strat <- svyglm(f_strat, design = sub_design, family = quasibinomial())
    
    tidy(m_strat, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term != "(Intercept)") %>%
      mutate(subgroup_variable = subgroup_var, subgroup_level = lvl, n = nrow(sub_design))
  })
  
  # NEW: pairwise ratio-of-ORs across all levels, from the joint model
  pairwise_ror <- get_all_pairwise_ror(m_int, subgroup_var, levels_sub)
  
  list(
    m_int = m_int,                     # <-- keep the fitted model for reuse
    interaction_terms = int_terms,
    joint_interaction_pvalue = joint_test$p,
    stratum_estimates = strata_results,
    pairwise_ror = pairwise_ror         # <-- NEW # This is exp(β4) # differ from ref level
  )
}

# -----------------------------------------------------------
# STEP 2: run_full_forest now also returns the ROR tables
# -----------------------------------------------------------
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
      or_ci = sprintf("%.2f (%.2f\u2013%.2f)", estimate, conf.low, conf.high)
    )
  
  list(
    forest_data  = forest_data,
    p_syndrome   = res_syndrome$joint_interaction_pvalue,
    p_pathogen   = res_pathogen$joint_interaction_pvalue,
    ror_syndrome = res_syndrome$pairwise_ror,   # e.g. BSI vs VAP row lives here
    ror_pathogen = res_pathogen$pairwise_ror
  )
}

itt <- run_full_forest(design_itt)
pp  <- run_full_forest(design_pp)

# -----------------------------------------------------------
itt$ror_syndrome
itt$ror_pathogen

# -----------------------------------------------------------



### within-group treatment effect
### "Among CRE patients only, what is the odds ratio of death for combo vs mono therapy?"
# -----------------------------------------------------------
get_group_or <- function(model, subgroup_var, level, treat_term = "treat") {
  b  <- coef(model)
  V  <- vcov(model)
  nm <- names(b)
  
  treat_main <- nm[grepl(paste0("^", treat_term, "[A-Za-z0-9_]*$"), nm) & !grepl(":", nm)]
  int_term   <- nm[grepl(paste0(treat_term, "[A-Za-z0-9_]*:", subgroup_var, level, "$"), nm)]
  
  cvec <- setNames(rep(0, length(b)), nm)
  cvec[treat_main] <- 1
  if (length(int_term) > 0) cvec[int_term] <- cvec[int_term] + 1  # 0 if level is reference
  
  est <- sum(cvec * b)
  se  <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
  
  tibble(
    subgroup_level = level,
    or      = exp(est),
    ci_low  = exp(est - 1.96 * se),
    ci_high = exp(est + 1.96 * se)
  )
}

# CRE-specific OR of combo vs mono, from the joint model:
get_group_or(res_pathogen$m_int, "org", "CRE")

# Cross-check against the subset-model version:
res_pathogen$stratum_estimates %>% filter(subgroup_level == "CRE")

#-------------------------------------#----------------------------------------#



### Treatment effect and 95% CI for the non-reference subgroup?
# -----------------------------------------------------------
### The mechanism: subset, then refit from scratch
### 
sub_design <- subset(design, design$variables[[subgroup_var]] == lvl)
f_strat <- as.formula(paste0(outcome, " ~ treat"))
m_strat <- svyglm(f_strat, design = sub_design, family = quasibinomial())

# OR_CI
tidy(m_strat, conf.int = TRUE, exponentiate = TRUE)




###-----------------------------------------------------------------------------
### ^.^ ###