
### IPW for Initiated Treatment Strategy with Grace period
### Treatment strategy: Initiate Combination Therapy or Monotherapy on or within 1 day after the AST report
### Sustain 3 days: The 3-day landmark analysis
### 20260922

# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(openxlsx,
                 magrittr,
                 dplyr,
                 purrr,
                 labelled,
                 gtsummary,
                 gt,
                 htmltools,
                 extrafont,
                 tibble,
                 tidyr,
                 writexl,
                 survival,
                 survminer,
                 ggplot2,
                 boot,
                 Hmisc,
                 Cairo,
                 MatchIt,
                 cobalt,
                 splitstackshape,
                 speedglm,
                 broom,
                 patchwork,
                 ggtext,
                 survey,
                 tableone,
                 WeightIt,
                 marginaleffects,
                 sandwich,
                 lmtest,
                 car,
                 scales,
                 ggbreak,
                 grid)
})

# Define working directory
wd <- "./"
setwd(wd)

### Load data
#
df <- readRDS("data/clean data/clean/clean/combine_group/baseline_all_clean_combine.RData")
df_bl <- readRDS("data/clean data/clean/clean/combine_group/ipw/baseline_ipw_clean_combine.RData")
df_lm <- readRDS("data/clean data/clean/clean/combine_group/ipw/landmark_ipw_clean_combine.RData")

#
df_event <- readRDS("data/clean data/clean/clean/combine_group/event_KM.RData")
df_event_bl <- readRDS("data/clean data/clean/clean/combine_group/ipw/ipw_bl.RData")
df_event_lm <- readRDS("data/clean data/clean/clean/combine_group/ipw/ipw_lm.RData")
df_event_lm_all <- readRDS("data/clean data/clean/clean/combine_group/ipw/ipw_lm_all.RData")


### Load fonts
loadfonts()

#-------------------------------------
### IPW1 Baseline
### Extract relevant covariates
### Baseline covariates ~ at randomization (Time Zero)
### Same for 2 IPWs
### Wide dataset
#
vars_bl <- c("recordid", 
            "age_new", "sex", "country_income", 
            "comorbidities_CCI",  
            "cmb_rhe", 
            "sofa_score", "vaso_use", "severity_score",
            "icu_at_onset", "mv_at_onset", 
            "length_before_onset", 
            "org_r_combined", "infection_types",
            "d21_death", "d21_random_days",
            "d07_death", "d07_random_days",
            "d14_death", "d14_random_days")

#
covariates_outcome_bl <- df_bl %>%
  dplyr::select(all_of(vars_bl))

#-------------------------------------
# Shorten Variable Name
covariates_outcome_bl <- rename(covariates_outcome_bl,
                  id = recordid,
                  age = age_new,
                  income = country_income,
                  cci = comorbidities_CCI,
                  rhe = cmb_rhe,
                  sofa = sofa_score,
                  vaso = vaso_use,
                  severity = severity_score,
                  icu = icu_at_onset,
                  mv = mv_at_onset,
                  los = length_before_onset,
                  org = org_r_combined,
                  types = infection_types,
                  death_21 = d21_death,
                  time_to_death_21 = d21_random_days,
                  death_07 = d07_death,
                  time_to_death_07 = d07_random_days,
                  death_14 = d14_death,
                  time_to_death_14 = d14_random_days)

#-------------------------------------
# Standardize Variable Type
covariates_outcome_bl$age <- as.numeric(covariates_outcome_bl$age)
covariates_outcome_bl$cci <- as.numeric(covariates_outcome_bl$cci)

# F ~ 0; M ~ 1
covariates_outcome_bl <- covariates_outcome_bl %>%
  mutate(sex = case_when(
    sex == "F" ~ 0,
    sex == "M" ~ 1
  ))

#
covariates_outcome_bl <- covariates_outcome_bl %>%
  mutate(income = case_when(
    income == "High income" ~ 0,
    income == "Upper middle income" ~ 1,
    income == "Lower middle income" ~ 2
  ))

#
covariates_outcome_bl <- covariates_outcome_bl %>%
  mutate(types = case_when(
    types == "VAP" ~ 0,
    types == "BSI" ~ 1
  ))


#
covariates_outcome_bl <- covariates_outcome_bl %>%
  mutate(org = case_when(
    org == "CRA" ~ 0,
    org == "CRE" ~ 1,
    org == "CRP" ~ 2
  ))

#-------------------------------------
### Select variables
df_bl_used <- df_event_bl %>%
  dplyr::select(recordid, treatment_strategy, matched_day) %>%
  mutate(treatment_strategy = dplyr::recode(treatment_strategy, "Monotherapy" = "mono", 
                                            "Combination therapy" = "combo"))


### Rename variables
df_bl_used <- rename(df_bl_used,
                  id = recordid,
                  init_day = matched_day,
                  treat = treatment_strategy)


### Combine covariates_outcome
###
final_data_bl <- df_bl_used %>%
  left_join(covariates_outcome_bl, by = "id")

### Check NA 
colSums(is.na(final_data_bl))

#-------------------------------------#----------------------------------------#



## =============================================================
## 1. Baseline treatment weights (IPW) -- stabilized.
## =============================================================
### 
# Keep one baseline row per person
# Use complete cases before fitting
baseline_for_ps <- final_data_bl %>%
  distinct(id, treat, types,
           age, sex, income,
           cci, rhe,
           severity, vaso,
           icu, mv,
           los) %>%
  filter(complete.cases(treat, 
                        age, sex, income,
                        cci, rhe,
                        severity, vaso,
                        icu, mv,
                        los))

### Fit the propensity score model
### Treatment ~ binary: combo vs mono
### The prob that a person gets combo treatment
### given their baseline covariates
### P(Combo=1|baseline covariates L)
###
ps_bl_model <- glm(I(treat == "combo") ~ age + sex + income 
                  + cci + rhe
                  + severity + vaso 
                  + icu + mv 
                  + los,
                  data = baseline_for_ps, family = binomial())


###
baseline_for_ps <- baseline_for_ps %>%
  mutate(
    ps_combo = predict(ps_bl_model, type = "response"), # individual predicted prob of receiving combo
    # stabilized numerator: marginal probability of observed treat
    p_treat_marginal = mean(treat == "combo"), # compute the marginal treatment prob ~ overall proportion of people received combo in the data (stab NUMerator)
    sw_treat1 = case_when(  # stab treatment weights
      treat == "combo" ~ p_treat_marginal / ps_combo, # assign a weight based on assigned strategy
      treat == "mono"  ~ (1 - p_treat_marginal) / (1 - ps_combo)
    )
  ) %>%
  dplyr::select(id, sw_treat1)

###
summary(baseline_for_ps$sw_treat1)
sd(baseline_for_ps$sw_treat1)

#-------------------------------------
### Truncate stabilized weights at the 99th percentile ###; 
# Truncate at the 99th percentile
### quantile()
threshold_bl_99 <- quantile(baseline_for_ps$sw_treat1, 0.99) # assign 99th percentile of Weights
baseline_for_ps$sw_treat1_99 <- baseline_for_ps$sw_treat1
baseline_for_ps$sw_treat1_99[baseline_for_ps$sw_treat1_99 > threshold_bl_99] <- threshold_bl_99


###  Min, 25th percentile, median, mean, SD, 75th percentile, and max: 
### stabilized weights truncated at 99th percentile ###
summary(baseline_for_ps$sw_treat1_99)
sd(baseline_for_ps$sw_treat1_99)

#-------------------------------------#----------------------------------------#



### Checking (Luca's code for TTE)
#=========================================================
# Propensity-score model for the combined cohort
#=========================================================
###
baseline_for_ps_c <- final_data_bl %>%
  distinct(id, treat, types,
           age, sex, income,
           cci, rhe,
           severity, vaso,
           icu, mv,
           los) %>%
  filter(complete.cases(treat, 
                        age, sex, income,
                        cci, rhe,
                        severity, vaso,
                        icu, mv,
                        los))


###
m_ps_bl <- glm(I(treat == "combo") ~ age + sex + income 
                                      + cci + rhe
                                      + severity + vaso 
                                      + icu + mv 
                                      + los,
                                    data = baseline_for_ps_c,
                                    family = binomial()
               )


#=========================================================
# Predict propensity scores
#=========================================================
###
baseline_for_ps_c <- baseline_for_ps_c %>%
  mutate(
    ps_bl = predict(
      m_ps_bl,
      newdata = baseline_for_ps_c,
      type = "response"
    )
  )

#
summary(baseline_for_ps_c$ps_bl)


#=========================================================
# Calculate stabilized IPTW
#=========================================================
# Marginal probability of combination therapy
p_treated <- mean(
  baseline_for_ps_c$treat == "combo",
  na.rm = TRUE
)

#
baseline_for_ps_c <- baseline_for_ps_c %>%
  mutate(
    iptw_bl = case_when(
      treat == "combo" ~ p_treated / ps_bl,
      treat == "mono" ~ (1 - p_treated) / (1 - ps_bl),
      TRUE ~ NA_real_
    )
  )

#
summary(baseline_for_ps_c$iptw_bl)

#=========================================================
# Truncate weights at the 1st and 99th percentiles
#=========================================================
###
weight_limits <- quantile(
  baseline_for_ps_c$iptw_bl,
  probs = c(0.01, 0.99),
  na.rm = TRUE
)

#
baseline_for_ps_c <- baseline_for_ps_c %>%
  mutate(
    iptw_overall_trunc = pmin(
      pmax(iptw_bl, weight_limits[[1]]),
      weight_limits[[2]]
    )
  )

#
summary(baseline_for_ps_c$iptw_overall_trunc)

#-------------------------------------#----------------------------------------#
#-------------------------------------#----------------------------------------#





#-------------------------------------
### Filter IDs
final_data_bl2 <- final_data_bl[which(final_data_bl$id %in% baseline_for_ps$id),]

###
final_data_bl2 <- final_data_bl2 %>%
  left_join(baseline_for_ps, by = "id")

#-------------------------------------
### At the beginning of Results
### CRA ~ 0, CRE ~ 1, CRP ~ 2
table(final_data_bl2$org)

###
table(final_data_bl2$treat)


#-------------------------------------
### Save Baseline dataset for Tables
saveRDS(final_data_bl2, file = "data/clean data/clean/clean/final/tables_bl.RData")

#-------------------------------------#----------------------------------------#
### Checking
###
table(final_data_bl2$types)

###
covariates_outcome_bl_check <- covariates_outcome_bl[which(covariates_outcome_bl$id %in% baseline_for_ps$id),]
table(covariates_outcome_bl_check$types)




## =============================================================
## 1.1. Compute SMD for each covariate - after IPW.
## =============================================================
### Define covariates
covars_bl <- c("age", "sex", "income",
              "cci", "rhe",
              "severity", "vaso",
              "icu", "mv", 
              "los")

### Unweighted balance check
tab_bl_unw <- CreateTableOne(vars = covars_bl, strata = "treat",
                          data = final_data_bl2, test = FALSE)

print(tab_bl_unw, smd = TRUE)

### Weighted balance check using IPW
tab_bl_w <- svydesign(ids = ~1, data = final_data_bl2, weights = ~sw_treat1_99)

bal_bl <- svyCreateTableOne(vars = covars_bl, strata = "treat", data = tab_bl_w, test = FALSE)

print(bal_bl, smd = TRUE)

### Data for Plot
###
final_data_bl_plot <- final_data_bl2

###
final_data_bl_plot <- rename(final_data_bl_plot,
                          Age = age,
                          Sex = sex,
                          Country_income_level = income,
                          Charlson_comorbidity_index = cci,
                          Recent_healthcare_exposure = rhe,
                          qSOFA = severity,
                          Vasopressor_use = vaso, 
                          ICU_admission = icu,
                          MV_status = mv,
                          Hospital_LOS_before = los)


### Love plot
lp_bl <- bal.tab(
  treat ~ Age + Sex + Country_income_level
  + Charlson_comorbidity_index + Recent_healthcare_exposure
  + qSOFA + Vasopressor_use
  + ICU_admission + MV_status 
  + Hospital_LOS_before,
  data = final_data_bl_plot,
  weights = final_data_bl_plot$sw_treat1_99,
  method = "weighting",
  s.d.denom = "pooled",
  binary = "std"      # <-- forces true SMD for binary covariates too
) |>
  love.plot(
    bal_bl,
    stats = "mean.diffs",
    abs = TRUE,
    thresholds = c(m = 0.1),
    var.order = "unadjusted",
    line = TRUE,
    size = 3,
    alpha = 0.9,
    colors = c("#5B2D90", "#E11D8C"),
    shapes = c(16, 17),
    sample.names = c("Unweighted", "Weighted"),
    xlab = "Standardized Mean Difference",
    title = "Covariate balance (Intention-to-treat population)",
    theme = theme_minimal(base_size = 12) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
  ) +
  geom_vline(xintercept = c(0, 0.1), linetype = c("dashed", "dashed"),
             color = "grey50") + 
  labs(x = "Standardized Mean Difference")

###
lp_bl

###-----------------------------------------------------------------------------

#colors = c("#F8766D", "#00BFC4")




### IPW2 Landmark
#-------------------------------------
#-------------------------------------
### Sustained 3d
### dev = 0 ~ sustain = 1
df_lm$sustain_3d <- ifelse(df_lm$dev_since_matched_day == 1, 0, 1)

#
table(df_lm$treatment_strategy, df_lm$sustain_3d)

#
df_lm_st <- df_lm %>%
  filter(sustain_3d == 1)

### Extract relevant covariates
### Baseline covariates ~ at randomization (Time Zero)
### Wide dataset
#
vars_lm <- c("recordid", 
            "age_new", "sex", "country_income", 
            "comorbidities_CCI",
            "cmb_rhe", 
            "sofa_score", "vaso_use", "severity_score",
            "icu_at_onset", "mv_at_onset", 
            "length_before_onset", 
            "org_r_combined", "infection_types",
            "d21_death", "d21_random_days",
            "d07_death", "d07_random_days",
            "d14_death", "d14_random_days",
            "sustain_3d")

###
covariates_outcome_lm <- df_lm_st %>%
  dplyr::select(all_of(vars_lm))

#-------------------------------------
# Shorten Variable Name
covariates_outcome_lm <- rename(covariates_outcome_lm,
                                id = recordid,
                                age = age_new,
                                income = country_income,
                                cci = comorbidities_CCI,
                                rhe = cmb_rhe,
                                sofa = sofa_score,
                                vaso = vaso_use,
                                severity = severity_score,
                                icu = icu_at_onset,
                                mv = mv_at_onset,
                                los = length_before_onset,
                                org = org_r_combined,
                                types = infection_types,
                                death_21 = d21_death,
                                time_to_death_21 = d21_random_days,
                                death_07 = d07_death,
                                time_to_death_07 = d07_random_days,
                                death_14 = d14_death,
                                time_to_death_14 = d14_random_days)

#-------------------------------------
# Standardize Variable Type
covariates_outcome_lm$age <- as.numeric(covariates_outcome_lm$age)
covariates_outcome_lm$cci <- as.numeric(covariates_outcome_lm$cci)

# F ~ 0; M ~ 1
covariates_outcome_lm <- covariates_outcome_lm %>%
  mutate(sex = case_when(
    sex == "F" ~ 0,
    sex == "M" ~ 1
  ))

#
covariates_outcome_lm <- covariates_outcome_lm %>%
  mutate(income = case_when(
    income == "High income" ~ 0,
    income == "Upper middle income" ~ 1,
    income == "Lower middle income" ~ 2
  ))

#
covariates_outcome_lm <- covariates_outcome_lm %>%
  mutate(types = case_when(
    types == "VAP" ~ 0,
    types == "BSI" ~ 1
  ))


#
covariates_outcome_lm <- covariates_outcome_lm %>%
  mutate(org = case_when(
    org == "CRA" ~ 0,
    org == "CRE" ~ 1,
    org == "CRP" ~ 2
  ))

#-------------------------------------
###
df_event_lm_st <- df_event_lm[which(df_event_lm$recordid %in% df_lm_st$recordid),]

### Select variables
df_event_lm_st_used <- df_event_lm_st %>%
  dplyr::select(recordid, treatment_strategy, matched_day) %>%
  mutate(treatment_strategy = dplyr::recode(treatment_strategy, "Monotherapy" = "mono", 
                                     "Combination therapy" = "combo"))


### Rename variables
df_event_lm_st_used <- rename(df_event_lm_st_used,
                     id = recordid,
                     init_day = matched_day,
                     treat = treatment_strategy)


### Combine covariates_outcome
###
final_data_lm_st <- df_event_lm_st_used %>%
  left_join(covariates_outcome_lm, by = "id")

### Check NA 
colSums(is.na(final_data_lm_st))

#-------------------------------------#----------------------------------------#



## =============================================================
## 2. Landmark treatment weights (IPW) -- stabilized.
## =============================================================
### 
# Keep one baseline row per person
# Use complete cases before fitting
landmark_for_ps <- final_data_lm_st %>%
  distinct(id, treat, 
           age, sex, income,
           cci, rhe,
           severity, vaso,
           icu, mv,
           los) %>%
  filter(complete.cases(treat, 
                        age, sex, income,
                        cci, rhe,
                        severity, vaso,
                        icu, mv,
                        los))

### Fit the propensity score model
### The prob that a person gets combo treatment
### given their baseline covariates
### P(Combo=1|baseline covariates L)
###
ps_lm_model <- glm(I(treat == "combo") ~ age + sex + income 
                   + cci + rhe
                   + severity + vaso 
                   + icu + mv 
                   + los,
                   data = landmark_for_ps, family = binomial())


###
landmark_for_ps <- landmark_for_ps %>%
  mutate(
    ps_combo = predict(ps_lm_model, type = "response"), # individual predicted prob of receiving combo
    # stabilized numerator: marginal probability of observed treat
    p_treat_marginal = mean(treat == "combo"), # compute the marginal treatment prob ~ overall proportion of people received combo in the data (stab NUMerator)
    sw_treat2 = case_when(  # stab treatment weights
      treat == "combo" ~ p_treat_marginal / ps_combo, # assign a weight based on assigned strategy
      treat == "mono"  ~ (1 - p_treat_marginal) / (1 - ps_combo)
    )
  ) %>%
  dplyr::select(id, sw_treat2)

###
summary(landmark_for_ps$sw_treat2)
sd(landmark_for_ps$sw_treat2)

#-------------------------------------
### Truncate stabilized weights at the 99th percentile ###; 
# Truncate at the 99th percentile
### quantile()
threshold_lm_99 <- quantile(landmark_for_ps$sw_treat2, 0.99) # assign 99th percentile of Weights
landmark_for_ps$sw_treat2_99 <- landmark_for_ps$sw_treat2
landmark_for_ps$sw_treat2_99[landmark_for_ps$sw_treat2_99 > threshold_lm_99] <- threshold_lm_99


###  Min, 25th percentile, median, mean, SD, 75th percentile, and max: 
### stabilized weights truncated at 99th percentile ###
summary(landmark_for_ps$sw_treat2_99)
sd(landmark_for_ps$sw_treat2_99)


#-------------------------------------
### Filter IDs
final_data_lm2 <- final_data_lm_st[which(final_data_lm_st$id %in% landmark_for_ps$id),]

###
final_data_lm2 <- final_data_lm2 %>%
  left_join(landmark_for_ps, by = "id")

#-------------------------------------
### Save Baseline dataset for Tables
saveRDS(final_data_lm2, file = "data/clean data/clean/clean/final/tables_lm.RData")

#-------------------------------------#----------------------------------------#



## =============================================================
## 1.2. Compute SMD for each covariate - after IPW.
## =============================================================
### Define covariates
covars_lm <- c("age", "sex", "income",
               "cci", "rhe",
               "severity", "vaso",
               "icu", "mv", 
               "los")

### Unweighted balance check
tab_lm_unw <- CreateTableOne(vars = covars_lm, strata = "treat",
                             data = final_data_lm2, test = FALSE)

print(tab_lm_unw, smd = TRUE)

### Weighted balance check using IPW
tab_lm_w <- svydesign(ids = ~1, data = final_data_lm2, weights = ~sw_treat2_99)

bal_lm <- svyCreateTableOne(vars = covars_lm, strata = "treat", data = tab_lm_w, test = FALSE)

print(bal_lm, smd = TRUE)

### Data for Plot
###
final_data_lm_plot <- final_data_lm2

###
final_data_lm_plot <- rename(final_data_lm_plot,
                             Age = age,
                             Sex = sex,
                             Country_income_level = income,
                             Charlson_comorbidity_index = cci,
                             Recent_healthcare_exposure = rhe,
                             qSOFA = severity,
                             Vasopressor_use = vaso, 
                             ICU_admission = icu,
                             MV_status = mv,
                             Hospital_LOS_before = los)


### Love plot
lp_lm <- bal.tab(
  treat ~ Age + Sex + Country_income_level
  + Charlson_comorbidity_index + Recent_healthcare_exposure
  + qSOFA + Vasopressor_use
  + ICU_admission + MV_status 
  + Hospital_LOS_before,
  data = final_data_lm_plot,
  weights = final_data_lm_plot$sw_treat2_99,
  method = "weighting",
  s.d.denom = "pooled",  # standardize differences by pooled SD
  binary = "std"      # <-- forces true SMD for binary covariates too
) |>
  love.plot(
    bal_lm,
    stats = "mean.diffs",  # mean differences
    abs = TRUE,  # absolute value
    thresholds = c(m = 0.1),
    var.order = "unadjusted",
    line = TRUE,
    size = 3,
    alpha = 0.9,
    colors = c("#5B2D90", "#E11D8C"),
    shapes = c(16, 17),
    sample.names = c("Unweighted", "Weighted"),
    xlab = "Standardized Mean Difference",
    title = "Covariate balance (Per-protocol population)",
    theme = theme_minimal(base_size = 12) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
  ) +
  scale_x_continuous(breaks = seq(0, 0.7, by = 0.1)) +
  geom_vline(xintercept = c(0, 0.1), linetype = c("dashed", "dashed"),
             color = "grey50") + 
  labs(x = "Standardized Mean Difference")

###
lp_lm

#-------------------------------------
### Combine 2 love plots
library(gridExtra)
lp <- grid.arrange(lp_bl, lp_lm, nrow = 2)

# ---- 10. Save ---------------------------------------------------------------
ggsave("sfig1_love_plot_combo_vs_mono.png", lp, width = 12.5, height = 6, dpi = 320, bg = "white")

###-----------------------------------------------------------------------------


### NA covariates
table(is.na(df$severity_score))
table(is.na(df$length_before_onset))
table(is.na(df$vaso_use))

#-------------------------------------
### CONSORT diagram

### Filter IDs 
df_event_lm_all <- df_event_lm_all[which(df_event_lm_all$recordid %in% final_data_bl2$id),]

### survivors 3 days
table(df_event_lm_all$treatment_strategy, df_event_lm_all$event_type_within_3day)

### sustained 3 days
table(df_lm_st$treatment_strategy)

###
table(final_data_lm2$treat)

###-----------------------------------------------------------------------------



# =============================================================================
# 3. COMBINE WEIGHTS  W_final = SW_T × SW_S
#   Multiply baseline IPW by IPCW for each patient.
#   Only SUSTAINERS enter the outcome model; non-sustainers were censored.
# =============================================================================
### Combine Weights for Sustainers
#final_data_lm2 <- left_join(final_data_lm2, final_data_bl2[, c("id", "sw_treat1_99")], 
                            #by = "id")

###
#final_data_lm2 <- final_data_lm2 %>%
  #mutate(w_final = sw_treat1_99 * sw_treat2_99)

###  Min, 25th percentile, median, mean, SD, 75th percentile, and max: 
### stabilized weights truncated at 99th percentile ###
#summary(final_data_lm2$w_final)
#sd(final_data_lm2$w_final)

# =============================================================================
# 4.  WEIGHTED OUTCOME MODEL  — Marginal Structural Model ~ Weighted LR (GLM)
# Binary treatment ~ combo vs mono
# 4.1. Binary outcome ~ 21-day ACM
# =============================================================================
### Ref level ~ Mono, RD as Combo minus Mono
final_data_bl2$treat <- relevel(factor(final_data_bl2$treat), ref = "mono")

final_data_lm2$treat <- relevel(factor(final_data_lm2$treat), ref = "mono")

### Baseline cohort
#-------------------------------------
#-------------------------------------
### Risk Difference — identity link
msm_rd_bl_21 <- glm(
  death_21 ~ treat,
  data    = final_data_bl2,
  family  = binomial(link = "identity"),
  weights = sw_treat1_99,
  start   = c(mean(final_data_bl2$death_21), 0)   # helps convergence
)

#-------------------------------------
###
coeftest(msm_rd_bl_21, vcov = vcovHC(msm_rd_bl_21, type = "HC3"))
# coef["arm"] = RD directly

### Getting the 95% CI directly
# Extract from your existing model
rd_bl_21    <- coef(msm_rd_bl_21)["treatcombo"]
se_rd_bl_21 <- sqrt(vcovHC(msm_rd_bl_21, type = "HC3")["treatcombo", "treatcombo"])

lower_bl_21 <- rd_bl_21 - 1.96 * se_rd_bl_21
upper_bl_21 <- rd_bl_21 + 1.96 * se_rd_bl_21

cat("RD:        ", round(rd_bl_21*100, 1), "\n")
cat("95% CI:    ", round(lower_bl_21*100, 1), "to", round(upper_bl_21*100, 1), "\n")

#-------------------------------------
### Monotherapy risk

risk_mono_bl_21 <- coef(msm_rd_bl_21)["(Intercept)"]
se_mono_bl_21   <- sqrt(vcovHC(msm_rd_bl_21, type = "HC3")["(Intercept)", "(Intercept)"])

lower_mono_bl_21 <- risk_mono_bl_21 - 1.96 * se_mono_bl_21
upper_mono_bl_21 <- risk_mono_bl_21 + 1.96 * se_mono_bl_21

cat("Monotherapy risk: ", round(risk_mono_bl_21*100, 1), "\n")
cat("95% CI: ", round(lower_mono_bl_21*100, 1), "to", round(upper_mono_bl_21*100, 1), "\n")

### Combination risk
#
result_combo_bl_21 <- deltaMethod(msm_rd_bl_21, "`(Intercept)` + treatcombo", vcov. = vcovHC(msm_rd_bl_21, type = "HC3"))

# %
result_combo_bl_21_scaled <- round(result_combo_bl_21 * 100, 1)
print(result_combo_bl_21_scaled)

#-------------------------------------
### Risk Ratio
### Combination ÷ Monotherapy
rr_result_bl_21 <- deltaMethod(msm_rd_bl_21, 
                         "(`(Intercept)` + treatcombo) / `(Intercept)`", 
                         vcov. = vcovHC(msm_rd_bl_21, type = "HC3"))

print(rr_result_bl_21)

cat("RR: ", round(rr_result_bl_21$Estimate, 1), "\n")
cat("95% CI: ", round(rr_result_bl_21$`2.5 %`, 1), "to", round(rr_result_bl_21$`97.5 %`, 1), "\n")


#-------------------------------------#----------------------------------------#




### CHECKING (Luca's code for TTE)
#---------------------------------------------------------
# Function for IPTW risk difference
#---------------------------------------------------------
### combine event
baseline_for_ps_c <- left_join(baseline_for_ps_c, final_data_bl2[, c("id", "death_21")], 
                               by = "id")

## RISK DIFFERENCES - primary outcome
###
df_c <- baseline_for_ps_c %>% 
{ stopifnot(nrow(.) > 0); . }
risk_A1 <- with(df_c,
                sum(iptw_overall_trunc * death_21 * (treat == "combo")) / sum(iptw_overall_trunc * (treat == "combo"))
)
risk_A0 <- with(df_c,
                sum(iptw_overall_trunc * death_21 * (treat == "mono")) / sum(iptw_overall_trunc * (treat == "mono"))
)
RD_21 <- risk_A1 - risk_A0
set.seed(123)
B <- 1000
n <- nrow(df_c)
risk_A1_boot <- numeric(B)
risk_A0_boot <- numeric(B)
rd_boot <- numeric(B)
for (b in seq_len(B)) {
  idx <- sample.int(n, n, replace = TRUE)
  d_b <- df_c[idx, , drop = FALSE]
  
  r1 <- with(d_b,
             sum(iptw_overall_trunc * death_21 * (treat == "combo")) / sum(iptw_overall_trunc * (treat == "combo"))
  )
  
  r0 <- with(d_b,
             sum(iptw_overall_trunc * death_21 * (treat == "mono")) / sum(iptw_overall_trunc * (treat == "mono"))
  )
  risk_A1_boot[b] <- r1
  risk_A0_boot[b] <- r0
  rd_boot[b] <- r1 - r0
}
CI_A1 <- quantile(risk_A1_boot, c(0.025, 0.975))
CI_A0 <- quantile(risk_A0_boot, c(0.025, 0.975))
CI_RD <- quantile(rd_boot, c(0.025, 0.975))
list(
  risk_A1_21 = sprintf("%.1f%% (%.1f–%.1f)",
                       risk_A1*100, CI_A1[1]*100, CI_A1[2]*100),
  risk_A0_21 = sprintf("%.1f%% (%.1f–%.1f)",
                       risk_A0*100, CI_A0[1]*100, CI_A0[2]*100),
  RD_21 = sprintf("%.1f%% (%.1f–%.1f)",
                  RD_21*100, CI_RD[1]*100, CI_RD[2]*100)
)

#-------------------------------------#----------------------------------------#
#-------------------------------------#----------------------------------------#



### Landmark cohort
#-------------------------------------
#-------------------------------------
### Risk Difference — identity link
msm_rd_lm_21 <- glm(
  death_21 ~ treat,
  data    = final_data_lm2,
  family  = binomial(link = "identity"),
  weights = sw_treat2_99,
  start   = c(mean(final_data_lm2$death_21), 0)   # helps convergence
)

###
coeftest(msm_rd_lm_21, vcov = vcovHC(msm_rd_lm_21, type = "HC3"))
# coef["arm"] = RD directly

### Getting the 95% CI directly
# Extract from your existing model
rd_lm_21    <- coef(msm_rd_lm_21)["treatcombo"]
se_rd_lm_21 <- sqrt(vcovHC(msm_rd_lm_21, type = "HC3")["treatcombo", "treatcombo"])

lower_lm_21 <- rd_lm_21 - 1.96 * se_rd_lm_21
upper_lm_21 <- rd_lm_21 + 1.96 * se_rd_lm_21

cat("RD:        ", round(rd_lm_21*100, 1), "\n")
cat("95% CI:    ", round(lower_lm_21*100, 1), "to", round(upper_lm_21*100, 1), "\n")

#-------------------------------------
### Monotherapy risk

risk_mono_lm_21 <- coef(msm_rd_lm_21)["(Intercept)"]
se_mono_lm_21   <- sqrt(vcovHC(msm_rd_lm_21, type = "HC3")["(Intercept)", "(Intercept)"])

lower_mono_lm_21 <- risk_mono_lm_21 - 1.96 * se_mono_lm_21
upper_mono_lm_21 <- risk_mono_lm_21 + 1.96 * se_mono_lm_21

cat("Monotherapy risk: ", round(risk_mono_lm_21*100, 1), "\n")
cat("95% CI: ", round(lower_mono_lm_21*100, 1), "to", round(upper_mono_lm_21*100, 1), "\n")

### Combination risk
#
result_combo_lm_21 <- deltaMethod(msm_rd_lm_21, "`(Intercept)` + treatcombo", vcov. = vcovHC(msm_rd_lm_21, type = "HC3"))

# %
result_combo_lm_21_scaled <- round(result_combo_lm_21 * 100, 1)
print(result_combo_lm_21_scaled)

#-------------------------------------
### Risk Ratio
### Combination ÷ Monotherapy
rr_result_lm_21 <- deltaMethod(msm_rd_lm_21, 
                               "(`(Intercept)` + treatcombo) / `(Intercept)`", 
                               vcov. = vcovHC(msm_rd_lm_21, type = "HC3"))

print(rr_result_lm_21)

cat("RR: ", round(rr_result_lm_21$Estimate, 1), "\n")
cat("95% CI: ", round(rr_result_lm_21$`2.5 %`, 1), "to", round(rr_result_lm_21$`97.5 %`, 1), "\n")


#-------------------------------------#-----------------------------------------



# =============================================================================
# 4. WEIGHTED OUTCOME MODEL  — Marginal Structural Model ~ Weighted LR (GLM)
# Binary treatment ~ combo vs mono
# 4.2. Binary outcome ~ 07-day ACM
# =============================================================================
### Baseline cohort
#-------------------------------------
#-------------------------------------
### Risk Difference — identity link
msm_rd_bl_07 <- glm(
  death_07 ~ treat,
  data    = final_data_bl2,
  family  = binomial(link = "identity"),
  weights = sw_treat1_99,
  start   = c(mean(final_data_bl2$death_07), 0)   # helps convergence
)

###
coeftest(msm_rd_bl_07, vcov = vcovHC(msm_rd_bl_07, type = "HC3"))
# coef["arm"] = RD directly

### Getting the 95% CI directly
# Extract from your existing model
rd_bl_07    <- coef(msm_rd_bl_07)["treatcombo"]
se_rd_bl_07 <- sqrt(vcovHC(msm_rd_bl_07, type = "HC3")["treatcombo", "treatcombo"])

lower_bl_07 <- rd_bl_07 - 1.96 * se_rd_bl_07
upper_bl_07 <- rd_bl_07 + 1.96 * se_rd_bl_07

cat("RD:        ", round(rd_bl_07*100, 1), "\n")
cat("95% CI:    ", round(lower_bl_07*100, 1), "to", round(upper_bl_07*100, 1), "\n")

#-------------------------------------
### Monotherapy risk

risk_mono_bl_07 <- coef(msm_rd_bl_07)["(Intercept)"]
se_mono_bl_07   <- sqrt(vcovHC(msm_rd_bl_07, type = "HC3")["(Intercept)", "(Intercept)"])

lower_mono_bl_07 <- risk_mono_bl_07 - 1.96 * se_mono_bl_07
upper_mono_bl_07 <- risk_mono_bl_07 + 1.96 * se_mono_bl_07

cat("Monotherapy risk: ", round(risk_mono_bl_07*100, 1), "\n")
cat("95% CI: ", round(lower_mono_bl_07*100, 1), "to", round(upper_mono_bl_07*100, 1), "\n")

### Combination risk
#
result_combo_bl_07 <- deltaMethod(msm_rd_bl_07, "`(Intercept)` + treatcombo", vcov. = vcovHC(msm_rd_bl_07, type = "HC3"))

# %
result_combo_bl_07_scaled <- round(result_combo_bl_07 * 100, 1)
print(result_combo_bl_07_scaled)

#-------------------------------------
### Risk Ratio
### Combination ÷ Monotherapy
rr_result_bl_07 <- deltaMethod(msm_rd_bl_07, 
                               "(`(Intercept)` + treatcombo) / `(Intercept)`", 
                               vcov. = vcovHC(msm_rd_bl_07, type = "HC3"))

print(rr_result_bl_07)

cat("RR: ", round(rr_result_bl_07$Estimate, 1), "\n")
cat("95% CI: ", round(rr_result_bl_07$`2.5 %`, 1), "to", round(rr_result_bl_07$`97.5 %`, 1), "\n")


#-------------------------------------#-----------------------------------------



### Landmark cohort
#-------------------------------------
#-------------------------------------
### Risk Difference — identity link
msm_rd_lm_07 <- glm(
  death_07 ~ treat,
  data    = final_data_lm2,
  family  = binomial(link = "identity"),
  weights = sw_treat2_99,
  start   = c(mean(final_data_lm2$death_07), 0)   # helps convergence
)

###
coeftest(msm_rd_lm_07, vcov = vcovHC(msm_rd_lm_07, type = "HC3"))
# coef["arm"] = RD directly

### Getting the 95% CI directly
# Extract from your existing model
rd_lm_07    <- coef(msm_rd_lm_07)["treatcombo"]
se_rd_lm_07 <- sqrt(vcovHC(msm_rd_lm_07, type = "HC3")["treatcombo", "treatcombo"])

lower_lm_07 <- rd_lm_07 - 1.96 * se_rd_lm_07
upper_lm_07 <- rd_lm_07 + 1.96 * se_rd_lm_07

cat("RD:        ", round(rd_lm_07*100, 1), "\n")
cat("95% CI:    ", round(lower_lm_07*100, 1), "to", round(upper_lm_07*100, 1), "\n")

#-------------------------------------
### Monotherapy risk

risk_mono_lm_07 <- coef(msm_rd_lm_07)["(Intercept)"]
se_mono_lm_07   <- sqrt(vcovHC(msm_rd_lm_07, type = "HC3")["(Intercept)", "(Intercept)"])

lower_mono_lm_07 <- risk_mono_lm_07 - 1.96 * se_mono_lm_07
upper_mono_lm_07 <- risk_mono_lm_07 + 1.96 * se_mono_lm_07

cat("Monotherapy risk: ", round(risk_mono_lm_07*100, 1), "\n")
cat("95% CI: ", round(lower_mono_lm_07*100, 1), "to", round(upper_mono_lm_07*100, 1), "\n")

### Combination risk
#
result_combo_lm_07 <- deltaMethod(msm_rd_lm_07, "`(Intercept)` + treatcombo", vcov. = vcovHC(msm_rd_lm_07, type = "HC3"))

# %
result_combo_lm_07_scaled <- round(result_combo_lm_07 * 100, 1)
print(result_combo_lm_07_scaled)

#-------------------------------------
### Risk Ratio
### Combination ÷ Monotherapy
rr_result_lm_07 <- deltaMethod(msm_rd_lm_07, 
                               "(`(Intercept)` + treatcombo) / `(Intercept)`", 
                               vcov. = vcovHC(msm_rd_lm_07, type = "HC3"))

print(rr_result_lm_07)

cat("RR: ", round(rr_result_lm_07$Estimate, 1), "\n")
cat("95% CI: ", round(rr_result_lm_07$`2.5 %`, 1), "to", round(rr_result_lm_07$`97.5 %`, 1), "\n")


#-------------------------------------#----------------------------------------#



# =============================================================================
# 4. WEIGHTED OUTCOME MODEL  — Marginal Structural Model ~ Weighted LR (GLM)
# Binary treatment ~ combo vs mono
# 4.3. Binary outcome ~ 14-day ACM
# =============================================================================
### Baseline cohort
#-------------------------------------
#-------------------------------------
### Risk Difference — identity link
msm_rd_bl_14 <- glm(
  death_14 ~ treat,
  data    = final_data_bl2,
  family  = binomial(link = "identity"),
  weights = sw_treat1_99,
  start   = c(mean(final_data_bl2$death_14), 0)   # helps convergence
)

###
coeftest(msm_rd_bl_14, vcov = vcovHC(msm_rd_bl_14, type = "HC3"))
# coef["arm"] = RD directly

### Getting the 95% CI directly
# Extract from your existing model
rd_bl_14    <- coef(msm_rd_bl_14)["treatcombo"]
se_rd_bl_14 <- sqrt(vcovHC(msm_rd_bl_14, type = "HC3")["treatcombo", "treatcombo"])

lower_bl_14 <- rd_bl_14 - 1.96 * se_rd_bl_14
upper_bl_14 <- rd_bl_14 + 1.96 * se_rd_bl_14

cat("RD:        ", round(rd_bl_14*100, 1), "\n")
cat("95% CI:    ", round(lower_bl_14*100, 1), "to", round(upper_bl_14*100, 1), "\n")

#-------------------------------------
### Monotherapy risk

risk_mono_bl_14 <- coef(msm_rd_bl_14)["(Intercept)"]
se_mono_bl_14   <- sqrt(vcovHC(msm_rd_bl_14, type = "HC3")["(Intercept)", "(Intercept)"])

lower_mono_bl_14 <- risk_mono_bl_14 - 1.96 * se_mono_bl_14
upper_mono_bl_14 <- risk_mono_bl_14 + 1.96 * se_mono_bl_14

cat("Monotherapy risk: ", round(risk_mono_bl_14*100, 1), "\n")
cat("95% CI: ", round(lower_mono_bl_14*100, 1), "to", round(upper_mono_bl_14*100, 1), "\n")

### Combination risk
#
result_combo_bl_14 <- deltaMethod(msm_rd_bl_14, "`(Intercept)` + treatcombo", vcov. = vcovHC(msm_rd_bl_14, type = "HC3"))

# %
result_combo_bl_14_scaled <- round(result_combo_bl_14 * 100, 1)
print(result_combo_bl_14_scaled)

#-------------------------------------
### Risk Ratio
### Combination ÷ Monotherapy
rr_result_bl_14 <- deltaMethod(msm_rd_bl_14, 
                               "(`(Intercept)` + treatcombo) / `(Intercept)`", 
                               vcov. = vcovHC(msm_rd_bl_14, type = "HC3"))

print(rr_result_bl_14)

cat("RR: ", round(rr_result_bl_14$Estimate, 1), "\n")
cat("95% CI: ", round(rr_result_bl_14$`2.5 %`, 1), "to", round(rr_result_bl_14$`97.5 %`, 1), "\n")


#-------------------------------------#----------------------------------------#



### Landmark cohort
#-------------------------------------
#-------------------------------------
### Risk Difference — identity link
msm_rd_lm_14 <- glm(
  death_14 ~ treat,
  data    = final_data_lm2,
  family  = binomial(link = "identity"),
  weights = sw_treat2_99,
  start   = c(mean(final_data_lm2$death_14), 0)   # helps convergence
)

###
coeftest(msm_rd_lm_14, vcov = vcovHC(msm_rd_lm_14, type = "HC3"))
# coef["arm"] = RD directly

### Getting the 95% CI directly
# Extract from your existing model
rd_lm_14    <- coef(msm_rd_lm_14)["treatcombo"]
se_rd_lm_14 <- sqrt(vcovHC(msm_rd_lm_14, type = "HC3")["treatcombo", "treatcombo"])

lower_lm_14 <- rd_lm_14 - 1.96 * se_rd_lm_14
upper_lm_14 <- rd_lm_14 + 1.96 * se_rd_lm_14

cat("RD:        ", round(rd_lm_14*100, 1), "\n")
cat("95% CI:    ", round(lower_lm_14*100, 1), "to", round(upper_lm_14*100, 1), "\n")

#-------------------------------------
### Monotherapy risk

risk_mono_lm_14 <- coef(msm_rd_lm_14)["(Intercept)"]
se_mono_lm_14   <- sqrt(vcovHC(msm_rd_lm_14, type = "HC3")["(Intercept)", "(Intercept)"])

lower_mono_lm_14 <- risk_mono_lm_14 - 1.96 * se_mono_lm_14
upper_mono_lm_14 <- risk_mono_lm_14 + 1.96 * se_mono_lm_14

cat("Monotherapy risk: ", round(risk_mono_lm_14*100, 1), "\n")
cat("95% CI: ", round(lower_mono_lm_14*100, 1), "to", round(upper_mono_lm_14*100, 1), "\n")

### Combination risk
#
result_combo_lm_14 <- deltaMethod(msm_rd_lm_14, "`(Intercept)` + treatcombo", vcov. = vcovHC(msm_rd_lm_14, type = "HC3"))

# %
result_combo_lm_14_scaled <- round(result_combo_lm_14 * 100, 1)
print(result_combo_lm_14_scaled)

#-------------------------------------
### Risk Ratio
### Combination ÷ Monotherapy
rr_result_lm_14 <- deltaMethod(msm_rd_lm_14, 
                               "(`(Intercept)` + treatcombo) / `(Intercept)`", 
                               vcov. = vcovHC(msm_rd_lm_14, type = "HC3"))

print(rr_result_lm_14)

cat("RR: ", round(rr_result_lm_14$Estimate, 1), "\n")
cat("95% CI: ", round(rr_result_lm_14$`2.5 %`, 1), "to", round(rr_result_lm_14$`97.5 %`, 1), "\n")


#-------------------------------------#----------------------------------------#

#-------------------------------------#----------------------------------------#
#-------------------------------------#----------------------------------------#

### TD Plot for Publication
### Use 3_2_TD_plot Script

## -------------------------------------------------------------
## 4.2. Plot marginal KM-style survival curves
## Cumulative mortality
## -------------------------------------------------------------
### Step 1 — Fit weighted Kaplan-Meier for each arm
#-------------------------------------
### Weighted log-rank test
# Baseline cohort — ITT-analogue (IPW1 only)
km_bl <- survfit(
  Surv(time_to_death_21, death_21) ~ treat,
  data    = final_data_bl2,
  weights = sw_treat1_99
)

# Landmark cohort — PP-analogue (IPW2)
km_lm <- survfit(
  Surv(time_to_death_21, death_21) ~ treat,
  data    = final_data_lm2,
  weights = sw_treat2_99
)

#-------------------------------------#----------------------------------------#
### Step 2 — Extract KM data into a tidy dataframe
#-------------------------------------
extract_km <- function(km_fit, cohort_label) {
  
  # Extract summary
  km_summary <- summary(km_fit, times = sort(unique(km_fit$time)))
  
  data.frame(
    time   = km_summary$time,
    surv   = km_summary$surv,
    lower  = km_summary$lower,
    upper  = km_summary$upper,
    arm    = km_summary$strata,
    cohort = cohort_label
  ) %>%
    mutate(
      # Clean up strata label
      arm = case_when(
        grepl("combo", arm, ignore.case = TRUE) ~ "Combination therapy",
        grepl("mono",  arm, ignore.case = TRUE) ~ "Monotherapy",
        TRUE ~ as.character(arm)
      ),
      # Convert survival to mortality (1 - S(t))
      mortality = 1 - surv,
      mort_lo   = 1 - upper,   # note: CI bounds flip when inverting
      mort_hi   = 1 - lower
    )
}

df_km_bl <- extract_km(km_bl, "Baseline cohort (ITT-analogue)")
df_km_lm <- extract_km(km_lm, "Landmark cohort (PP-analogue)")

# Add time zero rows so curves start at 0
add_time_zero <- function(df) {
  arms <- unique(df$arm)
  t0 <- data.frame(
    time      = 0,
    surv      = 1,
    lower     = 1,
    upper     = 1,
    arm       = arms,
    cohort    = df$cohort[1],
    mortality = 0,
    mort_lo   = 0,
    mort_hi   = 0
  )
  bind_rows(t0, df) %>% arrange(arm, time)
}

df_km_bl <- add_time_zero(df_km_bl)
df_km_lm <- add_time_zero(df_km_lm)

# Combine both cohorts
df_km_all <- bind_rows(df_km_bl, df_km_lm)

#-------------------------------------#----------------------------------------#


#-------------------------------------
### Step 3 — Plot function (decorative 100% break, data untouched)
#-------------------------------------
plot_km_broken <- function(data, title_label, max_day = 21,
                           break_at = 0.5,     # top of the "real" linear scale (50%)
                           top_gap  = 0.08) {  # how far above break_at the "100%" tick sits
  
  top_pos <- break_at + top_gap
  
  ggplot(data, aes(x = time, y = mortality, colour = arm, fill = arm)) +
    geom_ribbon(aes(ymin = mort_lo, ymax = mort_hi), alpha = 0.15, colour = NA) +
    geom_step(linewidth = 0.9) +
    
    scale_colour_manual(
      values = c("Combination therapy" = "#d95f02", "Monotherapy" = "#2c7fb8"),
      name = NULL
    ) +
    scale_fill_manual(
      values = c("Combination therapy" = "#d95f02", "Monotherapy" = "#2c7fb8"),
      name = NULL
    ) +
    
    scale_x_continuous(
      limits = c(0, max_day),
      breaks = seq(0, max_day, by = 7),
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      limits = c(0, top_pos),
      breaks = c(seq(0, break_at, by = 0.1), top_pos),
      labels = c(
        scales::percent(seq(0, break_at, by = 0.1), accuracy = 1),
        "100%"
      ),
      expand = c(0, 0)
    ) +
    
    labs(
      title = title_label,
      x = "Days from treatment assignment",
      y = "Cumulative mortality"
    ) +
    
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "top",
      legend.justification = "left",
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.8, "lines"),
      legend.margin = margin(b = 2),
      panel.grid.major.y = element_line(
        colour = "grey92",
        linewidth = 0.3
      ),
      
      # Extra room on the left for the // symbol
      plot.margin = margin(15, 20, 15, 40)
    ) +
    
    #-------------------------------------
  # "//" break mark between 50% and 100%
  #-------------------------------------
  annotate(
    "text",
    x = 0,
    y = break_at + top_gap * 0.50,
    label = "//",
    angle = -45,
    hjust = 0.5,
    vjust = 0.5,
    size = 5,
    fontface = "bold",
    colour = "black"
  ) +
    
    coord_cartesian(
      clip = "off",
      ylim = c(0, top_pos)
    )
}


#-------------------------------------
### Step 4 — Generate plots
#-------------------------------------
p_bl <- plot_km_broken(df_km_bl, "Intention-to-treat population")

p_lm <- plot_km_broken(df_km_lm, "Per-protocol population") +
  theme(legend.position = "none")

combined <- p_bl / p_lm +
  plot_annotation(
    title = "Marginal Kaplan-Meier curves",
    caption = "Shaded areas represent 95% confidence intervals.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.caption = element_text(size = 8, colour = "grey50", hjust = 0)
    )
  )

combined


#-------------------------------------
# Save
# Without risk table (patchwork version)
ggsave("fig3_km_marginal_curves_21.png",
       plot   = combined,
       width  = 12.5,
       height = 6,
       dpi    = 320,
       bg = "white")

#-------------------------------------#----------------------------------------#
#-------------------------------------#----------------------------------------#



###--------------------------------------------------------------------------###
### ^.^ ###