
### Treatment strategy: Initiate Combination Therapy or Monotherapy on or within 1 day after the AST report
### 20260921

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
                 writexl)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean data/clean/baseline_all.RData")
df_antib <- readRDS("data/clean data/general/anti_treat_prior_01.RData")
df_ast <- readRDS("data/clean data/general/ast_index_01.RData")

### Load fonts
loadfonts()

### combined BSI types
df$infection_types <- ifelse(df$infection_types %in% c("Hospital-acquired BSI", "Healthcare-associated BSI"),
                             "BSI",
                             df$infection_types)

#-------------------------------------#----------------------------------------#



#-------------------------------------
#-------------------------------------
### Create day_by_day dataset (Thank Luca!)
# Select related variables
anti_treat_day <- df_antib %>%
  dplyr::select("recordid", "anti_start", "anti_end", "anti_names", 
         "inf_onset", "spec_date") 

# We assumed that AST date ~ spec_date + 3
anti_treat_day$ast_zero <- anti_treat_day$spec_date + 3

#-------------------------------------
# No NA in anti_names

### Check for duplicates with anti_names, anti_start in the same recordid (typo)
#
dup_check <- anti_treat_day %>%
  group_by(recordid, anti_names, anti_start) %>%
  filter(n() > 1)

# Keep only the row with the latest (max) anti_end when duplicates exist
anti_treat_day <- anti_treat_day %>%
  group_by(recordid, anti_names, anti_start) %>%
  slice_max(order_by = anti_end, with_ties = FALSE) %>%
  ungroup()

# anti_start <= anti_end (No NA in anti_start, anti_end)
anti_treat_day <- anti_treat_day %>%
  filter(anti_end >= anti_start)

#-------------------------------------
# Expand each treatment episode into one row per day of exposure
anti_treat_day <- anti_treat_day %>%
  rowwise() %>% # group df by row
  mutate(days_seq = list(seq.Date(anti_start, anti_end, by = "day"))) %>%
  unnest(days_seq) %>%
  mutate(day_relative_to_ast = as.integer(days_seq - ast_zero)) %>%
  ungroup()

#-------------------------------------
### Antibiotics Window
### Initiate within 2 days of ast & sustain >= 3 days since the initiation
anti_treat_day_used <- anti_treat_day %>%
  filter(day_relative_to_ast >= -3 & day_relative_to_ast <= 4)

#-------------------------------------#-------------------------------------####



#-------------------------------------
# Mono or Combo by day
treat_combo <- anti_treat_day_used %>%
  group_by(recordid, day_relative_to_ast, days_seq) %>%
  summarise(
    treat_combo = paste(sort(unique(anti_names)), collapse = "+"),
    .groups = "drop"
  ) %>%
  as.data.frame()

#-------------------------------------
# Define a function to sort alphabetically and standardize combinations
standardize_combo <- function(combo) {
  elements <- unlist(strsplit(combo, ", "))
  standardized <- paste(sort(elements), collapse = ", ")
  return(standardized)
}

# Apply the function to the pathogen_group_combined column
treat_combo$treat_combo <- sapply(treat_combo$treat_combo , standardize_combo)

#-------------------------------------
# Combine pathogens
treat_combo <- left_join(treat_combo, df[, c("recordid", "org_r_combined")], by = "recordid")

# Only monomicrobial
treat_combo_used <- treat_combo %>%
  filter(org_r_combined %in% c("CRA", "CRE", "CRP"))

#-------------------------------------
### Save clean data for next
saveRDS(treat_combo_used, file = "data/clean data/clean/clean/combine_group/treat_combo.RData")


#-------------------------------------#----------------------------------------#



# ── 1. Define drug components per organism ───────────────────────────────────
# Monotherapy: exact match (single drug only)
# Combination: substring-based — treat_combo must contain at least one drug
#              from each required component group

# CRA -------------------------------------------------------------------------
###
mono_CRA <- c("Colistin", "Polymyxin B", "Ampicillin/sulbactam", "Cefoperazone/sulbactam", "Tigecycline", "Eravacycline")

# Combination requires: (Colistin OR Polymyxin B) + (Ampicillin/sulbactam OR Cefoperazone/sulbactam)
combo_CRA_backbone <- c("Colistin", "Polymyxin B")
combo_CRA_partner  <- c("Meropenem", "Imipenem", "Doripenem", "Ampicillin/sulbactam", "Cefoperazone/sulbactam", "Tigecycline", "Eravacycline")

# CRE -------------------------------------------------------------------------
###
mono_CRE <- c("Colistin", "Polymyxin B", "Ceftazidime/avibactam", "Tigecycline", "Eravacycline")

# Combination requires: 2 pairs of backbone and partner
# pair 1
combo_CRE_backbone_A <- c("Colistin", "Polymyxin B")
combo_CRE_partner_A  <- c("Meropenem", "Imipenem", "Doripenem", "Tigecycline", "Eravacycline")

# pair 2
combo_CRE_backbone_B <- c("Ceftazidime/avibactam")
combo_CRE_partner_B  <- c("Tigecycline", "Eravacycline", "Amikacin", "Gentamicin", "Tobramycin")

# CRP -------------------------------------------------------------------------
###
mono_CRP <- c("Colistin", "Polymyxin B", "Ceftazidime/avibactam", "Ceftolozane/tazobactam")

# Combination requires: 3 pairs of backbone and partner
# pair 1
combo_CRP_backbone_A <- c("Colistin", "Polymyxin B")
combo_CRP_partner_A  <- c("Meropenem", "Imipenem", "Doripenem")

# pair 2
combo_CRP_backbone_B <- c("Ceftazidime/avibactam")
combo_CRP_partner_B  <- c("Amikacin", "Gentamicin", "Tobramycin")

# pair 3
combo_CRP_backbone_C <- c("Ceftolozane/tazobactam")
combo_CRP_partner_C  <- c("Amikacin", "Gentamicin", "Tobramycin")

#-------------------------------------#----------------------------------------#



# ── 2. Helper: check if a string contains any item from a list ───────────────
###
contains_any <- function(string, items) {
  if (is.na(string)) return(FALSE)
  any(sapply(items, function(x) grepl(x, string, fixed = TRUE)))
}

#-------------------------------------#----------------------------------------#



# ── 3. Classification function ───────────────────────────────────────────────
### Included matched day_relative_to_ast & matched treat_combo
###
classify_one <- function(t, org_up) {
  if (org_up == "CRA") {
    if (t %in% mono_CRA) return("Monotherapy")
    if (contains_any(t, combo_CRA_backbone) && contains_any(t, combo_CRA_partner)) return("Combination therapy")
  } else if (org_up == "CRE") {
    if (t %in% mono_CRE) return("Monotherapy")
    if ((contains_any(t, combo_CRE_backbone_A) && contains_any(t, combo_CRE_partner_A)) ||
        (contains_any(t, combo_CRE_backbone_B) && contains_any(t, combo_CRE_partner_B))) return("Combination therapy")
  } else if (org_up == "CRP") {
    if (t %in% mono_CRP) return("Monotherapy")
    if ((contains_any(t, combo_CRP_backbone_A) && contains_any(t, combo_CRP_partner_A)) ||
        (contains_any(t, combo_CRP_backbone_B) && contains_any(t, combo_CRP_partner_B)) ||
        (contains_any(t, combo_CRP_backbone_C) && contains_any(t, combo_CRP_partner_C))) return("Combination therapy")
  }
  return(NA_character_)
}

### Returns a one-row data.frame: strategy + the day/treat_combo that triggered it
classify_strategy_full <- function(days, treats, org) {
  org_up <- toupper(trimws(org))
  
  if (!org_up %in% c("CRA", "CRE", "CRP")) {
    return(data.frame(treatment_strategy = "Unclassified",
                      matched_day = NA_real_,
                      matched_treat_combo = NA_character_))
  }
  
  for (i in seq_along(treats)) {
    res <- classify_one(treats[i], org_up)
    if (!is.na(res)) {
      return(data.frame(treatment_strategy = res,
                        matched_day = days[i],
                        matched_treat_combo = treats[i]))
    }
  }
  
  data.frame(treatment_strategy = "Unclassified",
             matched_day = NA_real_,
             matched_treat_combo = NA_character_)
}

#-------------------------------------#----------------------------------------#



# ── 4. Apply to dataset ──────────────────────────────────────────────────────
###
assign_treatment_strategy <- function(df) {
  result <- df %>%
    filter(day_relative_to_ast %in% c(0, 1)) %>%
    arrange(recordid, org_r_combined, day_relative_to_ast) %>%
    group_by(recordid, org_r_combined) %>%
    group_modify(~ classify_strategy_full(.x$day_relative_to_ast, .x$treat_combo, .y$org_r_combined)) %>%
    ungroup() %>%
    as.data.frame()
  
  return(result)
}

#-------------------------------------#----------------------------------------#
# Extract df
df_patient_init <- assign_treatment_strategy(treat_combo_used)

#
table(df_patient_init$treatment_strategy)

#-------------------------------------#----------------------------------------#



### Any use of matched treat_combo for 3 days before ast?
#-------------------------------------
#-------------------------------------
df_patient_prior_3d <- df_patient_init %>%
  left_join(treat_combo_used, by = c("recordid", "org_r_combined")) %>%
  group_by(recordid, org_r_combined, treatment_strategy, matched_day, matched_treat_combo) %>%
  summarise(
    used_all_3_days_before = as.integer(
      all(c(-3, -2, -1) %in% day_relative_to_ast[treat_combo == matched_treat_combo])
    ),
    .groups = "drop"
  ) %>%
  mutate(
    used_all_3_days_before = ifelse(used_all_3_days_before, 1, 0)
  )

#
table(df_patient_prior_3d$used_all_3_days_before)

#
df_patient_excl <- df_patient_prior_3d %>%
  filter(used_all_3_days_before == 1)

#-------------------------------------
### Exclude IDs with prior use of matched treat_combo >= 3 days consecutive
### Final df for Initiated strategy ### ***
df_filtered_init <- df_patient_init %>%
  anti_join(df_patient_excl, by = "recordid")

#
table(df_filtered_init$treatment_strategy)

#-------------------------------------#----------------------------------------#
#-------------------------------------#----------------------------------------#



### Time-to-event
### Primary Outcome as 14-day ACM
### Secondary Outcomes as 7-day ACM & 21-day ACM
#-------------------------------------
#-------------------------------------
### filter IDs
df_event <- df_filtered_init %>%
  filter(
    treatment_strategy %in% c("Combination therapy", "Monotherapy")
  ) 

#-------------------------------------
### matched_day ~ time zero or randomized_date
###
treat_combo_event <- treat_combo_used %>%
  filter(day_relative_to_ast %in% c(0,1))


### Matching to find the time zero
df_event2 <- df_event %>%
  left_join(
    treat_combo_event %>%
      dplyr::select(recordid, day_relative_to_ast, treat_combo, days_seq),
    by = c(
      "recordid" = "recordid",
      "matched_day" = "day_relative_to_ast",
      "matched_treat_combo" = "treat_combo"
    )
  ) %>%
  rename(time_zero = days_seq)

#-------------------------------------
### day 21 new
df_event2$d21_date_new <- df_event2$time_zero + 21

### day 7 new 
df_event2$d07_date_new <- df_event2$time_zero + 7

### day 14 new 
df_event2$d14_date_new <- df_event2$time_zero + 14

#-------------------------------------
# combine mortality_date
df_event2 <- left_join(df_event2, df[, c("recordid", "infection_types", "mortality_date")], by = "recordid")

# relocate variable
df_event2 <- df_event2 %>%
  relocate(infection_types, .after = org_r_combined)

#-------------------------------------
### check 21-day mortality status within new day 21
#
df_event2 <- df_event2 %>%
  mutate(
    d21_observed_date = case_when(
      is.na(mortality_date) ~ d21_date_new,
      TRUE ~ pmin(d21_date_new, mortality_date)
    ),
    d21_death = case_when(
      d21_observed_date == mortality_date ~ 1,
      !is.na(d21_observed_date) ~ 0,
      TRUE ~ NA_integer_
    ),
    d21_random_days = as.numeric(d21_observed_date - df_event2$time_zero)
  )

#-------------------------------------
### check 7-day mortality status
#
df_event2 <- df_event2 %>%
  mutate(
    d07_observed_date = case_when(
      is.na(mortality_date) ~ d07_date_new,
      TRUE ~ pmin(d07_date_new, mortality_date)
    ),
    d07_death = case_when(
      d07_observed_date == mortality_date ~ 1,
      !is.na(d07_observed_date) ~ 0,
      TRUE ~ NA_integer_
    ),
    d07_random_days = as.numeric(d07_observed_date - df_event2$time_zero)
  )

#-------------------------------------
### check 14-day mortality status
#
df_event2 <- df_event2 %>%
  mutate(
    d14_observed_date = case_when(
      is.na(mortality_date) ~ d14_date_new,
      TRUE ~ pmin(d14_date_new, mortality_date)
    ),
    d14_death = case_when(
      d14_observed_date == mortality_date ~ 1,
      !is.na(d14_observed_date) ~ 0,
      TRUE ~ NA_integer_
    ),
    d14_random_days = as.numeric(d14_observed_date - df_event2$time_zero)
  )

#-------------------------------------
#-------------------------------------
### first_discharge_date = min(ho_discharge_date, mortality_date)
df_event2 <- left_join(df_event2, df[, c("recordid", "first_discharge_date")], by = "recordid")

#-------------------------------------#----------------------------------------#
#-------------------------------------#----------------------------------------#





### Prepare dataset for IPW1 and IPW2
## =============================================================
## 1.1. Baseline cohort
## Keep only patients who are alive and not discharged by day 1 after time zero
## =============================================================
### Discharged and Died separately
df_ipw_bl <- df_event2 %>%
  mutate(
    time_zero = as.Date(time_zero),
    mortality_date = as.Date(mortality_date),
    first_discharge_date = as.Date(first_discharge_date),
    day1_date = time_zero + 1,
    
    death_within_1day = !is.na(mortality_date) &
      mortality_date >= time_zero &
      mortality_date <= day1_date,
    
    discharge_within_1day = !is.na(first_discharge_date) &
      first_discharge_date >= time_zero &
      first_discharge_date <= day1_date &
      !(!is.na(mortality_date) & first_discharge_date == mortality_date),
    
    event_type_within_1day = case_when(
      death_within_1day ~ "death",
      discharge_within_1day ~ "discharge",
      TRUE ~ "none"
    )
  )

#-------------------------------------
### relocate mortality_date
df_ipw_bl <- df_ipw_bl %>%
  relocate(mortality_date, .before = first_discharge_date)

###
df_ipw_bl <- df_ipw_bl %>%
  relocate(time_zero, .before = mortality_date)

###
df_ipw_bl <- df_ipw_bl %>%
  relocate(day1_date, .before = mortality_date)

###
df_ipw_bl <- df_ipw_bl %>%
  relocate(d21_date_new, .before = mortality_date)


#-------------------------------------
###
table(df_ipw_bl$treatment_strategy, df_ipw_bl$event_type_within_1day)


###
df_ipw_bl_used <- df_ipw_bl %>%
  filter(event_type_within_1day == "none")

#-------------------------------------
###
df_event3 <- df_event2[which(df_event2$recordid %in% df_ipw_bl_used$recordid),]


#-------------------------------------#----------------------------------------#


### Continue from the ITT population (639)
## =============================================================
## 1.2. Day-3 landmark cohort
## Keep only patients who are alive and not discharged by day 3 after time zero, 
## then define follow-up starting at day 3.
## =============================================================
###
df_ipw_lm <- df_event3 %>%
  mutate(
    time_zero = as.Date(time_zero),
    mortality_date = as.Date(mortality_date),
    first_discharge_date = as.Date(first_discharge_date),
    day3_date = time_zero + 3,
    
    death_within_3day = !is.na(mortality_date) &
      mortality_date >= time_zero &
      mortality_date <= day3_date,
    
    discharge_within_3day = !is.na(first_discharge_date) &
      first_discharge_date >= time_zero &
      first_discharge_date <= day3_date &
      !(!is.na(mortality_date) & first_discharge_date == mortality_date),
    
    event_type_within_3day = case_when(
      death_within_3day ~ "death",
      discharge_within_3day ~ "discharge",
      TRUE ~ "none"
    )
  )

#-------------------------------------
### relocate mortality_date
df_ipw_lm <- df_ipw_lm %>%
  relocate(mortality_date, .before = first_discharge_date)

###
df_ipw_lm <- df_ipw_lm %>%
  relocate(time_zero, .before = mortality_date)

###
df_ipw_lm <- df_ipw_lm %>%
  relocate(day3_date, .before = mortality_date)

###
df_ipw_lm <- df_ipw_lm %>%
  relocate(d21_date_new, .before = mortality_date)


#-------------------------------------
###
table(df_ipw_lm$treatment_strategy, df_ipw_lm$event_type_within_3day)


### Survivors from Day 4
df_ipw_lm_used <- df_ipw_lm %>%
  filter(event_type_within_3day == "none")

#-------------------------------------#----------------------------------------#
#-------------------------------------#----------------------------------------#



### COMBINE-CR
### Initiate Combination Therapy vs Monotherapy on day 0-1
### df_filtered_init
#-------------------------------------
###
id_used <- df_filtered_init %>%
  filter(
    treatment_strategy %in% c("Combination therapy", "Monotherapy")
  ) %>%
  distinct(recordid) %>%
  pull(recordid)

### filter df_all
df_used <- df[which(df$recordid %in% id_used),]

#-------------------------------------
### combine treatment arms
df_used <- df_used %>%
  left_join(
    df_filtered_init %>%
      dplyr::select(
        recordid,
        treatment_strategy,
        matched_day,
        matched_treat_combo
      ),
    by = "recordid"
  )

# relocate org_r_combined
df_used <- df_used %>%
  relocate(org_r_combined, .after = WHO_group)


### combine time_zero
df_used <- df_used %>%
  left_join(
    df_event2 %>%
      dplyr::select(
        recordid,
        time_zero, d21_date_new, d21_death, d21_random_days,
        d07_date_new, d07_death, d07_random_days,
        d14_date_new, d14_death, d14_random_days,
      ),
    by = "recordid"
  )

#-------------------------------------#----------------------------------------#
#-------------------------------------#----------------------------------------#



### 3-day Landmark: Sustained
#-------------------------------------
#-------------------------------------
# ── 5. SUSTAIN function ───────────────────────────────────────────────
###
get_sustained_info <- function(matched_day, matched_treat_combo, org, days, treats) {
  org_up <- toupper(trimws(org))
  
  day_lookup <- function(d) {
    idx <- which(days == d)
    if (length(idx) == 0) return(c(NA_character_, NA_character_))
    tc <- treats[idx[1]]
    strat <- classify_one(tc, org_up)
    if (is.na(strat)) strat <- "Unclassified"
    c(tc, strat)
  }
  
  d0 <- day_lookup(matched_day)
  d1 <- day_lookup(matched_day + 1)
  d2 <- day_lookup(matched_day + 2)
  d3 <- day_lookup(matched_day + 3)
  
  data.frame(
    treat_combo_day0 = d0[1], strategy_day0 = d0[2],
    treat_combo_day1 = d1[1], strategy_day1 = d1[2],
    treat_combo_day2 = d2[1], strategy_day2 = d2[2],
    treat_combo_day3 = d3[1], strategy_day3 = d3[2],
    used_3d_since_matched_day = as.integer(
      all(c(matched_day, matched_day + 1, matched_day + 2, matched_day + 3) %in%
            days[treats == matched_treat_combo])
    )
  )
}

# ── 6. Apply to dataset ──────────────────────────────────────────────────────
### Survivors from Day 4
df_lm_init <- df_filtered_init[which(df_filtered_init$recordid %in% df_ipw_lm_used$recordid),]

###
df_patient_sust_3d <- df_lm_init %>%
  left_join(treat_combo_used, by = c("recordid", "org_r_combined")) %>% # combine treat_combo_used ~ antib dataset by day
  group_by(recordid, org_r_combined, treatment_strategy, matched_day, matched_treat_combo) %>%
  group_modify(~ get_sustained_info(
    matched_day         = .y$matched_day,
    matched_treat_combo = .y$matched_treat_combo,
    org                  = .y$org_r_combined,
    days                 = .x$day_relative_to_ast,
    treats               = .x$treat_combo
  )) %>%
  ungroup()

#-------------------------------------#----------------------------------------#

## -------------------------------------------------------------
## 1. Reshape to long: one row per id x day_since_matched
## -------------------------------------------------------------
#
wide_data <- df_patient_sust_3d %>%
  rename(id = recordid)

#
wide_data <- wide_data %>%
  mutate(treat = case_when(
    treatment_strategy == "Combination therapy" ~ "combo",
    treatment_strategy == "Monotherapy" ~ "mono",
    TRUE ~ NA_character_
  ))

#
day_cols <- wide_data %>%
  dplyr::select(matches("^(strategy|treat_combo)_day\\d+$")) %>%
  names()

#
long_data <- wide_data %>%
  dplyr::select(id, all_of(day_cols)) %>%
  pivot_longer(
    cols = -id,
    names_to = c(".value", "day_since_matched"),
    names_pattern = "^(strategy|treat_combo)_day(\\d+)$"
  ) %>%
  mutate(day_since_matched = as.integer(day_since_matched)) %>%
  arrange(id, day_since_matched)

###
length(unique(long_data$id))

## sanity check: every id should have exactly (max day + 1) rows,
## and exactly one row per id where day_since_matched == 0
stopifnot(
  long_data %>% count(id) %>% pull(n) %>% unique() %>% length() == 1,
  long_data %>% filter(day_since_matched == 0) %>% count(id) %>% pull(n) %>% max() == 1
)

## -------------------------------------------------------------
## 2. Within each id, compare every day to the DAY-0 (baseline)
##    strategy and regimen, and flag the type of deviation.
##
##    KEY CHANGE: the underlying trigger is now a single exact-
##    match test (treat_combo today == treat_combo day0). The
##    switched/added/unclassified/discont flags are derived FROM
##    that same trigger, purely to LABEL what kind of deviation
##    it was -- they do not independently drive censoring.
## -------------------------------------------------------------
### Excellent code!
long_flagged <- long_data %>%
  group_by(id) %>%
  mutate(
    strategy_day0_ref     = strategy[day_since_matched == 0],
    treat_combo_day0_ref  = treat_combo[day_since_matched == 0],
    
    ## --- raw per-day deviation flags (only meaningful for day >= 1) ---
    
    # discontinued: regimen goes NA today, but was not NA on day 0
    is_discont_today = day_since_matched >= 1 &
      is.na(treat_combo) & !is.na(treat_combo_day0_ref),
    
    # unclassified: strategy becomes "Unclassified" today (its own
    # deviation category, regardless of day-0 strategy)
    is_unclassified_today = day_since_matched >= 1 &
      !is.na(strategy) & strategy == "Unclassified",
    
    # switched: strategy flips Monotherapy <-> Combination therapy
    # relative to day-0 strategy (excludes Unclassified, labeled above)
    is_switch_today = day_since_matched >= 1 &
      !is.na(strategy) & !is.na(strategy_day0_ref) &
      strategy %in% c("Monotherapy", "Combination therapy") &
      strategy_day0_ref %in% c("Monotherapy", "Combination therapy") &
      strategy != strategy_day0_ref,
    
    # added: stayed "Combination therapy" both at day0 and today,
    # but the regimen string changed (drug added/changed) vs day 0
    is_added_today = day_since_matched >= 1 &
      !is.na(strategy) & strategy == "Combination therapy" &
      strategy_day0_ref == "Combination therapy" &
      !is.na(treat_combo) & !is.na(treat_combo_day0_ref) &
      treat_combo != treat_combo_day0_ref,
    
    # CORE TRIGGER: any departure at all from the exact day-0
    # regimen string (this is what actually drives censoring;
    # the four flags above just explain WHY on a given day)
    is_any_deviation_today = is_discont_today | is_unclassified_today |
      is_switch_today | is_added_today |
      (day_since_matched >= 1 &
         !is.na(treat_combo) & !is.na(treat_combo_day0_ref) &
         treat_combo != treat_combo_day0_ref &
         strategy == "Monotherapy" & strategy_day0_ref == "Monotherapy")
    # ^ catches mono->different-drug-but-still-mono, e.g.
    #   Colistin -> Amikacin, both labeled "Monotherapy"
  ) %>%
  ungroup()


## -------------------------------------------------------------
## 3. Collapse to "first deviation wins, then sticks" per id,
##    restricted to the SUSTAIN WINDOW (day 1 to day 3 since
##    matched_day). Adjust `sustain_window_end` if your window
##    differs.
## -------------------------------------------------------------

sustain_window_end <- 3

### Type of deviations ###
deviation_summary <- long_flagged %>%
  filter(day_since_matched >= 1, day_since_matched <= sustain_window_end) %>%
  group_by(id) %>%
  summarise(
    first_discont_day      = suppressWarnings(min(day_since_matched[is_discont_today],      na.rm = TRUE)),
    first_switch_day       = suppressWarnings(min(day_since_matched[is_switch_today],        na.rm = TRUE)),
    first_added_day        = suppressWarnings(min(day_since_matched[is_added_today],         na.rm = TRUE)),
    first_unclassified_day = suppressWarnings(min(day_since_matched[is_unclassified_today],  na.rm = TRUE)),
    first_any_dev_day      = suppressWarnings(min(day_since_matched[is_any_deviation_today], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    across(starts_with("first_"), ~ ifelse(is.infinite(.x), NA_integer_, .x)),
    
    switched_since_matched_day      = as.integer(!is.na(first_switch_day)),
    added_since_matched_day         = as.integer(!is.na(first_added_day)),
    unclassified_since_matched_day  = as.integer(!is.na(first_unclassified_day)),
    discont_since_matched_day       = as.integer(!is.na(first_discont_day)),
    dev_since_matched_day           = as.integer(!is.na(first_any_dev_day)),
  ) %>%
  dplyr::select(id, switched_since_matched_day, added_since_matched_day,
         unclassified_since_matched_day, discont_since_matched_day, dev_since_matched_day, 
         first_switch_day, first_added_day, first_unclassified_day, 
         first_discont_day, first_any_dev_day)

## -------------------------------------------------------------
## 4. Join back to wide_data
## -------------------------------------------------------------
###
wide_data <- wide_data %>%
  dplyr::select(id, treat) %>%             
  left_join(deviation_summary %>%
              dplyr::select(
                id,
                dev_since_matched_day,
                first_any_dev_day
              ), by = "id") 

###
table(wide_data$treat, wide_data$dev_since_matched_day)


###
wide_data <- wide_data %>%
  rename(recordid = id)

## -------------------------------------------------------------
### Sustained for landmark dataset
df_lm_used <- df_used[which(df_used$recordid %in% df_ipw_lm_used$recordid),]

###
df_lm_used <- left_join(df_lm_used, wide_data[, c("recordid", "dev_since_matched_day", "first_any_dev_day")], 
                     by = "recordid")

###
table(df_lm_used$treatment_strategy, df_lm_used$dev_since_matched_day)


### relocate
df_lm_used <- df_lm_used %>% 
  relocate(dev_since_matched_day, .after = matched_treat_combo)

df_lm_used <- df_lm_used %>% 
  relocate(first_any_dev_day, .after = dev_since_matched_day)

### combine relevant variables
df_lm_used <- left_join(df_lm_used, df_ipw_lm[, c("recordid", "day3_date", 
                                                        "event_type_within_3day")], by = "recordid")

## -------------------------------------------------------------
### Baseline cohort (Overall dataset)
df_bl_used <- df_used[which(df_used$recordid %in% df_ipw_bl_used$recordid),]

#-------------------------------------#----------------------------------------#



### Duration of assigned treatment strategy for Baseline cohort
#-------------------------------------
#-------------------------------------
#
df_01 <- df_antib %>%
  dplyr::select("recordid", 
         "anti_start", "anti_end",
         "anti_names"
         )

# 691
df_01 <- df_01[which(df_01$recordid %in% df_ipw_bl_used$recordid),]

#
df_01 <- left_join(
  df_01,
  df_ipw_bl_used %>% dplyr::select(recordid, time_zero, treatment_strategy, matched_treat_combo, org_r_combined),
  by = "recordid"
)

#-------------------------------------#----------------------------------------#



#-------------------------------------
#-------------------------------------
### Save clean data for next
saveRDS(df_event2, file = "data/clean data/clean/clean/combine_group/event_KM.RData")
saveRDS(df_used, file = "data/clean data/clean/clean/combine_group/baseline_all_combine.RData")
saveRDS(df_bl_used, file = "data/clean data/clean/clean/combine_group/baseline_all_combine_bl.RData")
saveRDS(df_lm_used, file = "data/clean data/clean/clean/combine_group/baseline_all_combine_lm.RData")
saveRDS(df_ipw_bl_used, file = "data/clean data/clean/clean/combine_group/ipw/ipw_bl.RData")
saveRDS(df_ipw_lm_used, file = "data/clean data/clean/clean/combine_group/ipw/ipw_lm.RData")
saveRDS(df_ipw_lm, file = "data/clean data/clean/clean/combine_group/ipw/ipw_lm_all.RData")

# Treatment duration 
saveRDS(df_01, file = "data/clean data/clean/clean/combine_group/treat_dur.RData")

#-------------------------------------#----------------------------------------#





###-----------------------------------------------------------------------------
### ^.^ ###