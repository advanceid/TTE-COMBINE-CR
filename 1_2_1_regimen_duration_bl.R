# =========================================================
# Calculate duration of matched treatment regimen using
# day-by-day deviation flags from the DAY-0 (baseline)
# observed regimen.
#
# Pipeline:
#   1. Load raw start/end data (recordid, anti_start, anti_end,
#      anti_names, time_zero, matched_treat_combo)
#   2. Expand into one row per id per calendar day, classify
#      each day's strategy (Monotherapy / Combination therapy /
#      Unclassified) and treat_combo (sorted drug-name string)
#   3. Flag deviations from day-0 regimen (discontinued,
#      unclassified, switched, added, or mono->different mono)
#   4. is_any_deviation_today is the SINGLE censoring trigger;
#      the other flags only label WHY a day deviated
#   5. Duration = last day_since_matched before first deviation
#
# Strategy classification rule (edit if yours differs):
# Following treatment groups
# =========================================================
# Clear
rm(list = ls())

###
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# ---- 1. Load raw data ----
###
df <- readRDS("data/clean data/clean/clean/combine_group/treat_dur.RData")
df_bl <- readRDS("data/clean data/clean/clean/final/tables_bl.RData")

### Filter IDs
df <- df[which(df$recordid %in% df_bl$id),]

# ---- 2. Expand to one row per id per calendar day, then classify ----
###
build_long <- function(sub) {
  tz       <- unique(sub$time_zero)
  last_day <- max(sub$anti_end)
  days     <- seq(tz, last_day, by = "day")
  map_dfr(days, function(d) {
    active_rows <- sub %>%
      filter(anti_start <= d & anti_end >= d)
    
    active <- active_rows %>%
      pull(anti_names) %>%
      unique() %>%
      sort()
    n_drugs <- length(active)
    
    # Pull treatment_strategy directly from df instead of recomputing from n_drugs
    treatment_strategy <- if (n_drugs == 0) {
      NA_character_
    } else {
      active_rows %>%
        pull(treatment_strategy) %>%
        unique() %>%
        paste(collapse = "; ")   # flags if rows disagree on treatment_strategy for the same day
    }
    
    treat_combo <- if (n_drugs == 0) NA_character_ else paste(active, collapse = "+")
    
    tibble(
      date               = d,
      day_since_matched  = as.numeric(d - tz),
      strategy = treatment_strategy,
      treat_combo        = treat_combo
    )
  })
}
###
long_data <- df %>%
  group_by(recordid) %>%
  group_modify(~ build_long(.x)) %>%
  ungroup()


# ---- Optional sanity check: does observed day-0 regimen match the ----
# ---- originally assigned matched_treat_combo column?               ----
day0_check <- long_data %>%
  filter(day_since_matched == 0) %>%
  dplyr::select(recordid, treat_combo) %>%
  left_join(df %>% distinct(recordid, matched_treat_combo), by = "recordid") %>%
  mutate(
    combo_day0_sorted = map_chr(str_split(treat_combo, "\\+"), ~ paste(sort(str_trim(.x)), collapse = "+")),
    combo_target_sorted = map_chr(str_split(matched_treat_combo, "\\+"), ~ paste(sort(str_trim(.x)), collapse = "+")),
    mismatch = combo_day0_sorted != combo_target_sorted
  ) %>%
  filter(mismatch)

if (nrow(day0_check) > 0) {
  warning(nrow(day0_check), " id(s) have day-0 observed regimen != matched_treat_combo. See 'day0_check'.")
}

# ---- 3. Flag deviations from day-0 (baseline) strategy and regimen ----
long_flagged <- long_data %>%
  group_by(recordid) %>%
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

# ---- 4. Label the reason for each deviating day ----
long_flagged <- long_flagged %>%
  mutate(
    deviation_reason = case_when(
      !is_any_deviation_today ~ NA_character_,
      is_discont_today        ~ "Discontinued",
      is_unclassified_today   ~ "Unclassified",
      is_switch_today         ~ "Switched (mono<->combo)",
      is_added_today          ~ "Added/changed drug (combo->combo)",
      TRUE                    ~ "Changed drug (mono->different mono)"
    )
  )

# ---- 5. Censor each id at the FIRST deviation; compute duration ----
### include end_date in the total duration
###
result <- long_flagged %>%
  group_by(recordid) %>%
  summarise(

    matched_treat_combo = first(treat_combo_day0_ref),
    
    treatment_strategy = first(strategy),

    first_deviation_day = if (any(is_any_deviation_today, na.rm = TRUE)) {
      min(day_since_matched[is_any_deviation_today], na.rm = TRUE)
    } else {
      NA_real_
    },

    # duration = last day the patient followed the assigned regimen
    duration_days = if (!is.na(first_deviation_day)) {
      first_deviation_day  # also calculate the time zero
    } else {
      max(day_since_matched) + 1  # regimen held through all available data, plus 1 as starting from D0
    },

    deviation_reason = if (!is.na(first_deviation_day)) {
      deviation_reason[day_since_matched == first_deviation_day]
    } else {
      NA_character_            # no deviation observed -> not censored
    },

    .groups = "drop"
  )

###
table(is.na(result$duration_days))

# ---- Summary of duration_days by treatment_strategy ----
result_01 <- result %>%
  group_by(treatment_strategy) %>%
  summarise(
    n      = sum(!is.na(duration_days)),
    min    = min(duration_days, na.rm = TRUE),
    q1     = quantile(duration_days, 0.25, na.rm = TRUE),
    median = median(duration_days, na.rm = TRUE),
    mean   = mean(duration_days, na.rm = TRUE),
    q3     = quantile(duration_days, 0.75, na.rm = TRUE),
    max    = max(duration_days, na.rm = TRUE),
    sd     = sd(duration_days, na.rm = TRUE),
    .groups = "drop"
  )


# ---- 6. Inspect / export ----
#print(result)
# write.csv(long_flagged, "daily_flagged_data.csv", row.names = FALSE)
# write.csv(result, "matched_regimen_duration_with_reasons.csv", row.names = FALSE)
