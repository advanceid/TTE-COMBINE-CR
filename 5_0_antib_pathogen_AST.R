
# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(openxlsx,
                 data.table,
                 tidyverse,
                 devtools,
                 magrittr,
                 dplyr,
                 countrycode,
                 lubridate,
                 purrr,
                 stringr,
                 rlang,
                 ggpubr,
                 plotly,
                 janitor,
                 snakecase)
}) 

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean data/clean/clean/combine_group/baseline_all_clean_combine.RData")
id_bl <- readRDS("data/clean data/clean/clean/final/tables_bl.RData")
id_lm <- readRDS("data/clean data/clean/clean/final/tables_lm.RData")
ast <- readRDS("data/clean data/ast_all_index.RData")

#-------------------------------------
### Filter IDs
ast_bl <- ast[which(ast$recordid %in% id_bl$id),]
ast_lm <- ast[which(ast$recordid %in% id_lm$id),]
df_bl <- df[which(df$recordid %in% id_bl$id),]
df_lm <- df[which(df$recordid %in% id_lm$id),]

#-------------------------------------
### 
#
ast_bl <- ast_bl %>% clean_names()
names(ast_bl)

# Combine to AST data
ast_bl <- left_join(ast_bl, df[, c("recordid", "matched_treat_combo")], by = "recordid")


# Specify column with antib
antib_cols <- c("colistin", "polymyxin_b",
                "ampicillin_sulbactam", "cefoperazone_sulbactam",
                "tigecycline",
                "doripenem", "imipenem", "meropenem",
                "ceftazidime_avibactam", "ceftolozane_tazobactam",
                "amikacin", "gentamicin")

#-------------------------------------
###
ast_bl_long <- ast_bl %>%
  separate_rows(matched_treat_combo, sep = "\\+") %>%
  mutate(matched_treat_combo = trimws(matched_treat_combo)) %>%
  mutate(matched_treat_combo_clean = to_snake_case(matched_treat_combo))

# check for anything that still won't match
setdiff(unique(ast_bl_long$matched_treat_combo_clean), antib_cols)

#-------------------------------------
### Matching
antib_mat <- as.matrix(ast_bl_long[antib_cols])
col_idx <- match(ast_bl_long$matched_treat_combo_clean, antib_cols)
ast_bl_long$antib_ast <- antib_mat[cbind(seq_len(nrow(ast_bl_long)), col_idx)]

#-------------------------------------#----------------------------------------#



###
ast_bl_long <- ast_bl_long %>%
  group_by(recordid) %>%
  mutate(
    therapy_type = if (n() > 1) "combination" else "monotherapy",
    antib_act = case_when(
      therapy_type == "monotherapy" & antib_ast == 1            ~ 0, #inactive
      therapy_type == "monotherapy" & antib_ast %in% c(2, 3)    ~ 1, # active
      therapy_type == "monotherapy" & antib_ast %in% c(4, 5)    ~ 2, # unknown
      therapy_type == "combination" & any(antib_ast == 1)       ~ 0,
      therapy_type == "combination" & any(antib_ast %in% c(2,3))~ 1,
      therapy_type == "combination" & all(antib_ast %in% c(4,5))~ 2,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

#
table(ast_bl_long$therapy_type, ast_bl_long$antib_act)

#-------------------------------------
###
table(ast_bl_long$therapy_type)

#
table(ast_bl_long$therapy_type, ast_bl_long$antib_ast)

#-------------------------------------
###
ast_bl_dis <- ast_bl_long %>%
  group_by(therapy_type, antib_ast) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(therapy_type) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

#-------------------------------------#----------------------------------------#



#-------------------------------------
#-------------------------------------
### First spec to calculate AST time 
ast_bl_first <- ast_bl %>%
  group_by(recordid) %>%
  filter(spec_date == min(spec_date, na.rm = TRUE)) %>%
  ungroup()

### Time to AST from onset
ast_bl_first$onset_to_ast <- as.numeric(ast_bl_first$ast_date - ast_bl_first$inf_onset)

#
table(is.na(ast_bl_first$onset_to_ast))

table(ast_bl_first$onset_to_ast)

#-------------------------------------
### Example threshold: negatives down to -1 day become 0; less than -1 become NA
#
thr <- 1  # days
uhr <- 7

# negatives less than -1 become NA
ast_bl_first$onset_to_ast[ast_bl_first$onset_to_ast < -thr] <- NA

ast_bl_first$onset_to_ast[ast_bl_first$onset_to_ast > uhr] <- NA

# negatives down to -1 day become 0
ast_bl_first$onset_to_ast[ast_bl_first$onset_to_ast < 0] <- 0

## Change 0 -> 1 (only work for descriptive analysis)
ast_bl_first$onset_to_ast[ast_bl_first$onset_to_ast == 0] <- 1

#-------------------------------------
#
table(is.na(ast_bl_first$onset_to_ast))

table(ast_bl_first$onset_to_ast)

#-------------------------------------
# Combine to df
df_bl <- left_join(df_bl, ast_bl_first[, c("recordid", "onset_to_ast")], by = "recordid")

#-------------------------------------
#-------------------------------------
### Summary
summary_by_stg_bl <- df_bl %>%
  group_by(treatment_strategy) %>%
  summarise(
    n = sum(!is.na(onset_to_ast)),
    median = median(onset_to_ast, na.rm = TRUE),
    Q1 = quantile(onset_to_ast, 0.25, na.rm = TRUE),
    Q3 = quantile(onset_to_ast, 0.75, na.rm = TRUE),
    IQR = IQR(onset_to_ast, na.rm = TRUE)
  )

#-------------------------------------#----------------------------------------#



### Median AST from onset - PP population
#-------------------------------------
#-------------------------------------
### First spec to calculate AST time 
ast_lm_first <- ast_lm %>%
  group_by(recordid) %>%
  filter(spec_date == min(spec_date, na.rm = TRUE)) %>%
  ungroup()

### Time to AST from onset
ast_lm_first$onset_to_ast <- as.numeric(ast_lm_first$ast_date - ast_lm_first$inf_onset)

#
table(is.na(ast_lm_first$onset_to_ast))

table(ast_lm_first$onset_to_ast)

#-------------------------------------
### Example threshold: negatives down to -1 day become 0; less than -1 become NA
#
thr <- 1  # days
uhr <- 7

# negatives less than -1 become NA
ast_lm_first$onset_to_ast[ast_lm_first$onset_to_ast < -thr] <- NA

ast_lm_first$onset_to_ast[ast_lm_first$onset_to_ast > uhr] <- NA

# negatives down to -1 day become 0
ast_lm_first$onset_to_ast[ast_lm_first$onset_to_ast < 0] <- 0

## Change 0 -> 1 (only work for descriptive analysis)
ast_lm_first$onset_to_ast[ast_lm_first$onset_to_ast == 0] <- 1

#-------------------------------------
#
table(is.na(ast_lm_first$onset_to_ast))

table(ast_lm_first$onset_to_ast)

#-------------------------------------
# Combine to df
df_lm <- left_join(df_lm, ast_lm_first[, c("recordid", "onset_to_ast")], by = "recordid")

#-------------------------------------
#-------------------------------------
### Summary
summary_by_stg_lm <- df_lm %>%
  group_by(treatment_strategy) %>%
  summarise(
    n = sum(!is.na(onset_to_ast)),
    median = median(onset_to_ast, na.rm = TRUE),
    Q1 = quantile(onset_to_ast, 0.25, na.rm = TRUE),
    Q3 = quantile(onset_to_ast, 0.75, na.rm = TRUE),
    IQR = IQR(onset_to_ast, na.rm = TRUE)
  )

#-------------------------------------#----------------------------------------#

### Enterobacterales species
#
ast_bl_ent <- left_join(ast_bl, df_bl[, c("recordid", "org_r_combined", "treatment_strategy")], by = "recordid")

#
ast_bl_ent <- ast_bl_ent %>%
  filter(org_r_combined == "CRE")

#
table(ast_bl_ent$org_combined)

# NS or S to Ertapenem
table(ast_bl_ent$treatment_strategy, ast_bl_ent$ertapenem)





### CRE, NS to Erta vs NS to Mero/Imi
#-------------------------------------
#-------------------------------------
#


#-------------------------------------#----------------------------------------#



###-----------------------------------------------------------------------------
### ^.^ ###
