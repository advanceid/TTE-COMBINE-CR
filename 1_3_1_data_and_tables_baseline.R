
### Treatment strategy: Initiate Combination Therapy or Monotherapy on or within 1 day after the AST report
### 20260922
### Primary Outcome: 14-day ACM
### Secondary Outcome: 7-day & 21-day ACM


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
                 smd)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean data/clean/clean/combine_group/baseline_all_clean_combine.RData")
id_bl <- readRDS("data/clean data/clean/clean/final/tables_bl.RData")

### Load fonts
loadfonts()

#-------------------------------------
### Filter IDs
df <- df[which(df$recordid %in% id_bl$id),]

#-------------------------------------
### Discharged or died from D0 to D21?
###
df <- df %>%
  mutate(
    time_zero = as.Date(time_zero),
    mortality_date = as.Date(mortality_date),
    first_discharge_date = as.Date(first_discharge_date),
    d21_date_new = as.Date(d21_date_new),
    
    death_d21 = !is.na(mortality_date) &
      mortality_date >= time_zero &
      mortality_date <= d21_date_new,
    
    discharge_d21 = !is.na(first_discharge_date) &
      first_discharge_date >= time_zero &
      first_discharge_date <= d21_date_new &
      !(!is.na(mortality_date) & first_discharge_date == mortality_date),
    
    status_d21 = case_when(
      death_d21 ~ "death",
      discharge_d21 ~ "discharge",
      TRUE ~ "none"
    )
  )

### relocate variables
df <- df %>% relocate(status_d21, .before = d07_date_new)

df <- df %>% relocate(first_discharge_date, .before = time_zero)

df <- df %>% relocate(mortality_date, .before = time_zero)

#-------------------------------------



###-----------------------------------------------------------------------------
###
# Set variables
df_factor <- c("sex",
               "country_region", "country_income",
               "hpd_admtype", "adm_reason_01", 
               "cmb_rhe",
               "cmb_comorbidities___none",
               "icu_at_onset", 
               "mv_at_onset",
               "abnormal_temp",
               "vaso_use",
               "infection_types",
               "org_r_combined",
               "status_d21",
               "d21_death", "d07_death", "d14_death",
               "treatment_strategy")

df_num <- c("age_new", "comorbidities_CCI",
            "length_before_onset", 
            "sofa_score", "sofa_map", "cre_mg",
            "severity_score",
            "d21_random_days")


df_factor_index <- which(names(df) %in% df_factor)
df_num_index <- which(names(df) %in% df_num)

#
df <- as.data.frame(df)
for (i in df_factor_index){df[,i] = as.factor(df[,i])}
for (i in df_num_index){df[,i] = as.numeric(df[,i])}

#-------------------------------------
#-------------------------------------
###
df <- within(df, {c(
###
  treatment_strategy <- factor(treatment_strategy, 
                            levels = c("Combination therapy", "Monotherapy")),

###
  org_r_combined <- factor(org_r_combined, 
                           levels = c("CRA", "CRE", "CRP")),
  
    
###  
  sex <- factor(sex, levels = c("F", "M"), labels = c("Female", "Male")),

###  
  country_region <- factor(country_region, 
                           levels = c("Western Pacific Region", "South-East Asian Region",
                                      "Eastern Mediterranean Region")),

###  
  country_income <- factor(country_income, 
                           levels = c("High income", "Upper middle income",
                                      "Lower middle income")),
  
###  
  cmb_comorbidities___none <- factor(cmb_comorbidities___none, 
                                   levels = c(0,1), labels = c("Yes", "No")), # number of comorbidities
  
  cmb_comorbidities___aids <- factor(cmb_comorbidities___aids, 
                                     levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___onc <- factor(cmb_comorbidities___onc, 
                                    levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___cpd <- factor(cmb_comorbidities___cpd, 
                                    levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___cog <- factor(cmb_comorbidities___cog, 
                                    levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___rheu <- factor(cmb_comorbidities___rheu, 
                                     levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___dem <- factor(cmb_comorbidities___dem, 
                                    levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___diab <- factor(cmb_comorbidities___diab, 
                                     levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___diad <- factor(cmb_comorbidities___diad, 
                                     levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___hop <- factor(cmb_comorbidities___hop, 
                                    levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___hivwa <- factor(cmb_comorbidities___hivwa, 
                                      levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___hivna <- factor(cmb_comorbidities___hivna, 
                                      levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___mlr <- factor(cmb_comorbidities___mlr, 
                                    levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___mal <- factor(cmb_comorbidities___mal, 
                                    levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___mst <- factor(cmb_comorbidities___mst, 
                                    levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___mld <- factor(cmb_comorbidities___mld, 
                                    levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___liv <- factor(cmb_comorbidities___liv, 
                                    levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___pep <- factor(cmb_comorbidities___pep, 
                                    levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___renal <- factor(cmb_comorbidities___renal, 
                                      levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___tub <- factor(cmb_comorbidities___tub, 
                                    levels = c(1,0), labels = c("Yes", "No")),
  
  cmb_comorbidities___oth <- factor(cmb_comorbidities___oth, 
                                    levels = c(1,0), labels = c("Yes", "No")),
  
###  
  cmb_rhe <- factor(cmb_rhe, 
                    levels = c(1,0), labels = c("Yes", "No")),

###
  hpd_admtype <- factor(hpd_admtype, levels = c("EMR", "ELT"), labels = c("Emergency", "Elective")),
  

###  
  adm_reason_01 <- factor(adm_reason_01,
                          levels = c("Internal diseases", "Infectious disease", "Surgical", 
                                     "Hae_Onc", "Others")),
  
###
  abnormal_temp <- factor(abnormal_temp, 
                          levels = c(1,0), labels = c("Yes", "No")),

###
  vaso_use <- factor(vaso_use, 
                        levels = c(1,0), labels = c("Yes", "No")),


###  
  icu_at_onset <- factor(icu_at_onset,
                          levels = c(1,0), labels = c("Yes", "No")),

  mv_at_onset <- factor(mv_at_onset,
                        levels = c(1,0), labels = c("Yes", "No")),

###
  infection_types <- factor(infection_types,
                            levels = c("BSI", 
                                       "VAP")),

###
  status_d21 <- factor(status_d21, 
               levels = c("death", "discharge", "none")),


###
  d21_death <- factor(d21_death, 
                      levels = c(1,0), labels = c("Yes", "No")),

  d07_death <- factor(d07_death, 
                    levels = c(1,0), labels = c("Yes", "No")),

  d14_death <- factor(d14_death, 
                    levels = c(1,0), labels = c("Yes", "No"))
  
)})

#-------------------------------------
#-------------------------------------
# Revise labels
#
set_labels <- function(df, labels) {
  for (col in names(labels)) {
    df <- labelled::set_variable_labels(df, !!sym(col) := labels[[col]])
  }
  return(df)
}

#
labels <- list(
  age_new = "Age (years)",
  sex = "Sex",
  country_region = "WHO regions",
  country_income = "World Bank income classification",
  
  cmb_comorbidities___none = "Comorbidities",
  cmb_comorbidities___aids = "AIDS",
  cmb_comorbidities___onc = "Cancer/leukaemia",
  cmb_comorbidities___cpd = "Chronic pulmonary disease",
  cmb_comorbidities___cog = "Congestive heart failure",
  cmb_comorbidities___rheu = "Connective tissue/rheumatologic disease",
  cmb_comorbidities___dem = "Dementia",
  cmb_comorbidities___diab = "Diabetes",
  cmb_comorbidities___diad = "Diabetes with end organ damage",
  cmb_comorbidities___hop = "Hemi- or paraplegia",
  cmb_comorbidities___hivwa = "HIV on ART",
  cmb_comorbidities___hivna = "HIV without ART",
  cmb_comorbidities___mlr = "Malaria",
  cmb_comorbidities___mal = "Malnutrition",
  cmb_comorbidities___mst = "Metastatic solid tumour",
  cmb_comorbidities___mld = "Mild liver disease",
  cmb_comorbidities___liv = "Moderate or severe liver disease",
  cmb_comorbidities___pep = "Peptic ulcer",
  cmb_comorbidities___renal = "Renal disease",
  cmb_comorbidities___tub = "Tuberculosis",
  cmb_comorbidities___oth = "Others",
  
  comorbidities_CCI = "Charlson comorbidity index",
  cmb_rhe = "Recent healthcare exposure within 3 months",
  
  hpd_admtype = "Admission type",
  adm_reason_01 = "Primary admission reason",
  
  length_before_onset = "Hospital LOS before onset (days)",
  sofa_score = "SOFA score",
  severity_score = "qSOFA score",
  
  vaso_use = "Vasopressure use",
  abnormal_temp = "Abnormal Temperature",
  
  sofa_map = "MAP",
  cre_mg = "Creatinine (mg/dL)",
  
  icu_at_onset = "ICU admission at onset",
  mv_at_onset = "MV status at onset",

  
  infection_types = "Infection types",
  
  status_d21 = "Follow-up 21",

  d21_death = "21-day ACM",
  d07_death = "07-day ACM",
  d14_death = "14-day ACM",

  
  d21_random_days = "Length of Follow-up (days)",
  
  org_r_combined = "CR-GNB pathogens",

  treatment_strategy = "Treatment"
)

#
df <- set_labels(df, labels)

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------


### Table 1 - CR cohort
#-------------------------------------
#-------------------------------------
#
table_combine_baseline <- df %>%
  dplyr::select(age_new, #SMD on mean age
         sex,
         country_region, country_income, 
         cmb_comorbidities___none,
         cmb_comorbidities___aids, cmb_comorbidities___onc,
         cmb_comorbidities___cpd, cmb_comorbidities___cog,
         cmb_comorbidities___rheu, cmb_comorbidities___dem,
         cmb_comorbidities___diab, cmb_comorbidities___diad,
         cmb_comorbidities___hop, cmb_comorbidities___hivwa,
         cmb_comorbidities___hivna, cmb_comorbidities___mlr,
         cmb_comorbidities___mal, cmb_comorbidities___mst,
         cmb_comorbidities___mld, cmb_comorbidities___liv,
         cmb_comorbidities___pep, cmb_comorbidities___renal,
         cmb_comorbidities___tub, cmb_comorbidities___oth,
         comorbidities_CCI, 
         cmb_rhe,
         hpd_admtype, adm_reason_01, 
         length_before_onset, 
         icu_at_onset, mv_at_onset,
         abnormal_temp, vaso_use, 
         sofa_map, cre_mg, 
         sofa_score, severity_score,
         infection_types,
         org_r_combined,
         treatment_strategy) %>%
  tbl_summary(
    by = treatment_strategy,
    missing = "ifany", # showed only if having missing values
    missing_text = "Missing",
    type = list(
      all_continuous() ~ "continuous2",
      comorbidities_CCI ~ "continuous2",
      cmb_comorbidities___none ~ "dichotomous", 
      cmb_comorbidities___aids ~ "dichotomous", 
      cmb_comorbidities___onc ~ "dichotomous",
      cmb_comorbidities___cpd ~ "dichotomous", 
      cmb_comorbidities___cog ~ "dichotomous",
      cmb_comorbidities___rheu ~ "dichotomous", 
      cmb_comorbidities___dem ~ "dichotomous",
      cmb_comorbidities___diab ~ "dichotomous", 
      cmb_comorbidities___diad ~ "dichotomous",
      cmb_comorbidities___hop ~ "dichotomous", 
      cmb_comorbidities___hivwa ~ "dichotomous",
      cmb_comorbidities___hivna ~ "dichotomous", 
      cmb_comorbidities___mlr ~ "dichotomous",
      cmb_comorbidities___mal ~ "dichotomous", 
      cmb_comorbidities___mst ~ "dichotomous",
      cmb_comorbidities___mld ~ "dichotomous", 
      cmb_comorbidities___liv ~ "dichotomous",
      cmb_comorbidities___pep ~ "dichotomous", 
      cmb_comorbidities___renal ~ "dichotomous",
      cmb_comorbidities___tub ~ "dichotomous", 
      cmb_comorbidities___oth ~ "dichotomous",
      
      length_before_onset ~ "continuous2",

      sofa_score ~ "continuous2",
      severity_score ~ "continuous2"
    ),
    statistic = list(
      age_new ~ "{median} ({p25}-{p75})",
      comorbidities_CCI ~ "{median} ({p25}-{p75})",
      cmb_comorbidities___none ~ "{n} ({p}%)",
      cmb_comorbidities___aids ~ "{n} ({p}%)", 
      cmb_comorbidities___onc ~ "{n} ({p}%)",
      cmb_comorbidities___cpd ~ "{n} ({p}%)",
      cmb_comorbidities___cog ~ "{n} ({p}%)",
      cmb_comorbidities___rheu ~ "{n} ({p}%)", 
      cmb_comorbidities___dem ~ "{n} ({p}%)",
      cmb_comorbidities___diab ~ "{n} ({p}%)", 
      cmb_comorbidities___diad ~ "{n} ({p}%)",
      cmb_comorbidities___hop ~ "{n} ({p}%)",
      cmb_comorbidities___hivwa ~ "{n} ({p}%)",
      cmb_comorbidities___hivna ~ "{n} ({p}%)",
      cmb_comorbidities___mlr ~ "{n} ({p}%)",
      cmb_comorbidities___mal ~ "{n} ({p}%)", 
      cmb_comorbidities___mst ~ "{n} ({p}%)",
      cmb_comorbidities___mld ~ "{n} ({p}%)", 
      cmb_comorbidities___liv ~ "{n} ({p}%)",
      cmb_comorbidities___pep ~ "{n} ({p}%)", 
      cmb_comorbidities___renal ~ "{n} ({p}%)",
      cmb_comorbidities___tub ~ "{n} ({p}%)",
      cmb_comorbidities___oth ~ "{n} ({p}%)",
      length_before_onset  ~ "{median} ({p25}-{p75})",
      
      sofa_map ~ "{mean} ({sd})",
      cre_mg ~ "{median} ({p25}-{p75})",

      sofa_score ~ "{median} ({p25}-{p75})",
      severity_score ~ "{mean} ({sd})"
    ),
    digits = list(all_continuous() ~ 0, # round all continuous variables to 0 decimals
                  cre_mg ~ 2,
                  all_categorical() ~ c(0, 0)), # counts 0 decimals, percentage 0 decimals
    percent = "column",
    sort = adm_reason_01 ~ "frequency",
    label = adm_reason_01 ~ "Primary admission reason"
  ) %>%
  add_overall() %>%
  add_difference(test = everything() ~ "smd") %>%  #SMD in Table 1
  modify_column_hide(columns = ci) %>%
  modify_header(
    update = list(
      label ~ "**Characteristics**",
      all_stat_cols() ~ gt::html("**{level}**<br>N = {format(n, big.mark = ',')}")
    )
  ) %>%
  bold_labels() %>%
  remove_bold(columns = label, rows = label == "AIDS") %>%
  remove_bold(columns = label, rows = label == "Cancer/leukaemia") %>%
  remove_bold(columns = label, rows = label == "Chronic pulmonary disease") %>%
  remove_bold(columns = label, rows = label == "Congestive heart failure") %>%
  remove_bold(columns = label, rows = label == "Connective tissue/rheumatologic disease") %>%
  remove_bold(columns = label, rows = label == "Dementia") %>%
  remove_bold(columns = label, rows = label == "Diabetes") %>%
  remove_bold(columns = label, rows = label == "Diabetes with end organ damage") %>%
  remove_bold(columns = label, rows = label == "Hemi- or paraplegia") %>%
  remove_bold(columns = label, rows = label == "HIV on ART") %>%
  remove_bold(columns = label, rows = label == "HIV without ART") %>%
  remove_bold(columns = label, rows = label == "Malaria") %>%
  remove_bold(columns = label, rows = label == "Malnutrition") %>%
  remove_bold(columns = label, rows = label == "Metastatic solid tumour") %>%
  remove_bold(columns = label, rows = label == "Mild liver disease") %>%
  remove_bold(columns = label, rows = label == "Moderate or severe liver disease") %>%
  remove_bold(columns = label, rows = label == "Peptic ulcer") %>%
  remove_bold(columns = label, rows = label == "Renal disease") %>%
  remove_bold(columns = label, rows = label == "Tuberculosis") %>%
  remove_bold(columns = label, rows = label == "Others") %>%
  remove_bold(columns = label, rows = label == "Cardiovascular failure") %>%
  remove_bold(columns = label, rows = label == "Trauma") %>%
  remove_bold(columns = label, rows = label == "Sepsis") %>%
  remove_bold(columns = label, rows = label == "Respiratory failure") %>%
  remove_bold(columns = label, rows = label == "Neurological failure") %>%
  remove_bold(columns = label, rows = label == "Post-operative care") %>%
  remove_bold(columns = label, rows = label == "Metabolic acidosis") %>%
  remove_bold(columns = label, rows = label == "COVID-19") %>%
  modify_footnote(update = everything() ~ NA) %>%
  as_gt() %>%
  tab_style(
    style = cell_text(align = "left", v_align = "middle"),
    locations = cells_column_labels(columns = 1:5)
  ) 

#-------------------------------------
#-------------------------------------
# Customize the headers to split lines 
gt_table_combine_baseline <- table_combine_baseline %>%
  tab_options(
    column_labels.border.top.color = "black",
    column_labels.border.top.width = px(2.5),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(2.5),
    table_body.hlines.color = "white",
    table_body.hlines.width = px(0),
    table.border.bottom.color = "white",
    table.border.bottom.width = px(0),
    data_row.padding = px(0),
    footnotes.padding = px(0)
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "black",
      weight = px(2.5)
    ),
    locations = cells_body(
      rows = nrow(table_combine_baseline$`_data`)
    )
  ) %>%
  cols_width(
    1 ~ px(600),
    2:6 ~ px(300)
  ) %>%
  tab_style(
    style = cell_text(
      font = c("Poppins"), #"Times New Roman"
      size = px(14)
    ),
    locations = list(
      cells_title(groups = "title"),
      cells_title(groups = "subtitle"),
      cells_column_labels(),
      cells_body()
    )
  ) %>%
  tab_style(
    style = cell_text(
      font = c("Poppins"),
      size = px(12)
    ),
    locations = cells_footnotes()
  ) %>%
  tab_style(
    style = cell_text(
      font = c("Poppins"),
      size = px(12)
    ),
    locations = cells_source_notes()
  ) 

# Print the table
gt_table_combine_baseline

# Save
gtsave(table_combine_baseline, filename = "output/tables/table_baseline_bl.html")

###------------------------------

###------------------------------
table_combine_outcome <- df %>% 
  dplyr::select(d21_death, d07_death, d14_death,
         d21_random_days,
         status_d21,
         treatment_strategy) %>%
  tbl_summary(
    by = treatment_strategy,
    missing = "ifany", # showed only if having missing values
    missing_text = "Missing",
    type = list(
      all_continuous() ~ "continuous2",
      d21_death ~ "categorical",
      d07_death ~ "categorical",
      d14_death ~ "categorical",
      d21_random_days ~ "continuous2",
      status_d21 ~ "categorical"
    ),
    statistic = list(
      d21_death ~ "{n} ({p}%)",
      d07_death ~ "{n} ({p}%)",
      d14_death ~ "{n} ({p}%)",
      d21_random_days ~ "{median} ({p25}-{p75})",
      status_d21 ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 0, # round all continuous variables to 0 decimals
                  all_categorical() ~ c(0, 0)), # counts 0 decimals, percentage 0 decimals
    percent = "column"
  ) %>%
  add_overall() %>%
  modify_header(
    update = list(
      label ~ "**Characteristics**",
      all_stat_cols() ~ gt::html("**{level}**<br>N = {format(n, big.mark = ',')}")
    )
  ) %>%
  bold_labels() %>%
  modify_footnote(update = everything() ~ NA) %>%
  as_gt() %>%
  tab_style(
    style = cell_text(align = "left", v_align = "middle"),
    locations = cells_column_labels(columns = 1:5)
  ) 

# Customize the headers to split lines 
gt_table_combine_outcome <- table_combine_outcome %>%
  tab_options(
    column_labels.border.top.color = "black",
    column_labels.border.top.width = px(2.5),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(2.5),
    table_body.hlines.color = "white",
    table_body.hlines.width = px(0),
    table.border.bottom.color = "white",
    table.border.bottom.width = px(0),
    data_row.padding = px(0),
    footnotes.padding = px(0)
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "black",
      weight = px(2.5)
    ),
    locations = cells_body(
      rows = nrow(table_combine_outcome$`_data`)
    )
  ) %>%
  cols_width(
    1 ~ px(600),
    2:6 ~ px(300)
  ) %>%
  tab_style(
    style = cell_text(
      font = c("Poppins"), #"Times New Roman"
      size = px(14)
    ),
    locations = list(
      cells_title(groups = "title"),
      cells_title(groups = "subtitle"),
      cells_column_labels(),
      cells_body()
    )
  ) %>%
  tab_style(
    style = cell_text(
      font = c("Poppins"),
      size = px(12)
    ),
    locations = cells_footnotes()
  ) %>%
  tab_style(
    style = cell_text(
      font = c("Poppins"),
      size = px(12)
    ),
    locations = cells_source_notes()
  ) 

# Print the table
gt_table_combine_outcome 

# Save
gtsave(table_combine_outcome, filename = "output/tables/table_outcome_bl.html")

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------





###-----------------------------------------------------------------------------
### ^.^ ###