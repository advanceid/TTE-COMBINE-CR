
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
                 mice)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean data/baseline_outcomes_index.RData")
ast <- readRDS("data/clean data/ast_all_index.RData")
antib <- readRDS("data/clean data/anti_treat_index_poly.RData")

###
# 
df[df == ""] <- NA

#-------------------------------------
#-------------------------------------
# adult VAP (>= 18)
df_01 <- df[which(df$age_group == 1),]

# 
table(df_01$infection_types)

# filter ast with adult patients
ast_01 <- ast[which(ast$recordid %in% df_01$recordid), ]

###-------------------------------
### Combine with AST
#
df_used <- left_join(df_01, ast_01[, c("recordid", "org_names_all", "org_combined",
                                                         "pathogen_combined_types_ori", "pathogen_combined_types")], by = "recordid")

# distinct by IDs and orgs
df_used <- df_used %>%
  distinct(recordid, org_names_all, .keep_all = TRUE)

###-------------------------------
### NA ~ Other
#
df_used$pathogen_combined_types_ori[is.na(df_used$pathogen_combined_types_ori)] <- "Other bacteria"

#
df_used$pathogen_combined_types <- as.character(df_used$pathogen_combined_types)

#
df_used$pathogen_combined_types[is.na(df_used$pathogen_combined_types)] <- "Other"

###-------------------------------
### GNB or GPB with Fungi ~ GNB or GPB, not Polymicrobial !!!
#
df_used$pathogen_combined_types[df_used$pathogen_combined_types_ori == "GNB with fungi"] <- "Gram-negative bacteria"

#
df_used$pathogen_combined_types[df_used$pathogen_combined_types_ori == "GPB with fungi"] <- "Gram-positive bacteria"

###-------------------------------
### delete Fungi
#
df_used <- df_used %>%
  filter(pathogen_combined_types != "Fungi")

# distinct by pathogen types
df_baseline <- df_used %>%
  distinct(recordid, pathogen_combined_types, .keep_all = TRUE)

#-------------------------------------
#-------------------------------------
### Save clean data for Descriptive analysis
#
saveRDS(df_baseline, file = "data/clean data/clean/baseline_adult.RData") 

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------



###-------------------------------
### Specific pathogens
# Delete patients with NA pathogens
df_used <- df_used %>%
  filter(pathogen_combined_types != "Other")

# Delete GPB Mono
df_used <- df_used %>%
  filter(pathogen_combined_types != "Gram-positive bacteria")

###-------------------------------
### ENT
### Salmonella as exclusion criterion
ent_list <- c("E. coli", "K. pneumoniae", "Klebsiella", 
              "Enterobacter", "Serratia", "Proteus", "Morganella",
              "Citrobacter", "Providencia", "Raoultella", "Kluyvera",
              "Cronobacter", "Leclercia", "Pluralibacter", 
              "Edwardsiella", "Pantoea", "Escherichia",
              "Kalamiella") #Kalamiella ~ Pantoea

###-------------------------------
### 
df_used$org_names_all[is.na(df_used$org_names_all) & df_used$org_combined == "Acinetobacter spp."] <- "Acinetobacter spp."

###-------------------------------
#
org_r <- matrix(NA, ncol = 1, nrow = nrow(df_used))

#
for (i in 1:nrow(df_used)) {
  if (df_used$org_names_all[i] == "Acinetobacter spp." & df_used$aci_car[i] == 1) {
    org_r[i,1] <- "CRA"
  } else if (df_used$org_names_all[i] %in% ent_list & df_used$ent_car[i] == 1) {
    org_r[i,1] <- "CRE"
  } else if (df_used$org_names_all[i] == "Pseudomonas spp." & df_used$pse_car[i] == 1) {
    org_r[i,1] <- "CRP"
  } else {
    org_r[i,1] <- "Other"
  }
}

#
df_used <- df_used %>%
  mutate(org_r)

#
df_used$org_r <- as.character(df_used$org_r)

###-------------------------------
### Filter with CRA, CRE, CRP
df_used <- df_used %>%
  filter(org_r != "Other")

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------





### Combine with baseline
#-------------------------------------
#-------------------------------------
### combine org_r
df_used <- df_used %>%
  group_by(recordid) %>%
  mutate(
    org_r_combined = paste(unique(org_r), collapse = ", ")
  ) %>%
  ungroup()

# Define a function to sort alphabetically and standardize combinations
standardize_combination <- function(combination) {
  elements <- unlist(strsplit(combination, ", "))
  standardized <- paste(sort(elements), collapse = ", ")
  return(standardized)
}

# Apply the function to the pathogen_group_combined column
df_used$org_r_combined <- sapply(df_used$org_r_combined, standardize_combination)

#-------------------------------------
# distinct by org_r_combined
df_cr <- df_used %>%
  distinct(recordid, org_r_combined, .keep_all = TRUE)

#-------------------------------------#----------------------------------------#



#-------------------------------------
#-------------------------------------
### Save clean data for Descriptive analysis
#
saveRDS(df_cr, file = "data/clean data/clean/baseline_cr.RData") 

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
### ^.^ ###