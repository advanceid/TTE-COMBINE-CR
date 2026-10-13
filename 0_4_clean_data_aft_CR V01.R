
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
df <- readRDS("data/clean data/general/baseline_all.RData")
hos_reason <- read.xlsx("data/hos_reasons.xlsx")


###
# 
df[df == ""] <- NA

#-------------------------------------
#-------------------------------------
### Age 
# Distribution
table(df$age_new, useNA = "ifany")

prop.table(table(df$age_new))

par(mfrow=c(1,1), mar=c(4, 4, 2, 1))
plot(1:10, main="Test Plot")
# barplot() expects frequencies, not the raw data, so use table inside barplot()
barplot(table(df$age_new),
        ylab = "Frequency", xlab = "Age (years)")

#-------------------------------------
# No NA in age_new



#-------------------------------------
#-------------------------------------
## length_before_onset
#
table(df$length_before_onset, useNA = "ifany")

prop.table(table(df$length_before_onset))

par(mfrow=c(1,1), mar=c(4, 4, 2, 1))
plot(1:10, main="Test Plot")
# barplot() expects frequencies, not the raw data, so use table inside barplot()
barplot(table(df$length_before_onset),
        ylab = "Frequency", xlab = "length_before_onset (days)")

#-------------------------------------
## exclude neg days
#df <- df[which(df$length_before_onset >= 0), ]

### Example threshold: negatives down to -1 day become 0; less than -1 become NA
#
thr <- 1  # days

# negatives less than -1 become NA
df$length_before_onset[df$length_before_onset < -thr] <- NA

# negatives down to -1 day become 0
df$length_before_onset[df$length_before_onset < 0] <- 0

## Change 0 -> 1 (only work for descriptive analysis)
df$length_before_onset[df$length_before_onset == 0] <- 1



#-------------------------------------
#-------------------------------------
## icu_to_onset
#
table(df$icu_to_onset, useNA = "ifany")

prop.table(table(df$icu_to_onset))

par(mfrow=c(1,1), mar=c(4, 4, 2, 1))
plot(1:10, main="Test Plot")
# barplot() expects frequencies, not the raw data, so use table inside barplot()
barplot(table(df$icu_to_onset),
        ylab = "Frequency", xlab = "icu_to_onset (days)")

#-------------------------------------
## exclude neg days
#df <- df[which(df$icu_to_onset >= 0), ]

# negatives less than -1 become NA
df$icu_to_onset[df$icu_to_onset < -thr] <- NA

# negatives down to -1 day become 0
df$icu_to_onset[df$icu_to_onset < 0] <- 0

## Change 0 -> 1 (only work for descriptive analysis)
df$icu_to_onset[df$icu_to_onset == 0] <- 1



#-------------------------------------
#-------------------------------------
## mv_to_onset
#
table(df$mv_to_onset, useNA = "ifany")

prop.table(table(df$mv_to_onset))

par(mfrow=c(1,1), mar=c(4, 4, 2, 1))
plot(1:10, main="Test Plot")
# barplot() expects frequencies, not the raw data, so use table inside barplot()
barplot(table(df$mv_to_onset),
        ylab = "Frequency", xlab = "mv_to_onset (days)")

#-------------------------------------
## exclude neg days
#df <- df[which(df$mv_to_onset >= 0), ]

# negatives less than -1 become NA
df$mv_to_onset[df$mv_to_onset < -thr] <- NA

# negatives down to -1 day become 0
df$mv_to_onset[df$mv_to_onset < 0] <- 0


## Change 0 -> 1 (only work for descriptive analysis)
df$mv_to_onset[df$mv_to_onset == 0] <- 1



#-------------------------------------
#-------------------------------------
### Hospital LOS
#
table(df$first28_hosp, useNA = "ifany")

prop.table(table(df$first28_hosp))

par(mfrow=c(1,1), mar=c(4, 4, 2, 1))
plot(1:10, main="Test Plot")
# barplot() expects frequencies, not the raw data, so use table inside barplot()
barplot(table(df$first28_hosp),
        ylab = "Frequency", xlab = "Hospital LOS (days)")

#-------------------------------------
# exclude neg days
#df <- df[which(df$first28_hosp >= 0), ]

# negatives less than -1 become NA
df$first28_hosp[df$first28_hosp < -thr] <- NA

# negatives down to -1 day become 0
df$first28_hosp[df$first28_hosp < 0] <- 0


## Change 0 -> 1 (only work for descriptive analysis)
df$first28_hosp[df$first28_hosp == 0] <- 1



#-------------------------------------
#-------------------------------------
## ICU LOS
#
table(df$first28_icu, useNA = "ifany")

prop.table(table(df$first28_icu))

par(mfrow=c(1,1), mar=c(4, 4, 2, 1))
plot(1:10, main="Test Plot")
# barplot() expects frequencies, not the raw data, so use table inside barplot()
barplot(table(df$first28_icu),
        ylab = "Frequency", xlab = "ICU LOS (days)")

#-------------------------------------
# exclude neg days
#df <- df[which(df$first28_icu >= 0), ]

# negatives less than -1 become NA
df$first28_icu[df$first28_icu < -thr] <- NA

# negatives down to -1 day become 0
df$first28_icu[df$first28_icu < 0] <- 0


## Change 0 -> 1 (only work for descriptive analysis)
df$first28_icu[df$first28_icu == 0] <- 1



#-------------------------------------
#-------------------------------------
## MV duration
#
table(df$first28_mv, useNA = "ifany")

prop.table(table(df$first28_mv))

par(mfrow=c(1,1), mar=c(4, 4, 2, 1))
plot(1:10, main="Test Plot")
# barplot() expects frequencies, not the raw data, so use table inside barplot()
barplot(table(df$first28_mv),
        ylab = "Frequency", xlab = "MV duration (days)")

#-------------------------------------
# exclude neg days
#df <- df[which(df$first28_mv >= 0), ]

# negatives less than -1 become NA
df$first28_mv[df$first28_mv < -thr] <- NA

# negatives down to -1 day become 0
df$first28_mv[df$first28_mv < 0] <- 0


## Change 0 -> 1 (only work for descriptive analysis)
df$first28_mv[df$first28_mv == 0] <- 1

#-------------------------------------#-----------------------------------------
#-------------------------------------#-----------------------------------------



#-------------------------------------
#-------------------------------------
### Name primary admission reason
###
# Split Primary admission reason
# 
split_condition <- function(condition) {
  # Find the first " "
  pos <- regexpr(" ", condition)
  # 
  code <- substr(condition, 1, pos - 1)
  description <- substr(condition, pos + 1, nchar(condition))
  return(c(code, description))
}

#
split_hos <- t(sapply(hos_reason$Primary.admission.reason, split_condition))

# 
split_hos_reason <- data.frame(abb = split_hos[, 1], 
                               names = split_hos[, 2], 
                               stringsAsFactors = FALSE)

## Set new levels for hpd_admreason
new_levels <- c(setdiff(levels(df$hpd_admreason), "OTH"), "OTH")
df$hpd_admreason <- as.factor(df$hpd_admreason)
level_admreason <- levels(df$hpd_admreason)

# Match new names
match_indices <- match(level_admreason, split_hos_reason$abb)
#
new_level_admreason <- split_hos_reason$names[match_indices]
#
df$hpd_admreason <- factor(df$hpd_admreason, levels = levels(df$hpd_admreason), 
                           labels = seq_along(levels(df$hpd_admreason)))

#
df$hpd_admreason  <- factor(df$hpd_admreason, 
                            levels = c(1:17),
                            labels = new_level_admreason)
###-----------------------------------------------------------------------------



#-------------------------------------
#-------------------------------------
### Save clean data for Descriptive & Regression analysis
#
saveRDS(df, file = "data/clean data/clean/baseline_all.RData") 


###-----------------------------------------------------------------------------
### ^.^ ###
