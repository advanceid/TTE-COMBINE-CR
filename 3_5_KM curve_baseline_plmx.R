
### Treatment strategy: Initiate Combination Therapy or Monotherapy on or within 1 day after the AST report
###20260904
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
                 writexl,
                 survival,
                 survminer,
                 ggplot2,
                 boot,
                 Hmisc,
                 Cairo,
                 patchwork,
                 cowplot)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean data/clean/clean/combine_group/sens/event_KM.RData")
id_bl <- readRDS("data/clean data/clean/clean/final/sens/tables_bl.RData")

### Load fonts
loadfonts()

#-------------------------------------
### Filter IDs
df <- df[which(df$recordid %in% id_bl$id),]

#-------------------------------------
###
df_used <- df %>%
  dplyr::select("recordid", "treatment_strategy",
         "d21_random_days", "d21_death",
         "org_r_combined", "infection_types")

### change column names
df_used <- df_used %>%
  rename(time = d21_random_days,
         event = d21_death,
         group = treatment_strategy
         )

### Ref ~ Combination therapy
df_used$group <- factor(df_used$group)

df_used$group <- relevel(df_used$group, ref = "Combination therapy")

#-------------------------------------#----------------------------------------#



### KM
### Nonparametric estimator 
### Cumulative incidence (risk) curves for D28 ACM
### Estimate the ITT effect of randomization antib groups
#-------------------------------------
#-------------------------------------
### Create KM curves by randomization group ###
### Cumulative Mortality
fit <- survfit(Surv(time, event) ~ group, data = df_used)

plot_02 <- ggsurvplot(
  fit,
  fun = "event",
  data = df_used,
  conf.int = TRUE,
  conf.int.alpha = 0.15,                 # softer, less visually heavy bands
  risk.table = TRUE,
  risk.table.height = 0.25,
  pval = TRUE,
  pval.size = 4.5,
  pval.coord = c(0, 0.12),               # reposition away from curve crossing
  censor.shape = "|",
  censor.size = 3,
  xlab = "Days from treatment assignment",
  xlim = c(0, 22),                       # <- add a little headroom past 21
  break.x.by = 7,
  ylab = "Cumulative mortality",
  title = "Unadjusted KM Survival Curves",
  subtitle = "Initiate Combination Therapy vs Monotherapy",
  legend.labs = c("Polymyxin-based Combination therapy", "Polymyxin Monotherapy"),
  legend.title = "Treatment strategy",
  palette = c("#d95f02", "#2c7fb8"),
  risk.table.y.text.col = TRUE,
  risk.table.y.text = TRUE,
  risk.table.title = "Number at risk",
  fontsize = 4,
  ggtheme = theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(color = "grey92", linewidth = 0.35),
      axis.line          = element_line(color = "grey20", linewidth = 0.4),
      axis.ticks         = element_line(color = "grey20", linewidth = 0.4),
      plot.title         = element_text(face = "bold", size = 16, margin = margin(b = 2)),
      plot.subtitle      = element_text(size = 12, color = "grey35", margin = margin(b = 10)),
      legend.title       = element_text(face = "bold", size = 11),
      legend.text        = element_text(size = 10.5),
      legend.position    = "top",
      plot.margin        = margin(t = 10, r = 25, b = 5, l = 10)   # <- right padding
    ),
  tables.theme = theme_minimal(base_size = 12) +
    theme(
      panel.grid    = element_blank(),
      axis.line.x   = element_line(color = "grey20", linewidth = 0.4),
      axis.ticks.x  = element_line(color = "grey20", linewidth = 0.4),
      axis.line.y   = element_line(color = "grey20", linewidth = 0.4),
      axis.ticks.y  = element_line(color = "grey20", linewidth = 0.4),
      axis.text.x   = element_text(size = 10, color = "grey30"),
      plot.title    = element_text(size = 12, face = "plain", color = "grey30", hjust = 0),
      plot.margin   = margin(t = 5, r = 25, b = 5, l = 10)         # <- match plot's right padding
    )
)

###
plot_02$plot <- plot_02$plot +
  scale_y_continuous(
    limits = c(0, 1.0),
    breaks = seq(0, 1.0, by = 0.1),
    labels = scales::percent_format(accuracy = 1)
  )

###
plot_02

# ---- 10. Save ---------------------------------------------------------------
png("sens_unadj_KM_plmx_21.png", width = 12, height = 6, units = "in", res = 320, bg = "white")
print(plot_02)
dev.off()

#-------------------------------------#----------------------------------------#



### Log-rank test
#-------------------------------------
#-------------------------------------
logrank <- survdiff(Surv(time, event) ~ group, data = df_used)

print(logrank)

#
p_value <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)
cat("Log-rank p-value =", p_value, "\n")

#-------------------------------------
#-------------------------------------
### Cox regression
#
cox_model <- coxph(Surv(time, event) ~ group, data = df_used)

#
summary(cox_model)

#
cbind(
  HR = exp(coef(cox_model)),
  exp(confint(cox_model))
)

#-------------------------------------#-----------------------------------------





###-----------------------------------------------------------------------------
### ^.^ ###