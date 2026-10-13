
### 20260904
### Antibiotic Regimens


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
                 stringr)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean data/clean/clean/combine_group/baseline_all_clean_combine.RData")
treat_combo <- readRDS("data/clean data/clean/clean/combine_group/treat_combo.RData")
id_bl <- readRDS("data/clean data/clean/clean/final/tables_bl.RData")

### Load fonts
loadfonts()

#-------------------------------------
### Filter IDs
df <- df[which(df$recordid %in% id_bl$id),]
treat_combo <- treat_combo[which(treat_combo$recordid %in% id_bl$id),]

#-------------------------------------

## ── 1. Define drug components per organism ──────────────────────────────
# Monotherapy: exact match (single drug only)
# Combination: POSITIONAL — take antibiotics in the order listed (split on "+"),
#              slide a 2-drug window left-to-right, and use the FIRST window
#              that satisfies one backbone+partner rule for that organism.
#              Any remaining drugs in the string beyond that pair are ignored.
# CRA -----------------------------------------------------------------
mono_CRA <- c("Colistin", "Polymyxin B", "Ampicillin/sulbactam", "Cefoperazone/sulbactam", "Tigecycline", "Eravacycline")
combo_CRA_backbone <- c("Colistin", "Polymyxin B")
combo_CRA_partner  <- c("Meropenem", "Imipenem", "Doripenem", "Ampicillin/sulbactam", "Cefoperazone/sulbactam", "Tigecycline", "Eravacycline")
# CRE -----------------------------------------------------------------
mono_CRE <- c("Colistin", "Polymyxin B", "Ceftazidime/avibactam", "Tigecycline", "Eravacycline")
combo_CRE_backbone_A <- c("Colistin", "Polymyxin B")
combo_CRE_partner_A  <- c("Meropenem", "Imipenem", "Doripenem", "Tigecycline", "Eravacycline")
combo_CRE_backbone_B <- c("Ceftazidime/avibactam")
combo_CRE_partner_B  <- c("Tigecycline", "Eravacycline", "Amikacin", "Gentamicin", "Tobramycin")
# CRP -----------------------------------------------------------------
mono_CRP <- c("Colistin", "Polymyxin B", "Ceftazidime/avibactam", "Ceftolozane/tazobactam")
combo_CRP_backbone_A <- c("Colistin", "Polymyxin B")
combo_CRP_partner_A  <- c("Meropenem", "Imipenem", "Doripenem")
combo_CRP_backbone_B <- c("Ceftazidime/avibactam")
combo_CRP_partner_B  <- c("Amikacin", "Gentamicin", "Tobramycin")
combo_CRP_backbone_C <- c("Ceftolozane/tazobactam")
combo_CRP_partner_C  <- c("Amikacin", "Gentamicin", "Tobramycin")

## ── 2. Sub-groups used only to build a readable label for each combo ────
partner_subgroups_CRA <- list(
  "Meropenem/Imipenem/Doripenem"          = c("Meropenem", "Imipenem", "Doripenem"),
  "Ampicillin/sulbactam or Cefoperazone/sulbactam" = c("Ampicillin/sulbactam", "Cefoperazone/sulbactam"),
  "Tigecycline/Eravacycline"              = c("Tigecycline", "Eravacycline")
)
partner_subgroups_CRE_A <- list(
  "Meropenem/Imipenem/Doripenem" = c("Meropenem", "Imipenem", "Doripenem"),
  "Tigecycline/Eravacycline"     = c("Tigecycline", "Eravacycline")
)
partner_subgroups_CRE_B <- list(
  "Tigecycline/Eravacycline"        = c("Tigecycline", "Eravacycline"),
  "Amikacin/Gentamicin/Tobramycin"  = c("Amikacin", "Gentamicin", "Tobramycin")
)
partner_subgroups_CRP_A <- list("Meropenem/Imipenem/Doripenem" = c("Meropenem", "Imipenem", "Doripenem"))
partner_subgroups_CRP_B <- list("Amikacin/Gentamicin/Tobramycin" = c("Amikacin", "Gentamicin", "Tobramycin"))
partner_subgroups_CRP_C <- list("Amikacin/Gentamicin/Tobramycin" = c("Amikacin", "Gentamicin", "Tobramycin"))

backbone_label <- function(v) paste(v, collapse = " or ")

## ── 3. Helpers: normalized, token-level drug matching ───────────────────
## Normalizes formatting differences like "Ampicillin-sulbactam" vs
## "Ampicillin/sulbactam" vs spacing/case, so tokens compare reliably.
normalize_drug <- function(x) {
  x %>% str_trim() %>% tolower() %>% str_replace_all("[-/]", "") %>% str_replace_all("\\s+", "")
}

## Is a single drug token an exact (normalized) match to any drug in drug_list?
is_drug_in <- function(token, drug_list) {
  any(normalize_drug(token) == normalize_drug(drug_list))
}

## Which named subgroup does a single drug token belong to?
which_subgroup_token <- function(token, subgroup_list) {
  ntoken <- normalize_drug(token)
  hit <- names(subgroup_list)[map_lgl(subgroup_list, ~ any(ntoken == normalize_drug(.x)))]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

## Check whether an unordered pair of tokens satisfies one backbone+partner rule.
## Returns list(match = TRUE/FALSE, label = ...)
check_pair_generic <- function(tok1, tok2, backbone_list, partner_list, subgroup_list) {
  if (is_drug_in(tok1, backbone_list) && is_drug_in(tok2, partner_list)) {
    sg <- which_subgroup_token(tok2, subgroup_list)
    return(list(match = TRUE, label = paste(backbone_label(backbone_list), "plus", sg)))
  }
  if (is_drug_in(tok2, backbone_list) && is_drug_in(tok1, partner_list)) {
    sg <- which_subgroup_token(tok1, subgroup_list)
    return(list(match = TRUE, label = paste(backbone_label(backbone_list), "plus", sg)))
  }
  list(match = FALSE)
}

## ── 4. Per-organism pair checkers (try each rule set in priority order) ─
check_pair_CRA <- function(tok1, tok2) {
  check_pair_generic(tok1, tok2, combo_CRA_backbone, combo_CRA_partner, partner_subgroups_CRA)
}

check_pair_CRE <- function(tok1, tok2) {
  r <- check_pair_generic(tok1, tok2, combo_CRE_backbone_A, combo_CRE_partner_A, partner_subgroups_CRE_A)
  if (r$match) return(r)
  check_pair_generic(tok1, tok2, combo_CRE_backbone_B, combo_CRE_partner_B, partner_subgroups_CRE_B)
}

check_pair_CRP <- function(tok1, tok2) {
  r <- check_pair_generic(tok1, tok2, combo_CRP_backbone_A, combo_CRP_partner_A, partner_subgroups_CRP_A)
  if (r$match) return(r)
  r <- check_pair_generic(tok1, tok2, combo_CRP_backbone_B, combo_CRP_partner_B, partner_subgroups_CRP_B)
  if (r$match) return(r)
  check_pair_generic(tok1, tok2, combo_CRP_backbone_C, combo_CRP_partner_C, partner_subgroups_CRP_C)
}

## ── 5. Core classifier for one treat_combo string given its organism ────
classify_combo <- function(org, text) {
  # Split on "+" preserving order, e.g. "Ampicillin-sulbactam+Colistin+Tigecycline"
  # -> c("Ampicillin-sulbactam", "Colistin", "Tigecycline")
  tokens <- str_split(text, "\\+")[[1]] %>% str_trim()
  tokens <- tokens[tokens != ""]

  if (length(tokens) == 0) {
    return(tibble(category = "Unclassified", label = str_trim(text)))
  }

  # ---- Monotherapy: single drug only ----
  if (length(tokens) == 1) {
    mono_list <- switch(org, "CRA" = mono_CRA, "CRE" = mono_CRE, "CRP" = mono_CRP, character(0))
    if (is_drug_in(tokens[1], mono_list)) {
      idx <- which(normalize_drug(mono_list) == normalize_drug(tokens[1]))[1]
      return(tibble(category = "Monotherapy", label = mono_list[idx]))
    }
    return(tibble(category = "Unclassified", label = tokens[1]))
  }

  # ---- Combination: anchor on each drug in order, scan forward for its ----
  # ---- first valid partner (may skip over unrelated drugs in between)  ----
  # e.g. "Ceftazidime/avibactam+Meropenem+Tigecycline" -> anchor on
  # Ceftazidime/avibactam fails with Meropenem, but matches Tigecycline,
  # so Meropenem is skipped and that pair is used.
  n <- length(tokens)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tok1 <- tokens[i]
      tok2 <- tokens[j]
      result <- switch(org,
        "CRA" = check_pair_CRA(tok1, tok2),
        "CRE" = check_pair_CRE(tok1, tok2),
        "CRP" = check_pair_CRP(tok1, tok2),
        list(match = FALSE)
      )
      if (isTRUE(result$match)) {
        return(tibble(category = "Combination", label = result$label))
      }
    }
  }

  # ---- No valid consecutive pair found anywhere in the string ----
  tibble(category = "Unclassified", label = str_trim(text))
}

## ── 5.5. Grouped display labels for monotherapy (mirrors combo subgroups) ──
## Drugs listed together here will be combined into ONE row in the summary
## tables, e.g. "Ampicillin-sulbactam or Cefoperazone-sulbactam".
## Any monotherapy drug NOT listed here keeps its own individual label.
mono_label_groups <- list(
  "Colistin or Polymyxin B"                        = c("Colistin", "Polymyxin B"),
  "Ampicillin/sulbactam or Cefoperazone/sulbactam" = c("Ampicillin/sulbactam", "Cefoperazone/sulbactam"),
  "Tigecycline or Eravacycline"                    = c("Tigecycline", "Eravacycline")
)

group_mono_label <- function(drug) {
  hit <- names(mono_label_groups)[map_lgl(mono_label_groups, ~ drug %in% .x)]
  if (length(hit) == 0) return(drug)
  hit[1]
}

## ── 6. Apply row by row ──────────────────────────────────────────────
classified <- map2_dfr(df$org_r_combined, df$matched_treat_combo, classify_combo)
df <- bind_cols(df, classified)

## Grouped label used for reporting/summarizing (keeps `label` as the
## original per-drug value in case you need it for anything else)
df <- df %>%
  mutate(label_grouped = if_else(
    category == "Monotherapy",
    map_chr(label, group_mono_label),
    label
  ))

## Inspect anything that didn't match a rule — fix drug lists above if needed
unclassified <- df %>% filter(category == "Unclassified") %>%
  count(org_r_combined, matched_treat_combo, sort = TRUE)
if (nrow(unclassified) > 0) {
  message("Unclassified combos - review and adjust component lists:")
  print(unclassified)
}

## ── 7. Summarize n and % within each organism ────────────────────────
## % is calculated within organism, using the GROUPED label
## (i.e., % of all CRA regimens that were "Colistin or Polymyxin B").
## matched_treat_combo is kept so you can see exactly which raw strings
## rolled up into each label.

# 7a. Top-level summary (one row per organism/category/grouped label) ---
summary_tbl <- df %>%
  count(org_r_combined, category, label_grouped, name = "n") %>%
  group_by(org_r_combined) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  arrange(org_r_combined, category, desc(n))

#print(summary_tbl, n = 100)

# 7b. Detailed summary: breaks each grouped label down by raw matched_treat_combo
summary_detail_tbl <- df %>%
  count(org_r_combined, category, label_grouped, matched_treat_combo, name = "n") %>%
  group_by(org_r_combined) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%   # % of all regimens for that organism
  ungroup() %>%
  arrange(org_r_combined, category, label_grouped, desc(n))

#print(summary_detail_tbl, n = 200)

report_tbl <- summary_detail_tbl %>%
  mutate(n_pct = sprintf("%d (%.1f%%)", n, pct)) %>%
  dplyr::select(org_r_combined, category, label_grouped, matched_treat_combo, n_pct)

#print(report_tbl, n = 200)

# write.csv(summary_tbl, "combo_therapy_summary.csv", row.names = FALSE)

###-------------------------------------
### Distribution
#
table(df$infection_types)

table(df$infection_types, df$org_r_combined)

###
regimen_summary <- df %>%
  count(infection_types, label_grouped, name = "n") %>%
  group_by(infection_types) %>%
  mutate(
    total = sum(n),
    percent = round(100 * n / total, 1)
  ) %>%
  ungroup() %>%
  arrange(infection_types, desc(n))


###-----------------------------------------------------------------------------
### ^.^ ###
