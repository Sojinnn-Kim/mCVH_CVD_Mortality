library(readr); library(dplyr); library(lubridate); library(stringr)
library(survival); library(survminer); library(rms)
library(VIM); library(mice); library(tidyr)

rm(list = ls())
data <- read_csv("mortality data_20251027_163020/mortality data_20251027_163020.csv")
data %>% distinct(new_patientid) %>% nrow()

# Remove participants with missing family history variables at baseline
data <- data %>%
  filter(visit == 1) %>%
  filter(!(is.na(family_stroke) | is.na(family_heart_dis)))
data %>% distinct(new_patientid) %>% nrow()

# Exclusion: baseline cancer
cancer_id <- data %>%
  filter(history_cancer == 1) %>%
  distinct(new_patientid) %>% pull(new_patientid)

# Exclusion: baseline CVD
cvd_id <- data %>%
  filter(history_heart_dis == 1 | history_stroke == 1 | history_coronary_dis == 1) %>%
  distinct(new_patientid) %>% pull(new_patientid)

# Exclusion: missing key covariates
var_list <- c('sbp', 'dbp', 'history_diabetes', 'hba1c', 'tchol', 'hdl', 'smk', 'bmi', 'met')
kcov_id <- data %>%
  filter(if_any(all_of(var_list), is.na)) %>%
  distinct(new_patientid) %>% pull(new_patientid)

length(unique(c(cancer_id, cvd_id, kcov_id)))

df <- data.frame(data) %>%
  filter(!new_patientid %in% c(cancer_id, cvd_id, kcov_id))
df %>% distinct(new_patientid) %>% nrow()

# Define family history, CVD death outcome (ICD-10 I-codes), and follow-up time
df <- df %>%
  mutate(
    family_CVD = ifelse(family_stroke == 1 | family_heart_dis == 1, 1, 0),
    death = case_when(
      str_detect(new_dth_code1, "^I[0-9]{2}") | str_detect(new_dth_code2, "^I[0-9]{2}") ~ 1L,
      TRUE ~ 0L
    ),
    death_date = as.Date(ISOdate(dth_date1, dth_date2, dth_date3))
  ) %>%
  group_by(new_patientid) %>%
  mutate(first_visit = min(dov)) %>%
  ungroup() %>%
  mutate(fu_years = case_when(
    death == 1 ~ as.numeric(difftime(death_date, first_visit, units = "days")) / 365.25,
    death == 0 ~ as.numeric(difftime(as.Date('2022-12-31'), first_visit, units = 'days')) / 365.25,
    TRUE ~ NA_real_
  ))
df %>% distinct(new_patientid) %>% nrow()
df %>% filter(death == 1) %>% distinct(new_patientid) %>% nrow()

df_base <- df %>%
  arrange(new_patientid) %>%
  filter(visit == 1)
df_base %>% filter(death == 1) %>% distinct(new_patientid) %>% nrow()


# ============================== mCVH Score ==============================
# Compute component scores following the Life's Essential 8 framework (sleep excluded → LE6)

qs <- quantile(df_base$psqi, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE)

df_base <- df_base %>%
  mutate(
    nonHDL  = tchol - hdl,
    dov     = as.Date(dov),
    year    = year(dov),
    smk_quit = year - smk_endyr,
    met_h   = met / 60,
    LE8_bloodpressure = case_when(
      sbp < 115 & dbp < 75                          ~ 100,
      sbp >= 115 & sbp < 125 & dbp < 75             ~ 75,
      (sbp >= 125 & sbp < 135) | (dbp >= 75 & dbp < 85) ~ 50,
      sbp >= 135 & sbp < 155  | (dbp >= 85 & dbp < 95)  ~ 25,
      sbp >= 155 | dbp >= 95                         ~ 0,
      TRUE ~ NA_real_
    ),
    LE8_hba1c = case_when(
      history_diabetes == 0 & hba1c < 5.7             ~ 100,
      history_diabetes == 0 & hba1c >= 5.7 & hba1c <= 6.4 ~ 60,
      hba1c < 7                                        ~ 40,
      hba1c >= 7  & hba1c < 8                          ~ 30,
      hba1c >= 8  & hba1c < 9                          ~ 20,
      hba1c >= 9  & hba1c < 10                         ~ 10,
      hba1c >= 10                                      ~ 0,
      TRUE ~ NA_real_
    ),
    LE8_bloodlipid = case_when(
      nonHDL < 130                       ~ 100,
      nonHDL >= 130 & nonHDL < 160       ~ 60,
      nonHDL >= 160 & nonHDL < 190       ~ 40,
      nonHDL >= 190 & nonHDL < 220       ~ 20,
      nonHDL >= 220                      ~ 0,
      TRUE ~ NA_real_
    ),
    LE8_nicotine = case_when(
      smk == 0 ~ 100,
      smk == 1 ~ 50,
      smk == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    LE8_bmi = case_when(
      bmi < 25                 ~ 100,
      bmi >= 25 & bmi < 30     ~ 70,
      bmi >= 30 & bmi < 35     ~ 30,
      bmi >= 35 & bmi < 40     ~ 15,
      bmi >= 40                ~ 0,
      TRUE ~ NA_real_
    ),
    LE8_physicalactivity = case_when(
      met_h >= 10              ~ 100,
      met_h >= 8  & met_h < 10 ~ 90,
      met_h >= 6  & met_h < 8  ~ 80,
      met_h >= 4  & met_h < 6  ~ 60,
      met_h >= 2  & met_h < 4  ~ 40,
      met_h >  0  & met_h < 2  ~ 20,
      met_h == 0               ~ 0,
      TRUE ~ NA_real_
    ),
    LE8_sleep = case_when(
      psqi <= qs[1] ~ 100,
      psqi <= qs[2] ~ 90,
      psqi <= qs[3] ~ 70,
      psqi <= qs[4] ~ 40,
      psqi <= qs[5] ~ 20,
      psqi >  qs[5] ~ 0,
      TRUE ~ NA_real_
    ),
    LE7 = rowMeans(across(c(LE8_bloodpressure, LE8_hba1c, LE8_bloodlipid,
                             LE8_nicotine, LE8_bmi, LE8_physicalactivity, LE8_sleep))),
    LE6 = rowMeans(across(c(LE8_bloodpressure, LE8_hba1c, LE8_bloodlipid,
                             LE8_nicotine, LE8_bmi, LE8_physicalactivity)), na.rm = TRUE),
    LE5 = rowMeans(across(c(LE8_bloodpressure, LE8_hba1c, LE8_bloodlipid,
                             LE8_nicotine, LE8_bmi)), na.rm = TRUE)
  )

df_whole <- df_base %>%
  filter(!is.na(family_CVD)) %>%
  mutate(
    family_CVD       = factor(family_CVD, levels = c('0', '1')),
    physical_activity = if_else(vigorous_exercise_freq >= 3, 1, 0),
    alc_intake        = if_else(
      (sex == 1 & alc_amount_grams >= 20) | (sex == 0 & alc_amount_grams >= 10), 1, 0
    )
  )

df_male   <- df_whole %>% filter(sex == 1)
df_female <- df_whole %>% filter(sex == 0)

df_analysis <- df_whole
df_analysis %>% filter(death == 1) %>% distinct(new_patientid) %>% nrow()


# ============================== Cox Models ==============================

# Extract HR (95% CI) for the four mCVH x family history groups from an interaction model
Group_HR_table <- function(model,
                           level_LE6_low = "LE6_groupLow",
                           level_FH_yes  = "family_CVD1",
                           interaction   = "LE6_groupLow:family_CVD1") {
  b <- coef(model)
  V <- vcov(model)

  beta_ref  <- 0
  beta_low  <- b[level_LE6_low]
  beta_fh   <- b[level_FH_yes]
  beta_both <- b[level_LE6_low] + b[level_FH_yes] + b[interaction]

  var_ref  <- 0
  var_low  <- V[level_LE6_low, level_LE6_low]
  var_fh   <- V[level_FH_yes, level_FH_yes]
  var_both <- (
    V[level_LE6_low, level_LE6_low] +
      V[level_FH_yes, level_FH_yes] +
      V[interaction, interaction] +
      2 * (V[level_LE6_low, level_FH_yes] +
             V[level_LE6_low, interaction] +
             V[level_FH_yes, interaction])
  )

  HR    <- exp(c(beta_ref, beta_fh, beta_low, beta_both))
  SE    <- sqrt(c(var_ref, var_fh, var_low, var_both))
  lower <- exp(c(beta_ref, beta_fh, beta_low, beta_both) - 1.96 * SE)
  upper <- exp(c(beta_ref, beta_fh, beta_low, beta_both) + 1.96 * SE)

  tibble::tibble(
    Group = c('High_LE6, No_FH', 'High_LE6, Yes_FH', 'Low_LE6, No_FH', 'Low_LE6, Yes_FH'),
    HR    = round(HR, 2),
    lower = round(lower, 2),
    upper = round(upper, 2)
  )
}


# Categorize mCVH score: High (top tertile) vs. Low (bottom two tertiles)
df_analysis <- df_analysis %>%
  mutate(
    LE6_group = ntile(LE6, 3),
    LE6_group = factor(if_else(LE6_group == 3, 'High', 'Low'), levels = c('High', 'Low'))
  )

table(df_analysis$LE6_group)

# Incidence rates by mCVH group and family history
incidence_table <- df_analysis %>%
  group_by(LE6_group, family_CVD) %>%
  summarise(
    N        = n(),
    events   = sum(death),
    duration = sum(fu_years),
    IR       = (events / duration) * 1000,
    .groups  = 'drop'
  )
as.data.frame(incidence_table)
summary(df_analysis$fu_years)

# Age-adjusted model
model_age <- coxph(Surv(fu_years, death) ~ LE6_group * family_CVD + age, df_analysis)
Group_HR_table(model_age)
tail(data.frame(summary(model_age)$coef), 1)

# Model 1: sociodemographic adjustment (age, sex, center, examination year, education)
adj_vars <- c('sex', 'age', 'center', 'year', 'education')
cox_formula <- as.formula(paste("Surv(fu_years, death) ~ LE6_group*family_CVD +",
                                paste(adj_vars, collapse = "+")))
cox_model1 <- coxph(cox_formula, data = df_analysis)
Group_HR_table(cox_model1)
tail(data.frame(summary(cox_model1)$coef), 1)

# Model 2: fully adjusted (Model 1 + LDL, hypertension, dyslipidemia)
adj_vars2 <- c(adj_vars, 'ldl', 'history_hypertension', 'history_dyslipidemia')
cox_formula <- as.formula(paste("Surv(fu_years, death) ~ LE6_group*family_CVD +",
                                paste(adj_vars2, collapse = "+")))
cox_model2 <- coxph(cox_formula, data = df_analysis)
summary(cox_model2)
Group_HR_table(cox_model2)
tail(data.frame(summary(cox_model2)$coef), 1)

# Proportional hazards assumption test
ph_test <- cox.zph(cox_model2)
ph_test


# ============================== Adjusted Survival Curves ==============================

df_curve <- df_analysis
df_curve$group <- interaction(df_curve$LE6_group, df_curve$family_CVD, drop = TRUE)

adj_vars  <- c('sex', 'age', 'center', 'year', 'education')
adj_vars2 <- c(adj_vars, 'ldl', 'history_hypertension', 'history_dyslipidemia')
cox_formula <- as.formula(paste("Surv(fu_years, death) ~ group +",
                                paste(adj_vars2, collapse = "+")))
cox_model2 <- coxph(cox_formula, data = df_curve)

set.seed(1886)
df_sample <- df_curve[sample(nrow(df_curve), 10000), ]
res <- ggadjustedcurves(cox_model2,
                        data     = as.data.frame(df_sample),
                        variable = 'group',
                        method   = 'average',
                        type     = 'survival')
df_survprob        <- res$data
df_survprob$event  <- 1 - df_survprob$surv
write.csv(df_survprob, 'df_survprob.csv', row.names = FALSE)
