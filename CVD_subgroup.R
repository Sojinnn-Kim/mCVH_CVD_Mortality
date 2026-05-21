# Note: this script assumes df_analysis is loaded from CVD_analysis.R
library(survminer)

df_subgroup <- df_analysis %>%
  mutate(
    smk_group1   = factor(if_else(smk %in% c(1, 2), 1, 0)),
    smk_group2   = factor(if_else(smk %in% c(0, 1), 0, 1)),
    bmi_group1   = factor(ifelse(bmi < 25, 0, 1)),
    bmi_group2   = factor(case_when(
      bmi >= 25              ~ 2,
      bmi >= 18 & bmi < 25  ~ 1,
      bmi < 18               ~ 0,
      TRUE ~ NA_real_
    )),
    met_h_group  = factor(if_else(met_h >= 6, 0, 1)),
    pa_group     = factor(physical_activity),
    psqi_group   = factor(if_else(psqi > 5, 0, 1)),
    bp_group1    = factor(if_else(sbp < 120 & dbp < 80, 0, 1)),
    bp_group2    = factor(if_else(sbp < 140 & dbp < 90, 0, 1)),
    hba1c_group1 = factor(if_else(hba1c < 6.5, 0, 1)),
    hba1c_group2 = factor(if_else(hba1c <= 5.6, 0, 1)),
    hba1c_group3 = factor(case_when(
      hba1c <= 5.6             ~ 0,
      hba1c > 5.6 & hba1c < 6.5 ~ 1,
      hba1c >= 6.5             ~ 2,
      TRUE ~ NA_real_
    )),
    nonHDL_group = factor(if_else(nonHDL < 130, 0, 1))
  )

# Extract HR (95% CI) table from a subgroup Cox model
cox_HR_table <- function(model, digits = 3) {
  res <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(
      HR    = round(estimate, digits),
      Lower = round(conf.low, digits),
      Upper = round(conf.high, digits)
    ) %>%
    select(term, HR, Lower, Upper)
  rbind(c('group0.0', 1, 1, 1), res[1:3, ])
}


# ============================== Subgroup Analysis ==============================
# Set var to run a specific subgroup: 'smk' | 'sbp' | 'hba1c' | 'hdl'

var <- 'smk'

if (var == 'smk') {
  var_interest <- 'smk_group1'
} else if (var == 'sbp') {
  var_interest <- 'bp_group2'
} else if (var == 'hba1c') {
  var_interest <- 'hba1c_group2'
} else if (var == 'hdl') {
  var_interest <- 'nonHDL_group'
}

imp_list  <- c('sbp', 'dbp', 'history_diabetes', 'hba1c', 'tchol', 'hdl', 'smk', 'bmi', 'met')
adj_vars  <- c('sex', 'age', 'center', 'year', 'education')
adj_vars2 <- c(adj_vars, 'ldl', 'history_hypertension', 'history_dyslipidemia')
sub_adj_vars <- setdiff(imp_list, var)

# Incidence rate by subgroup x family history
df_subgroup <- df_subgroup %>%
  mutate(group = interaction(!!sym(var_interest), family_CVD))
incidence_table <- df_subgroup %>%
  group_by(group) %>%
  summarise(
    N        = n(),
    events   = sum(death),
    duration = sum(fu_years),
    IR       = (events / duration) * 1000,
    .groups  = 'drop'
  )
as.data.frame(incidence_table)
summary(df_subgroup$fu_years)

# Cox model
cox_formula <- as.formula(paste("Surv(fu_years, death) ~ group +",
                                paste(c(adj_vars2, sub_adj_vars), collapse = "+")))
sub_fit <- coxph(cox_formula, data = df_subgroup)
summary(sub_fit)
cox_HR_table(sub_fit)

# p for interaction (subgroup variable x family history)
cox_formula <- as.formula(paste("Surv(fu_years, death) ~", var_interest, "* family_CVD +",
                                paste(c(adj_vars2, sub_adj_vars), collapse = "+")))
cox_int <- coxph(cox_formula, data = df_subgroup)
summary(cox_int)$coefficients

# Adjusted survival curves for the subgroup
set.seed(1886)
sub_sample <- df_subgroup[sample(nrow(df_subgroup), 5000), ]
new_df <- df_subgroup %>%
  group_by(group) %>%
  summarise(across(all_of(sub_adj_vars), ~mean(.x, na.rm = TRUE))) %>%
  mutate(group = unique(df_subgroup$group)) %>%
  as.data.frame()
res <- ggadjustedcurves(sub_fit,
                        data        = as.data.frame(sub_sample),
                        newdata     = new_df,
                        variable    = 'group',
                        return_data = TRUE)
df_curve       <- res$data
df_curve$event <- 1 - df_curve$surv
write.csv(df_curve, paste0('../analysis_results/CVD_subgroup_', var, '.csv'), row.names = FALSE)
