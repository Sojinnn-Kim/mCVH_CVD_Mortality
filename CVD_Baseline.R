# Note: this script assumes df_analysis is loaded from CVD_analysis.R

df_base <- df_analysis

# Overall N
df_base %>% distinct(new_patientid) %>% count()

# N by family history status
df_base %>% distinct(new_patientid, .keep_all = TRUE) %>% count(family_CVD)

# Helper functions for p-values
get_anova_p <- function(varname, data, group) {
  formula <- as.formula(paste(varname, "~", group))
  fit <- aov(formula, data = data)
  p <- round(summary(fit)[[1]][["Pr(>F)"]][1], 4)
  return(ifelse(p < 0.001, "p<0.001", p))
}
get_chisq_p <- function(varname, data, group) {
  tbl <- table(data[[varname]], data[[group]])
  p <- round(chisq.test(tbl)$p.value, 4)
  return(ifelse(p < 0.001, "p<0.001", p))
}

df_base <- df_base %>%
  mutate(t_group = ntile(LE6, 3))

cont_list <- c('age', 'sbp', 'dbp', 'hba1c', 'nonHDL', 'bmi',
               'met_h', 'homair', 'hscrp', 'ldl', 'LE6', 't_group')
cat_list  <- c('sex', 'education', 'history_diabetes', 'smk', 'med_hypertension', 't_group')


# ============================== Overall Summary ==============================

# Continuous variables
res_list <- list()
for (var in cont_list) {
  mu      <- mean(df_base[[var]], na.rm = TRUE)
  sd      <- sd(df_base[[var]], na.rm = TRUE)
  p_value <- get_anova_p(var, df_base, 'family_CVD')
  res_list[[var]] <- data.frame(Variable = var,
                                mu_sd    = paste0(round(mu, 2), "±", round(sd, 2)),
                                p_value)
}
result_df <- do.call(rbind, res_list)
rownames(result_df) <- NULL
result_df

# Categorical variables
res_list <- list()
for (var in cat_list) {
  tab     <- table(df_base[[var]])
  n       <- as.numeric(tab)
  pct     <- round(100 * n / sum(n), 2)
  p_value <- get_chisq_p(var, df_base, 'family_CVD')
  res_list[[var]] <- data.frame(Variable = var, level = names(tab),
                                n_pct    = paste0(formatC(n, format = 'd', big.mark = ","),
                                                  " (", pct, ")"),
                                p_value)
}
result_df <- do.call(rbind, res_list)
rownames(result_df) <- NULL
result_df <- result_df %>%
  filter(
    (Variable %in% c('history_diabetes', 'med_hypertension') & level == 1) |
      (Variable == 'education' & level == 5) |
      (Variable %in% c('sex', 'smk'))
  ) %>%
  select(-level)
result_df


# ============================== By Family History Status ==============================

# Continuous variables
res_list <- list()
for (var in cont_list) {
  tmp <- df_base %>%
    group_by(family_CVD) %>%
    summarise(
      mu = mean(!!sym(var), na.rm = TRUE),
      sd = sd(!!sym(var), na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(Variable = var,
           mu_sd    = paste0(round(mu, 2), "±", round(sd, 2)))
  res_list[[var]] <- tmp
}
result_wide <- bind_rows(res_list) %>%
  select(Variable, family_CVD, mu_sd) %>%
  pivot_wider(names_from = family_CVD, values_from = mu_sd)
result_wide

# Categorical variables
res_list <- list()
for (var in cat_list) {
  tmp <- df_base %>%
    group_by(family_CVD, !!sym(var)) %>%
    summarise(n = n(), .groups = 'drop') %>%
    group_by(family_CVD) %>%
    mutate(pct = round(100 * n / sum(n), 2)) %>%
    ungroup() %>%
    rename(level = !!sym(var)) %>%
    mutate(Variable = var,
           n_pct    = paste0(formatC(n, format = 'd', big.mark = ","), " (", pct, ")"))
  res_list[[var]] <- tmp
}
result_wide <- do.call(rbind, res_list) %>%
  filter(
    (Variable %in% c('history_diabetes', 'med_hypertension') & level == 1) |
      (Variable == 'education' & level == 5) |
      (Variable %in% c('sex', 'smk'))
  ) %>%
  select(Variable, level, family_CVD, n_pct) %>%
  pivot_wider(names_from = family_CVD, values_from = n_pct)
result_wide


# ============================== mCVH Score by Tertile ==============================

# Overall mean ± SD per tertile
result <- df_base %>%
  group_by(t_group) %>%
  summarise(
    mu = mean(LE6, na.rm = TRUE),
    sd = sd(LE6, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(mu_sd = paste0(round(mu, 2), "±", round(sd, 2)))
result

# Mean ± SD per tertile by family history
result <- df_base %>%
  group_by(t_group, family_CVD) %>%
  summarise(
    mu = mean(LE6, na.rm = TRUE),
    sd = sd(LE6, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(mu_sd = paste0(round(mu, 2), "±", round(sd, 2))) %>%
  select(t_group, family_CVD, mu_sd) %>%
  pivot_wider(names_from = family_CVD, values_from = mu_sd)
result

# t-test for LE6 by family history within each tertile
df_base %>%
  group_by(t_group) %>%
  do(broom::tidy(t.test(LE6 ~ family_CVD, data = .)))
