df_base = df_analysis

# 1
df_base %>% distinct(new_patientid) %>% count() 

# 2
df_base %>% distinct(new_patientid, .keep_all = TRUE) %>% count(family_CVD)

# 3
get_anova_p = function(varname, data, group){
  formula = as.formula(paste(varname, "~", group))
  fit = aov(formula, data=data)
  p = round(summary(fit)[[1]][["Pr(>F)"]][1],4)
  return (ifelse (p<0.001, "p<0.001",p))
}
get_chisq_p = function(varname, data, group){
  tbl = table(data[[varname]], data[[group]])
  p = round(chisq.test(tbl)$p.value, 4)
  return (ifelse (p<0.001, "p<0.001",p))
}

df_base = df_base %>% 
  mutate(t_group = ntile(LE6, 3))

# 3--cont
cont_list = c('age', 'sbp', 'dbp', 'hba1c', 'nonHDL', 'bmi',
              'met_h', 'homair', 'hscrp', 'ldl', 'LE6', 't_group')

res_list = list()
for (var in cont_list){
  mu = mean(df_base[[var]], na.rm=TRUE)
  sd = sd(df_base[[var]], na.rm=TRUE)
  p_value = get_anova_p(var, df_base, 'family_CVD')
  mu_sd = paste0(round(mu,2), "±", round(sd,2))
  res_list[[var]] = data.frame(Variable=var, mu_sd, p_value)
}

result_df = do.call(rbind, res_list)
rownames(result_df) = NULL
result_df
write.table(result_df, file='clipboard', sep='\t', row.names = FALSE, col.names = FALSE)

# 3--cat

cat_list = c('sex', 'education', 'history_diabetes', 'smk', 'med_hypertension', 't_group')

res_list = list()
for (var in cat_list){
  tab = table(df_base[[var]])
  n = as.numeric(tab)
  pct = round(100*n/sum(n),2)
  n_pct = paste0(formatC(n, format='d', big.mark = ",")," (",pct,")")
  p_value = get_chisq_p(var, df_base, 'family_CVD')
  res_list[[var]] = data.frame(Variable=var, level = names(tab), 
                               n_pct, p_value)
}

result_df = do.call(rbind, res_list)
rownames(result_df) = NULL
result_df = result_df %>%
  filter((Variable %in% c('history_diabetes', 'med_hypertension') & level ==1) |
           (Variable == 'education' & level==5) |
           (Variable %in% c('sex', 'smk'))
  ) %>%
  select(-level)
result_df
write.table(result_df, file='clipboard', sep='\t', row.names = FALSE, col.names = FALSE)

# 4
# --4-cont
res_list = list()
for (var in cont_list){
  tmp = df_base %>%
    group_by(family_CVD) %>%
    summarise(
      mu = mean(!!sym(var), na.rm=TRUE),
      sd = sd(!!sym(var), na.rm=TRUE),
      .groups = 'drop'
    ) %>%
    mutate(Variable = var,
           mu_sd = paste0(round(mu,2), "±", round(sd,2)))
  res_list[[var]] = tmp
}
result_df = bind_rows(res_list)
result_wide = result_df %>%
  select(Variable, family_CVD, mu_sd) %>%
  pivot_wider(names_from=family_CVD, values_from=mu_sd)
result_wide
write.table(result_wide[2:3], file='clipboard', sep='\t', row.names = FALSE, col.names = FALSE)

# --4-cat
res_list = list()
for (var in cat_list){
  tmp = df_base %>%
    group_by(family_CVD, !!sym(var)) %>%
    summarise(
      n = n(),.groups = 'drop') %>%
    group_by(family_CVD) %>%
    mutate(pct = round(100*n/sum(n),2)) %>%
    ungroup() %>%
    rename(level = !!sym(var)) %>%
    mutate(Variable = var,
           n_pct = paste0(formatC(n, format='d', big.mark = ",")," (",pct,")"))
  res_list[[var]] = tmp
}
result_df = do.call(rbind, res_list)
result_wide = result_df %>%
  filter((Variable %in% c('history_diabetes', 'med_hypertension') & level ==1) |
           (Variable == 'education' & level==5) |
           (Variable %in% c('sex', 'smk'))
  ) %>%
  select(Variable, level, family_CVD, n_pct)
result_wide = result_wide %>%
  pivot_wider(names_from=family_CVD, values_from=n_pct)
write.table(result_wide[,3:4], file='clipboard', sep='\t', row.names = FALSE, col.names = FALSE)


result = df_base %>%
  group_by(t_group) %>%
  summarise(
    mu = mean(LE6, na.rm=TRUE),
    sd = sd(LE6, na.rm=TRUE)
  ) %>%
  mutate(Variable = var,
         mu_sd = paste0(round(mu,2), "±", round(sd,2))) %>%
  select(Variable, mu_sd)
result[,2]
write.table(result[,2], file='clipboard', sep='\t', row.names = FALSE, col.names = FALSE)

result = df_base %>%
  group_by(t_group, family_CVD) %>%
  summarise(
    mu = mean(LE6, na.rm=TRUE),
    sd = sd(LE6, na.rm=TRUE),
    p_value = get_anova_p(var, df_base, 'family_CVD')
  ) %>%
  mutate(Variable = var,
         mu_sd = paste0(round(mu,2), "±", round(sd,2))) %>%
  select(Variable, family_CVD, mu_sd, p_value) %>%
  pivot_wider(names_from=family_CVD, values_from=mu_sd)
result[,4:5]
write.table(result[,4:5], file='clipboard', sep='\t', row.names = FALSE, col.names = FALSE)

df_base %>%
  group_by(t_group) %>%
  do(broom::tidy(t.test(LE6 ~ family_CVD, data=.)))
