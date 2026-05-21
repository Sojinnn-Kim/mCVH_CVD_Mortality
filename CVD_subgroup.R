df_subgroup = df_analysis %>%
  mutate(smk_group1 = if_else(smk %in% c(1,2), 1, 0),
         smk_group2 = if_else(smk %in% c(0,1), 0, 1),
         bmi_group1 = ifelse(bmi<25, 0, 1),
         bmi_group2 = case_when(
           bmi>=25 ~ 2,
           bmi>=18 & bmi<25 ~ 1,
           bmi<18 ~ 0,
           TRUE ~ NA_real_
         ),
         met_h_group = if_else(met_h>=6, 0, 1),
         pa_group = physical_activity,
         psqi_group = if_else(psqi>5, 0, 1),
         bp_group1 = if_else(sbp<120 & dbp<80, 0, 1),
         bp_group2 = if_else(sbp<140 & dbp<90, 0, 1),
         hba1c_group1 = if_else(hba1c<6.5, 0, 1),
         hba1c_group2 = if_else(hba1c<=5.6, 0, 1),
         hba1c_group3 = case_when(
           hba1c<=5.6 ~ 0,
           hba1c>5.6 & hba1c<6.5 ~ 1,
           hba1c>=6.5 ~ 2,
           TRUE ~ NA_real_
         ),
         nonHDL_group = if_else(nonHDL<130, 0, 1)
  ) %>%
  mutate(smk_group1 = factor(smk_group1),
         smk_group2 = factor(smk_group2),
         bmi_group1 = factor(bmi_group1),
         bmi_group2 = factor(bmi_group2),
         met_h_group = factor(met_h_group),
         pa_group = factor(pa_group),
         psqi_group = factor(psqi_group),
         bp_group1 = factor(bp_group1),
         bp_group2 = factor(bp_group2),
         hba1c_group1 = factor(hba1c_group1),
         hba1c_group2 = factor(hba1c_group2),
         hba1c_group3 = factor(hba1c_group3),
         nonHDL_group = factor(nonHDL_group))

# ----------
cox_HR_table = function(model, digits=3, copy_clipboard=TRUE){
  res = broom::tidy(model, exponentiate=TRUE, conf.int=TRUE) %>%
    mutate(HR = round(estimate, digits),
           Lower = round(conf.low, digits),
           Upper = round(conf.high, digits)) %>%
    select(term, HR, Lower, Upper)
  final = rbind(c('group0.0',1,1,1), res[1:3,])
  if (copy_clipboard) {write.table(final[,2:4], file='clipboard', sep='\t', row.names = FALSE, col.names = FALSE)}
  return(final)
}


#########################

# subgroup cox
imp_list = c('sbp', 'dbp', 'history_diabetes', 'hba1c', 'tchol', 'hdl', 'smk',  'bmi','met')

var = 'hdl' # smk sbp hba1c hdl

if (var == 'smk') {var_interest = 'smk_group1'
}else if (var == 'sbp') {var_interest = 'bp_group2'
}else if (var == 'hba1c') {var_interest = 'hba1c_group2'
}else if (var == 'hdl') {var_interest = 'nonHDL_group'}


# incidence_table
df_subgroup = df_subgroup %>%
  mutate(group = interaction(!!sym(var_interest), family_CVD))
incidence_table = df_subgroup  %>%  
  group_by(group) %>%
  summarise(N = n(),
            events = sum(death),
            duration = sum(fu_years),
            IR = (events/duration)*1000)
as.data.frame(incidence_table)
write.table(incidence_table[,3:5], file='clipboard', sep='\t', row.names = FALSE, col.names = FALSE)
summary(df_subgroup$fu_years)

# cox
sub_adj_vars = setdiff(imp_list, var)
adj_vars = c('sex', 'age','center', 'year', 'education')
adj_vars2 = c(adj_vars, 'ldl', 'history_hypertension', 'history_dyslipidemia')
sub_adj_vars = setdiff(imp_list, var)
cox_formula = as.formula(paste("Surv(fu_years, death) ~ group +", paste(c(adj_vars2, sub_adj_vars), collapse="+")))
sub_fit = coxph(cox_formula, data=df_subgroup)
summary(sub_fit)
cox_HR_table(sub_fit)

# p for interaction
cox_formula = as.formula(paste("Surv(fu_years, death) ~", var_interest, "*family_CVD +", paste(c(adj_vars2, sub_adj_vars), collapse="+")))
cox_model2 = coxph(cox_formula, data=df_subgroup)
summary(cox_model2)$coefficients

# data for KM curve
sub_sample = df_subgroup[sample(nrow(df_subgroup), 5000),]
new_df = df_subgroup %>%
  group_by(group) %>%
  summarise(across(all_of(sub_adj_vars), ~mean(.x, na.rm=TRUE))) %>%
  mutate(group = unique(df_subgroup$group)) %>%
  as.data.frame()
res = ggadjustedcurves(sub_fit,
                       data = as.data.frame(sub_sample),
                       newdata = new_df,
                       variable = 'group',
                       # fun = 'event',
                       return_data = TRUE)
df_curve = res$data
df_curve$event = 1 - df_curve$surv
write.csv(df_curve, paste('./CVD_Mortality/CVD_subgroup_', var, '.csv'), row.names = FALSE)

