library(readxl)
library(ggplot2)
library(dplyr)
library(grid)
library(patchwork)

file_path <- "../analysis_results"

# Sheet list (one sheet per subgroup)
sheet_list <- c("smk1", "sbp2", "hba1c", "nonhdl")

# Color mapping: 0/1 = family history status, .1 suffix = unfavorable subgroup
value_colors <- c(
  "0"   = "#ABD9E9",  # No FH, favorable
  "0.1" = "#2B83BA",  # No FH, unfavorable
  "1"   = "#FDAE61",  # FH, favorable
  "1.1" = "#D7191C"   # FH, unfavorable
)

# Plot adjusted cumulative mortality curves for a given subgroup sheet (Figure 3)
plot_by_sheet <- function(sheet_name) {

  df <- read_excel(paste0(file_path, "/CVD_mortality_results_v5.xlsx"),
                   sheet = sheet_name) %>%
    mutate(variable = as.character(variable))

  # Create group labels
  df <- df %>%
    mutate(Group = case_when(
      sheet_name == "smk1"   & variable == "0"   ~ "No FH, Non-smoker",
      sheet_name == "smk1"   & variable == "0.1" ~ "No FH, Former/Current smoker",
      sheet_name == "smk1"   & variable == "1"   ~ "FH, Non-smoker",
      sheet_name == "smk1"   & variable == "1.1" ~ "FH, Former/Current smoker",

      sheet_name == "bmi1"   & variable == "0"   ~ "No FH, BMI<25",
      sheet_name == "bmi1"   & variable == "0.1" ~ "No FH, BMI≥25",
      sheet_name == "bmi1"   & variable == "1"   ~ "FH, BMI<25",
      sheet_name == "bmi1"   & variable == "1.1" ~ "FH, BMI≥25",

      sheet_name == "PA"     & variable == "0"   ~ "No FH, ≥3 times/week",
      sheet_name == "PA"     & variable == "0.1" ~ "No FH, <3 times/week",
      sheet_name == "PA"     & variable == "1"   ~ "FH, ≥3 times/week",
      sheet_name == "PA"     & variable == "1.1" ~ "FH, <3 times/week",

      sheet_name == "bp1"    & variable == "0"   ~ "No FH, Normal BP (<120/<80)",
      sheet_name == "bp1"    & variable == "0.1" ~ "No FH, Elevated BP",
      sheet_name == "bp1"    & variable == "1"   ~ "FH, Normal BP (<120/<80)",
      sheet_name == "bp1"    & variable == "1.1" ~ "FH, Elevated BP",

      sheet_name == "sbp2"   & variable == "0"   ~ "No FH, Normal BP (<140/<90)",
      sheet_name == "sbp2"   & variable == "0.1" ~ "No FH, Elevated BP",
      sheet_name == "sbp2"   & variable == "1"   ~ "FH, Normal BP (<140/<90)",
      sheet_name == "sbp2"   & variable == "1.1" ~ "FH, Elevated BP",

      sheet_name == "hba1c"  & variable == "0"   ~ "No FH, HbA1c≤5.6",
      sheet_name == "hba1c"  & variable == "0.1" ~ "No FH, HbA1c>5.6",
      sheet_name == "hba1c"  & variable == "1"   ~ "FH, HbA1c≤5.6",
      sheet_name == "hba1c"  & variable == "1.1" ~ "FH, HbA1c>5.6",

      sheet_name == "nonhdl" & variable == "0"   ~ "No FH, Non-HDL<130",
      sheet_name == "nonhdl" & variable == "0.1" ~ "No FH, Non-HDL≥130",
      sheet_name == "nonhdl" & variable == "1"   ~ "FH, Non-HDL<130",
      sheet_name == "nonhdl" & variable == "1.1" ~ "FH, Non-HDL≥130"
    )) %>%
    mutate(
      time     = as.numeric(time),
      event    = as.numeric(event),
      variable = factor(variable, levels = c("0", "0.1", "1", "1.1"))
    )

  title_text <- switch(sheet_name,
    "smk1"   = "Smoking status",
    "bmi1"   = "BMI group",
    "PA"     = "Physical activity",
    "bp1"    = "Blood pressure (cutoff 120/80)",
    "sbp2"   = "Blood pressure",
    "hba1c"  = "HbA1c level",
    "nonhdl" = "Non-HDL cholesterol"
  )

  labels_vec <- df %>%
    distinct(variable, Group) %>%
    arrange(factor(variable, levels = c("0", "0.1", "1", "1.1"))) %>%
    pull(Group)

  ggplot(df, aes(time, event, color = variable)) +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    scale_color_manual(values = value_colors, labels = labels_vec) +
    labs(
      title = title_text,
      x     = "Follow-up time (years)",
      y     = "Cumulative mortality"
    ) +
    coord_cartesian(ylim = c(0, 0.0025)) +
    theme_classic(base_size = 14) +
    theme(
      plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title       = element_text(size = 14, face = "bold"),
      axis.text        = element_text(size = 12),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 12),
      legend.margin    = margin(t = -5, b = 0),
      plot.margin      = margin(t = 5, r = 11, b = 10, l = 10)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
}

plots <- lapply(sheet_list, plot_by_sheet)

subgroup_fig <- (plots[[1]] | plots[[2]]) /
  (plots[[3]] | plots[[4]])

ggsave(paste0(file_path, "/figure/subgroup_figure.png"),
       plot = subgroup_fig, width = 12, height = 8, dpi = 600)

ggsave(paste0(file_path, "/figure/subgroup_figure.tiff"),
       plot = subgroup_fig, width = 12, height = 8, dpi = 300)
