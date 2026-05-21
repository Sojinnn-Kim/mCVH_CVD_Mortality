library(readxl)
library(ggplot2)
library(dplyr)

file_path <- "../analysis_results"

# Plot adjusted cumulative mortality curves for the whole cohort (Figure 2)
plot_whole <- function() {

  df <- read_excel(paste0(file_path, "/CVD_mortality_results_v5.xlsx"), sheet = "whole") %>%
    mutate(
      time     = as.numeric(time),
      event    = as.numeric(event),
      variable = as.character(variable)
    )

  # Create group labels based on family history and mCVH score
  df <- df %>%
    mutate(
      FH        = ifelse(grepl("\\.1$", variable), "Yes FH", "No FH"),
      Lifestyle = ifelse(grepl("^High", variable), "High mCVH score", "Low mCVH score"),
      Group     = paste(FH, Lifestyle, sep = ", ")
    ) %>%
    mutate(
      Group = factor(Group, levels = c(
        "No FH, High mCVH score",
        "No FH, Low mCVH score",
        "Yes FH, High mCVH score",
        "Yes FH, Low mCVH score"
      ))
    )

  color_map <- c(
    "No FH, High mCVH score"  = "#ABD9E9",
    "No FH, Low mCVH score"   = "#2B83BA",
    "Yes FH, High mCVH score" = "#FDAE61",
    "Yes FH, Low mCVH score"  = "#D7191C"
  )

  ggplot(df, aes(x = time, y = event, color = Group)) +
    geom_line(linewidth = 1.2, alpha = 0.95) +
    scale_color_manual(values = color_map) +
    labs(
      title = "",
      x     = "Follow-up time (years)",
      y     = "Cumulative mortality",
      color = ""
    ) +
    coord_cartesian(ylim = c(0, max(df$event, na.rm = TRUE))) +
    theme_classic(base_size = 14) +
    theme(
      plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title      = element_text(size = 14, face = "bold"),
      axis.text       = element_text(size = 12),
      legend.position = "bottom",
      legend.title    = element_text(size = 13),
      legend.text     = element_text(size = 12),
      legend.margin   = margin(t = -5, b = 0)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
}

whole_plot <- plot_whole()
whole_plot

ggsave(paste0(file_path, "/figure/main_figure.png"),
       plot = whole_plot, width = 10, height = 6, dpi = 600)

ggsave(paste0(file_path, "/figure/main_figure.tiff"),
       plot = whole_plot, width = 10, height = 6, dpi = 300)
