# Load required libraries
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(cowplot)
library(boot)

# Set paths
out_dir <- "C:/Users/ms/Desktop/r/chem_data/metabo2/result/section3"
raw_data_path <- "C:/Users/ms/Desktop/r/chem_data/metabo2/Merged_VIP_hub_r_Path2.csv"
tissue_coord_path <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/initial_stat/initial_stat6e/Summary/tissue_comparison/tissue_coordination_summary.csv"

# Read data
raw_data <- read.csv(raw_data_path)
tissue_coord_data <- read.csv(tissue_coord_path)

# Nature theme
nature_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

# Panel A: Metabolic Trajectories
raw_data$Tissue <- ifelse(raw_data$Tissue.type == "L", "Leaf", "Root")
metabolic_trajectories <- raw_data %>%
  group_by(Genotype, Tissue, Day) %>%
  summarise(
    mean_value = mean(Metabolite_Value, na.rm = TRUE),
    se_value = sd(Metabolite_Value, na.rm = TRUE) / sqrt(n()),
    n = n(),
    ci95 = qt(0.975, df = n - 1) * se_value,
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = Day, y = mean_value, color = Genotype, 
             group = interaction(Genotype, Tissue))) +
  geom_ribbon(aes(ymin = mean_value - ci95, 
                  ymax = mean_value + ci95,
                  fill = Genotype), alpha = 0.2) +
  geom_line(size = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
  scale_fill_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
  facet_wrap(~Tissue, scales = "free_y") +
  labs(title = "Molecular Feature Trajectories",
       x = "Time (Days)",
       y = "Molecular Response (a.u.)") +
  scale_x_continuous(breaks = 1:3) +
  nature_theme

# Panel B: Temporal Coordination
temporal_corr_data <- data.frame(
  Day = c(1, 2, 3, 1, 2, 3),
  Genotype = rep(c("G1", "G2"), each = 3),
  Correlation = c(0.546, 0.448, 0.350,  # G1 values from manuscript
                  0.236, 0.262, 0.288),  # G2 values from manuscript
  CI_Width = c(0.05, 0.05, 0.05,        # Estimated CI widths
               0.04, 0.04, 0.04)        # Adjust these based on actual data
)

temporal_coordination <- ggplot(temporal_corr_data) +
  geom_line(aes(x = Day, y = Correlation, color = Genotype, 
                group = Genotype), size = 1) +
  geom_point(aes(x = Day, y = Correlation, color = Genotype), 
             size = 3) +
  geom_errorbar(aes(x = Day, 
                    ymin = Correlation - CI_Width, 
                    ymax = Correlation + CI_Width,
                    color = Genotype),
                width = 0.1, alpha = 0.5) +
  scale_color_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
  scale_x_continuous(breaks = 1:3) +
  scale_y_continuous(limits = c(0, 0.6),
                     breaks = seq(0, 0.6, by = 0.1)) +
  labs(title = "Temporal Coordination",
       x = "Time (Days)",
       y = "Cross-tissue Correlation (ρ)") +
  nature_theme

# Panel C: Tissue Coordination
tissue_coordination <- raw_data %>%
  ggplot(aes(x = Tissue, y = Metabolite_Value, fill = Genotype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2),
              aes(color = Genotype),
              size = 1,
              alpha = 0.4) +
  scale_fill_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
  scale_color_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
  labs(title = "Tissue Coordination",
       y = "Molecular Response",
       x = "Tissue") +
  nature_theme

# Panel D: Response Metrics
response_metrics <- raw_data %>%
  group_by(Tissue, Genotype) %>%
  summarise(
    mean_response = mean(abs(Metabolite_Value), na.rm = TRUE),
    sd_response = sd(abs(Metabolite_Value), na.rm = TRUE),
    n = n(),
    ci95 = qt(0.975, df = n - 1) * (sd_response/sqrt(n)),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = Genotype, y = mean_response, fill = Tissue)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_response - ci95, 
                    ymax = mean_response + ci95),
                position = position_dodge(width = 0.9),
                width = 0.25) +
  scale_fill_manual(values = c("Leaf" = "#27ae60", "Root" = "#33d6d3")) +
  labs(title = "Response Metrics",
       y = "Mean Molecular Response",
       x = "Genotype") +
  nature_theme

# Combine plots
combined_plot <- plot_grid(
  metabolic_trajectories, temporal_coordination,
  tissue_coordination, response_metrics,
  ncol = 2,
  align = 'hv',
  axis = 'tblr',
  rel_heights = c(1, 1)
)

# Add title
title <- ggdraw() + 
  draw_label(
    "Figure 2 | Temporal Dynamics & Cross-tissue Coordination in Response to Osmotic Stress",
    fontface = 'bold',
    x = 0,
    hjust = 0,
    size = 14
  ) +
  theme(
    plot.margin = margin(0, 0, 5, 7)
  )

# Create final plot with title
final_plot <- plot_grid(
  title, combined_plot,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

# Save plots in multiple formats
ggsave(file.path(out_dir, "Fig2_temporal_dynamics.pdf"),
       final_plot,
       width = 190,
       height = 160,
       units = "mm",
       dpi = 300)

ggsave(file.path(out_dir, "Fig2_temporal_dynamics.tiff"),
       final_plot,
       width = 190,
       height = 160,
       units = "mm",
       dpi = 300,
       compression = "lzw")

print("Plots saved successfully in PDF and TIFF formats")