

###########################################
#############Plot 5 #######################New code 
########### effect size  ##################
###########################################


library(ggplot2)
library(dplyr)
library(tidyr)
library(ggridges)
library(patchwork)  # for combining plots

data_dir <- "C:/Users/ms/Desktop/r/chem_data/metabo2/result/section3"
save_dir <- "C:/Users/ms/Desktop/r/chem_data/metabo2/result/section3"

# Read data
effects_data <- read.csv(file.path(data_dir, "effect_sizes.csv"))

# Define theme with larger font sizes
large_font_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "top"
  )

# Effect Size Distribution
p1 <- ggplot(effects_data, aes(x = effect_size, y = Tissue.type, fill = Genotype)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5) +
  geom_vline(xintercept = c(-0.474, -0.33, -0.147, 0, 0.147, 0.33, 0.474), 
             linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_fill_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
  large_font_theme +
  labs(
    x = "Effect Size (Cliff's Delta |d|)",
    y = "Tissue Type",
    title = "Effect Size Distribution"
  ) +
  scale_y_discrete(labels = c("L" = "Leaf", "R" = "Root"))

# Create data for leaf:root ratios
ratio_data <- data.frame(
  Genotype = c("G1", "G1", "G2", "G2"),
  Category = c("≥ 0.474", "0.33-0.474", "≥ 0.474", "0.33-0.474"),
  Ratio = c(10.2, 5.1, 0.4, 1.2)
)

# Leaf:Root Ratios
p2 <- ggplot(ratio_data, aes(x = Category, y = Ratio, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
  large_font_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14)
  ) +
  labs(
    x = "Effect Size Category",
    y = "Leaf:Root Ratio",
    title = "Tissue Response Ratio"
  )

# Combine plots
combined_plot <- p1 + p2 +
  plot_layout(ncol = 2, widths = c(3, 2))

# Save combined plot
ggsave(file.path(save_dir, "effect_size_and_ratio.pdf"), 
       combined_plot, width = 10, height = 4, dpi = 300)
ggsave(file.path(save_dir, "effect_size_and_ratio.png"), 
       combined_plot, width = 10, height = 4, dpi = 300)