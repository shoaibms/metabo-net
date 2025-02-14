


##############################
########## Hub Distribution and Hub Rank Decay Plot V2 # will use this version
#############################



library(ggplot2)
library(dplyr)

# Read files
input_dir <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E"
output_dir <- file.path(input_dir, "plots")
dir.create(output_dir, showWarnings = FALSE)

files <- list(
  leaf_g1 = read.csv(file.path(input_dir, "leaf_g1_hub_metabolites.csv")),
  leaf_g2 = read.csv(file.path(input_dir, "leaf_g2_hub_metabolites.csv")), 
  root_g1 = read.csv(file.path(input_dir, "root_g1_hub_metabolites.csv")),
  root_g2 = read.csv(file.path(input_dir, "root_g2_hub_metabolites.csv"))
)

# Prepare combined dataset
combined_df <- bind_rows(
  mutate(files$leaf_g1, Tissue = "Leaf", Genotype = "G1"),
  mutate(files$leaf_g2, Tissue = "Leaf", Genotype = "G2"),
  mutate(files$root_g1, Tissue = "Root", Genotype = "G1"),
  mutate(files$root_g2, Tissue = "Root", Genotype = "G2")
)

# Hub distribution plot
hub_dist <- ggplot(combined_df, aes(x = Tissue, y = Degree, fill = Genotype)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.7) +
  # Add individual points for better data transparency
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
  # Add statistical annotations
  stat_summary(fun.data = "mean_sdl", geom = "pointrange", position = position_dodge(0.9)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),  # Remove minor gridlines for clarity
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  ) +
  labs(title = "Hub Molecular Features  Distribution",
       subtitle = "Showing degree centrality across tissues and genotypes",
       y = "Degree Centrality", 
       x = "Tissue")

# Hub rank decay plot 
rank_decay <- ggplot(combined_df, aes(x = 1:nrow(combined_df), y = Degree, 
                                      color = interaction(Tissue, Genotype))) +
  geom_line(size = 1) +
  # Add confidence bands
  geom_smooth(se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("Leaf.G1" = "#2ecc71", 
                                "Leaf.G2" = "#27ae60",
                                "Root.G1" = "#33d6d3", 
                                "Root.G2" = "#1e597d")) +
  # Add break points for better readability
  scale_x_continuous(breaks = seq(0, nrow(combined_df), by = 200)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "right"
  ) +
  labs(title = "Hub Connectivity Decay",
       subtitle = "Showing degree distribution by rank order",
       x = "Hub Rank (ordered by degree)", 
       y = "Degree Centrality",
       color = "Tissue-Genotype")

# Save plots function
# Save plots function with adjusted width and height
save_plot <- function(plot, name, output_dir, width, height) {
  formats <- c("pdf", "png")
  
  for (format in formats) {
    filename <- file.path(output_dir, paste0(name, ".", format))
    ggsave(
      filename, 
      plot = plot, 
      width = width, 
      height = height, 
      dpi = 300, 
      device = format
    )
  }
}

# Adjust text sizes in the theme
hub_dist <- hub_dist +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

rank_decay <- rank_decay +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

# Save plots in reduced size with larger text
save_plot(hub_dist, "hub_distribution", output_dir, 6, 4)
save_plot(rank_decay, "hub_rank_decay", output_dir, 6, 4)





