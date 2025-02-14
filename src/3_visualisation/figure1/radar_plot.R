


##############################
########## Radar plot 
#############################
# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggradar)
library(scales)
library(gridExtra)

# Set paths
input_path <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/network_metrics_summary-Copy.csv"
output_dir <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/plots"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE)

# Read the data
network_metrics <- read.csv(input_path)

# Prepare data for radar plot with full labels
radar_data <- network_metrics %>%
  mutate(group = paste(Tissue.type, Genotype)) %>%
  select(group, density, transitivity, modularity, 
         temporal_coherence, path_length, directionality)

# Scale the data between 0 and 1
radar_data_scaled <- radar_data %>%
  mutate(across(-group, function(x) (x - min(x)) / (max(x) - min(x))))

# Format data for ggradar with full labels
radar_data_formatted <- radar_data_scaled %>%
  rename(
    "Density" = density,
    "Transitivity" = transitivity,
    "Modularity" = modularity,
    "Temporal Coherence" = temporal_coherence,
    "Path Length" = path_length,
    "Directionality" = directionality
  )

# Create radar plot with optimized settings
radar_plot <- ggradar(
  radar_data_formatted,
  grid.min = 0,
  grid.mid = 0.5,
  grid.max = 1,
  values.radar = c("0", "0.5", "1"),
  group.point.size = 3,
  group.line.width = 1,
  group.colours = c("#2ecc71", "#27ae60", "#33d6d3", "#1e597d"),
  background.circle.colour = "white",
  gridline.min.colour = "grey90",
  gridline.mid.colour = "grey85",
  gridline.max.colour = "grey80",
  legend.position = "top",
  legend.title = "Tissue",
  plot.title = "Network Metrics Comparison",
  fill = TRUE,
  fill.alpha = 0.25,
  axis.label.size = 7,     # Large label size
  grid.label.size = 8,     # Grid value size
  plot.extent.x.sf = 1.2,  # Extend plot width to accommodate labels
  plot.extent.y.sf = 1.2  # Extend plot height to accommodate labels
) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "top",
    legend.justification = "center",
    plot.margin = margin(5, 10, 5, 10),  # Increased margins to prevent cutoff
    legend.spacing.x = unit(0.3, 'cm'),
    legend.margin = margin(0, 0, 10, 0)
  )

# Save the plot with larger dimensions
ggsave(file.path(output_dir, "network_metrics_radar.pdf"), 
       radar_plot, 
       width = 7,     # Increased width
       height = 6,    # Increased height
       dpi = 300)

ggsave(file.path(output_dir, "network_metrics_radar.png"), 
       radar_plot, 
       width = 7,     # Increased width
       height = 6,    # Increased height
       dpi = 300)

