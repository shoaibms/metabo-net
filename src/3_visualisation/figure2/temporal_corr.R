# Load required libraries
library(tidyverse)
library(cowplot)
library(boot)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(grid)


# Define paths
leaf_data_path <- "C:/Users/ms/Desktop/data_chem_3_10/data/data/n_p_l.csv"
root_data_path <- "C:/Users/ms/Desktop/data_chem_3_10/data/data/n_p_r.csv"
vip_path <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/vip_bonferroni/VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv"
metrics_path <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/network_metrics_summary.csv"
output_dir <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/temporal_analysis"

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
leaf_data <- read.csv(leaf_data_path)
root_data <- read.csv(root_data_path)
vip_data <- read.csv(vip_path)
metrics_data <- read.csv(metrics_path)

# Nature style theme
nature_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 20, color = "gray30"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.text = element_text(size = 18),
    legend.position = "right"
  )

# Function for temporal correlation calculation
calculate_temporal_correlation <- function(x, y, n_bootstrap = 10) {
  correlation <- cor(x, y, method = "spearman", use = "pairwise.complete.obs")
  
  boot_cor <- function(data, indices) {
    x_boot <- data$x[indices]
    y_boot <- data$y[indices]
    return(cor(x_boot, y_boot, method = "spearman", use = "pairwise.complete.obs"))
  }
  
  boot_data <- data.frame(x = x, y = y)
  boot_results <- boot(boot_data, boot_cor, R = n_bootstrap)
  ci <- boot.ci(boot_results, type = "perc")
  
  return(list(
    correlation = correlation,
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5]
  ))
}

# Calculate temporal metrics
calculate_temporal_metrics <- function(leaf_data, root_data) {
  results <- list()
  
  get_molecular_features <- function(data) {
    molecular_cols <- grep("^(N_Cluster_|P_Cluster_)", names(data), value = TRUE)
    return(molecular_cols)
  }
  
  leaf_features <- get_molecular_features(leaf_data)
  root_features <- get_molecular_features(root_data)
  common_features <- intersect(leaf_features, root_features)
  
  for(genotype in c("G1", "G2")) {
    for(day in unique(leaf_data$Day)) {
      leaf_subset <- leaf_data %>%
        filter(Genotype == genotype, Day == day) %>%
        select(all_of(common_features)) %>%
        as.matrix() %>%
        as.vector()
      
      root_subset <- root_data %>%
        filter(Genotype == genotype, Day == day) %>%
        select(all_of(common_features)) %>%
        as.matrix() %>%
        as.vector()
      
      if(length(leaf_subset) == length(root_subset)) {
        cor_results <- calculate_temporal_correlation(leaf_subset, root_subset)
        
        results[[paste(genotype, day)]] <- c(
          Genotype = genotype,
          Day = day,
          Correlation = cor_results$correlation,
          CI_Lower = cor_results$ci_lower,
          CI_Upper = cor_results$ci_upper
        )
      }
    }
  }
  
  return(bind_rows(results))
}

# Calculate response magnitudes
calculate_response_magnitudes <- function(leaf_data, root_data) {
  leaf_molecular_cols <- grep("^(N_Cluster_|P_Cluster_)", names(leaf_data), value = TRUE)
  root_molecular_cols <- grep("^(N_Cluster_|P_Cluster_)", names(root_data), value = TRUE)
  
  leaf_magnitude <- leaf_data %>%
    select(Genotype, Day, any_of(leaf_molecular_cols)) %>%
    pivot_longer(cols = any_of(leaf_molecular_cols),
                 names_to = "Feature",
                 values_to = "Value") %>%
    group_by(Genotype, Day) %>%
    summarise(
      mean_response = mean(abs(Value), na.rm = TRUE),
      se_response = sd(abs(Value), na.rm = TRUE)/sqrt(n()),
      .groups = 'drop'
    ) %>%
    mutate(Tissue = "Leaf")
  
  root_magnitude <- root_data %>%
    select(Genotype, Day, any_of(root_molecular_cols)) %>%
    pivot_longer(cols = any_of(root_molecular_cols),
                 names_to = "Feature",
                 values_to = "Value") %>%
    group_by(Genotype, Day) %>%
    summarise(
      mean_response = mean(abs(Value), na.rm = TRUE),
      se_response = sd(abs(Value), na.rm = TRUE)/sqrt(n()),
      .groups = 'drop'
    ) %>%
    mutate(Tissue = "Root")
  
  return(list(leaf = leaf_magnitude, root = root_magnitude))
}

# Panel A: Temporal Correlation Plot
create_temporal_correlation_plot <- function(temporal_metrics) {
  temporal_metrics <- temporal_metrics %>%
    mutate(
      Day = factor(Day),
      Correlation = as.numeric(Correlation),
      CI_Lower = as.numeric(CI_Lower),
      CI_Upper = as.numeric(CI_Upper)
    )
  
  p <- ggplot(temporal_metrics, 
              aes(x = Day, y = Correlation, color = Genotype, group = Genotype)) +
    geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper, fill = Genotype), 
                alpha = 0.2) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
    scale_fill_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
    scale_y_continuous(
      limits = c(0.2, 0.4),
      breaks = seq(0.2, 0.4, by = 0.05)
    ) +
    nature_theme +
    labs(
      title = "Cross-tissue Correlation",
      x = "Time (Days)",
      y = expression(paste("Cross-tissue Correlation (", rho, ")"))
    )
  
  return(p)
}

# Panel B: Response Magnitude Plot
create_magnitude_comparison_plot <- function(response_magnitudes) {
  leaf_magnitude <- response_magnitudes$leaf
  root_magnitude <- response_magnitudes$root
  
  p <- ggplot() +
    geom_line(data = leaf_magnitude, 
              aes(x = Day, y = mean_response, color = Genotype, linetype = "Leaf"),
              size = 1) +
    geom_line(data = root_magnitude,
              aes(x = Day, y = mean_response, color = Genotype, linetype = "Root"),
              size = 1) +
    geom_point(data = leaf_magnitude,
               aes(x = Day, y = mean_response, color = Genotype),
               size = 3) +
    geom_point(data = root_magnitude,
               aes(x = Day, y = mean_response, color = Genotype),
               size = 3) +
    geom_errorbar(data = leaf_magnitude,
                  aes(x = Day, 
                      ymin = mean_response - se_response,
                      ymax = mean_response + se_response,
                      color = Genotype),
                  width = 0.2) +
    geom_errorbar(data = root_magnitude,
                  aes(x = Day, 
                      ymin = mean_response - se_response,
                      ymax = mean_response + se_response,
                      color = Genotype),
                  width = 0.2) +
    scale_color_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
    nature_theme +
    labs(
      title = "Response Magnitude Distribution",
      x = "Time (Days)",
      y = "Response Magnitude",
      linetype = "Tissue",
      color = "Genotype"
    )
  
  return(p)
}

library(grid)  # For custom line ends

# Panel D: Feature Distribution Plot
library(grid)  # For custom line ends

# Panel D: Feature Distribution Plot
create_feature_distribution_plot <- function(leaf_data, root_data) {
  feature_dist <- data.frame(
    Group = factor(c("G1", "G1", "G2", "G2")),
    Tissue = factor(c("Leaf", "Root", "Leaf", "Root")),
    Resilient_Features = c(7.80, 13.00, 16.36, 16.95)
  ) %>%
    group_by(Group) %>%
    mutate(
      Differential = abs(diff(Resilient_Features)),
      Position = ifelse(Group == "G1", 1, 2)
    )
  
  # Calculate the maximum y-value for positioning the lines and labels
  max_value <- max(feature_dist$Resilient_Features)
  
  # Create the base plot
  p <- ggplot() +
    geom_bar(data = feature_dist,
             aes(x = Tissue, y = Resilient_Features, fill = Group),
             stat = "identity",
             position = position_dodge(0.8),
             width = 0.7) +
    
    # Add black lines with round heads for differential
    geom_segment(data = unique(feature_dist[c("Group", "Differential")]),
                 aes(x = 1, xend = 2,
                     y = max_value * 1.05,
                     yend = max_value * 1.05),
                 color = "black",
                 size = 0.8,
                 lineend = "round") +
    
    # Add black text for differential labels
    geom_text(data = unique(feature_dist[c("Group", "Differential")]),
              aes(x = 1.5, y = max_value * 1.08),
              label = paste0("Δ ", format(round(unique(feature_dist$Differential), 2), nsmall = 2), "%"),
              color = "black",
              size = 5) +
    
    scale_fill_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3"),
                      name = "Genotype") +
    
    facet_grid(~ Group, scales = "free_x", space = "free") +
    
    nature_theme +
    labs(
      title = "Feature Redistribution",
      x = NULL,
      y = "Resilient Features (%)"
    ) +
    theme(
      strip.text = element_blank(),
      panel.spacing = unit(1.5, "lines"),
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      panel.grid.major.x = element_blank(),
      legend.position = "right"
    )
  
  return(p)
}


# Panel D: Feature Distribution Plot
library(tidyverse)
library(grid) # for unit()

# Sample data
feature_dist <- data.frame(
  Group = factor(c("G1", "G1", "G2", "G2")),
  Tissue = factor(c("Leaf", "Root", "Leaf", "Root")),
  Resilient_Features = c(7.80, 13.00, 16.36, 16.95)
) %>%
  group_by(Group) %>%
  mutate(
    Differential = abs(diff(Resilient_Features)),
    Position = ifelse(Group == "G1", 1, 2)
  ) %>%
  ungroup()

# Create a parseable label as a string:
# "Delta~5.20~'%'"
feature_dist <- feature_dist %>%
  group_by(Group) %>%
  mutate(
    Label = {
      val <- format(round(Differential[1], 2), nsmall = 2)
      paste0("Delta~", val, "~'%'")
    }
  ) %>%
  ungroup()

# Extract annotation data frame (one row per group)
df_annot <- distinct(feature_dist, Group, Label)

max_value <- max(feature_dist$Resilient_Features)

p4 <- ggplot(feature_dist, aes(x = Tissue, y = Resilient_Features, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  
  # Draw arrows for each group using the annotation data
  geom_segment(
    data = df_annot,
    aes(
      x = 1,
      xend = 2,
      y = max_value * 1.05,
      yend = max_value * 1.05,
      group = Group
    ),
    arrow = arrow(ends = "both", type = "closed", length = unit(0.1, "cm")),
    color = "black",
    size = 0.5,
    inherit.aes = FALSE
  ) +
  
  # Add text annotations for each group using parse = TRUE
  geom_text(
    data = df_annot,
    aes(
      x = 1.5,
      y = max_value * 1.08,
      label = Label,
      group = Group
    ),
    parse = TRUE, # Interpret Label as a plotmath expression
    color = "black",
    size = 5,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3"), name = "Genotype") +
  
  facet_grid(~ Group, scales = "free_x", space = "free") +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.text = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  ) +
  
  labs(
    title = "Feature Redistribution",
    x = NULL,
    y = "Resilient Features (%)"
  )

p4



# Main execution function
main <- function() {
  temporal_metrics <- calculate_temporal_metrics(leaf_data, root_data)
  response_magnitudes <- calculate_response_magnitudes(leaf_data, root_data)
  
  p1 <- create_temporal_correlation_plot(temporal_metrics)
  p2 <- create_magnitude_comparison_plot(response_magnitudes)
  p3 <- create_temporal_trends_plot()
  p4 <- create_feature_distribution_plot(leaf_data, root_data)
  
  combined_plot <- plot_grid(
    p1, p2, p3, p4,
    ncol = 2,
    align = 'h',
    axis = 'tblr',
    rel_widths = c(1, 1),
    rel_heights = c(1, 1)
  )
  
  title <- ggdraw() + 
    draw_label(
      "Figure 2 | Temporal dynamics and coordination of tissue-specific molecular responses under osmotic stress",
      fontface = 'bold',
      x = 0,
      hjust = 0,
      size = 16
    ) +
    theme(plot.margin = margin(0, 0, 10, 7))
  
  final_plot <- plot_grid(
    title, combined_plot,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )
  
  ggsave(file.path(output_dir, "Fig2_temporal_dynamics.pdf"),
         final_plot, width = 12, height = 10, dpi = 300)
  ggsave(file.path(output_dir, "Fig2_temporal_dynamics.png"),
         final_plot, width = 12, height = 10, dpi = 300)
  
  return(final_plot)
}

# Execute
plot <- main()