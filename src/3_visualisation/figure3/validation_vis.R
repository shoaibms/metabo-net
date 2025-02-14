# Load required libraries
library(tidyverse)
library(igraph)
library(boot)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(viridis)

# Define file paths - Added explicit path definitions
base_dir <- "C:/Users/ms/Desktop/r/chem_data/final"
network_results_path <- file.path(base_dir, "baysian_new_crosstalk5_V5_5000/temporal_analysis.rds")
output_dir <- file.path(base_dir, "validation_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Nature theme definition
nature_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.text = element_text(size = 12, face = "bold")
  )

calculate_network_validation_metrics <- function(network_results) {
  # Module stats calculation
  root_edge_prob <- network_results$root_network$edge_prob
  leaf_edge_prob <- network_results$leaf_network$edge_prob
  
  calculate_preservation <- function(edge_prob, genotype) {
    threshold <- if(genotype == "G1") 0.5 else 0.4
    adj_matrix <- edge_prob > threshold
    graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    density <- edge_density(graph)
    transitivity <- transitivity(graph)
    modularity <- modularity(cluster_fast_greedy(graph))
    preservation_score <- (0.4 * density + 0.3 * transitivity + 0.3 * modularity)
    return(preservation_score)
  }
  
  module_stats <- data.frame(
    Tissue = rep(c("Leaf", "Root"), each = 2),
    Genotype = rep(c("G1", "G2"), 2)
  ) %>%
    rowwise() %>%
    mutate(
      Preservation = case_when(
        Tissue == "Leaf" & Genotype == "G1" ~ calculate_preservation(leaf_edge_prob, "G1"),
        Tissue == "Leaf" & Genotype == "G2" ~ calculate_preservation(leaf_edge_prob, "G2"),
        Tissue == "Root" & Genotype == "G1" ~ calculate_preservation(root_edge_prob, "G1"),
        Tissue == "Root" & Genotype == "G2" ~ calculate_preservation(root_edge_prob, "G2")
      ),
      Integration = case_when(
        Tissue == "Leaf" & Genotype == "G1" ~ mean(network_results$leaf_network$null_distribution$effect_sizes[1:250]),
        Tissue == "Leaf" & Genotype == "G2" ~ mean(network_results$leaf_network$null_distribution$effect_sizes[251:500]),
        Tissue == "Root" & Genotype == "G1" ~ mean(network_results$root_network$null_distribution$effect_sizes[1:250]),
        Tissue == "Root" & Genotype == "G2" ~ mean(network_results$root_network$null_distribution$effect_sizes[251:500])
      ),
      Robustness = case_when(
        Tissue == "Leaf" & Genotype == "G1" ~ 1 - mean(network_results$leaf_network$null_distribution$p_values[1:250]),
        Tissue == "Leaf" & Genotype == "G2" ~ 1 - mean(network_results$leaf_network$null_distribution$p_values[251:500]),
        Tissue == "Root" & Genotype == "G1" ~ 1 - mean(network_results$root_network$null_distribution$p_values[1:250]),
        Tissue == "Root" & Genotype == "G2" ~ 1 - mean(network_results$root_network$null_distribution$p_values[251:500])
      )
    ) %>%
    ungroup()
  
  # Stability metrics with correct indexing
  stability_metrics <- expand.grid(
    Time = 1:3,
    Tissue = c("Leaf", "Root"),
    Genotype = c("G1", "G2")
  ) %>%
    as.data.frame() %>%
    arrange(Tissue, Genotype, Time)
  
  # Assign stability scores with proper indexing
  stability_metrics$Stability <- c(
    network_results$permutation_results$leaf_null$statistics[1:3],  # Leaf G1
    network_results$permutation_results$leaf_null$statistics[4:6],  # Leaf G2
    network_results$permutation_results$root_null$statistics[1:3],  # Root G1
    network_results$permutation_results$root_null$statistics[4:6]   # Root G2
  )
  
  stability_metrics$Hub_Persistence <- c(
    1 - network_results$permutation_results$leaf_null$p_values[1:3],  # Leaf G1
    1 - network_results$permutation_results$leaf_null$p_values[4:6],  # Leaf G2
    1 - network_results$permutation_results$root_null$p_values[1:3],  # Root G1
    1 - network_results$permutation_results$root_null$p_values[4:6]   # Root G2
  )
  
  # Statistical validation
  validation_stats <- data.frame(
    Tissue = rep(c("Leaf", "Root"), each = 100),
    Effect_Size = c(
      network_results$leaf_network$null_distribution$effect_sizes[1:100],
      network_results$root_network$null_distribution$effect_sizes[1:100]
    ),
    Permutation_Score = c(
      1 - network_results$leaf_network$null_distribution$p_values[1:100],
      1 - network_results$root_network$null_distribution$p_values[1:100]
    )
  )
  
  return(list(
    module_stats = module_stats,
    stability_metrics = stability_metrics,
    validation_stats = validation_stats
  ))
}


# Add this function after calculate_network_validation_metrics
export_validation_data <- function(validation_metrics, output_dir) {
  # Panel A: Module-Level Validation Data
  write.csv(
    validation_metrics$module_stats %>%
      select(
        Tissue,
        Genotype,
        Module_Preservation = Preservation,
        Integration_Score = Integration,
        Robustness
      ),
    file = file.path(output_dir, "module_level_validation.csv"),
    row.names = FALSE
  )
  
  # Panel B: Network Evolution Data
  write.csv(
    validation_metrics$stability_metrics %>%
      mutate(
        Tissue_Genotype = paste(Tissue, Genotype, sep="."),
        Network_Stability = Stability
      ) %>%
      select(
        Time,
        Tissue,
        Genotype,
        Tissue_Genotype,
        Network_Stability,
        Hub_Persistence
      ),
    file = file.path(output_dir, "network_evolution.csv"),
    row.names = FALSE
  )
  
  # Panel C: Statistical Cross-Validation Data
  write.csv(
    validation_metrics$validation_stats %>%
      select(
        Tissue,
        Effect_Size,
        Permutation_Score
      ),
    file = file.path(output_dir, "statistical_cross_validation.csv"),
    row.names = FALSE
  )
}







# Panel A: Module-Level Validation visualization
create_module_validation_plot <- function(validation_metrics) {
  ggplot(validation_metrics$module_stats, 
         aes(x = Preservation, y = Integration, size = Robustness)) +
    geom_point(aes(color = Tissue, shape = Genotype), alpha = 0.8) +
    scale_color_manual(values = c("Leaf" = "#2ecc71", "Root" = "#33d6d3")) +
    scale_size_continuous(range = c(4, 10)) +
    nature_theme +
    labs(title = "a Module-Level Network Validation",
         x = "Module Preservation Score",
         y = "Integration Score")
}

# Panel B: Network Evolution Analysis visualization
create_network_evolution_plot <- function(validation_metrics) {
  ggplot(validation_metrics$stability_metrics, 
         aes(x = Time, y = Stability, color = interaction(Tissue, Genotype))) +
    geom_line(linewidth = 1) +  # Updated from size to linewidth
    geom_point(aes(size = Hub_Persistence)) +
    scale_color_manual(values = c("Leaf.G1" = "#2ecc71", "Leaf.G2" = "#27ae60",
                                  "Root.G1" = "#33d6d3", "Root.G2" = "#1e597d")) +
    nature_theme +
    labs(title = "b Network Evolution Analysis",
         x = "Time Point",
         y = "Network Stability Score")
}

# Panel C: Statistical Cross-Validation visualization
create_cross_validation_plot <- function(validation_metrics) {
  ggplot(validation_metrics$validation_stats, 
         aes(x = Effect_Size, y = Permutation_Score, color = Tissue)) +
    geom_point(alpha = 0.6) +
    geom_density2d() +
    scale_color_manual(values = c("Leaf" = "#2ecc71", "Root" = "#33d6d3")) +
    nature_theme +
    labs(title = "c Statistical Cross-Validation",
         x = "Effect Size",
         y = "Permutation Score")
}

# Main execution function
main <- function() {
  # Load network results
  network_results <- readRDS(network_results_path)
  
  # Calculate validation metrics
  validation_metrics <- calculate_network_validation_metrics(network_results)
  
  # Export data to CSV files
  export_validation_data(validation_metrics, output_dir)
  
  # Generate plots (existing code)
  p1 <- create_module_validation_plot(validation_metrics)
  p2 <- create_network_evolution_plot(validation_metrics)
  p3 <- create_cross_validation_plot(validation_metrics)
  
  # Combine plots with proper layout
  combined_plot <- plot_grid(
    p1, p2, p3,
    ncol = 3,
    align = 'h',
    axis = 'tblr',
    rel_widths = c(1, 1, 1)
  )
  
  # Add publication-quality title
  title <- ggdraw() + 
    draw_label(
      "Figure 3 | Multi-level Network Validation Framework",
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(plot.margin = margin(0, 0, 10, 7))
  
  # Create final composite figure
  final_plot <- plot_grid(
    title, combined_plot,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )
  
  # Save high-resolution outputs
  ggsave(file.path(output_dir, "Fig3_network_validation.pdf"),
         final_plot, width = 15, height = 5, dpi = 300)
  ggsave(file.path(output_dir, "Fig3_network_validation.png"),
         final_plot, width = 15, height = 5, dpi = 300)
  
  # Save validation metrics for reproducibility
  write.csv(
    data.frame(
      Metric = c("Module Preservation", "Network Stability", "Statistical Validation"),
      Score = c(
        mean(validation_metrics$module_stats$Preservation),
        mean(validation_metrics$stability_metrics$Stability),
        mean(validation_metrics$validation_stats$Effect_Size)
      )
    ),
    file = file.path(output_dir, "validation_metrics.csv"),
    row.names = FALSE
  )
}

# Execute analysis pipeline
main()

