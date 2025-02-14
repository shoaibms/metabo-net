



################################################
################################################
########### Module−Level Organisation V2 #######
################################################
################################################
# Load required libraries
library(igraph)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Enable debug mode
DEBUG <- TRUE

debug_print <- function(message, data = NULL) {
  if (DEBUG) {
    cat("\nDEBUG:", message, "\n")
    if (!is.null(data)) {
      print(head(data))
      cat("Data dimensions:", dim(data), "\n")
      if (is.data.frame(data)) {
        cat("Column names:", paste(colnames(data), collapse = ", "), "\n")
      }
    }
  }
}

# File paths
base_dir <- "C:/Users/ms/Desktop/r/chem_data/final/baysian_new_crosstalk5_V5_5000"
output_dir <- file.path(base_dir, "module_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read network results
network_results <- readRDS(file.path(base_dir, "network_results.rds"))
debug_print("Network results loaded", network_results)

# Analyze modules
analyze_modules <- function(network_data, tissue_type, genotype) {
  debug_print(sprintf("Analyzing modules for %s %s", tissue_type, genotype))
  
  adj_matrix <- network_data$edge_prob > 0.5
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
  
  louvain_comm <- cluster_louvain(g)
  module_sizes <- sizes(louvain_comm)
  
  module_assignments <- data.frame(
    node = V(g)$name,
    module_id = as.numeric(membership(louvain_comm)),
    module_size = as.numeric(module_sizes[membership(louvain_comm)]),
    tissue = tissue_type,
    genotype = genotype
  )
  
  debug_print("Module assignments created", module_assignments)
  
  return(module_assignments)
}

# Run analysis and collect module assignments
leaf_g1 <- analyze_modules(network_results$leaf_network, "Leaf", "G1")
leaf_g2 <- analyze_modules(network_results$leaf_network, "Leaf", "G2")
root_g1 <- analyze_modules(network_results$root_network, "Root", "G1")
root_g2 <- analyze_modules(network_results$root_network, "Root", "G2")

all_assignments <- bind_rows(leaf_g1, leaf_g2, root_g1, root_g2)

# Function to save plots
save_publication_plots <- function(plot, filename_base, output_dir, width = 8, height = 6) {
  ggsave(file.path(output_dir, paste0(filename_base, ".pdf")), plot, width = width, height = height, dpi = 300)
  ggsave(file.path(output_dir, paste0(filename_base, ".tiff")), plot, width = width, height = height, dpi = 300, compression = "lzw")
  ggsave(file.path(output_dir, paste0(filename_base, ".png")), plot, width = width, height = height, dpi = 600)
}

# Option to manually set colors
custom_colors <- NULL  # Replace NULL with your custom color vector, e.g., c("#FF5733", "#33FF57", "#3357FF")

# Create module organization plot
if (!is.null(all_assignments) && nrow(all_assignments) > 0) {
  # Use custom colors if provided; otherwise, use a default gradient palette
  module_colors <- if (!is.null(custom_colors)) {
    custom_colors
  } else {
    colorRampPalette(c("#4F6D7A", "#4A86B4", "#5D9C59", "#78BE8A", "#8FC741", "#A9C66E", "#BCCB56", "#D4D156", "#E6E04B", "#F0F06F"))(length(unique(all_assignments$module_id)))
  }
  
  module_plot <- ggplot(all_assignments, 
                        aes(x = as.numeric(module_size), 
                            y = tissue)) +
    geom_point(aes(color = factor(module_id), 
                   size = as.numeric(module_size)), 
               alpha = 0.8) +
    facet_wrap(~genotype) +
    # Apply the color palette
    scale_color_manual(values = module_colors) +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() +
    theme(
      text = element_text(size = 18),  # Increased font size
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 16, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      legend.key.size = unit(1.5, "lines"),  # Increase legend dot size
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      strip.text = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 20, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(
      title = "Module-Level Organisation in Root & Leaf Networks",
      x = "Module Size",
      y = "Tissue Type",
      color = "Module ID",
      size = "Module Size"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))  # Adjust the size here
  
  # Save the updated module plot
  save_publication_plots(module_plot, "module_organization", output_dir, width = 8, height = 6)
}

# Summary Statistics
cat("\nSummary Statistics:\n")
if (!is.null(all_assignments)) {
  print(summary(all_assignments))
}

