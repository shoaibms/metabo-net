################################################
################################################
########### Network Temporal Stability V3 ###### 
################################################
################################################


# Load required libraries
library(tidyverse)
library(igraph)
library(RColorBrewer)

# File paths
base_dir <- "C:/Users/ms/Desktop/r/chem_data/final/baysian_new_crosstalk5_V5_5000"
output_dir <- file.path(base_dir, "temporal_stability")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read network results
network_results <- readRDS(file.path(base_dir, "processed_data.rds"))
pivot_data <- network_results$pivot_data

# Function to calculate stability metrics for a single timepoint
calculate_timepoint_stability <- function(data_subset, tissue, genotype) {
  tryCatch({
    # Create correlation matrix
    cor_matrix <- cor(data_subset, method = "spearman", use = "pairwise.complete.obs")
    
    # Create graph from correlation matrix
    g <- graph_from_adjacency_matrix(
      cor_matrix > 0.7, 
      mode = "undirected", 
      weighted = TRUE
    )
    
    # Calculate metrics
    edge_density_val <- edge_density(g)
    mod_val <- tryCatch(
      modularity(cluster_louvain(g)),
      error = function(e) NA
    )
    hub_val <- mean(sort(degree(g), decreasing = TRUE)[1:min(10, vcount(g))])
    
    metrics <- data.frame(
      tissue = tissue,
      genotype = genotype,
      edge_consistency = edge_density_val,
      module_preservation = mod_val,
      hub_conservation = hub_val
    )
    
    return(metrics)
  }, error = function(e) {
    warning(paste("Error in calculation for", tissue, genotype, ":", e$message))
    return(NULL)
  })
}

# Calculate stability metrics for each timepoint
stability_results <- list()

# Process data by tissue and genotype
tissues <- c("L", "R")
genotypes <- c("G1", "G2")
days <- sort(unique(pivot_data$Day))

for(tissue in tissues) {
  for(genotype in genotypes) {
    # Get metabolite columns for this tissue
    metabolite_cols <- colnames(pivot_data)[
      grep(paste0("^", tissue, "_"), colnames(pivot_data))
    ]
    
    for(day in days) {
      # Subset data for this condition
      subset_data <- pivot_data %>%
        filter(Day == day, Genotype == genotype) %>%
        select(all_of(metabolite_cols))
      
      if(ncol(subset_data) > 0 && nrow(subset_data) > 0) {
        # Calculate stability metrics
        metrics <- calculate_timepoint_stability(subset_data, tissue, genotype)
        if(!is.null(metrics)) {
          metrics$day <- day
          stability_results[[length(stability_results) + 1]] <- metrics
        }
      }
    }
  }
}

# Combine results
stability_df <- bind_rows(stability_results)

# Create group label for plotting
stability_df$group <- paste(stability_df$tissue, stability_df$genotype, sep="_")

# Create integrated stability plot with Nature-style theme

integrated_stability <- ggplot(stability_df, 
                               aes(x = edge_consistency, 
                                   y = module_preservation,
                                   size = hub_conservation,
                                   color = group)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(name = expression(bold("Tissue-Genotype")),
                     values = c("L_G1" = "#34eb81", 
                                "L_G2" = "#1e9651",
                                "R_G1" = "#33d6d3", 
                                "R_G2" = "#1e597d"),
                     labels = c("L_G1" = "Leaf G1",
                                "L_G2" = "Leaf G2",
                                "R_G1" = "Root G1",
                                "R_G2" = "Root G2")) +
  scale_size_continuous(range = c(3, 10)) +
  theme_minimal() +
  theme(
    # Text sizes
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14, color = "black"),
    
    # Legend formatting
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    
    # Grid and panel customization
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    
    # Axis lines
    axis.line = element_line(color = "black", linewidth = 0.5),
    
    # Remove plot background
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    x = "Network Density",
    y = "Module Preservation",
    size = "Hub Conservation"
  )

# Save plot in PDF, TIFF, and PNG formats
ggsave(file.path(output_dir, "integrated_stability.pdf"),
       integrated_stability, width = 7, height = 5)

ggsave(file.path(output_dir, "integrated_stability.tiff"),
       integrated_stability, width = 7, height = 5, dpi = 300, device = "tiff")

ggsave(file.path(output_dir, "integrated_stability.png"),
       integrated_stability, width = 7, height = 5, dpi = 300, device = "png")




