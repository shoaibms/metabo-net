
# Load required packages
suppressPackageStartupMessages({
  library(ggraph)
  library(igraph)
  library(plotly)
  library(reshape2)
  library(bnlearn)
})

# Function to check if directory exists and create if needed
create_dir_safe <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

# File paths
filtered_metabolites_file <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/group/merge2/Merged_VIP_hub_name_only.csv"
merged_data_file <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/group/merge2/Merged_VIP_hub.csv"
output_dir <- "C:/Users/ms/Desktop/r/chem_data/final/baysian_root_shoot_plot2"
output_2d_plot <- file.path(output_dir, "root_shoot_interaction_network_2d.png")
output_3d_plot <- file.path(output_dir, "root_shoot_interaction_network_3d.html")

# Create output directory
create_dir_safe(output_dir)

# Check if input files exist
if (!file.exists(filtered_metabolites_file) || !file.exists(merged_data_file)) {
  stop("Input files not found. Please check file paths.")
}

# Load and process data
filtered_metabolites <- read.csv(filtered_metabolites_file)
merged_data <- read.csv(merged_data_file)

# Filter metabolites for root and shoot tissues
root_metabolites <- filtered_metabolites[filtered_metabolites$Tissue.type == "R", ]
leaf_metabolites <- filtered_metabolites[filtered_metabolites$Tissue.type == "L", ]

# Identify crosstalk metabolites
crosstalk_metabolites <- merge(root_metabolites, leaf_metabolites, by = "Metabolite", suffixes = c("_root", "_leaf"))

# Filter and prepare data for network construction
filtered_merged_data <- merged_data[merged_data$Metabolite %in% crosstalk_metabolites$Metabolite, ]
filtered_merged_data$Metabolite <- paste0(substr(filtered_merged_data$Tissue.type, 1, 1), "_", filtered_merged_data$Metabolite)

# Pivot data
pivot_data <- dcast(filtered_merged_data, Vac_id + Genotype + Entry + Batch + Treatment + Replication + Day ~ Metabolite, 
                    value.var = "Metabolite_Value")

# Prepare data for network construction
data_for_bn <- pivot_data[, -(1:7)]
root_data <- data_for_bn[, grepl("^R_", colnames(data_for_bn))]
leaf_data <- data_for_bn[, grepl("^L_", colnames(data_for_bn))]

# Learn network structure
bn_structure_root <- hc(root_data)
bn_structure_leaf <- hc(leaf_data)

# Create network graphs
root_arcs <- arcs(bn_structure_root)
leaf_arcs <- arcs(bn_structure_leaf)

if (nrow(root_arcs) == 0 || nrow(leaf_arcs) == 0) {
  stop("No edges found in network structure")
}

graph_root <- graph_from_data_frame(as.data.frame(root_arcs), directed = TRUE)
graph_leaf <- graph_from_data_frame(as.data.frame(leaf_arcs), directed = TRUE)

# Create combined graph manually
all_vertices <- unique(c(V(graph_root)$name, V(graph_leaf)$name))
combined_graph <- make_empty_graph(n = length(all_vertices), directed = TRUE)
V(combined_graph)$name <- all_vertices

# Add edges from both graphs
root_edges <- as_edgelist(graph_root)
leaf_edges <- as_edgelist(graph_leaf)
all_edges <- rbind(root_edges, leaf_edges)
combined_graph <- add_edges(combined_graph, edges = as.vector(t(all_edges)))

# Add cross-tissue connections
common_metabolites <- intersect(gsub("^R_", "", V(graph_root)$name), 
                                gsub("^L_", "", V(graph_leaf)$name))
cross_edges <- data.frame(
  from = paste0("R_", common_metabolites),
  to = paste0("L_", common_metabolites)
)

# Add cross edges
combined_graph <- add_edges(combined_graph, edges = as.vector(t(cross_edges)))

# Calculate node sizes using a more stable approach
node_degrees <- sapply(V(combined_graph), function(v) {
  length(incident(combined_graph, v, mode="all"))
})

# Normalize scores (1-10 range)
if (length(unique(node_degrees)) > 1) {
  node_scores_norm <- 1 + 9 * (node_degrees - min(node_degrees)) / (max(node_degrees) - min(node_degrees))
} else {
  node_scores_norm <- rep(5, length(node_degrees))  # Default size if all nodes have same degree
}

# Set node sizes
V(combined_graph)$size <- ifelse(V(combined_graph)$name %in% unlist(cross_edges), 
                                 node_scores_norm * 5,
                                 node_scores_norm)

# Set colors
V(combined_graph)$color <- ifelse(grepl("^R_", V(combined_graph)$name), "#1e597d",
                                  ifelse(grepl("^L_", V(combined_graph)$name), "#1b941b", "black"))

# Set edge colors
edge_list <- as_edgelist(combined_graph)
E(combined_graph)$color <- sapply(1:nrow(edge_list), function(i) {
  from <- edge_list[i, 1]
  to <- edge_list[i, 2]
  if (grepl("^R_", from) && grepl("^R_", to)) return("#246b24")      # Root to Root
  if (grepl("^L_", from) && grepl("^L_", to)) return("#1e597d")      # Leaf to Leaf
  return("orange")                                                    # Cross connections
})

# 2D Visualization with error handling
tryCatch({
  g_2d <- ggraph(combined_graph, layout = "fr") +
    geom_edge_link(aes(color = color),
                   arrow = arrow(length = unit(2, 'mm')),
                   end_cap = circle(2, 'mm'),
                   alpha = 0.7,
                   edge_width = 0.5) +
    geom_node_point(aes(size = size, color = I(color)), alpha = 0.7) +
    scale_edge_color_identity() +
    scale_size_continuous(range = c(1, 20)) +
    theme_void() +
    labs(title = "Root-Shoot Interaction Network - Crosstalk Molecular Features") +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 26),
          text = element_text(size = 12)) +
    annotate("text", x = Inf, y = -Inf,
             label = "Blue = Root, Green = Leaf \nEdges: Green = Root-Root, Blue = Leaf-Leaf, Orange = Root-Leaf",
             hjust = 1.1, vjust = -1, size = 10, color = "black")##############7
  
  ggsave(output_2d_plot, g_2d, width = 15, height = 15, dpi = 300)
}, error = function(e) {
  message("Error in 2D plot creation: ", e$message)
})

# 3D Visualization
tryCatch({
  set.seed(42)  # for reproducibility
  layout_3d <- layout_with_fr(combined_graph, dim = 3)
  
  # Prepare node data for 3D
  node_data <- data.frame(
    x = layout_3d[,1],
    y = layout_3d[,2],
    z = layout_3d[,3],
    name = V(combined_graph)$name,
    size = V(combined_graph)$size,
    type = ifelse(grepl("^R_", V(combined_graph)$name), "Root", "Leaf")
  )
  node_data$color <- ifelse(node_data$type == "Root", "#1b941b", "#1e597d")
  
  # Prepare edge data for 3D
  edge_list <- as_edgelist(combined_graph)
  edge_data <- data.frame(
    x = rep(NA, nrow(edge_list) * 3),
    y = rep(NA, nrow(edge_list) * 3),
    z = rep(NA, nrow(edge_list) * 3),
    type = character(nrow(edge_list) * 3)
  )
  
  for (i in 1:nrow(edge_list)) {
    from <- which(V(combined_graph)$name == edge_list[i, 1])
    to <- which(V(combined_graph)$name == edge_list[i, 2])
    idx <- (i-1)*3 + 1
    edge_data[idx:(idx+2), c("x", "y", "z")] <- rbind(
      layout_3d[from, ],
      layout_3d[to, ],
      c(NA, NA, NA)
    )
    edge_data$type[idx:(idx+2)] <- if (grepl("^R_", edge_list[i, 1]) && grepl("^R_", edge_list[i, 2])) {
      "Root-Root"
    } else if (grepl("^L_", edge_list[i, 1]) && grepl("^L_", edge_list[i, 2])) {
      "Leaf-Leaf"
    } else {
      "Root-Leaf"
    }
  }
  
  # Create 3D plot
  p_3d <- plot_ly() %>%
    add_trace(
      data = edge_data,
      x = ~x, y = ~y, z = ~z,
      type = 'scatter3d', mode = 'lines',
      line = list(width = 1),
      color = ~type,
      colors = c("Root-Root" = "#246b24", "Leaf-Leaf" = "#1e597d", "Root-Leaf" = "orange"),
      name = ~type,
      legendgroup = ~type,
      hoverinfo = 'none'
    ) %>%
    add_trace(
      data = node_data,
      x = ~x, y = ~y, z = ~z,
      type = 'scatter3d', mode = 'markers',
      marker = list(size = ~size, color = ~color),
      text = ~name,
      hoverinfo = 'text',
      name = ~type,
      legendgroup = ~type
    ) %>%
    layout(
      scene = list(
        xaxis = list(title = '', showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
        yaxis = list(title = '', showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
        zaxis = list(title = '', showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE)
      ),
      title = "3D Root-Shoot Interaction Network",
      legend = list(title = list(text = 'Node and Edge Types'))
    )
  
  # Save 3D plot
  htmlwidgets::saveWidget(p_3d, output_3d_plot)
}, error = function(e) {
  message("Error in 3D plot creation: ", e$message)
})

print(paste("2D network visualization saved to:", output_2d_plot))
print(paste("3D network visualization saved to:", output_3d_plot))
