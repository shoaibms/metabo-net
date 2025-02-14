


################################################
################################################
########### Module Stability V2 ########new code 
################################################
################################################
# Set seed for reproducibility
set.seed(42)

# File paths
leaf_data_path <- "C:/Users/ms/Desktop/data_chem_3_10/data/data/n_p_l.csv"
root_data_path <- "C:/Users/ms/Desktop/data_chem_3_10/data/data/n_p_r.csv"
vip_path <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/vip_bonferroni/VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv"
output_dir <- "C:/Users/ms/Desktop/r/chem_data/final/baysian_new_crosstalk5_V5_5000/module_stability"

# Set number of permutations
n_permutations <- 10############1000
# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read data
leaf_data <- read.csv(leaf_data_path)
root_data <- read.csv(root_data_path)
vip_data <- read.csv(vip_path)

# Initialize empty data frame with correct structure
stability_metrics <- data.frame(
  Tissue = character(),
  Genotype = character(),
  Module_Preservation = numeric(),
  Module_Coherence = numeric(),
  P_Preservation = numeric(),
  P_Coherence = numeric(),
  stringsAsFactors = FALSE
)

# Function to calculate module metrics
calculate_module_metrics <- function(data, genotype, n_permutations = 10) {###########1000
  set.seed(42)
  
  feature_cols <- grep("^(N_Cluster_|P_Cluster_)", names(data), value = TRUE)
  genotype_data <- data[data$Genotype == genotype, feature_cols]
  
  # Observed metrics
  correlation_matrix <- cor(genotype_data, use = "pairwise.complete.obs")
  observed_preservation <- mean(abs(correlation_matrix), na.rm = TRUE)
  observed_coherence <- mean(correlation_matrix > 0.7, na.rm = TRUE)
  
  # Permutation test
  null_preservation <- numeric(n_permutations)
  null_coherence <- numeric(n_permutations)
  
  for(i in 1:n_permutations) {
    permuted_data <- apply(genotype_data, 2, sample)
    perm_cor <- cor(permuted_data, use = "pairwise.complete.obs")
    
    null_preservation[i] <- mean(abs(perm_cor), na.rm = TRUE)
    null_coherence[i] <- mean(perm_cor > 0.7, na.rm = TRUE)
  }
  
  p_preservation <- mean(null_preservation >= observed_preservation)
  p_coherence <- mean(null_coherence >= observed_coherence)
  
  return(c(
    observed_preservation,
    observed_coherence,
    p_preservation,
    p_coherence
  ))
}

# Calculate metrics
for(genotype in c("G1", "G2")) {
  leaf_metrics <- calculate_module_metrics(leaf_data, genotype, n_permutations)
  root_metrics <- calculate_module_metrics(root_data, genotype, n_permutations)
  
  stability_metrics <- rbind(stability_metrics, data.frame(
    Tissue = "Leaf",
    Genotype = genotype,
    Module_Preservation = leaf_metrics[1],
    Module_Coherence = leaf_metrics[2],
    P_Preservation = leaf_metrics[3],
    P_Coherence = leaf_metrics[4]
  ))
  
  stability_metrics <- rbind(stability_metrics, data.frame(
    Tissue = "Root",
    Genotype = genotype,
    Module_Preservation = root_metrics[1],
    Module_Coherence = root_metrics[2],
    P_Preservation = root_metrics[3],
    P_Coherence = root_metrics[4]
  ))
}

stability_metrics$Group <- paste(stability_metrics$Tissue, stability_metrics$Genotype)

# Create visualization
library(ggplot2)
set.seed(42)

p <- ggplot(stability_metrics, 
            aes(x = Module_Preservation, y = Module_Coherence, color = Group)) +
  geom_point(size = 6, alpha = 0.7) +
  scale_color_manual(name = "Tissue-Genotype",
                     values = c("Leaf G1" = "#41e085", 
                                "Leaf G2" = "#209150",
                                "Root G1" = "#33d6d3", 
                                "Root G2" = "#1e597d")) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 18, ),
    axis.text = element_text(size = 18, color = "black"),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 18),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA)
  ) +
  labs(
    title = "Module Stability Analysis",
    subtitle = paste("Permutation test:", formatC(n_permutations, big.mark=","), "iterations"),
    x = "Module Preservation Score",
    y = "Module Coherence"
  )




# Save outputs
ggsave(file.path(output_dir, "module_stability.pdf"), p, width = 7, height = 5)
ggsave(file.path(output_dir, "module_stability.png"), p, width = 7, height = 5, dpi = 300)
write.csv(stability_metrics, file.path(output_dir, "module_stability_metrics.csv"), row.names = FALSE)