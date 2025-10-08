setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/UMAP")

#############################
# Load libraries
library(umap)
library(ggplot2)

# Example gene expression data (assuming you are reading from CSV files)
gene_expression_data <- read.csv("13_T_10_RNA_UMAP_Count.csv", row.names = 1)

# Example metadata (assuming you are reading from CSV files)
metadata <- read.csv("13_RNA_Meta.csv")

# Convert Condition to factor
metadata$group <- factor(metadata$group)

# Perform UMAP dimensionality reduction
umap_result <- umap(gene_expression_data, n_neighbors = 5, min_dist = 0.5)

# Extract UMAP coordinates and combine with metadata
umap_df <- data.frame(
  X1 = umap_result$layout[, 1],  # UMAP component 1
  X2 = umap_result$layout[, 2],  # UMAP component 2
  metadata
)

# Plot using ggplot2
ggplot(umap_df, aes(x = X1, y = X2, color = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "red"), labels = c("0", "1")) +  # Define colors and labels for conditions
  labs(title = "UMAP_Top_10_RNA",
       x = "UMAP Component 1", y = "UMAP Component 2") +
  theme_minimal()



############################################# seed ###################################

# Load necessary libraries
library(umap)
library(ggplot2)

# Set random seed for reproducibility
set.seed(123)

# Example gene expression data (assuming you are reading from CSV files)
gene_expression_data <- read.csv("13_T_10_RNA_UMAP_Count.csv", row.names = 1)

# Example metadata (assuming you are reading from CSV files)
metadata <- read.csv("13_RNA_Meta.csv")

# Convert Condition to factor
metadata$group <- factor(metadata$group)

# Perform UMAP dimensionality reduction
umap_result <- umap(gene_expression_data, n_neighbors = 5, min_dist = 0.5)

# Extract UMAP coordinates and combine with metadata
umap_df <- data.frame(
  X1 = umap_result$layout[, 1],  # UMAP component 1
  X2 = umap_result$layout[, 2],  # UMAP component 2
  metadata
)

# Perform k-means clustering on the UMAP coordinates
kmeans_result <- kmeans(umap_df[, c("X1", "X2")], centers = 2)  # Assuming 2 clusters for simplicity

# Add the cluster assignments to the umap_df dataframe
umap_df$cluster <- as.factor(kmeans_result$cluster)

# Plot using ggplot2
ggplot(umap_df, aes(x = X1, y = X2, color = group, shape = cluster)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "red"), labels = c("0", "1")) +  # Define colors and labels for conditions
  labs(title = "UMAP_Top_10_RNA",
       x = "UMAP Component 1", y = "UMAP Component 2") +
  theme_minimal()
