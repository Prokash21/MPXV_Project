# Example data matrix, saved as 'gene_signature_data.csv'
write.csv(data.frame(
  sample_id = c("sample_1", "sample_2", "sample_3", "sample_4", "sample_5", "sample_6"),
  gene1 = c(2.3, 2.1, 2.4, 4.3, 4.1, 4.4),
  gene2 = c(1.2, 1.5, 1.1, 3.2, 3.5, 3.1),
  gene3 = c(3.4, 3.2, 3.5, 5.4, 5.2, 5.5),
  gene4 = c(2.5, 2.7, 2.6, 4.5, 4.7, 4.6),
  gene5 = c(1.9, 1.8, 1.7, 3.9, 3.8, 3.7),
  condition = c("control", "control", "control", "disease", "disease", "disease")
), "gene_signature_data.csv", row.names = FALSE)

# Load the data
data <- read.csv("gene_signature_data.csv")

# Duplicate the sample data
augmented_data <- rbind(data, data)

# Convert 'condition' to a factor
augmented_data$condition <- as.factor(augmented_data$condition)

# Separate features and target variable
features <- augmented_data[, -ncol(augmented_data)]
target <- augmented_data$condition


# Define parameter grid for RF
rf_grid <- expand.grid(
  mtry = c(2, 3, 4),
  ntree = c(100, 200, 500),
  nodesize = c(1, 5, 10)
)

# Function to perform grid search for RF
rf_grid_search <- function(features, target, grid) {
  best_model <- NULL
  best_accuracy <- 0
  
  for (i in 1:nrow(grid)) {
    model <- randomForest(
      x = features,
      y = target,
      mtry = grid$mtry[i],
      ntree = grid$ntree[i],
      nodesize = grid$nodesize[i]
    )
    
    predictions <- predict(model, features)
    accuracy <- mean(predictions == target)
    
    if (accuracy > best_accuracy) {
      best_model <- model
      best_accuracy <- accuracy
    }
  }
  return(best_model)
}

library("randomForest")
# Run grid search
best_rf_model <- rf_grid_search(features, target, rf_grid)

###################################################################################

BiocManager::install("Rtsne")
library("Rtsne")

# t-SNE
tsne_results <- Rtsne(as.matrix(features), dims = 2, perplexity = 3)

BiocManager::install("umap")
library("umap")
# UMAP
umap_results <- umap(as.matrix(features), n_neighbors = 15, min_dist = 0.1)

# Visualization of t-SNE results
plot(tsne_results$Y, col = as.numeric(target), main = "t-SNE Clustering", pch = 19)

# Visualization of UMAP results
plot(umap_results$layout, col = as.numeric(target), main = "UMAP Clustering", pch = 19)

# Feature importance from the best RF model
rf_importance <- importance(best_rf_model)
top_features <- names(sort(rf_importance, decreasing = TRUE)[1:5])  # Select top 5 features for simplicity

# t-SNE on top features
tsne_top_features <- Rtsne(as.matrix(features[, top_features]), dims = 2, perplexity = 1)
plot(tsne_top_features$Y, col = as.numeric(target), main = "t-SNE on Top Features", pch = 19)

# UMAP on top features
umap_top_features <- umap(as.matrix(features[, top_features]), n_neighbors = 2, min_dist = 0.1)
plot(umap_top_features$layout, col = as.numeric(target), main = "UMAP on Top Features", pch = 19)
