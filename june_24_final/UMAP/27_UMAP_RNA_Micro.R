
setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/UMAP")

#############################
# Load libraries
library(umap)
library(ggplot2)

# Example gene expression data (assuming you are reading from CSV files)
gene_expression_data <- read.csv("27_T_10_RNA_UMAP_Count.csv", row.names = 1)

# Example metadata (assuming you are reading from CSV files)
metadata <- read.csv("27_T_10_RNA_UMAP_Meta.csv")

# Convert Condition to factor
metadata$Condition <- factor(metadata$Condition)

# Perform UMAP dimensionality reduction
umap_result <- umap(gene_expression_data, n_neighbors = 5, min_dist = 0.5)

# Extract UMAP coordinates and combine with metadata
umap_df <- data.frame(
  X1 = umap_result$layout[, 1],  # UMAP component 1
  X2 = umap_result$layout[, 2],  # UMAP component 2
  metadata
)

# Plot using ggplot2
ggplot(umap_df, aes(x = X1, y = X2, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "red"), labels = c("0", "1")) +  # Define colors and labels for conditions
  labs(title = "UMAP_Top_10_RNA",
       x = "UMAP Component 1", y = "UMAP Component 2") +
  theme_minimal()


############################ Micro ###############################

# Load libraries
library(umap)
library(ggplot2)

# Example gene expression data (assuming you are reading from CSV files)
gene_expression_data <- read.csv("27_T_10_Micro_UMAP_Count.csv", row.names = 1)

# Example metadata (assuming you are reading from CSV files)
metadata <- read.csv("27_T_10_Micro_UMAP_Meta_114.csv")

# Convert Condition to factor
metadata$Condition <- factor(metadata$Condition)

#########
any(is.na(gene_expression_data))

gene_expression_data <- gene_expression_data[complete.cases(gene_expression_data), ]


###########3
# Perform UMAP dimensionality reduction
umap_result <- umap(gene_expression_data, n_neighbors = 5, min_dist = 0.5)

# Check for missing values in umap_result$layout (optional step)
if (any(is.na(umap_result$layout))) {
  umap_result$layout <- umap_result$layout[complete.cases(umap_result$layout), ]
}



# Create umap_df by combining UMAP coordinates with metadata
umap_df <- data.frame(
  X1 = umap_result$layout[, 1],  # UMAP component 1
  X2 = umap_result$layout[, 2],  # UMAP component 2
  metadata
)

# Plot using ggplot2
library(ggplot2)
ggplot(umap_df, aes(x = X1, y = X2, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "red"), labels = c("0", "1")) +  # Define colors and labels for conditions
  labs(title = "UMAP_Top_10_Micro",
       x = "UMAP Component 1", y = "UMAP Component 2") +
  theme_minimal()


########


# Perform UMAP dimensionality reduction
umap_result <- umap(gene_expression_data, n_neighbors = 5, min_dist = 0.5)

# Extract UMAP coordinates and combine with metadata
umap_df <- data.frame(
  X1 = umap_result$layout[, 1],  # UMAP component 1
  X2 = umap_result$layout[, 2],  # UMAP component 2
  metadata
)

# Plot using ggplot2
ggplot(umap_df, aes(x = X1, y = X2, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "red"), labels = c("0", "1")) +  # Define colors and labels for conditions
  labs(title = "UMAP_Top_10_Micro",
       x = "UMAP Component 1", y = "UMAP Component 2") +
  theme_minimal()

###################################################################


# Load necessary libraries
library(umap)
library(ggplot2)
library(missForest)  # for imputation

# Example gene expression data (assuming you are reading from CSV files)
gene_expression_data <- read.csv("27_T_10_Micro_UMAP_Count.csv", row.names = 1)

# Example metadata (assuming you are reading from CSV files)
metadata <- read.csv("27_T_10_Micro_UMAP_Meta.csv")

# Convert Condition to factor
metadata$Condition <- factor(metadata$Condition)

# Check for missing values
if (any(is.na(gene_expression_data))) {
  # Impute missing values using missForest
  imputed_data <- missForest(gene_expression_data)
  gene_expression_data <- imputed_data$ximp
}

# Perform UMAP dimensionality reduction
umap_result <- umap(gene_expression_data, n_neighbors = 5, min_dist = 0.5)

# Check for missing values in umap_result$layout (optional step)
if (any(is.na(umap_result$layout))) {
  umap_result$layout <- umap_result$layout[complete.cases(umap_result$layout), ]
}

# Create umap_df by combining UMAP coordinates with metadata
umap_df <- data.frame(
  X1 = umap_result$layout[, 1],  # UMAP component 1
  X2 = umap_result$layout[, 2],  # UMAP component 2
  metadata
)

# Plot using ggplot2
ggplot(umap_df, aes(x = X1, y = X2, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "red"), labels = c("0", "1")) +  # Define colors and labels for conditions
  labs(title = "UMAP_Top_10_Micro",
       x = "UMAP Component 1", y = "UMAP Component 2") +
  theme_minimal()


#******************** seed ***************************#######

# Load necessary libraries
library(umap)
library(ggplot2)
library(missForest)  # for imputation

# Example gene expression data (assuming you are reading from CSV files)
gene_expression_data <- read.csv("27_T_10_Micro_UMAP_Count.csv", row.names = 1)

# Example metadata (assuming you are reading from CSV files)
metadata <- read.csv("27_T_10_Micro_UMAP_Meta.csv")

# Convert Condition to factor
metadata$Condition <- factor(metadata$Condition)

# Check for missing values
if (any(is.na(gene_expression_data))) {
  # Impute missing values using missForest
  imputed_data <- missForest(gene_expression_data)
  gene_expression_data <- imputed_data$ximp
}

# Set random seed for reproducibility
set.seed(123)  # You can choose any number

# Perform UMAP dimensionality reduction
umap_result <- umap(gene_expression_data, n_neighbors = 5, min_dist = 0.5)

# Check for missing values in umap_result$layout (optional step)
if (any(is.na(umap_result$layout))) {
  umap_result$layout <- umap_result$layout[complete.cases(umap_result$layout), ]
}

# Create umap_df by combining UMAP coordinates with metadata
umap_df <- data.frame(
  X1 = umap_result$layout[, 1],  # UMAP component 1
  X2 = umap_result$layout[, 2],  # UMAP component 2
  metadata
)

# Plot using ggplot2
ggplot(umap_df, aes(x = X1, y = X2, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "red"), labels = c("0", "1")) +  # Define colors and labels for conditions
  labs(title = "UMAP_Top_10_Micro",
       x = "UMAP Component 1", y = "UMAP Component 2") +
  theme_minimal()

