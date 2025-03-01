
########################################################################################
# Install required packages if not already installed
#install.packages("Rtsne")
#install.packages("ggplot2")
packageVersion("Rtsne")
# Load the libraries
library(Rtsne)
library(ggplot2)
setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/PCA")

# Load the gene expression dataset
data_matrix <- read.csv("27_Rna_count.csv", header = TRUE, row.names = 1)
data_matrix <- read.csv("27_Rna_all_tsne.csv", header = TRUE, row.names = 1)

# Load the metadata
metadata <- read.csv("27_Rna_meta.csv", header = TRUE, row.names = 1)

# Ensure that the metadata has a column 'ID' for grouping
grouping_variable <- metadata$group
names(grouping_variable) <- rownames(metadata)

# Check that the metadata matches the samples in the data_matrix
if (!all(colnames(data_matrix) %in% names(grouping_variable))) {
  stop("Mismatch between sample names in the data matrix and the metadata")
}

# Reorder the grouping variable to match the sample order in the data matrix
grouping_variable <- grouping_variable[colnames(data_matrix)]

# Convert the grouping variable to a factor
grouping_variable <- as.factor(grouping_variable)

set.seed(42)
# Perform t-SNE
# Note: It's important to scale/normalize your data before applying t-SNE
tsne_result <- Rtsne(t(data_matrix), dims=2, perplexity=3, verbose=TRUE, max_iter=500)

# Create a data frame for plotting
tsne_data <- data.frame(
  X = tsne_result$Y[,1],
  Y = tsne_result$Y[,2],
  Sample = colnames(data_matrix),
  Group = grouping_variable
)

# Plot the t-SNE results using ggplot2
ggplot(tsne_data, aes(x = X, y = Y, color = Group, label = Sample)) +
  geom_point(size = 3) +
  
  theme_minimal() +
  ggtitle("t-SNE Plot of Gene Expression Data") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  scale_color_manual(values = c("0" = "blue", "1" = "red"))  # Adjust colors as necessary


ggplot(tsne_data, aes(x = X, y = Y, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Group), type = "norm", level = 0.95) +  # Add ellipses around groups
  theme_minimal() +
  ggtitle("t-SNE Plot of Gene Expression Data") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  scale_color_manual(values = c("0" = "blue", "1" = "red"))  # Adjust colors as necessary


############################################  micro ############################################
# Install required packages if not already installed
install.packages("Rtsne")
install.packages("ggplot2")

# Load the libraries
library(Rtsne)
library(ggplot2)

# Load the gene expression dataset
data_matrix <- read.csv("27_micro_top_10_tSNE.csv", header = TRUE, row.names = 1)
#data_matrix <- read.csv("27_Rna_all_tsne.csv", header = TRUE, row.names = 1)

# Remove columns with missing values
data_matrix <- data_matrix[, colSums(is.na(data_matrix)) == 0]

# Load the metadata
metadata <- read.csv("Meta_Biggest_GSE11234.csv", header = TRUE, row.names = 1) 

# Ensure that the metadata has a column 'ID' for grouping
grouping_variable <- metadata$group
names(grouping_variable) <- rownames(metadata)

# Check that the metadata matches the samples in the data_matrix
if (!all(colnames(data_matrix) %in% names(grouping_variable))) {
  stop("Mismatch between sample names in the data matrix and the metadata")
}

# Reorder the grouping variable to match the sample order in the data matrix
grouping_variable <- grouping_variable[colnames(data_matrix)]

# Convert the grouping variable to a factor
grouping_variable <- as.factor(grouping_variable)

# Perform t-SNE
# Note: It's important to scale/normalize your data before applying t-SNE
tsne_result <- Rtsne(t(data_matrix), dims=2, perplexity=3, verbose=TRUE, max_iter=500)

# Create a data frame for plotting
tsne_data <- data.frame(
  X = tsne_result$Y[,1],
  Y = tsne_result$Y[,2],
  Sample = colnames(data_matrix),
  Group = grouping_variable
)

# Plot the t-SNE results using ggplot2
ggplot(tsne_data, aes(x = X, y = Y, color = Group, label = Sample)) +
  geom_point(size = 3) +
  
  theme_minimal() +
  ggtitle("t-SNE Plot of Gene Expression Data") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  scale_color_manual(values = c("0" = "blue", "1" = "red"))  # Adjust colors as necessary


ggplot(tsne_data, aes(x = X, y = Y, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Group), type = "norm", level = 0.95) +  # Add ellipses around groups
  theme_minimal() +
  ggtitle("t-SNE Plot of Gene Expression Data") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  scale_color_manual(values = c("0" = "blue", "1" = "red"))  # Adjust colors as necessary

############################################### 13 rna ###############################################


########################################################################################
# Install required packages if not already installed
install.packages("Rtsne")
install.packages("ggplot2")

# Load the libraries
library(Rtsne)
library(ggplot2)

# Load the gene expression dataset
data_matrix <- read.csv("13_rna_count2.csv", header = TRUE, row.names = 1)
#data_matrix <- read.csv("27_Rna_all_tsne.csv", header = TRUE, row.names = 1)
data_matrix <- read.csv("27_RNA_PLAGL1.csv", header = TRUE, row.names = 1)
data <- read.csv("13_RNA_Count.csv")
data <- data[-1]
data <- data[c(6, 9), ]
data_matrix <- data

data_matrix <- read.csv("6_AUC_RNA_Count.csv", header = TRUE, row.names = 1)
data_matrix <- data_matrix[c(3,5), ]

# Load the metadata
metadata <- read.csv("13_RNA_Meta.csv", header = TRUE, row.names = 1)

# Ensure that the metadata has a column 'ID' for grouping
grouping_variable <- metadata$group
names(grouping_variable) <- rownames(metadata)

# Check that the metadata matches the samples in the data_matrix
if (!all(colnames(data_matrix) %in% names(grouping_variable))) {
  stop("Mismatch between sample names in the data matrix and the metadata")
}

# Reorder the grouping variable to match the sample order in the data matrix
grouping_variable <- grouping_variable[colnames(data_matrix)]

# Convert the grouping variable to a factor
grouping_variable <- as.factor(grouping_variable)

# Perform t-SNE
# Note: It's important to scale/normalize your data before applying t-SNE
tsne_result <- Rtsne(t(data_matrix), dims=2, perplexity=5, verbose=TRUE, max_iter=500)

# Create a data frame for plotting
tsne_data <- data.frame(
  X = tsne_result$Y[,1],
  Y = tsne_result$Y[,2],
  Sample = colnames(data_matrix),
  Group = grouping_variable
)

# Plot the t-SNE results using ggplot2
ggplot(tsne_data, aes(x = X, y = Y, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 1.5) +
  theme_minimal() +
  ggtitle("t-SNE Plot of Gene Expression Data") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  scale_color_manual(values = c("0" = "blue", "1" = "red"))  # Adjust colors as necessary

library(ggplot2)

# Assuming tsne_data is a data frame with columns X, Y, Group, and Sample

ggplot(tsne_data, aes(x = X, y = Y, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Group), type = "norm", level = 0.95) +  # Add ellipses around groups
  theme_minimal() +
  ggtitle("t-SNE Plot of Gene Expression Data") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  scale_color_manual(values = c("0" = "blue", "1" = "red"))  # Adjust colors as necessary
