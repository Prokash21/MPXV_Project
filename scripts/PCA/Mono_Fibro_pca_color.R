library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

# Set working directory and load data
setwd("E:/DWCT/DEG/GSE219036_6_june/P.PCA/GSE219036_TPM_for_GEO_upload_skin.txt")

data <- read.csv("11234_expression_Mono.csv", header = TRUE, row.names = 1)
metadata <- read.csv("11234_expression_Mono_Meta.csv", header = TRUE, row.names = 1) 

# Quality Control: detect outlier genes
gsg <- goodSamplesGenes(t(data))
summary(gsg)

# Remove outlier genes
data <- data[gsg$goodGenes == TRUE,]

# Detect outlier samples using hierarchical clustering
htree <- hclust(dist(t(data)), method = "average")

# Plot hierarchical clustering with colors based on metadata groups
plot(htree, labels = metadata$group, main = "Hierarchical Clustering Dendrogram")
rect.hclust(htree, k = length(unique(metadata$group)), border = 2:5)

# Check for NA or Infinite values
summary(data)
is.na(data)
is.infinite(data)

# Replace NA and Infinite values with zero
data[is.na(data)] <- 0
data[is.infinite(data)] <- 0

# Verify no NA or Infinite values remain
summary(data)

# Remove non-numeric columns for PCA
data_numeric <- data[, sapply(data, is.numeric)]

# Perform PCA
pca <- prcomp(t(data_numeric))

# View the PCA results
summary(pca)

# Prepare PCA data for plotting
pca.dat <- as.data.frame(pca$x)
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var / sum(pca.var) * 100, digits = 2)

# Merge PCA data with metadata
pca.dat <- cbind(pca.dat, metadata)

# Plot PCA with metadata groups
ggplot(pca.dat, aes(PC1, PC2, color = group)) +
  geom_point() +
  geom_text(aes(label = rownames(pca.dat)), hjust = 0, vjust = 1) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %')) +
  theme_minimal() +
  theme(legend.title = element_blank())


###########################################################################################


setwd("E:/DWCT/DEG/GSE219036_6_june/P.PCA/GSE219036_TPM_for_GEO_upload_skin.txt")

data <- read.csv("Exprs_GSE11234_8vs8_normalized_Fibro.csv", header = TRUE, row.names = 1)
metadata <- read.csv("Metadata_GSE11234_8vs8_normalized_Fibro.csv", header = TRUE, row.names = 1) 

# Quality Control: detect outlier genes
gsg <- goodSamplesGenes(t(data))
summary(gsg)

# Remove outlier genes
data <- data[gsg$goodGenes == TRUE,]

# Detect outlier samples using hierarchical clustering
htree <- hclust(dist(t(data)), method = "average")

# Plot hierarchical clustering with colors based on metadata groups
plot(htree, labels = metadata$group, main = "Hierarchical Clustering Dendrogram")
rect.hclust(htree, k = length(unique(metadata$group)), border = 2:5)

# Check for NA or Infinite values
summary(data)
is.na(data)
is.infinite(data)

# Replace NA and Infinite values with zero
data[is.na(data)] <- 0
data[is.infinite(data)] <- 0

# Verify no NA or Infinite values remain
summary(data)

# Remove non-numeric columns for PCA
data_numeric <- data[, sapply(data, is.numeric)]

# Perform PCA
pca <- prcomp(t(data_numeric))

# View the PCA results
summary(pca)

# Prepare PCA data for plotting
pca.dat <- as.data.frame(pca$x)
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var / sum(pca.var) * 100, digits = 2)

# Merge PCA data with metadata
pca.dat <- cbind(pca.dat, metadata)

# Plot PCA with metadata groups
ggplot(pca.dat, aes(PC1, PC2, color = group)) +
  geom_point() +
  geom_text(aes(label = rownames(pca.dat)), hjust = 0, vjust = 1) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %')) +
  theme_minimal() +
  theme(legend.title = element_blank())


#######################################################################################

setwd("E:/DWCT/DEG/GSE219036_6_june/P.PCA/GSE219036_TPM_for_GEO_upload_skin.txt")

data <- read.csv("Exprs_GSE11234_8vs8_normalized_Fibro.csv", header = TRUE, row.names = 1)
metadata <- read.csv("Metadata_GSE11234_8vs8_normalized_Fibro.csv", header = TRUE, row.names = 1) 

# Quality Control: detect outlier genes
gsg <- goodSamplesGenes(t(data))
summary(gsg)

# Remove outlier genes
data <- data[gsg$goodGenes == TRUE,]

# Detect outlier samples using hierarchical clustering
htree <- hclust(dist(t(data)), method = "average")

# Plot hierarchical clustering with colors based on metadata groups
plot(htree, labels = metadata$group, main = "Hierarchical Clustering Dendrogram")
rect.hclust(htree, k = length(unique(metadata$group)), border = 2:5)

# Check for NA or Infinite values
summary(data)
is.na(data)
is.infinite(data)

# Replace NA and Infinite values with zero
data[is.na(data)] <- 0
data[is.infinite(data)] <- 0

# Verify no NA or Infinite values remain
summary(data)

# Remove non-numeric columns for PCA
data_numeric <- data[, sapply(data, is.numeric)]

# Perform PCA
pca <- prcomp(t(data_numeric))

# View the PCA results
summary(pca)

# Prepare PCA data for plotting
pca.dat <- as.data.frame(pca$x)
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var / sum(pca.var) * 100, digits = 2)

# Merge PCA data with metadata
pca.dat <- cbind(pca.dat, metadata)

# Plot PCA with metadata groups
ggplot(pca.dat, aes(PC1, PC2, color = group)) +
  geom_point() +
  geom_text(aes(label = rownames(pca.dat)), hjust = 0, vjust = 1) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %')) +
  theme_minimal() +
  theme(legend.title = element_blank())







