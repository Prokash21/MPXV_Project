


library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

setwd("E:/DWCT/DEG/GSE219036_6_june/P.PCA/GSE219036_TPM_for_GEO_upload_skin.txt")

data <- read.csv("11234_expression_Mono.csv", header = TRUE, row.names = 1)
metadata <- read.csv("11234_expression_Mono_Meta.csv", header = TRUE, row.names = 1) 

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

# Plot hierarchical clustering with colors based on metadata groups
plot(htree, labels = metadata$group, main = "Hierarchical Clustering Dendrogram")
rect.hclust(htree, k = length(unique(metadata$group)), border = 2:5)

# Check for NA or Infinite values
summary(data)
is.na(data)
is.infinite(data)

# Replace NA with the column mean (or any other imputation method)
data[is.na(data)] <- 0
data[is.infinite(data)] <- 0

# Alternatively, you can remove rows with NA or Infinite values
# data <- na.omit(data)
# data <- data[!rowSums(!is.finite(data)), ]

# Verify no NA or Infinite values remain
summary(data)

# Remove non-numeric columns for PCA (e.g., gene_id)
data_numeric <- data[, sapply(data, is.numeric)]

# Perform PCA
pca <- prcomp(t(data_numeric))

# View the PCA results
summary(pca)
plot(pca)

#pca - method 2

#pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))





######################################################################################


setwd("E:/DWCT/DEG/GSE219036_6_june/P.PCA/GSE219036_TPM_for_GEO_upload_skin.txt")

data <- read.csv("11234_expression_Fibro.csv", header = TRUE, row.names = 1)


# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

# Check for NA or Infinite values
summary(data)
is.na(data)
is.infinite(data)

# Replace NA with the column mean (or any other imputation method)
data[is.na(data)] <- 0
data[is.infinite(data)] <- 0

# Alternatively, you can remove rows with NA or Infinite values
# data <- na.omit(data)
# data <- data[!rowSums(!is.finite(data)), ]

# Verify no NA or Infinite values remain
summary(data)

# Remove non-numeric columns for PCA (e.g., gene_id)
data_numeric <- data[, sapply(data, is.numeric)]

# Perform PCA
pca <- prcomp(t(data_numeric))

# View the PCA results
summary(pca)
plot(pca)

# pca - method 2

#pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
