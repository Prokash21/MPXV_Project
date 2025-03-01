# transpose the data because in variance_stabilised_counts the rows are the variables and the columns correspond to the samples
t_variance_stabilised_counts <- t(count_data)
install.packages("mypca")
library("myPCA")
BiocManager::install("myPCA")
# before computing the PCA, check that samples are in rows and genes in columns

X<-as.data.frame(dds)


rld_ds <- rlog(dds, blind=FALSE, normalized = FALSE) #regularized log-transformation

# PCA plot 2 using ggplot2
(pcadata <- plotPCA(rld_ds, intgroup = c( "Cell_type", "Treatment"), returnData=TRUE))


  
# Extract percent variance
percentVar <- round(100 * attr(pcadata, "percentVar"))

# Plot PCA using ggplot2
ggplot(pcadata, aes(PC1, PC2, color = Treatment, label = name)) +
  geom_point(size = 3) +
  geom_text(hjust = 0.5, vjust = -1, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))



# Step 1: Create a DESeqDataSet object from your raw count data
dds <- DESeqDataSetFromMatrix(countData = your_raw_counts_matrix,
                              colData = your_sample_metadata,
                              design = ~ condition)

# Step 2: Filter low count genes if needed
dds <- dds[rowSums(counts(dds)) > 1,]

# Step 3: Transform the count data (e.g., using log transformation)
log_counts <- log2(counts(dds) + 1) # adding 1 to avoid log(0)

# Step 4: Perform PCA on the log-transformed counts
pca_result <- prcomp(t(log_counts))

# Step 5: Plot PCA
plot(pca_result$x[,1], pca_result$x[,2], 
     xlab = paste0("PC1: ", round(summary(pca_result)$importance[2,1] * 100, 2), "% variance"),
     ylab = paste0("PC2: ", round(summary(pca_result)$importance[2,2] * 100, 2), "% variance"),
     main = "PCA Plot of Raw Count Data",
     pch = 19)
pca_result$x
# Optional: Customize the plot as needed

  

# Step 1: Create a DESeqDataSet object from your raw count data
dds <- DESeqDataSetFromMatrix(countData = your_raw_counts_matrix,
                              colData = your_sample_metadata,
                              design = ~ condition)

# Step 2: Perform normalization and transformation (e.g., using DESeq2's default normalization method)
dds <- DESeq(dds)

# Step 3: Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts1 <- counts(dds, normalized = FALSE)

# Step 4: Perform PCA on the normalized counts
pca_result <- prcomp(t(normalized_counts))
pca_result1 <- prcomp(t(normalized_counts1))
# Step 5: Plot PCA
plot(pca_result$x[,1], pca_result$x[,2], 
     xlab = paste0("PC1: ", round(summary(pca_result)$importance[2,1] * 100, 2), "% variance"),
     ylab = paste0("PC2: ", round(summary(pca_result)$importance[2,2] * 100, 2), "% variance"),
     main = "PCA Plot of Raw Count Data",
     pch = 19)
# Adding labels to points
text(pca_result$x[,1], pca_result$x[,2], labels = colnames(sample_data), pos = 1)
# Optional: Customize the plot as needed
