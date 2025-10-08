# Load required libraries
library(DESeq2)
library(ggplot2)

# Assuming 'dds' is your DESeqDataSet object

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Get raw counts
raw_counts <- counts(dds, normalized = FALSE)

# Perform PCA on normalized counts
pca_normalized <- prcomp(t(normalized_counts))

# Perform PCA on raw counts
pca_raw <- prcomp(t(raw_counts))

# Function to plot PCA with sample labels
plot_pca <- function(pca_data, title, labels) {
  pca_df <- as.data.frame(pca_data$x)
  pca_df$sample <- labels
  ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
    geom_point() +
    geom_text(nudge_x = 0.1, nudge_y = 0.1) +
    labs(title = title,
         x = paste0("PC1 (", round(100 * summary(pca_data)$importance[2,1], 1), "%)"),
         y = paste0("PC2 (", round(100 * summary(pca_data)$importance[2,2], 1), "%)"))
}

# Extract sample names
sample_names <- colnames(normalized_counts)

# Plot PCA for normalized counts with sample labels
plot_normalized <- plot_pca(pca_normalized, "PCA Plot: Normalized Counts", sample_names)

# Plot PCA for raw counts with sample labels
plot_raw <- plot_pca(pca_raw, "PCA Plot: Raw Counts", sample_names)

# Display the plots
plot_normalized
plot_raw



































# Load required libraries
library(DESeq2)
library(ggplot2)

# Assuming 'dds' is your DESeqDataSet object

# Get VSD (Variance Stabilizing Transformation)
vsd <- vst(dds, blind = FALSE)

# Perform PCA on VSD
pca_vsd <- prcomp(assay(vsd))

# Function to plot PCA with sample labels
plot_pca_vsd <- function(pca_data, title, labels) {
  pca_df <- as.data.frame(pca_data$x)
  pca_df$sample <- labels
  ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
    geom_point() +
    geom_text(nudge_x = 0.1, nudge_y = 0.1) +
    labs(title = title,
         x = paste0("PC1 (", round(100 * summary(pca_data)$importance[2,1], 1), "%)"),
         y = paste0("PC2 (", round(100 * summary(pca_data)$importance[2,2], 1), "%)"))
}

# Extract sample names
sample_names <- colnames(assay(vsd))

# Plot PCA for VSD with sample labels
plot_vsd <- plot_pca_vsd(pca_vsd, "PCA Plot: VSD", sample_names)

# Display the plot
plot_vsd




# Load required libraries
library(DESeq2)

# Assuming 'dds' is your DESeqDataSet object

# Compute VSD for normalized counts
vsd_normalized <- varianceStabilizingTransformation(dds, blind=TRUE)

# Compute VSD for raw counts
vsd_raw <- varianceStabilizingTransformation(dds, blind=TRUE, fitType="parametric")

# Plot PCA for normalized counts
plotPCA(vsd_normalized, intgroup = NULL, main = "PCA Plot: Normalized Data")

# Plot PCA for raw counts
plotPCA(vsd_raw, intgroup = NULL, main = "PCA Plot: Raw Data")

plotPCA(vsd_normalized,intgroup=c("Cell_type", "Treatment"))

plotPCA(vsd_raw,intgroup=c("Cell_type", "Treatment"))















