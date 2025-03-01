# Load necessary libraries
library(ggplot2)

# Step 1: Generate sample gene expression data
set.seed(123)  # for reproducibility

# Generate data for 5 different samples, each with 30 observations
sample_data <- data.frame(
  Sample = rep(paste("Sample", 1:5), each = 30),
  Expression = c(rnorm(30, mean = 10, sd = 2),  # Sample 1
                 rnorm(30, mean = 12, sd = 2.5),  # Sample 2
                 rnorm(30, mean = 11, sd = 1.5),  # Sample 3
                 rnorm(30, mean = 13, sd = 3),  # Sample 4
                 rnorm(30, mean = 9, sd = 2))  # Sample 5
)

# Introduce some outliers
sample_data$Expression[c(5, 35, 65, 95, 125)] <- c(20, 5, 25, 2, 22)

# Step 2: Create a boxplot with ggplot2
p <- ggplot(sample_data, aes(x = Sample, y = Expression)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +
  labs(title = "Gene Expression Data",
       x = "Sample",
       y = "Expression Level") +
  theme_minimal()

# Display the plot
print(p)

###########################

# Example with random gene expression data
set.seed(123)
gene_data <- matrix(rnorm(100), nrow = 10)  # Example data matrix
colnames(gene_data) <- paste0("Sample_", 1:10)  # Sample names

# Create boxplot
boxplot(gene_data,
        main = "Boxplot of Gene Expression Data",
        xlab = "Samples",
        ylab = "Expression Level",
        col = "lightblue",
        notch = TRUE,
        outline = TRUE
)

# Example with random gene expression data
set.seed(123)
gene_data <- matrix(rnorm(100), nrow = 10)  # Example data matrix
colnames(gene_data) <- paste0("Sample_", 1:10)  # Sample names

# Perform PCA
pca_result <- prcomp(gene_data, scale = TRUE)

# Biplot of PCA
biplot(pca_result, scale = 0, cex = 0.7)
