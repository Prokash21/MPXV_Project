
#14_2022_MPXV

setwd("E:/DWCT/DEG/GSE219036_6_june/Co-relation")

# Load necessary libraries
#install.packages("corrplot")
library(corrplot)
packageVersion("corrplot")

# Load your gene expression data
gene_expression_data_14_genes <- read.csv("14.1_clade_IIb_pROC.csv" , header = TRUE,row.names = 1)

head(gene_expression_data_14_genes)

marker_genes <- colnames(gene_expression_data_14_genes[-14])


# Compute Pearson correlation between the marker genes
gene_correlation_matrix <- cor(gene_expression_data_14_genes, method = "pearson")

# Visualize the correlation matrix
corrplot(gene_correlation_matrix, method = "color", 
         type = "upper", tl.col = "black", tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.7, 
         col = colorRampPalette(c("red", "white", "blue"))(200))


gene_correlation_matrix <- as.data.frame(gene_correlation_matrix)
write.csv(gene_correlation_matrix, file = "13_gene_correlation_matrix.csv")

# 11_all_MPXV

setwd("E:/DWCT/DEG/GSE219036_6_june/Co-relation")

# Load necessary libraries
#install.packages("corrplot")
library(corrplot)
packageVersion("corrplot")

# Load your gene expression data
gene_expression_data_14_genes <- read.csv("11_your_target.csv" , header = TRUE,row.names = 1)

head(gene_expression_data_14_genes)

marker_genes <- colnames(gene_expression_data_14_genes)


# Compute Pearson correlation between the marker genes
gene_correlation_matrix <- cor(gene_expression_data_14_genes, method = "pearson")


# Visualize the correlation matrix
corrplot(gene_correlation_matrix, method = "color", 
         type = "upper", tl.col = "black", tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.7, 
         col = colorRampPalette(c("red", "white", "blue"))(200))

gene_correlation_matrix <- as.data.frame(gene_correlation_matrix)
write.csv(gene_correlation_matrix, file = "11_gene_correlation_matrix.csv")


# 6_old_MPXV

setwd("E:/DWCT/DEG/GSE219036_6_june/Co-relation")

# Load necessary libraries
install.packages("corrplot")
library(corrplot)
packageVersion("corrplot")

# Load your gene expression data
gene_expression_data_14_genes <- read.csv("6_your_target.csv" , header = TRUE,row.names = 1)

head(gene_expression_data_14_genes)

marker_genes <- colnames(gene_expression_data_14_genes)


# Compute Pearson correlation between the marker genes
gene_correlation_matrix <- cor(gene_expression_data_14_genes, method = "pearson")

# Visualize the correlation matrix
corrplot(gene_correlation_matrix, method = "color", 
         type = "upper", tl.col = "black", tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.7, 
         col = colorRampPalette(c("red", "white", "blue"))(200))

gene_correlation_matrix <- as.data.frame(gene_correlation_matrix)
write.csv(gene_correlation_matrix, file = "6_gene_correlation_matrix.csv")



