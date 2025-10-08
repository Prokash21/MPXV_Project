

###################################################### RNA #############################

setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Correlation")

# Load necessary libraries
#install.packages("corrplot")
library(corrplot)
packageVersion("corrplot")

# Load your gene expression data
gene_expression_data_14_genes <- read.csv("Fibro_27.csv" , header = TRUE,row.names = 1)

head(gene_expression_data_14_genes)

marker_genes <- colnames(gene_expression_data_14_genes)


# Compute Pearson correlation between the marker genes
gene_correlation_matrix <- cor(gene_expression_data_14_genes, method = "pearson")

# Visualize the correlation matrix
corrplot(gene_correlation_matrix, method = "color", 
         type = "upper", tl.col = "black", tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.7, 
         col = colorRampPalette(c("red", "white", "blue"))(200))
getwd()

gene_correlation_matrix <- as.data.frame(gene_correlation_matrix)
write.csv(gene_correlation_matrix, file = "27_fibro_correlation_matrix.csv")

# 11_all_MPXV

#setwd("E:/DWCT/DEG/GSE219036_6_june/Co-relation")

# Load necessary libraries
#install.packages("corrplot")
library(corrplot)
packageVersion("corrplot")

# Load your gene expression data
gene_expression_data_14_genes <- read.csv("RNA_27.csv" , header = TRUE,row.names = 1)

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
write.csv(gene_correlation_matrix, file = "27_Rna_correlation_matrix.csv")


# 6_old_MPXV

setwd("E:/DWCT/DEG/GSE219036_6_june/Co-relation")

# Load necessary libraries
install.packages("corrplot")
library(corrplot)
packageVersion("corrplot")

# Load your gene expression data
gene_expression_data_14_genes <- read.csv("27_LFC_Cor.csv" , header = TRUE,row.names = 1)

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
write.csv(gene_correlation_matrix, file = "27_LFC_correlation_matrix.csv")



