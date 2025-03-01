library(pROC)
library(randomForest)
getwd()

setwd("E:/DWCT/DEG/GSE219036_6_june/POMA/14-marker")
#BiocManager::install("pROC")
#package_version("pROC")
library(pROC)

X <- read.csv("6_pROC.csv")
# Define the outcome variable and the hub gene expression
obese <- X$condition
# Assuming 'obese' is your outcome variable
obese <- X$condition

# List of gene expression variables (adjust as per your actual column names)
gene_list <- c("PLEKHM2", "CCN1", "SPATA2", "C1S", "PCBP2", "MLX")

# Set up the plot layout
par(mfrow=c(4, 4))  # 4 rows and 4 columns

# Loop through each gene
for (i in seq_along(gene_list)) {
  gene <- gene_list[i]
  weight <- X[[gene]]  # Select the gene expression variable
  
  # Logistic regression model
  glm.fit <- glm(obese ~ weight, family=binomial)
  roc_glm <- roc(obese, glm.fit$fitted.values)
  
  # Random forest model
  rf.model <- randomForest(factor(obese) ~ weight)
  roc_rf <- roc(obese, rf.model$votes[,1])
  
  # Plot ROC curve for logistic regression
  plot(roc_glm, col="#377eb8", lwd=2, percent=TRUE, 
       xlab="False Positive Percentage", ylab="True Positive Percentage",
       main=paste( gene))
  
  # Add AUC value of logistic regression
  text(0.5, 0.4, paste("AUC =", round(auc(roc_glm), 3)),col="#377eb8", adj=0, cex=0.8)
  
  # Add ROC curve for random forest
  plot.roc(roc_rf, col="#4daf4a", lwd=2, percent=TRUE, add=TRUE)
  
  # Add AUC value of random forest below the plot
  text(0.5, 0.3, paste("AUC =", round(auc(roc_rf), 3)), col="#4daf4a", adj=0, cex=0.8)
  
  # Add legend
  #legend("bottomright", legend=c("Logistic Regression", "Random Forest"),
  #      col=c("#377eb8", "#4daf4a"), lwd=2)
}


# Select the 14 genes
genes <- X[, c("PLEKHM2", "CCN1", "SPATA2", "C1S", "PCBP2", "MLX")]

# Convert genes to data frame
genes <- as.data.frame(genes)

# Check the length of the obese vector
length(obese)
factor(obese)
# Check the number of rows in the genes data frame
nrow(genes)


# Random Forest with all genes
rf.model = randomForest(factor(obese) ~ ., data = genes)


# Plotting model 
plot(rf.model) 

# Importance plot 
importance(rf.model) 

# Variable importance plot 
varImpPlot(rf.model) 




# Set a random seed for reproducibility
set.seed(42)

# Random Forest with all genes
rf.model = randomForest(factor(obese) ~ ., data = genes, ntree = 500)

# Plotting model
plot(rf.model)

# Importance plot
importance(rf.model)

# Variable importance plot
varImpPlot(rf.model)



