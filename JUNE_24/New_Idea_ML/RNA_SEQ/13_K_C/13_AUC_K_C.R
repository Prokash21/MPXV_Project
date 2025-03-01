library(pROC)
library(randomForest)
getwd()

setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/New_Idea_ML/RNA_SEQ/13_K_C")

# List the files in the directory
files_in_directory <- list.files()

# Print the list of files
print(files_in_directory)

#BiocManager::install("pROC")
#package_version("pROC")
library(pROC)

X <- read.csv("13_K_C_ex_run.csv")
#colon all 0
#X <- read.csv("13_3X1_K_C_ex_run.csv")


# Define the outcome variable and the hub gene expression
obese <- X$condition

# List of gene expression variables (adjust as per your actual column names)
gene_list <- c(
  "AURKA",
  "CDC20",
  "CCNB1",
  "KIF11",
  "HSPG2",
  "CCNA2",
  "BUB1B",
  "CCNB2",
  "BUB1",
  "RRM2",
  "CTBP2",
  "UBE2C",
  "TXNRD1"
)

# Set up the plot layout
par(mfrow=c(3, 5))  # 4 rows and 4 columns

# Loop through each gene
for (i in seq_along(gene_list)) {
  gene <- gene_list[i]
  weight <- X[[gene]]  # Select the gene expression variable
  
  # Logistic regression model
  glm.fit <- glm(obese ~ weight, family=binomial)
  roc_glm <- roc(obese, glm.fit$fitted.values)
  
  # Random forest model
  set.seed(123) 
  rf.model <- randomForest(factor(obese) ~ weight, ntree=500)
  roc_rf <- roc(obese, rf.model$votes[,1])
  
  # Plot ROC curve for logistic regression
  plot(roc_glm, col="#377eb8", lwd=2, percent=TRUE, 
       xlab="Specificity ", ylab="Sensitivity",
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
par(mfrow=c(1, 1))

# Select the 14 genes
genes <- X[, c(
  "AURKA",
  "CDC20",
  "CCNB1",
  "KIF11",
  "HSPG2",
  "CCNA2",
  "BUB1B",
  "CCNB2",
  "BUB1",
  "RRM2",
  "CTBP2",
  "UBE2C",
  "TXNRD1"
)]

# Convert genes to data frame
genes <- as.data.frame(genes)



################################################################################


# Set a random seed for reproducibility
set.seed(123)
#set.seed(123)

# Random Forest with all genes
rf.model = randomForest(factor(obese) ~ ., data = genes, ntree = 500)

# Plotting model
plot(rf.model)

# Add a color legend
legend("topright", legend=c("Error Rate", "Class 1", "Class 2"), 
       col=c("black", "red", "green"), lty=1, cex=0.8)
# Importance plot
importance(rf.model)

# Variable importance plot
varImpPlot(rf.model)

# Add a color legend
legend("top", legend=c("13_RNA"), 
       )
