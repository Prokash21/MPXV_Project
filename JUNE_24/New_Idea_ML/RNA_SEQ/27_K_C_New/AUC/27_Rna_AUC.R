
install.packages("missForest")
library("missForest")
library(pROC)
library(randomForest)
getwd()
setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/New_Idea_ML/RNA_SEQ/27_K_C_New/AUC")
#BiocManager::install("pROC")
#package_version("pROC")
library(pROC)

X <- read.csv("27_K_C_run.csv")
#X <- read.csv("27_K_run.csv")

head(X)

#gene_data_clean <- na.omit(X)

# Define the outcome variable and the hub gene expression
obese <- X$condition


# List of gene expression variables (adjust as per your actual column names)
gene_list <- c(
  "ATP1B3", "ATP5F1B", "BSG", "GATA6", "GM2A", "IDH2", "INHBA", "ITGB6", "NEFM", "CNOT4",
  "OVGP1", "PLAGL1", "SLC8A2", "ICAM5", "XRCC4", "ZNF212", "ELL", "DHX38", "TUBA1B", "ATP5MG",
  "NFAT5", "OGFR", "WDR37", "CBX6", "CLDN15", "UBXN7", "FGFR1OP2"
)


# Set up the plot layout
par(mfrow=c(4, 7))  # 4 rows and 4 columns

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

par(mfrow=c(1,1))  # 4 rows and 4 columns

####################################################
# Select the 14 genes
genes <- X[, c(
  "ATP1B3", "ATP5F1B", "BSG", "GATA6", "GM2A", "IDH2", "INHBA", "ITGB6", "NEFM", "CNOT4",
  "OVGP1", "PLAGL1", "SLC8A2", "ICAM5", "XRCC4", "ZNF212", "ELL", "DHX38", "TUBA1B", "ATP5MG",
  "NFAT5", "OGFR", "WDR37", "CBX6", "CLDN15", "UBXN7", "FGFR1OP2"
)]


# Convert genes to data frame
genes <- as.data.frame(genes)



################################################################################


# Set a random seed for reproducibility
set.seed(123)

# Random Forest with all genes
rf.model = randomForest(factor(obese) ~ ., data = genes, ntree = 500)

# Plotting model
plot(rf.model)

# Importance plot
importance(rf.model)

# Variable importance plot
varImpPlot(rf.model)

