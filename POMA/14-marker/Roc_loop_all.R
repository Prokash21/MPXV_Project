getwd()
setwd("E:/DWCT/DEG/GSE219036_6_june/pROC")
#BiocManager::install("pROC")
#package_version("pROC")
library(pROC)

X <- read.csv("14_clade_IIb_pROC.csv")
# Define the outcome variable and the hub gene expression
obese <- X$condition

################################# 14 gene ############################
#1
weight <- X$ENSG00000141646
#2
weight <- X$ENSG00000175029
#3
weight <- X$ENSG00000198431
#4
weight <- X$ENSG00000142798
#5
weight <- X$ENSG00000145386
#6
weight <- X$ENSG00000134057
#7auc 79.6
weight <- X$ENSG00000157456
#8 96.3,77.8
weight <- X$ENSG00000156970
#9
weight <- X$ENSG00000169679
#10
weight <- X$ENSG00000087586
#11
weight <- X$ENSG00000138160
#12
weight <- X$ENSG00000117399
#13
weight <- X$ENSG00000171848
#14
weight <- X$ENSG00000175063


############################################################

plot(x=weight, y=obese)

y <- as.data.frame(weight)
glm.fit=glm(obese ~ weight, family=binomial)

par(pty = "s") ## pty sets the aspect ratio of the plot region. Two options:
##                "s" - creates a square plotting region
##                "m" - (the default) creates a maximal plotting region
roc(obese, glm.fit$fitted.values, plot=TRUE)

## We can also change the color of the ROC line, and make it wider...
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4)
## We can calculate the area under the curve...
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)


library("randomForest")
rf.model <- randomForest(factor(obese) ~ weight)

## ROC for random forest
roc(obese, rf.model$votes[,1], plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#4daf4a", lwd=4, print.auc=TRUE)


roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)

plot.roc(obese, rf.model$votes[,1], percent=TRUE, col="#4daf4a", lwd=4, print.auc=TRUE, add=TRUE, print.auc.y=40)
legend("bottomright", legend=c("Logisitic Regression", "Random Forest"), col=c("#377eb8", "#4daf4a"), lwd=4)


#######################################################

obese <- X$condition  # Assuming this is your outcome variable

# List of gene expression variables (adjust as per your actual column names)
gene_list <- c("ENSG00000141646", "ENSG00000175029", "ENSG00000198431", "ENSG00000142798",
               "ENSG00000145386", "ENSG00000134057", "ENSG00000157456", "ENSG00000156970",
               "ENSG00000169679", "ENSG00000087586", "ENSG00000138160", "ENSG00000117399",
               "ENSG00000171848", "ENSG00000175063")

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
       main=paste("ROC Curve for", gene,
                  print.auc=TRUE))
  
  # Add ROC curve for random forest
  plot.roc(roc_rf, col="#4daf4a", lwd=2, percent=TRUE, add=TRUE,
           print.auc=TRUE)
  
  # Add legend
  #legend("bottomright", legend=c("Logistic Regression", "Random Forest"),
   #      col=c("#377eb8", "#4daf4a"), lwd=2)
}

##########################################################################################################

library(pROC)
library(randomForest)

# Assuming 'obese' is your outcome variable
obese <- X$condition

# List of gene expression variables (adjust as per your actual column names)
gene_list <- c("ENSG00000141646", "ENSG00000175029", "ENSG00000198431", "ENSG00000142798",
               "ENSG00000145386", "ENSG00000134057", "ENSG00000157456", "ENSG00000156970",
               "ENSG00000169679", "ENSG00000087586", "ENSG00000138160", "ENSG00000117399",
               "ENSG00000171848", "ENSG00000175063")

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
       main=paste("ROC Curve for", gene, "- Logistic Regression"),
       print.auc=TRUE)
  
  # Add ROC curve for random forest
  plot.roc(roc_rf, col="#4daf4a", lwd=2, percent=TRUE, add=TRUE,
           print.auc.y=40)
  
  # Add legend
  #legend("bottomright", legend=c("Logistic Regression", "Random Forest"),
   #      col=c("#377eb8", "#4daf4a"), lwd=2)
}


############################################################################################3

library(pROC)
library(randomForest)

# Assuming 'obese' is your outcome variable
obese <- X$condition

# List of gene expression variables (adjust as per your actual column names)
gene_list <- c("ENSG00000141646", "ENSG00000175029", "ENSG00000198431", "ENSG00000142798",
               "ENSG00000145386", "ENSG00000134057", "ENSG00000157456", "ENSG00000156970",
               "ENSG00000169679", "ENSG00000087586", "ENSG00000138160", "ENSG00000117399",
               "ENSG00000171848", "ENSG00000175063")

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
       main=paste("ROC Curve for", gene, "- Logistic Regression"))
  
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
