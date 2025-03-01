setwd("E:/DWCT/DEG/GSE219036_6_june/POMA/14-marker")

# Load necessary libraries
library(pROC)
library(randomForest)

# Read the dataset
X <- read.csv("14_clade_IIb_pROC.csv", header = TRUE,row.names = 1)
#X <- rbind(data, data)

# Define the outcome variable and the gene expression data
obese <- X$condition

# Select the 14 genes
genes <- X[, c("ENSG00000141646", "ENSG00000175029", "ENSG00000198431", "ENSG00000142798", 
               "ENSG00000145386", "ENSG00000134057", "ENSG00000157456", "ENSG00000156970", 
               "ENSG00000169679", "ENSG00000087586", "ENSG00000138160", "ENSG00000117399", 
               "ENSG00000171848", "ENSG00000175063")]

# Convert genes to data frame
genes <- as.data.frame(genes)

# Logistic Regression with all genes
glm.fit = glm(obese ~ ., data = genes, family = binomial)

# Predictions and ROC for Logistic Regression
glm.probs = predict(glm.fit, type = "response")
roc(obese, glm.probs, plot = TRUE, legacy.axes = TRUE, percent = TRUE,
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#377eb8", lwd = 4, print.auc = TRUE)

# Random Forest with all genes
rf.model = randomForest(factor(obese) ~ ., data = genes)

# Predictions and ROC for Random Forest
rf.probs = rf.model$votes[, 2]  # Probability for the positive class
roc(obese, rf.probs, plot = TRUE, legacy.axes = TRUE, percent = TRUE,
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#4daf4a", lwd = 4, print.auc = TRUE)

# Layer Logistic Regression and Random Forest ROC graphs
roc(obese, glm.probs, plot = TRUE, legacy.axes = TRUE, percent = TRUE,
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#377eb8", lwd = 4, print.auc = TRUE)
plot.roc(obese, rf.probs, percent = TRUE, col = "#4daf4a", lwd = 4, print.auc = TRUE, add = TRUE, print.auc.y = 40)
legend("bottomright", legend = c("Logistic Regression", "Random Forest"), col = c("#377eb8", "#4daf4a"), lwd = 4)

# Plot the Random Forest model
plot(rf.model)

# Variable importance plot for Random Forest
varImpPlot(rf.model)

# Reset par() variables
par(pty = "m")
# Importance plot
importance(rf.model)
varImpPlot(rf.model)


#######################################################################

# Fit the logistic regression model
glm.fit <- glm(obese ~ ., data=genes, family=binomial)

# Get predicted probabilities
pred.probs <- predict(glm.fit, type="response")

# Plot ROC curve for the multivariate logistic regression model
roc(obese, pred.probs, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Positive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)

# If you want to find the optimal threshold
roc.info <- roc(obese, pred.probs, legacy.axes=TRUE)
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive percentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)
roc.df[roc.df$tpp > 60 & roc.df$tpp < 80,]


# Install and load the randomForest library if not already installed
# install.packages("randomForest")
library(randomForest)

# Fit the random forest model
rf.model <- randomForest(factor(obese) ~ ., data=genes)

# Plot ROC curve for the random forest model
roc(obese, rf.model$votes[,1], plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Positive Percentage", col="#4daf4a", lwd=4, print.auc=TRUE)

# Plotting model
plot(rf.model)

# Importance plot
importance(rf.model)
varImpPlot(rf.model)

# Overlay logistic regression and random forest ROC curves
roc(obese, pred.probs, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Positive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)
plot.roc(obese, rf.model$votes[,1], percent=TRUE, col="#4daf4a", lwd=4, print.auc=TRUE, add=TRUE, print.auc.y=40)
legend("bottomright", legend=c("Logistic Regression", "Random Forest"), col=c("#377eb8", "#4daf4a"), lwd=4)


# Install and load the randomForest library if not already installed
# install.packages("randomForest")
library(randomForest)

# Fit the random forest model
rf.model <- randomForest(factor(obese) ~ ., data=genes)

# Plot ROC curve for the random forest model
roc(obese, rf.model$votes[,1], plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Positive Percentage", col="#4daf4a", lwd=4, print.auc=TRUE)

# Plotting model
plot(rf.model)

# Importance plot
importance(rf.model)
varImpPlot(rf.model)

# Overlay logistic regression and random forest ROC curves
roc(obese, pred.probs, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Positive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)
plot.roc(obese, rf.model$votes[,1], percent=TRUE, col="#4daf4a", lwd=4, print.auc=TRUE, add=TRUE, print.auc.y=40)
legend("bottomright", legend=c("Logistic Regression", "Random Forest"), col=c("#377eb8", "#4daf4a"), lwd=4)


# AUC comparison
logistic_auc <- roc(obese, pred.probs)$auc
random_forest_auc <- roc(obese, rf.model$votes[,1])$auc
print(paste("Logistic Regression AUC:", logistic_auc))
print(paste("Random Forest AUC:", random_forest_auc))

# Feature importance
print(importance(rf.model))

