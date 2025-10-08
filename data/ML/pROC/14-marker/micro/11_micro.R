getwd()
setwd("E:/DWCT/DEG/GSE219036_6_june/POMA/14-marker/micro")
#BiocManager::install("pROC")
#package_version("pROC")
library(pROC)

X <- read.csv("11_micro_marker.csv" , header = TRUE,row.names = 1)
# Define the outcome variable and the hub gene expression
obese <- X$group

################################# 14 gene ############################

#1
weight <- X$WDR37
#2
weight <- X$NFAT5
#3
weight <- X$ICAM5
#4
weight <- X$ELL
#5
weight <- X$PLAGL1
#6
weight <- X$DHX38
#7
weight <- X$GATA6
#8 
weight <- X$ZNF212
#9
weight <- X$CBX6
#10
weight <- X$ATP5MG
#11
weight <- X$GM2A



############################################################

plot(x=weight, y=obese)

y <- as.data.frame(weight)
glm.fit=glm(obese ~ weight, family=binomial)

H <- as.data.frame(glm.fit$fitted.values)


lines(weight, glm.fit$fitted.values)



roc(obese, glm.fit$fitted.values, plot=TRUE)



## Now let's configure R so that it prints the graph as a square.
##
par(pty = "s") ## pty sets the aspect ratio of the plot region. Two options:
##                "s" - creates a square plotting region
##                "m" - (the default) creates a maximal plotting region
roc(obese, glm.fit$fitted.values, plot=TRUE)

## NOTE: By default, roc() uses specificity on the x-axis and the values range
## from 1 to 0. This makes the graph look like what we would expect, but the
## x-axis itself might induce a headache. To use 1-specificity (i.e. the 
## False Positive Rate) on the x-axis, set "legacy.axes" to TRUE.
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE)

## If you want to rename the x and y axes...
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage")

## We can also change the color of the ROC line, and make it wider...
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4)

## If we want to find out the optimal threshold we can store the 
## data used to make the ROC graph in a variable...
roc.info <- roc(obese, glm.fit$fitted.values, legacy.axes=TRUE)
str(roc.info)

## and then extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df) ## head() will show us the values for the upper right-hand corner
## of the ROC graph, when the threshold is so low 
## (negative infinity) that every single sample is called "obese".
## Thus TPP = 100% and FPP = 100%

tail(roc.df) ## tail() will show us the values for the lower left-hand corner
## of the ROC graph, when the threshold is so high (infinity) 
## that every single sample is called "not obese". 
## Thus, TPP = 0% and FPP = 0%

## now let's look at the thresholds between TPP 60% and 80%...
roc.df[roc.df$tpp > 60 & roc.df$tpp < 80,]

## We can calculate the area under the curve...
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)

## ...and the partial area under the curve.
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, print.auc.x=45, partial.auc=c(100, 90), auc.polygon = TRUE, auc.polygon.col = "#377eb822")


#######################################
##
## Now let's fit the data with a random forest...
##
#######################################
#install.packages("randomForest")
#BiocManager::install("randomForest")
#install.packages("randomForest") 
library("randomForest")
rf.model <- randomForest(factor(obese) ~ weight)

## ROC for random forest
roc(obese, rf.model$votes[,1], plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#4daf4a", lwd=4, print.auc=TRUE)

# Plotting model 
plot(rf.model) 

# Importance plot 
importance(rf.model) 

# Variable importance plot 
varImpPlot(rf.model) 
#######################################
##
## Now layer logistic regression and random forest ROC graphs..
##
#######################################
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)

plot.roc(obese, rf.model$votes[,1], percent=TRUE, col="#4daf4a", lwd=4, print.auc=TRUE, add=TRUE, print.auc.y=40)
legend("bottomright", legend=c("Logisitic Regression", "Random Forest"), col=c("#377eb8", "#4daf4a"), lwd=4)


#######################################
##
## Now that we're done with our ROC fun, let's reset the par() variables.
## There are two ways to do it...
##
#######################################
par(pty = "m")


# Select the 11 genes
genes <- X[, c("WDR37", "NFAT5", "ICAM5", "ELL", "PLAGL1", "DHX38", "GATA6", "ZNF212", "CBX6", "ATP5MG", "GM2A")]

# Convert genes to data frame
genes <- as.data.frame(genes)
head(genes)
# Random Forest with all genes
rf.model = randomForest(factor(obese) ~ ., data = genes)


# Plotting model 
plot(rf.model) 

# Importance plot 
importance(rf.model) 

# Variable importance plot 
varImpPlot(rf.model) 



#########################################################################################################

# Load necessary libraries
library(randomForest)
library(pROC)
library(ggplot2)
getwd()
setwd("E:/DWCT/DEG/GSE219036_6_june/POMA/14-marker/micro")
#BiocManager::install("pROC")
#package_version("pROC")
library(pROC)

X <- read.csv("11_micro_marker.csv" , header = TRUE,row.names = 1)

# Define the outcome variable and the gene expressions
obese <- X$group

genes <- X[, c("WDR37", "NFAT5", "ICAM5", "ELL", "PLAGL1", "DHX38", "GATA6", "ZNF212", "CBX6", "ATP5MG", "GM2A")]


# Convert genes to data frame
genes <- as.data.frame(genes)

# Fit Random Forest model with all genes
rf.model <- randomForest(factor(obese) ~ ., data = genes)

# Plotting Random Forest model
plot(rf.model)

# Importance plot
importance(rf.model)

# Variable importance plot
varImpPlot(rf.model)

# Extract and print importance scores
importance_scores <- importance(rf.model)
print(importance_scores)

# Fit logistic regression model with all genes
glm.fit <- glm(obese ~ ., data = genes, family = binomial)

# Extract coefficients
coefficients <- summary(glm.fit)$coefficients
print(coefficients)

# Create a data frame with coefficients for better visualization
coeff_df <- data.frame(Gene = rownames(coefficients), Coefficient = coefficients[, "Estimate"])
print(coeff_df)

# Plot the coefficients
ggplot(coeff_df, aes(x = Gene, y = Coefficient)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Logistic Regression Coefficients")

# Layering logistic regression and random forest ROC graphs
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, 
    xlab="False Positive Percentage", ylab="True Positive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)

# Adding ROC for Random Forest model
plot.roc(obese, rf.model$votes[,1], percent=TRUE, col="#4daf4a", lwd=4, print.auc=TRUE, add=TRUE, print.auc.y=40)

# Adding legend to the ROC plot
legend("bottomright", legend=c("Logistic Regression", "Random Forest"), col=c("#377eb8", "#4daf4a"), lwd=4)










# Extract and prepare importance scores
importance_scores <- data.frame(
  Gene = rownames(importance(rf.model)),
  Importance = importance(rf.model)[, 1]
)

# Plot the importance scores
ggplot(importance_scores, aes(x = reorder(Gene, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "orange") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Random Forest Feature Importance",
       x = "Gene",
       y = "Importance Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))








