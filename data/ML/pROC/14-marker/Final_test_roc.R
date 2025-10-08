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

H <- as.data.frame(glm.fit$fitted.values)


lines(weight, glm.fit$fitted.values)



#roc(obese, glm.fit$fitted.values, plot=TRUE)



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
#roc.info <- roc(obese, glm.fit$fitted.values, legacy.axes=TRUE)
#str(roc.info)

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


# Select the 14 genes
genes <- X[, c("ENSG00000141646", "ENSG00000175029", "ENSG00000198431", "ENSG00000142798", 
               "ENSG00000145386", "ENSG00000134057", "ENSG00000157456", "ENSG00000156970", 
               "ENSG00000169679", "ENSG00000087586", "ENSG00000138160", "ENSG00000117399", 
               "ENSG00000171848", "ENSG00000175063")]

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







