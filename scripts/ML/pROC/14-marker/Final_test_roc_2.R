# Load necessary libraries
library(randomForest)
library(pROC)
library(ggplot2)

# Set working directory (modify this as needed)
setwd("E:/DWCT/DEG/GSE219036_6_june/pROC")

# Read the data
X <- read.csv("14_clade_IIb_pROC.csv")

# Define the outcome variable and the gene expressions
obese <- X$condition

# Select the 14 genes
genes <- X[, c("ENSG00000141646", "ENSG00000175029", "ENSG00000198431", "ENSG00000142798", 
               "ENSG00000145386", "ENSG00000134057", "ENSG00000157456", "ENSG00000156970", 
               "ENSG00000169679", "ENSG00000087586", "ENSG00000138160", "ENSG00000117399", 
               "ENSG00000171848", "ENSG00000175063")]

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






