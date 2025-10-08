setwd("E:/DWCT/DEG/GSE219036_6_june/POMA/shuvo vai")


install.packages("randomForest")
install.packages("caret")
install.packages("e1071")
install.packages("Rtsne")
install.packages("umap")
install.packages("ggplot2")


data <- read.csv("14_clade_IIb_pROC.csv", header = TRUE,row.names = 1)
#data <- rbind(data, data)
library(caret)
set.seed(42)
index <- createDataPartition(data$condition, p=0.7, list=FALSE)
train_data <- data[index,]
test_data <- data[-index,]

logistic_model <- glm(condition ~ ., data = train_data, family = binomial())
summary(logistic_model)

library(randomForest)
set.seed(42)
rf_model <- randomForest(condition ~ ., data=train_data, ntree=100)
print(rf_model)

importance(rf_model)
varImpPlot(rf_model)

# Plotting model 
plot(rf_model) 

# Importance plot 
importance(rf_model) 

# Variable importance plot 
varImpPlot(rf_model) 



library(Rtsne)
set.seed(42)
tsne_data <- Rtsne(train_data[, -which(names(train_data) == "condition")], dims = 2, perplexity = 3)
tsne_df <- data.frame(tsne_data$Y, train_data$condition)
colnames(tsne_df) <- c("Comp1", "Comp2", "Condition")
ggplot(tsne_df, aes(x=Comp1, y=Comp2, color=factor(Condition))) + geom_point() + labs(color = 'Condition')


library(umap)
umap_data <- umap(train_data[, -which(names(train_data) == "condition")])
umap_df <- data.frame(umap_data$layout, train_data$condition)
colnames(umap_df) <- c("Comp1", "Comp2", "Condition")
ggplot(umap_df, aes(x=Comp1, y=Comp2, color=factor(Condition))) + geom_point() + labs(color = 'Condition')

pca <- prcomp(train_data[, -which(names(train_data) == "condition")], center = TRUE, scale. = TRUE)
pca_df <- data.frame(pca$x[,1:2], train_data$condition)
colnames(pca_df) <- c("PC1", "PC2", "Condition")
ggplot(pca_df, aes(x=PC1, y=PC2, color=factor(Condition))) + geom_point() + labs(color = 'Condition')



