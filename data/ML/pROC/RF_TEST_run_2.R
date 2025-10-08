setwd("E:/DWCT/DEG/GSE219036_6_june/POMA/14-marker")

library(SummarizedExperiment)
library(POMA)



# Example matrix
data_matrix <- read.csv("your_target.csv", header = TRUE,row.names = 1)

sample_metadata<-read.csv("your_features.csv", header = TRUE,row.names = 1)

# Ensure 'group' is a factor
sample_metadata$condition <- as.factor(sample_metadata$condition)

# Create a SummarizedExperiment object
se <- SummarizedExperiment(
  assays = list(counts = data_matrix),
  colData = sample_metadata
)
# Double the sample data
#doubled_data <- rbind(se, se)

library(SummarizedExperiment)

library(randomForest) 
# Running PomaRandForest
poma_rf <- PomaRandForest(se, ntest = 10, nvar = 14)
Y<-PomaRandForest(
  se,
  ntest = NULL,
  ntree = 500,
  mtry = floor(sqrt(ncol(t(SummarizedExperiment::assay(se))))),
  nodesize = 1,
  nvar = 20
)

poma_rf$MeanDecreaseGini


poma_rf$MeanDecreaseGini_plot
Y$MeanDecreaseGini_plot
poma_rf$error_tree

## -----------------------------------------------------------------------------
X<- poma_rf$confusionMatrix$table

## -----------------------------------------------------------------------------
poma_rf$MeanDecreaseGini_plot

plot(poma_rf)







