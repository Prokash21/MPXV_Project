setwd("E:/DWCT/DEG/GSE219036_6_june/POMA")
library(SummarizedExperiment)
library(POMA)

# Example matrix
data_matrix <- matrix(runif(100), nrow=10)
rownames(data_matrix) <- paste0("Metabolite", 1:10)
colnames(data_matrix) <- paste0("Sample", 1:10)

# Example sample metadata
sample_metadata <- data.frame(
  samples = paste0("Sample", 1:10),
  group = rep(c("Control", "Treatment"), each=5)
)
sample_metadata<-read.csv("sample_metadata.csv", header = TRUE,row.names = 1)

# Ensure 'group' is a factor
sample_metadata$group <- as.factor(sample_metadata$group)

# Create a SummarizedExperiment object
se <- SummarizedExperiment(
  assays = list(counts = data_matrix),
  colData = sample_metadata
)

library(SummarizedExperiment)

library(randomForest) 

# Running PomaRandForest
poma_rf <- PomaRandForest(se, ntest = 10, nvar = 10)


poma_rf$MeanDecreaseGini


poma_rf$MeanDecreaseGini_plot



