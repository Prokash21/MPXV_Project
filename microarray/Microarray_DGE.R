
getwd()
setwd("E:/DWCT/DEG/GSE219036_6_june/microarray")
# script to perform RMA normalization
# setwd("~/Desktop/demo/microarray_RMA")
install.packages("affy")
BiocManager::install("affy")
library(affy)
library(GEOquery)
library(tidyverse)

# get supplementary files
getGEOSuppFiles("GSE24125")

# untar files
gse<-untar("GSE24125/GSE24125_RAW.tar", exdir = 'data/')

# reading in .cel files
raw.data <- ReadAffy(celfile.path = "data/")

# performing RMA normalization
normalized.data <- rma(raw.data)

# get expression estimates
normalized.expr <- as.data.frame(exprs(normalized.data))

# map probe IDs to gene symbols
gse <- getGEO("GSE148537", GSEMatrix = TRUE)

# fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE148537_series_matrix.txt.gz@featureData@data
# subset
feature.data <- feature.data[,c(1,11)]

normalized.expr <- normalized.expr %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature.data, by = 'ID')