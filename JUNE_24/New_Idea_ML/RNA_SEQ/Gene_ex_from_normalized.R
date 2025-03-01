


setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/New_Idea_ML/RNA_SEQ")

# Load necessary libraries
library(tidyverse)
library(dplyr)

library(readr)
heatmap_fibro <- read_csv("Colon_Normalized.csv")
#View(heatmap_fibro)

data <- heatmap_fibro

common <- c(
  "ENSG00000047056", "ENSG00000060491", "ENSG00000069849", "ENSG00000080802", 
  "ENSG00000085465", "ENSG00000102908", "ENSG00000104722", "ENSG00000105376", 
  "ENSG00000105656", "ENSG00000106404", "ENSG00000110955", "ENSG00000111790", 
  "ENSG00000115221", "ENSG00000118160", "ENSG00000118495", "ENSG00000122641", 
  "ENSG00000123416", "ENSG00000140829", "ENSG00000141448", "ENSG00000152422", 
  "ENSG00000163960", "ENSG00000167283", "ENSG00000170260", "ENSG00000172270", 
  "ENSG00000182054", "ENSG00000183741", "ENSG00000196743"
)




X <- data[data$...1 %in% common,]

write.csv(X,file = "27_RNa_Correlation_expression.csv")

################################################# 13 ###################################################

setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/New_Idea_ML/RNA_SEQ")

# Load necessary libraries
library(tidyverse)
library(dplyr)

library(readr)
heatmap_fibro <- read_csv("Colon_Normalized.csv")
#View(heatmap_fibro)

data <- heatmap_fibro

common <- c(
  "ENSG00000175029",
  "ENSG00000198431",
  "ENSG00000142798",
  "ENSG00000145386",
  "ENSG00000134057",
  "ENSG00000157456",
  "ENSG00000156970",
  "ENSG00000169679",
  "ENSG00000087586",
  "ENSG00000138160",
  "ENSG00000117399",
  "ENSG00000171848",
  "ENSG00000175063"
)




X <- data[data$...1 %in% common,]

write.csv(X,file = "13_Colon_Correlation_expression.csv")

####################################################################################################



