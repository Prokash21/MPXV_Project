


setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/New_Idea_ML/RNA_SEQ/13_K_C")

files <- list.files()
print(files)
# Load necessary libraries
library(tidyverse)
library(dplyr)

library(readr)
heatmap_fibro <- read_csv("27_34_micro_ex_raw.csv")
#View(heatmap_fibro)

data <- heatmap_fibro

common <- c(
  "5325", "7087", "2627", "10725", "9785", "7988", "23466", "506", "3694", "5016",
  "22884", "26043", "3418", "11054", "2760", "3624", "483", "26127", "7518", "8178",
  "6543", "24146", "4741", "4850", "682", "10376", "10632"
)




X <- data[data$ID %in% common,]

write.csv(X,file = "27_34_micro_P_1_sort.csv")


############################################## kera #####################################################3




setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/New_Idea_ML/RNA_SEQ/13_K_C")

# Load necessary libraries
library(tidyverse)
library(dplyr)

library(readr)
heatmap_fibro <- read_csv("Kera_normalized_data.csv")
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




X <- data[data$ID %in% common,]

#write.csv(X,file = "13_Colon_Correlation_expression.csv")