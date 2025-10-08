
################################################ Zaire #####################################################
setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Heatmap")

# Load necessary libraries
library(tidyverse)
library(dplyr)

library(readr)
heatmap_fibro <- read_csv("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Heatmap/heatmap_fibro.csv")
View(heatmap_fibro)

data <- heatmap_fibro


common <- data$common[1:27]

# Filtering rows where gene_name is in common
filtered_data <- data[data$gene_name %in% common, ]
filtered_id <- filtered_data$gene_id
gene_id <- data[data$ID %in% filtered_id,]


write.csv(gene_id[1:2],file = "LFC_27.csv")


lfc <- gene_id[1:2]
X <- filtered_data[4:5]

merged_data <- inner_join(lfc, X, by = c("ID" = "gene_id"))

write.csv(merged_data, file = "LFC_27_Firbro_merged_data.csv")


####################################################### IIb ######################################################

setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Heatmap")

# Load necessary libraries
library(tidyverse)
library(dplyr)

library(readr)
heatmap_fibro <- read_csv("resLFC_I_p_cut_clean.csv")
View(heatmap_fibro)

data <- heatmap_fibro

common <- data$common[1:27]

# Filtering rows where gene_name is in common
top_geneid <- read.csv("ensemble_geneid.csv")

id <- top_geneid[top_geneid$gene_name %in% common,]

ID <- id$gene_id

LFC <- data[data$ID %in% ID,]

LFC <- LFC[,c(1,3)]

merged <- inner_join(LFC,id,by =c("ID" = "gene_id"))

write.csv(merged, file = "LFC_27_I_merged_data.csv")


#################################################### IIa ##############################################


setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Heatmap")

# Load necessary libraries
library(tidyverse)
library(dplyr)
####################################################################
heatmap_fibro <- read_csv("resLFC_I_p_cut_clean.csv")
View(heatmap_fibro)

data <- heatmap_fibro

common <- data$common[1:27]
############################################################
library(readr)
heatmap_fibro <- read_csv("resLFC_IIa_p_cut_clean.csv")
View(heatmap_fibro)

data <- heatmap_fibro

#common <- data$common[1:27]

# Filtering rows where gene_name is in common
top_geneid <- read.csv("ensemble_geneid.csv")

id <- top_geneid[top_geneid$gene_name %in% common,]

ID <- id$gene_id

LFC <- data[data$ID %in% ID,]

LFC <- LFC[,c(1,3)]

merged <- inner_join(LFC,id,by =c("ID" = "gene_id"))

write.csv(merged, file = "LFC_27_IIa_merged_data.csv")


############################################### IIb ###########################################


setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Heatmap")

# Load necessary libraries
library(tidyverse)
library(dplyr)
####################################################################
heatmap_fibro <- read_csv("resLFC_I_p_cut_clean.csv")
#View(heatmap_fibro)

data <- heatmap_fibro

common <- data$common[1:27]
############################################################
library(readr)
heatmap_fibro <- read_csv("resLFC_IIb_p_cut_clean.csv")
View(heatmap_fibro)

data <- heatmap_fibro

#common <- data$common[1:27]

# Filtering rows where gene_name is in common
top_geneid <- read.csv("ensemble_geneid.csv")

id <- top_geneid[top_geneid$gene_name %in% common,]

ID <- id$gene_id

LFC <- data[data$ID %in% ID,]

LFC <- LFC[,c(1,3)]

merged <- inner_join(LFC,id,by =c("ID" = "gene_id"))

write.csv(merged, file = "LFC_27_IIb_merged_data.csv")



###########################################################################################################
#*************************************************Old***********************************************************#
###########################################################################################################



################################################ Zaire #####################################################
setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Heatmap")

# Load necessary libraries
library(tidyverse)
library(dplyr)

library(readr)
heatmap_fibro <- read_csv("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Heatmap/heatmap_fibro.csv")
View(heatmap_fibro)

data <- heatmap_fibro


common <- data$gene_name_old[1:14]

# Filtering rows where gene_name is in common
filtered_data <- data[data$gene_name %in% common, ]
filtered_id <- filtered_data$gene_id
gene_id <- data[data$ID %in% filtered_id,]


#write.csv(gene_id[1:2],file = "LFC_27.csv")


lfc <- gene_id[1:2]
X <- filtered_data[4:5]

merged_data <- inner_join(lfc, X, by = c("ID" = "gene_id"))

write.csv(merged_data, file = "LFC_27_Firbro_merged_data.csv")


####################################################### I ######################################################

setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Heatmap")

# Load necessary libraries
library(tidyverse)
library(dplyr)

library(readr)
heatmap_fibro <- read_csv("resLFC_I_p_cut_clean.csv")
View(heatmap_fibro)

data <- heatmap_fibro

#common <- data$common[1:27]

# Filtering rows where gene_name is in common
top_geneid <- read.csv("ensemble_geneid.csv")

id <- top_geneid[top_geneid$gene_name %in% common,]

ID <- id$gene_id

LFC <- data[data$ID %in% ID,]

LFC <- LFC[,c(1,3)]

merged <- inner_join(LFC,id,by =c("ID" = "gene_id"))

write.csv(merged, file = "LFC_27_I_merged_data.csv")


#################################################### IIa ##############################################


setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Heatmap")

# Load necessary libraries
library(tidyverse)
library(dplyr)
####################################################################
heatmap_fibro <- read_csv("resLFC_I_p_cut_clean.csv")
View(heatmap_fibro)

data <- heatmap_fibro

common <- data$common[1:27]
############################################################
library(readr)
heatmap_fibro <- read_csv("resLFC_IIa_p_cut_clean.csv")
View(heatmap_fibro)

data <- heatmap_fibro

#common <- data$common[1:27]

# Filtering rows where gene_name is in common
top_geneid <- read.csv("ensemble_geneid.csv")

id <- top_geneid[top_geneid$gene_name %in% common,]

ID <- id$gene_id

LFC <- data[data$ID %in% ID,]

LFC <- LFC[,c(1,3)]

merged <- inner_join(LFC,id,by =c("ID" = "gene_id"))

write.csv(merged, file = "LFC_27_IIa_merged_data.csv")


############################################### IIb ###########################################


setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Heatmap")

# Load necessary libraries
library(tidyverse)
library(dplyr)
####################################################################
heatmap_fibro <- read_csv("resLFC_I_p_cut_clean.csv")
#View(heatmap_fibro)

data <- heatmap_fibro

common <- data$common[1:27]
############################################################
library(readr)
heatmap_fibro <- read_csv("resLFC_IIb_p_cut_clean.csv")
View(heatmap_fibro)

data <- heatmap_fibro

#common <- data$common[1:27]

# Filtering rows where gene_name is in common
top_geneid <- read.csv("ensemble_geneid.csv")

id <- top_geneid[top_geneid$gene_name %in% common,]

ID <- id$gene_id

LFC <- data[data$ID %in% ID,]

LFC <- LFC[,c(1,3)]

merged <- inner_join(LFC,id,by =c("ID" = "gene_id"))

write.csv(merged, file = "LFC_27_IIb_merged_data.csv")

########################***************************************########################################
########################*    27 gene expression
########################*Fibro




setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Correlation")

# Load necessary libraries
library(tidyverse)
library(dplyr)

library(readr)
heatmap_fibro <- read_csv("Expression_GSE11234_Fibroblast.csv")
View(heatmap_fibro)

data <- heatmap_fibro

common <- c(5325, 7087, 2627, 10725, 9785, 7988,
            23466, 506, 3694, 5016, 22884, 26043,
            3418, 11054, 2760, 3624, 483, 26127, 7518,
            8178, 6543, 24146, 4741, 4850, 682, 10376, 10632)



X <- data[data$...1 %in% common,]

write.csv(X,file = "27_Firbro_expression.csv")



########################****************** Kera *********************########################################





setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Correlation")

# Load necessary libraries
library(tidyverse)
library(dplyr)

library(readr)
heatmap_fibro <- read_csv("normalized_counts_keratinocyte_I_vs_mock.csv")
View(heatmap_fibro)

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

