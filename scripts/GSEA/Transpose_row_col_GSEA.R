

setwd("D:/DWCT/DEG/GSE219036_6_june/GSEA/IIB_vs_mock.Gsea.1710420709511/edb")


# Read the CSV file into a data frame
data <- read.csv("gene_sets.gmt.csv")

# Transpose the data frame
transposed_data <- t(data)

# Write the transposed data to a new CSV file
write.csv(transposed_data, file = "transposed_data.csv")
