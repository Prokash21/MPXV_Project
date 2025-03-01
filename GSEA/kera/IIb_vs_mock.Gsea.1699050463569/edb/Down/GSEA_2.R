setwd("D:/DWCT/DEG/GSE219036_6_june/GSEA/kera/IIb_vs_mock.Gsea.1699050463569/edb/Down")
X<- read.csv("gene_sets_GSEA_IIb.csv")
head(X)
Y <- read.csv("248_genes.csv")
head(Y)



######################################
# Initialize an empty list to store results
results <- list()

# Loop through each row of X
for (i in 1:nrow(X)) {
  # Get the genes present in the current row of X
  genes_in_row <- X[i, -1]  # Exclude the first column which contains row names
  
  # Check which genes in Y are present in the current row of X
  matching_genes <- Y$Gene_ID[Y$Gene_ID %in% genes_in_row]
  
  # Store the matching genes in the results list
  results[[rownames(X)[i]]] <- matching_genes
}

# Convert the results list into a data frame
results_df <- data.frame(
  Row_Names = names(results),
  Genes_Present = sapply(results, function(x) paste(x, collapse = ", "))
)

# Print the results
#print(results_df)
head(results_df)
write.csv(results_df,file = "Matched_248.csv")
