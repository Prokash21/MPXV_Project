# Install and load the clusterProfiler package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

library(clusterProfiler)
#

# Load your gene list
my_gene_list <- c("WDR37", "NFAT5", "ICAM5", "ICAM5",
                  "ICAM5", "ELL", "PLAGL1", "DHX38", "GATA6", "ZNF212", "CBX6")  # Replace with your gene symbols

# Load the gene-to-GO mapping data (use the appropriate organism database)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# Perform GO enrichment analysis (Biological Process)
go_enrichment <- enrichGO(
  gene          = my_gene_list,
  OrgDb         = org.Hs.eg.db,  # Replace with the appropriate organism database
  keyType       = "SYMBOL",
  ont           = "CC",  # BP for Biological Process, MF for Molecular Function, CC for Cellular Component
  pAdjustMethod = "fdr",  # Adjust p-values for multiple testing
  pvalueCutoff  = 0.05
)

# View the results
summary(go_enrichment)
# View the results using as.data.frame
result_df <- as.data.frame(go_enrichment)
# View the results
print(result_df)
resul

go_enrichment@result[["Description"]]
x<-go_enrichment@result
