library(DESeq2)
library(gProfileR)
library(knitr)
# extract differential expression results
DEresults <- results(dds, contrast = c("Treatment", "MPXV clade I infected", "mock"))

#remove genes with NA values 
DE <- DEresults[!is.na(DEresults$padj),]
#select genes with adjusted p-values below 0.1
DE <- DE[DE$padj < 0.1,]
#select genes with absolute log2 fold change above 1 (two-fold change)
DE <- DE[abs(DE$log2FoldChange) > 1,]

#get the list of genes of interest
genesOfInterest <- rownames(DE)

#calculate enriched GO terms
goResults <- gprofiler(query = genesOfInterest, 
                       organism = 'hsapiens', 
                       src_filter = 'GO', 
                       hier_filtering = 'moderate')



install.packages("gprofiler2")

# Load the gprofiler2 library
library(gprofiler2)

# Query GO terms for the list of genes
result <- gprofiler(organism = "hsapiens", 
                     query = genesOfInterest, 
                     #sources = "GO", 
                     significant = TRUE)


# View the results
print(result)