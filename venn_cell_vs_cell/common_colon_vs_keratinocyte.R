################################################################################

#set the working directory
setwd("D:/DEG/GSE219036_6_june/venn_cell_vs_cell")

# Load your DESeq2 result files
common_colon <- read.csv("common_genes_colon.csv")
common_kera <- read.csv("common_genes_keratinocyte.csv")


# Extract upregulated genes from each comparison
common_colon1 <- common_colon$Gene_Name
common_kera1 <- common_kera$Gene_Name



# Create a list of upregulated genes from each comparison
upregulated_gene_lists <- list(colon = common_colon1,
                               keratinocyte = common_kera1)

# Define colors for the Venn diagram
colors <- c("dodgerblue", "darkorange")

# Create a Venn diagram
venn.plot <- venn.diagram(
  x = upregulated_gene_lists,
  category.names = c("colon_organoids", "keratinocyte"),
  filename = NULL,
  
)


library(gplots)
venn(upregulated_gene_lists)

# Plot the Venn diagram
pdf("venn_diagram.pdf")
grid.draw(venn.plot)
dev.off()

# Find common genes among the three comparisons
common_genes <- Reduce(intersect, list(common_colon1, common_kera1))
# Create a data frame with common genes
common_genes_df <- data.frame(Gene_Name = common_genes)
# Write the common genes to a CSV file
write.csv(common_genes_df, "common_colon_vs_keratinocyte.csv", row.names = FALSE)



