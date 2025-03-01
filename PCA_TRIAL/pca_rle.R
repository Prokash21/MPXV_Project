# Boxplot
#boxplot(count_data, outline = FALSE,
#   main="Boxplot of count data")
XXX<-read.csv("ex_Fibro_34.csv", header = TRUE, row.names = 1)
XXX<-read.csv("ex1_Fibro_34.csv", header = TRUE, row.names = 1)
XXX<-read.csv("exx_Fibro_34.csv", header = TRUE, row.names = 1)
XXX<-read.csv("exx1_Fibro_34.csv", header = TRUE, row.names = 1)



# Log2 transformation for count data
count_matrix <- XXX  # Adding 1 to avoid log(0)
log2_count_matrix <- log2(count_matrix)
boxplot(XXX, outline = FALSE, main = "Boxplot of Log2-transformed Count Data",
       
        ylab = "Log2-transformed Counts")

# Adding sample names with 90-degree rotation
text(x = 1:ncol(log2_count_matrix), y = par("usr")[3] - 1, labels = colnames(log2_count_matrix),
     srt = 90, adj = 1, xpd = TRUE, cex = 0.7)

# Log2 transformation for normalized count data
normalized_counts <- counts(dds, normalized = TRUE)
log2_normalized_counts <- log2(normalized_counts + 1)  # Adding 1 to avoid log(0)
#colnames(log2_normalized_counts) <- custom_col_names
boxplot(log2_normalized_counts, outline = FALSE,
        main = "Boxplot of Log2-transformed Normalized Count Data",
        xlab = "Sample Name",
        ylab = "Log2-transformed Counts")


library(EDASeq)
#par(mfrow = c(2, 1))
plotRLE(CC), 
        outline=FALSE, ylim=c(-2, 2), 
        main = 'Raw counts',
        xlab = "Sample Name",
        ylab = "LFC")

plotRLE(DESeq2::counts(dds, normalized = TRUE), 
        outline=FALSE, ylim=c(-2, 2), 
        col = as.numeric(sample_info$Treatment), 
        main = 'Normalized Counts (DESeq2)',
        xlab = "Sample Name",
        ylab = "LFC")









install.packages("scater")
library(scater)

BiocManager::install("scater")

#example
example_sce <- mockSCE()
example_sce <- logNormCounts(example_sce)

plotRLE(example_sce, colour_by = "Mutation_Status", style = "minimal")

plotRLE(example_sce, colour_by = "Mutation_Status", style = "full",
        outlier.alpha = 0.1, outlier.shape = 3, outlier.size = 0)




























