library(DESeq2)
install.packages("pheatmap")
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
install.packages("ggrepel")
library(ggrepel)

#set the working directory
setwd("D:/DEG/GSE219036_6_june")

#load the count data

############## type 1
library(readr)
merged <- read_csv("merged.csv")
View(merged)
count_data <- read.csv("merged.csv", header=TRUE,row.names = 1)
#############

# type 2
#$$$$$$$$ it gives 25 col so dont use it
count_data <- read.csv("merged.csv", header = TRUE)
#$$$$$$$

row.names(count_data) <- NULL

colnames(count_data)
head(count_data)

#load the sample info

sample_info <- read.csv("mata_data_serial.csv",header = TRUE,row.names = 1)
colnames(sample_info)
head(sample_info)
#set actor levels

sample_info$Treatment <- factor(sample_info$Treatment)
sample_info$DATASTORE.region <- factor(sample_info$DATASTORE.region)


#create a DEseq object and import the count data and sample information

###############################################################
#****************# dont req for type 1 so skip it 
# Exclude the "gene_id" column 
numeric_columns <- colnames(count_data)[-1]

# Round the numeric columns
count_data[, numeric_columns] <- round(count_data[, numeric_columns])

#****************#

# Convert non-integer values to integers in count data
count_data <- round(count_data)

# Create a new count data object
new_count_data <- as.matrix(count_data)

# Generate the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = new_count_data, colData = sample_info, design = ~ DATASTORE.region + Treatment)

# Perform DESeq2 analysis
dds <- DESeq(dds)
##############################################################
#dds <- DESeqDataSetFromMatrix(countData = count_data,colData = sample_info,design = ~DATASTORE.region + Treatment)


#set the referene for the treatment factor

dds$Treatment <- factor(dds$Treatment, levels = c ("mock","MPXV clade I infected", "MPXV clade IIb infected", "MPXV clade IIa infected")) 


#filter the genes

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]


#perform the statistical tests to identify differentialy expressed genes
dds <- DESeq(dds)

#######################################################################
# Set the Treatment factor levels
dds$Treatment <- factor(dds$Treatment, levels = c("mock", "MPXV clade I infected", "MPXV clade IIb infected", "MPXV clade IIa infected"))

# Filter low-count rows
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep, ]

# Perform DESeq analysis
dds <- DESeq(dds, maxit = 1000)  # Increase maxit argument to allow for more iterations


#######################################################################

deseq_results <- results(dds)
deseq_results


#Change DEseq object to R object(dataframe)

deseq_results <- as.data.frame(deseq_results)
class(deseq_results)

head(deseq_results)
#order the result table by increaing p value

deseq_results_ordered <- deseq_results[order(deseq_results$pvalue),]
head(deseq_results_ordered)


#is gene differentially expressed?
deseq_results["ENSG00000000971",]


#Is the ENSG00000000971  downregulated by the RNAi treatment?


#Extract the most differentially expressed genes due to the treatment.
#Select genes with a significant change in gene expression (adjusted p-value below 0.05)
#And log2fold change <1 and >1



#Step 1: filter based on p adjusted value
filtered <- deseq_results %>% filter(deseq_results$padj < 0.05)



#Step 2: filter based on fold changes. here we will use a threshold of 1
filtered <- filtered %>% filter(abs(filtered$log2FoldChange)>1)
dim(deseq_results)
dim(filtered)



# make queries





#Save the deseg result. We will save the both the original data(res) and the filtered one(hits)

write.csv(deseq_results,"de_result.all.csv")
write.csv(filtered,"de_result.filtered.csv")


#save the normalized counts
normalize_counts <- counts(dds,normalized=TRUE)
head(normalize_counts)
write.csv(normalize_counts,"normalized_counts.csv")


#visualization


#dispersion plot

plotDispEsts(dds)




#PCA
#pca stands for principal component analysis.
#it is a dimensionality reduction technique and in gene expression analysis, 139 #it can be used to explain the variance in gene expression datasets.
#to generate the pa plot we will first perform a variance stabilizing transformation.
#we will use the vst function in deseg.
#after that we use the transformed values to plot the PCA




#variance transformed values to generate a pca plot 
vsd <- vst (dds,blind = FALSE)


#use transformed values to generate a pcs plot 
plotPCA(vsd,intgroup=c("DATASTORE.region", "Treatment"))



#Heatmaps

#R Package: pheatmap
#Heatmap of sample-to-sample distance matrix (with clustering) based on the normalized counts.

#generate the distance matrix
sampledists <- dist (t (assay(vsd)))
sampleDistMatrix <- as.matrix(sampledists)
colnames(sampleDistMatrix)

#set color scheme
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
#generate the heatmap
pheatmap(sampleDistMatrix,clustering_distance_rows=sampledists,
          clustering_distance_cols = sampledists, col=colors)
#Heatmap of log transformed normalized counts. We will use the top 10 genes.
#Heatmap of Z scores. We will use the top 10 genes



# Heatmap of log transformed normalized counts . we will use the top 10 genes

#top 10 genes
top_hits <- deseq_results[order(deseq_results$padj),][1:20,]
top_hits <-row.names(top_hits)
top_hits
'"ENSG00000197061" "ENSG00000260996" "ENSG00000237338" "ENSG00000278463" "ENSG00000266472" "ENSG00000270022"
"ENSG00000287001" "ENSG00000120885" "ENSG00000128591" "ENSG00000088832"'''

rld <- rlog(dds,blind = FALSE)

pheatmap(assay(rld)[top_hits,], cluster_rows = FALSE, show_rownames = TRUE,cluster_cols = FALSE)
pheatmap(assay(rld)[top_hits,])


annot_info <- as.data.frame(colData(dds)[,c("DATASTORE.region", "Treatment")])
pheatmap(assay(rld)[top_hits,],cluster_rows = FALSE,show_rownames = TRUE,cluster_cols = FALSE,
         annotation_col = annot_info)


#heatmap of Z scores . we will use the top 10 genes
cal_z_score <- function(x){(x-mean(x)) / sd(x)}


zscore_all <- t(apply(normalize_counts, 1, cal_z_score))
zscore_subset <- zscore_all [top_hits,]
pheatmap(zscore_subset)



#MA plot


plotMA(dds,ylim = c (-2,2))


#remove the noise 
resLFC <- lfcShrink(dds ,coef = "Treatment_mock_vs_MPXV clade I infected", type="apeglm")

##################
# Identify available coefficient names
coeff_names <- resultsNames(dds)

# Print the coefficient names
print(coeff_names)

# Remove the noise using correct coefficient name
BiocManager::install("apeglm")
library(apeglm)


"Intercept"                                                                        
"DATASTORE.region_gs.US.s3.us.east.1.ncbi.public_vs_gs.US.ncbi.public.s3.us.east.1"
"DATASTORE.region_ncbi.public.gs.US.s3.us.east.1_vs_gs.US.ncbi.public.s3.us.east.1"
"DATASTORE.region_ncbi.public.s3.us.east.1.gs.US_vs_gs.US.ncbi.public.s3.us.east.1"
"DATASTORE.region_s3.us.east.1.gs.US.ncbi.public_vs_gs.US.ncbi.public.s3.us.east.1"
"DATASTORE.region_s3.us.east.1.ncbi.public.gs.US_vs_gs.US.ncbi.public.s3.us.east.1"
"Treatment_MPXV.clade.I.infected_vs_mock"                                          
"Treatment_MPXV.c.infected_vs_mock"                                        
"Treatment_MPXV.clade.IIalade.IIb.infected_vs_mock"    


resLFC <- lfcShrink(dds, coef = "Treatment_MPXV.clade.I.infected_vs_mock", type = "apeglm")

plotMA(resLFC,ylim=c(-2,2))

########################


#Volcano Plot

#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)
#label the genes 


resLFC$diffexpressed <- "NO"

resLFC$diffexpressed[resLFC$log2FoldChange>0.1 & resLFC$padj<0.05]<- "UP"

resLFC$diffexpressed[resLFC$log2FoldChange<0.1 & resLFC$padj<0.05]<- "DOWN"

resLFC$delabel<-NA


ggplot(data=resLFC,aes(x=log2FoldChange,y =-log10(pvalue),col=diffexpressed,label=delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values=c("blue","black","red"))+
  theme(text = element_text(size = 20))



























































