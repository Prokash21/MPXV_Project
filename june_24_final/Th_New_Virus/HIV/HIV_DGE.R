# script to perform differential gene expression analysis using DESeq2 package
# setwd("~/Desktop/demo/DESeq2_tutorial/data")
setwd("E:/DWCT/DEG/GSE219036_6_june/JUNE_24/Th_New_Virus/HIV")
# load libraries
library(DESeq2)
library(tidyverse)
#library(airway)

# Step 1: preparing count data ----------------

# read in counts data
counts_data_m <- read.csv('Count_matrix_HIV.csv' )
head(counts_data_m)
counts_data <- counts_data_m[-1]
head(counts_data)
rownames(counts_data) <- counts_data [,1]
head(counts_data)
counts_data <- counts_data[-1]
head(counts_data)

# read in sample info
colData <- read.csv('Meta_data_HIV.csv', header = TRUE, row.names = 1)


# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ treatment)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 7
dds <- dds[keep,]

dds

# set the factor level
dds$dexamethasone <- relevel(dds$treatment, ref = "healthy")

# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res



# Explore Results ----------------

summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)

# e.g.: treated_4hrs, treated_8hrs, untreated

results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))

# MA plot
plotMA(res)

########################################################################################################


#set factor levels
colData$treatment <- factor(colData$treatment)

# Convert non-integer values to integers in count data
counts_data <- round(counts_data)
head(counts_data)

# Create a new count data object
new_count_data <- as.matrix(counts_data)
head(new_count_data)


unique(colData$treatment)

# Generate the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = new_count_data, colData = colData, design = ~ treatment)

# Perform DESeq2 analysis
dds <- DESeq(dds)
head(dds)

#set the factor level
dds$treatment <- factor(dds$treatment, levels = c ("healthy","patient")) 


#filter the genes
keep <- rowSums(counts(dds)) >= 7
dds <- dds[keep,]
dds

#set the referene for the treatment factor
dds$treatment <- relevel(dds$treatment , ref = "healthy")
dds$treatment

#perform the statistical tests to identify differentialy expressed genes
dds <- DESeq(dds)
head(dds)

#save the normalized counts
normalize_counts <- counts(dds,normalized=TRUE)
raw_counts <- counts(dds,normalized=FALSE)
head(normalize_counts)
dim(normalize_counts)
write.csv(normalize_counts,"normalized_counts_HIV.csv")


#Identify available coefficient names
coeff_names <- resultsNames(dds)

#Print the coefficient names
print(coeff_names)

#[3] "treatment_patient_vs_healthy"

resLFC <- lfcShrink(dds, coef ="treatment_patient_vs_healthy"  , type = "apeglm")

#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)

res <- as.data.frame(res)  

find <- resLFC[ENSG00000118495,] 

write.csv(resLFC,"resLFC_HIV.csv")
