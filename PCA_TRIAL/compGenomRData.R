library(EDASeq)
# remove 'width' column from counts
countData <- as.matrix(count_data(counts, select = c(-width)))
# create a seqExpressionSet object using EDASeq package 
set <- newSeqExpressionSet(counts = count_data1,
                           phenoData = sample_info1)


devtools::install_github("compgenomr/compGenomRData")

library(compGenomRData)

BiocManager::install("compGenomRData")

library(compGenomRData)
install.packages("compGenomRData")




counts_file <- system.file('extdata/rna-seq/SRP049988.raw_counts.tsv', 
                           package = 'compGenomRData')
colData_file <- system.file('extdata/rna-seq/SRP049988.colData.tsv', 
                            package = 'compGenomRData')

counts <- read.table(counts_file)
colData <- read.table(colData_file, header = T, 
                      sep = '\t', stringsAsFactors = TRUE)
# simplify condition descriptions
colData$source_name <- ifelse(colData$group == 'CASE', 
                              'EHF_overexpression', 'Empty_Vector')