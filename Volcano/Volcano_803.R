library(org.Hs.eg.db)
res.df <- as.data.frame(resLFC)
res.df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res.df),keytype = "ENSEMBL", column = "SYMBOL")
res.df

BiocManager::install('EnhancedVolcano')
install.packages("textshaping")
library(EnhancedVolcano)
selected <- XX$Gene
selected <- XX$ID
exclusive<- read.csv("Volcano_396_407_21_3_24.csv")
res.df1 <- as.data.frame(exclusive) 
EnhancedVolcano(res.df ,x= "log2FoldChange",y="padj", lab = row.names(res.df))

EnhancedVolcano(res.df ,x= "log2FoldChange",y="padj", lab = res.df$symbol,
                pCutoff = 0.05,FCcutoff = 1, selectLab = exclusive$X)



EnhancedVolcano(res.df1 ,x= "log2FoldChange",y="padj", lab = res.df1$X,
                pCutoff = 0.05,FCcutoff = 1)









EnhancedVolcano(res.df, x = "log2FoldChange", y = "padj", lab = res.df$symbol,
                pCutoff = 0.05, FCcutoff = 1, selectLab = exclusive$X,
                col = ifelse(res.df$log2FoldChange < -1 & res.df$padj < 0.05, "blue"))
