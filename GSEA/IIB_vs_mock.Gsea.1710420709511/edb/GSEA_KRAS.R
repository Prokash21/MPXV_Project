setwd("D:/DWCT/DEG/GSE219036_6_june/GSEA/IIB_vs_mock.Gsea.1710420709511/edb")
X<-read.csv("KRAS_list.csv")
Y<- read.csv("Gene_396.csv")
XX <- X$SYMBOL
vv<- Y [ Y$Gene_ID %in% XX,]
View(vv)



Y<- read.csv("407_down_Gene_names.csv")
XX <- X$SYMBOL
vv<- Y [ Y$ID %in% XX,]
View(vv)
