
setwd("D:/DWCT/DEG/GSE219036_6_june/Heatmap/Padj_LFC")

X<-read.csv("396_up.csv")
Y<-read.csv("LFC_Fibroblast.csv", header=TRUE,row.names = 1)


y<-Y[X$converted_alias,]
x<-X$initial_alias
rownames(y) <- x
# Ensure uniqueness of values in X$initial_alias
unique_x <- make.unique(as.character(X$initial_alias))

# Assign unique values as row names to y
rownames(y) <- unique_x

write.csv(y,file="UP_Zaire_396.csv")
#########################################################################################
X<-read.csv("470_down.csv")
Y<-read.csv("LFC_Fibroblast.csv", header=TRUE,row.names = 1)


y<-Y[X$converted_alias,]
x<-X$initial_alias
rownames(y) <- x
# Ensure uniqueness of values in X$initial_alias
unique_x <- make.unique(as.character(X$initial_alias))

# Assign unique values as row names to y
rownames(y) <- unique_x

write.csv(y,file="DOWN_Zaire_407.csv")
