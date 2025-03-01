library(GOplot)
install.packages("GOplot")

setwd("E:/DWCT/DEG/GSE219036_6_june/GO/407_down")
X<-read.csv("BP_RAW.csv")
XXX<-data(EC)
as.data.frame(XXX)

T<- EC$eset
V<- EC$genelist
P<-EC$david
B<- EC$genes

head(EC$david)
head(EC$genelist)
circ <- circle_dat(EC$david, EC$genelist)
GOBar(subset(circ, category == 'BP'))
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)  
GOBubble(circ, labels = 3)



circ1<-read.csv("Down_407_go_all_sorted1.csv")
GOBubble(circ1, labels = 3)
#GOBar(subset(circ1, Category == 'BP'))
GOBubble(circ1, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)  


# Colour the background according to the category
GOBubble(circ1, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)  

# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ1, overlap = 0.75)
# ...and plot it
GOBubble(reduced_circ, labels = 2.8)

###################################################################

XX<-read.csv("david_407.csv")
YY<-read.csv("genelist_Bubble.csv")
XX<-as.data.frame(XX)
YY<-as.data.frame(YY)

circ <- circle_dat(XX, YY)
write.csv(circ,file = "407_Down_david_file.csv")
GOBubble(circ, labels = 3)

# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# ...and plot it
GOBubble(reduced_circ, labels = 2.8)

####UP_396




XX<-read.csv("UP_396_david.csv")
YY<-read.csv("genelist_Bubble.csv")
XX<-as.data.frame(XX)
YY<-as.data.frame(YY)

circ <- circle_dat(XX, YY)
write.csv(circ,file = "396_Up_david_file.csv")

GOBubble(circ, labels = 3)

# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# ...and plot it
GOBubble(reduced_circ, labels = 2.8)


# Colour the background according to the category
GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)  


# Add a title, change the colour of the circles, facet the plot according to the categories and change the label threshold
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)  
