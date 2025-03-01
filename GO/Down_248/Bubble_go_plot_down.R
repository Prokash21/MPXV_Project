library(GOplot)
#install.packages("GOplot")

setwd("D:/DWCT/DEG/GSE219036_6_june/GO/407_down")
X<-read.csv("BP_RAW.csv")
data(EC)
head(EC$david)
head(EC$genelist)
circ <- circle_dat(EC$david, EC$genelist)
GOBar(subset(circ, category == 'BP'))
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)  



circ1<-read.csv("Down_407_go_all_sorted1.csv")
#GOBar(subset(circ1, Category == 'BP'))
GOBubble(circ1, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)  


# Colour the background according to the category
GOBubble(circ1, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)  

# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ1, overlap = 0.75)
# ...and plot it
GOBubble(reduced_circ, labels = 2.8)
