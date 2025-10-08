install.packages("stats")
install.packages("dplyr")


library(stats)
library(dplyr)


mydata = count_data
cor(mydata)
mean(cor(mydata))


pca = princomp(mydata)
pca$loadings



pc=pca$scores
view(pc)
print(pc)
cor(pc)
# Plot the principal component scores for the first two components
plot(pc[, 1], pc[, 2], xlab = "PC1", ylab = "PC2", main = "Principal Component Scores")
var1<-round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
var2<-round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)

plot(pca)
library("factoextra")
fviz_eig(pca)
biplot(pca)

install.packages("ggfortify")
library(ggfortify)

library(ggfortify)
autoplot(pca)


pc<-prcomp(count_data,center = TRUE,scale. = TRUE)
summary(pc)
attributes(pc)
print(pc)
var1<- round(pc$sdev[1]^2/sum(pc$sdev^2)*100,2)
var2<- round(pc$sdev[2]^2/sum(pc$sdev^2)*100,2)

library("ggfortify")
autoplot(pc)

autoplot(pc,scale = 0)




fviz_pca_var(pc,
             col.var="")
install.packages("factoextra")
library(factoextra)

fviz_pca_var(pc)
fvi
fviz_pca_var(pc)
fviz_pca_biplot(pc)
fviz_pca_contrib(pc)
fviz_contrib(pc)
fviz_eig(pca)
biplot(pc)
autoplot(pc)
