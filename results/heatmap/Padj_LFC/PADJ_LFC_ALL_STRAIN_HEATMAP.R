

setwd("D:/DWCT/DEG/GSE219036_6_june/Heatmap/Padj_LFC/DOWN_407")

################################################## Down ###########################################

## clade I
X <- read.csv("407_IIb_down_Heat_All_genes_Heatmap.csv", header = TRUE,row.names = 1)
Y<- read.csv("DOWN_Zaire_407_final.csv", header = TRUE,row.names = 1)



# Create a new data frame to store the results
result_df <- data.frame(RowNames = rownames(X), Y_Data = NA)

# Initialize a variable to keep track of the current row index in Y
y_row_index <- 1

# Iterate through the row names of X
for (i in 1:nrow(result_df)) {
  # Check if the row name of X matches the current row name of Y
  if (result_df$RowNames[i] == rownames(Y)[y_row_index]) {
    # If there's a match, extract the data from the second column of Y and move to the next row of Y
    result_df$Y_Data[i] <- Y[y_row_index, 2]
    y_row_index <- y_row_index + 1
  } else {
    # If there's no match, assign NA
    result_df$Y_Data[i] <- NA
  }
}


write.csv(result_df,file = "DOWN_ZAIRE_407_LFC.csv")




#########################################################

setwd("D:/DWCT/DEG/GSE219036_6_june/Heatmap/Padj_LFC/UP_396")


## clade I
X <- read.csv("396_Up_IIb_Heat_All_genes_Heatmap.csv", header = TRUE,row.names = 1)
Y<- read.csv("UP_Zaire_396_final.csv", header = TRUE,row.names = 1)



# Create a new data frame to store the results
result_df <- data.frame(RowNames = rownames(X), Y_Data = NA)

# Initialize a variable to keep track of the current row index in Y
y_row_index <- 1

# Iterate through the row names of X
for (i in 1:nrow(result_df)) {
  # Check if the row name of X matches the current row name of Y
  if (result_df$RowNames[i] == rownames(Y)[y_row_index]) {
    # If there's a match, extract the data from the second column of Y and move to the next row of Y
    result_df$Y_Data[i] <- Y[y_row_index, 2]
    y_row_index <- y_row_index + 1
  } else {
    # If there's no match, assign NA
    result_df$Y_Data[i] <- NA
  }
}


write.csv(result_df,file = "UP_IIa_396_LFC_Final.csv")


