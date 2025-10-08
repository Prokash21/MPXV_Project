


# Set working directory
setwd("E:/DWCT/DEG/GSE219036_6_june/Drug")

#######################################################################################
library(cowplot)
library(ggplot2)
library(dplyr)

# Read the DSigDB_table.csv file
dsigdb <- read.csv("DSigDB_table_U.csv", header = TRUE, row.names = 1)

# Filter data based on selection criteria (p < 0.05 and combined score > 100)
filtered <- dsigdb %>%
  filter(P.value < 0.05 & Combined.Score > 100)

# Optionally, sort the data for better visualization
filtered <- filtered %>%
  arrange(desc(Combined.Score))

# Add a rank column for bubble plot size differentiation
filtered$Rank <- rank(-filtered$Combined.Score)

# Select the top 10 significant results
top10 <- head(filtered, 10)

# Sort top10 by -log10(P.value)
top10 <- top10 %>%
  arrange(desc(-log10(P.value)))

# Create a new rank column based on the sorted order by -log10(P.value)
top10$Rank <- 1:nrow(top10)

# Add rank to the drug names for y-axis labels and reverse the order
top10$Drug <- factor(paste0(top10$Rank, ": ", rownames(top10)), levels = rev(paste0(top10$Rank, ": ", rownames(top10))))

# Plot 1: Bubble plot with rownames(top10) on Y axis and log(p-value) on X axis
bubble_plot <- ggplot(top10, aes(x = -log10(P.value), y = Drug, size = Combined.Score)) +
  geom_point(alpha = 0.7) +
  scale_size_area(max_size = 10) +
  theme_minimal() +
  labs(title = "Bubble Plot: Top 10 Results",
       x = "-log10(P-value)",
       y = "Drugs")

# Plot 2: Bar plot with rownames(top10) on Y axis and combined score on X axis
bar_plot <- ggplot(top10, aes(x = Combined.Score, y = Drug)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Bar Plot: Top 10 Results",
       x = "Combined Score",
       y = "Drugs")

# Combine the two plots into one figure
combined_plot <- plot_grid(bubble_plot, bar_plot, labels = "AUTO", ncol = 1)

# Print the combined plot
print(combined_plot)

#######################################################################################


# Load required libraries
library(ggplot2)
library(dplyr)
library(cowplot)

# Read the DSigDB_table.csv file
dsigdb <- read.csv("DSigDB_table_U.csv", header = TRUE, row.names = 1)

# Filter data based on selection criteria (p < 0.05 and combined score > 100)
filtered <- dsigdb %>%
  filter(P.value < 0.05 & Combined.Score > 100)

# Optionally, sort the data for better visualization
filtered <- filtered %>%
  arrange(desc(Combined.Score))

# Add a rank column for bubble plot size differentiation
filtered$Rank <- rank(-filtered$Combined.Score)

# Print the first few rows of filtered data to check
head(filtered)

# Select the top 10 significant results
top10 <- head(filtered, 10)

# Print the top 10 rows of filtered data to check
head(top10)

# Plot 1: Bubble plot with rownames(top10) on Y axis and log(p-value) on X axis
bubble_plot <- ggplot(top10, aes(x = -log10(P.value), y = rownames(top10), size = Combined.Score)) +
  geom_point(alpha = 0.7) +
  scale_size_area(max_size = 10) +
  theme_minimal() +
  labs(title = "Bubble Plot: Top 10 Results",
       x = "-log10(P-value)",
       y = "Drugs")

# Plot 2: Bar plot with rownames(top10) on Y axis and combined score on X axis
bar_plot <- ggplot(top10, aes(x = Combined.Score, y = rownames(top10))) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Bar Plot: Top 10 Results",
       x = "Combined Score",
       y = "Drugs")

# Combine the two plots into one figure
combined_plot <- plot_grid(bubble_plot, bar_plot, labels = "AUTO", ncol = 1)

# Print the combined plot
print(combined_plot)


#***********************************************************************************


setwd("E:/DWCT/DEG/GSE219036_6_june/Drug")

# Load required libraries
library(ggplot2)
library(dplyr)

# Read the DSigDB_table.txt file (replace with your actual file path)
dsigdb <- read.csv("DSigDB_table_U.csv",header = TRUE,row.names = 1)

# Filter data based on your selection criteria (p < 0.05 and combined score > 100)
filtered <- dsigdb %>%
  filter(P.value < 0.05 & Combined.Score > 100)

# Optionally, you can sort the data for better visualization
filtered <- filtered %>%
  arrange(desc(Combined.Score))

# Add a rank column for bubble plot size differentiation (optional but recommended)
filtered$Rank <- rank(-filtered$Combined.Score)

# Print the first few rows of filtered data to check
head(filtered)

# Assuming 'filtered' is already filtered DSigDB data
top10 <- head(filtered, 10)  # Selecting the top 10 significant results

# Print the top 10 rows of filtered data to check
head(top10)



# Assuming 'top10' is your dataframe with Drug names, Combined Score, and P-values

# Load necessary packages
library(ggplot2)
library(scales)



# Log10 transformation of P-values
top10$logP <- -log10(top10$P.value)

# Order top10 by logP in descending order
top10 <- top10[order(top10$logP, decreasing = TRUE), ]

# Assign sequential row names as factor in correct order
top10$Y <- factor(top10$row.names(top10), levels = rev(top10$row.names(top10)))  # Factor levels in reverse order



# Plotting
ggplot(top10, aes(x = Combined.Score, y = Y)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_point(aes(size = logP), color = "blue", alpha = 0.6) +
  scale_size_continuous(range = c(2, 10)) +  # Size range for bubble sizes
  labs(title = "Combined Scores and -log10(P-Values) for Drugs",
       x = "Combined Score", y = "Drug Names", size = "-log10(P-Value)") +
  theme_minimal()

# Log10 transformation of P-values and sorting
top10$logP <- -log10(top10$P.value)
top10 <- top10[order(top10$logP, decreasing = TRUE), ]

top10$Y <- row.names(top10)



# Assign sequential row names
top10$Y <- factor(row.names(top10), levels = row.names(top10))  # Ensure Y is a factor for correct ordering

# Plotting
ggplot(top10, aes(x = Combined.Score, y = Y)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_point(aes(size = logP), color = "blue", alpha = 0.6) +
  scale_size_continuous(range = c(2, 10)) +  # Size range for bubble sizes
  labs(title = "Combined Scores and -log10(P-Values) for Drugs",
       x = "Combined Score", y = "Drug Names", size = "-log10(P-Value)") +
  theme_minimal()
# Plotting
ggplot(top10, aes(x = Combined.Score, y = Y)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_point(aes(size = logP), color = "blue", alpha = 0.6) +
  scale_size_continuous(range = c(2, 10)) +  # Size range for bubble sizes
  labs(title = "Combined Scores and -log10(P-Values) for Drugs",
       x = "Combined Score", y = "Drug Names", size = "-log10(P-Value)") +
  theme_minimal()








# Bubble plot for top 10 significant results
ggplot(top10, aes(x = P.value, y = Combined.Score, size = Combined.Score, color = P.value)) +
  geom_point(alpha = 0.7) +  # Adjust alpha for transparency
  scale_size_continuous(range = c(3, 10)) +  # Adjust bubble size range as needed
  scale_color_gradient(low = "blue", high = "red") +  # Adjust color gradient
  labs(x = "P Value", y = "Combined Score", title = "Top 10 Significant Pharmacological Molecules") +
  theme_minimal() +  # Customize theme if needed
  theme(legend.position = "none")  # Remove legend if not needed

# Bar plot for top 10 significant results
ggplot(top10, aes(x = row.names(top10), y = Combined.Score, fill = P.value)) +
  geom_bar(stat = "identity") +  # Create bar plot
  scale_fill_gradient(low = "blue", high = "red") +  # Color gradient for PValue
  labs(x = "Drug Name", y = "Combined Score", title = "Top 10 Significant Pharmacological Molecules") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  theme_minimal()  # Customize theme if needed

# Column diagram (bar plot) for top 10 significant results
ggplot(top10, aes(x = row.names(top10), y = Combined.Score, fill = P.value)) +
  geom_col() +  # Create column diagram (bar plot)
  scale_fill_gradient(low = "blue", high = "red") +  # Color gradient for PValue
  labs(x = "Drug Name", y = "Combined Score", title = "Top 10 Significant Pharmacological Molecules") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  theme_minimal()  # Customize theme if needed

# Bar chart for top 10 significant results
ggplot(top10, aes(x = row.names(top10), y = Combined.Score, fill = P.value)) +
  geom_bar(stat = "identity", color = "black") +  # Create bar chart
  scale_fill_gradient(low = "blue", high = "red") +  # Color gradient for PValue
  labs(x = "Drug Name", y = "Combined Score", title = "Top 10 Significant Pharmacological Molecules") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  theme_minimal()  # Customize theme if needed


# Assuming 'drug_data' is your dataset

# Transform p-value for better visualization (using -log10)
top10$p_value_transformed <- -log10(top10$P.value)

# Define the plot
ggplot(top10, aes(x = row.names(top10), y = Combined.Score, fill = source)) +
  geom_col() +  # Use geom_col for bar chart
  labs(title = "Candidate Drug Prediction Analysis",
       x = "Candidate Drug",
       y = "Combined Score",
       fill = "Source") +
  theme_minimal()  # Adjust theme for better presentation (optional)

# Plotting
ggplot(top10, aes(x = Combined.Score, y = row.names(top10), size = P.value)) +
  geom_point(color = "blue", alpha = 0.6) +
  scale_size_continuous(range = c(2, 10)) +  # Adjust bubble size range as needed
  labs(title = "Combined Scores and P-Values for Drugs",
       x = "Drug Names", y = "Combined Score", size = "P-Value") +
  theme_minimal()

# Assuming 'top10' is your dataframe with Drug names, Combined Score, and P-values

# Load necessary packages
library(ggplot2)
library(scales)



# Plotting
ggplot(top10, aes(x = Combined.Score, y = row.names(top10))) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_point(aes(size = P.value), color = "blue", alpha = 0.6) +
  scale_size_continuous(range = c(2, 10), trans = "log") +  # Log scale for bubble sizes
  labs(title = "Combined Scores and P-Values for Drugs",
       x = "Combined Score", y = "Drug Names", size = "P-Value (log scale)") +
  theme_minimal()

