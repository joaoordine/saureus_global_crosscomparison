# Visualizing mlst results - frequency bar plot of most representatives STs in R
library(ggplot2)
data <- just_STs_W_header_R

## Pre-processing 
filtered_data <- subset(data, Count > 10) # Filter the data to include only counts greater than 10
others_data <- data.frame(Count = sum(data$Count[data$Count <= 10]),
                          Sequence_Type = "Others") # Filter the data to include only counts greater than 10

combined_data <- rbind(filtered_data, others_data) # Combine the filtered data and 'Others' data

## Create a custom color palette
my_colors <- c("Others" = "black", "Unclassified" = "#303030")

unique_types <- unique(filtered_data$Sequence_Type)
for (st in unique_types) {
    my_colors <- c(my_colors, setNames("gray", st))
} 
## Reorder Sequence_Type by Count in ascending order
combined_data$Sequence_Type <- reorder(combined_data$Sequence_Type, combined_data$Count)

## Create the barplot with different fill colors
ggplot(combined_data, aes(x = Count, y = Sequence_Type, fill = Sequence_Type)) +
    geom_bar(stat = "identity") +
    labs(x = "Absolute Frequency", y = "Sequence Type") +
    scale_fill_manual(values = my_colors)
