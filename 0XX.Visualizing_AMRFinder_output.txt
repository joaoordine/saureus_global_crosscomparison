# Visualize AMRFinder output - all eSRGs in S aureus genomes

# Heat map 
library(pheatmap)

## Remove the first column (Genome) as it's not needed for the heatmap
heatmap_stress_df <- amrfinder.eSRGs_corrected[, -1]

## Create a custom color palette (0 = white; 1 = green)
color_palette <- colorRampPalette(c("white", "green"))(2)  # 0 to 1 inclusive

## Create the heatmap
pheatmap(
  heatmap_stress_df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = color_palette,
  border_color = "grey",
  main = "Heatmap of SRGs Presence in Genomes",
  breaks = seq(0, 1, by = 0.5)  # Specify the breaks from 0 to 1
)


# Barplot 

## Calculate the sum of values in each column
stress_hit_count <- colSums(heatmap_stress_df, na.rm = TRUE)

## Create a dataframe for plotting
stress_hit_count_df <- data.frame(SRGs = names(stress_hit_count), Count = (stress_hit_count))

## Sort the dataframe by Count in descending order
stress_hit_count_df <- stress_hit_count_df[order(-stress_hit_count_df$Count), ]

library(ggplot2)

## Create the bar plot with sorted data
ggplot(stress_hit_count_df, aes(x = Count, y = reorder(SRGs, --Count))) +
    geom_bar(stat = "identity", fill = "green") +
    labs(title = "SRGs Hit Count in all S. aureus genomes", x = "Count", y = "SRGs")

