# Visualizing SCCmec types across S. aureus genomes in R

# Bar Plot to Visualize SCCmec Type, mecA and ccr gene Counts
## Load necessary libraries
library(ggplot2)

## Remove the 'ccr4' column
summary_SCCmectable <- summary_SCCmectable[, !names(summary_SCCmectable) %in% "ccr4"]

## Count 'True' hits in each group (columns 2 to 12)
true_hits_counts <- colSums(summary_SCCmectable[, 2:11] == "True")

## Create a data frame for plotting
hit_counts_df <- data.frame(Group = names(true_hits_counts), Count = true_hits_counts)

## Identify genomes with "FALSE" in all SCCmec type columns and mecA gene 
all_false_genomes <- rowSums(summary_SCCmectable[, 2:11] == "False") == ncol(summary_SCCmectable[, 2:11])

## Count genomes with "FALSE" in all SCCmec type columns
false_all_sccmec_count <- sum(all_false_genomes) # 835


## Create a data frame for plotting
hit_counts_df2 <- data.frame(Group = c(names(true_hits_counts), 'SCCmec negative'), 
                             Count = c(true_hits_counts, false_all_sccmec_count))

## Reorder Group by Count in ascending order
hit_counts_df2$Group <- reorder(hit_counts_df2$Group, hit_counts_df2$Count)

## Create a custom color palette
my_colors <- c("mecA147" = "#006400", "SCCmec negative" = "#8B0000", "TYPE_I" = "green", "TYPE_IVa" = "green", "TYPE_II" = "green", "TYPE_V" = "green", "TYPE_III" = "green", "TYPE_IVc" = "green", "TYPE_IVb" = "green", "TYPE_IVd" = "green", "TYPE_IVe" = "green")


## Create the barplot with different fill colors
ggplot(hit_counts_df2, aes(x = Count, y = Group, fill = Group)) +
    geom_bar(stat = "identity") +
    labs(title = "Number of genomes in each SCCmec group", x = "Count", y = "SCCmec Type") +
    scale_fill_manual(values = my_colors)

## Save the final SCCmec count table 
write.table(hit_counts_df2, "SCCmec_count_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)















