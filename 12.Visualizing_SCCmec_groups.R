# Visualizing SCCmec types across S. aureus genomes in R

## Imported sccmec.staphopia table into R

## Combining all subtypes into a single column representing all genomes included in a ST 
library(dplyr)
library(readr)

### Define the columns to check for each group of subtypes
subtype_columns <- list(
  I = c("I", "Ia"),
  II = c("II", "IIa", "IIb"),
  III = c("III", "IIIa"),
  IV = c("IV", "IVa", "IVb", "IVc", "IVd", "IVg", "IVh")
)
#### For columns V, VI, VII, VIII and IX, you can just sum and add the final value afterwards, since they dn't have any subtypes 

### Loop through each subtype group and add new columns based on conditions
for (subtype in names(subtype_columns)) {
  columns_to_check <- subtype_columns[[subtype]]
  new_column_name <- paste0(subtype, "_all")
  
  # Add a new column based on the condition
  sccmec.staphopia <- mutate(
    sccmec.staphopia,
    !!new_column_name := ifelse(
      rowSums(select(sccmec.staphopia, all_of(columns_to_check)) == "True") > 0,
      "True",
      "False"
    )
  )
}

View(sccmec.staphopia)

## Making a new table with summed results
sum_sccmec_staphopia <- sccmec.staphopia[, c(1, 22:25, 6:11)]
View(sum_sccmec_staphopia)

## Adding a new column with True for SCCmec negative genomes 
all_false_rows <- rowSums(sum_sccmec_staphopia[, 2:11] == "False") == ncol(sum_sccmec_staphopia[, 2:11]) # make a logical vector to store genomes with False in all columns 
transposed_data <- data.frame(`SCCmec_neg` = all_false_rows) # make a new df with the transposed data
sum_sccmec_staphopia <- cbind(sum_sccmec_staphopia, transposed_data) # append the transposed vector as a new column

write.table(sum_sccmec_staphopia, file = "sum_sccmec-staphopia", sep = "\t", row.names = FALSE, col.names = TRUE)

## Correct columns names 
new_column_names <- c("Genome", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "mecA", "SCCmec_neg")
colnames(sum_sccmec_staphopia) <- new_column_names
colnames(sum_sccmec_staphopia)

## Correct last column issue of caps lock - Replace FALSE with False and TRUE with True in the last column
sum_sccmec_staphopia <- sum_sccmec_staphopia %>%
  mutate(SCCmec_neg = ifelse(SCCmec_neg == "TRUE", "True", "False"))

## Count 'True' hits in each category 
true_hits_counts <- colSums(sum_sccmec_staphopia[, 2:12] == "True")

## Create a data frame for plotting
hit_counts_df <- data.frame(Group = names(true_hits_counts), Count = true_hits_counts)

## Reorder Group by Count in ascending order
hit_counts_df$Group <- reorder(hit_counts_df$Group, hit_counts_df$Count)

## Save the final SCCmec count table 
write.table(hit_counts_df, "SCCmec_count_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## Numbers aren't matching - check where's the error 

### Identify rows where there's a discrepancy between columns 2-10 and column 11 (identified a SCCmec type but not mecA gene)
inconsistent_rows <- sum_sccmec_staphopia %>%
  filter((I == "True" | II == "True" | III == "True" | IV == "True" | V == "True" |
          VI == "True" | VII == "True" | VIII == "True" | IX == "True") & mecA != "True")
View(inconsistent_rows) # 16

### Identify rows where there's more than one identified SCCmec type - ambiguous 
columns_to_check <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX")
multiple_sccmec_types <- sum_sccmec_staphopia %>%
  filter(rowSums(select(., all_of(columns_to_check)) == "True") > 1)
View(multiple_sccmec_types) # 184

### Find genomes with both errors
genomes_with_both_errors <- semi_join(inconsistent_rows, multiple_sccmec_types, by = "Genome")
View(genomes_with_both_errors) # no genome has both errors 

# Add a new column 'Unassigned' based on the condition (multiple SCCmec types)
sum_sccmec_staphopia_CORRECT <- sum_sccmec_staphopia %>%
  mutate(Unassigned = ifelse(rowSums(select(., all_of(columns_to_check)) == "True") > 1, "True", "False"))
View(sum_sccmec_staphopia_CORRECT)

write.table(sum_sccmec_staphopia_CORRECT, "sum_sccmec-staphopia-CORRECT.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## Count 'True' hits in each category 
true_hits_counts <- colSums(sum_sccmec_staphopia_CORRECT[, 2:13] == "True")

## Create a data frame for plotting
hit_counts_df <- data.frame(Group = names(true_hits_counts), Count = true_hits_counts)

## Reorder Group by Count in ascending order
hit_counts_df$Group <- reorder(hit_counts_df$Group, hit_counts_df$Count)

write.table(hit_counts_df, "SCCmec_count_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# I've edited it a little bit by hand and imported it back to R (added a category, in which all SCCmec types identified are in the SCCmec_pos category, I've removed mecA row and SCCmec_neg is in it's own category and the same with Unassigned)
plot_sccmec_count <- plot_sccmec_count[order(plot_sccmec_count$Category, -plot_sccmec_count$Count),]
ggplot(plot_sccmec_count, aes(x = Category, y = Count, fill = Group)) +
  geom_bar(stat = "identity") +
  labs(title = "Stacked Bar Plot of SCCmec Counts by Category", x = "Category", y = "Number of Genomes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")  # Use a color palette from RColorBrewer








