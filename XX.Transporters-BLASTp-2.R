# Carbohydrate Transporters exploratatory analyses - processing BLASTp in R

library(dplyr)

## Group rows by unique genomes (col 7)
allCARB_transp_grouped <- carb_transporters_processed %>%
  group_by(Annotated_Protein)

## Filter out hits with low identity and bad evalues 
allCARB_transp_grouped <- allCARB_transp_grouped %>%
  filter(E_Value <= 0.001 & Identity > 25 & Bitscore 50) # Pearson, William R. "An introduction to sequence similarity (“homology”) searching." Current protocols in bioinformatics 42.1 (2013): 3-1. doi:10.1002/0471250953.bi0301s42.

# Filter rows where the values in Transporter_DB column start with "tr"
CARB_transp_filt <- allCARB_transp_grouped %>%
  filter(substr(Transporter_DB, 1, 2) == "tr")


## Removing duplicated hits (i.e. a hit for more than one annotated protein (prokka) in the same genome)	
CARB_transp_filt <- CARB_transp_filt %>%
  distinct(Genome, .keep_all = TRUE)

View(CARB_transp_filt)
write.table(CARB_transp_filt, file = "CARB_transporters_filt.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

# MAKE CORRESPONDENCE BETWEEN TRANSPORTERS UNIPROT CODE AND CORRESPONDING SUGAR --> THEN REPEAT THOSE CODES BELOW

## Group the dataframe by Transporter and calculate summary statistics
transporter_summary <- CARB_transp_filt %>%
  group_by(Transporter_DB) %>%
  summarise(
    Total_Transporters = n_distinct(Annotated_Protein),  # Count unique transporters per genome
    Max_Bitscore = max(Bitscore),                         # Maximum score per genome
    Min_Bitscore = min(Bitscore),                         # Minimum score per genome
    Mean_Bitscore = mean(Bitscore),                       # Mean score per genome
    Median_Bitscore = median(Bitscore),                    # Median score per genome
    Max_E_Value = max(E_Value),                     # Maximum E_Value per genome
    Min_E_Value = min(E_Value),                     # Minimum E_Value per genome
    Mean_E_Value = mean(E_Value),                   # Mean E_Value per genome
    Median_E_Value = median(E_Value),               # Median E_Value per genome
    Max_Identity = max(Identity),                     # Maximum E_Value per genome
    Min_Identity = min(Identity),                     # Minimum E_Value per genome
    Mean_Identity = mean(Identity),                   # Mean E_Value per genome
    Median_Identity = median(Identity),               # Median E_Value per genome
    Max_Length = max(Length),                   # Maximum Coverage per genome
    Min_Length = min(Length),                   # Minimum Coverage per genome
    Mean_Length = mean(Length),                 # Mean Coverage per genome
    Median_Length = median(Length)              # Median Coverage per genome
  ) %>%
  ungroup()  # Remove grouping for further operations

View(transporter_summary)
write.table(transporter_summary, file = "carb_transporter_summary.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

## Group the dataframe by Genome and calculate summary statistics
genome_summary <- CARB_transp_filt %>%
  group_by(Genome) %>%
  summarise(
    Total_Transporters = n_distinct(Annotated_Protein),  # Count unique transporters per genome
    Max_Bitscore = max(Bitscore),                         # Maximum score per genome
    Min_Bitscore = min(Bitscore),                         # Minimum score per genome
    Mean_Bitscore = mean(Bitscore),                       # Mean score per genome
    Median_Bitscore = median(Bitscore),                    # Median score per genome
    Max_E_Value = max(E_Value),                     # Maximum E_Value per genome
    Min_E_Value = min(E_Value),                     # Minimum E_Value per genome
    Mean_E_Value = mean(E_Value),                   # Mean E_Value per genome
    Median_E_Value = median(E_Value),               # Median E_Value per genome
    Max_Identity = max(Identity),                     # Maximum E_Value per genome
    Min_Identity = min(Identity),                     # Minimum E_Value per genome
    Mean_Identity = mean(Identity),                   # Mean E_Value per genome
    Median_Identity = median(Identity),               # Median E_Value per genome
    Max_Length = max(Length),                   # Maximum Coverage per genome
    Min_Length = min(Length),                   # Minimum Coverage per genome
    Mean_Length = mean(Length),                 # Mean Coverage per genome
    Median_Length = median(Length)              # Median Coverage per genome
  ) %>%
  ungroup()  # Remove grouping for further operations

View(genome_summary)

write.table(genome_summary, file = "genome_summary_CARBtransp.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

## Combining b_lactams table with the transporters summary 

### Remove suffix and rename columns
abricate_blactams$FILE <- gsub("_genomic-ARGs.txt", "", abricate_blactams$FILE)
names(abricate_blactams)[1] <- "Genome"
names(abricate_blactams)[2] <- "Total_ARGs"

### Merge dfs 
carbs_ARGs_df <- merge(abricate_blactams, genome_summary, by = "Genome", all.x = TRUE)
View(carbs_ARGs_df)

### Replace coverage values for blaZ and mec genes with hits count
count_semicolons <- function(x) {
  ifelse(x == 0, 0, nchar(x) - nchar(gsub(";", "", x)) + 1)
} # Create a function to count semicolons in each element of a column
carbs_ARGs_df <- carbs_ARGs_df %>%
  mutate(
    blaZ = count_semicolons(blaZ),
    mecA = count_semicolons(mecA),
    mecC = count_semicolons(mecC)
  ) # Replace values in columns 3, 4, and 5 based on the condition
 
carbs_ARGs_df$Total_Transporters[is.na(carbs_ARGs_df$Total_Transporters)] <- 0 # replace NAs with 0

## Check the distribution of your data
library(ggplot2)
library(ggsignif)
library(dplyr)
 
### Histogram with Transporter copy numbers in genomes
transp_distribution <- ggplot(carbs_ARGs_df, aes(x = reorder(Genome, -Total_Transporters), y = Total_Transporters)) +
  geom_bar(stat = "identity", fill = "gray") +
  labs(x = "Genome", y = "Total Carbohydrate Transporters", title = "Carbohydrate Transporter Distribution across Genome") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  # Rotate x-axis labels for better readability
transp_distribution

### Scatter plot ARG x Transporters
ggplot(carbs_ARGs_df, aes(x = Total_Transporters, y = Total_ARGs, color = Category)) +
  geom_point() +
  labs(x = "Carbohydrate Transporter copies/genome", y = "ARG copies/genome", color = "Category") +
  theme_minimal()   
  
## Make a box plot by b-lactam resistance category
transporter_boxplot <- ggplot(carbs_ARGs_df, aes(x = Category, y = Total_Transporters, fill = Category)) +
  geom_boxplot() +
  labs(x = "Blactam Category", y = "Carbohydrate Transporter copy number") +
  theme_minimal() 
transporter_boxplot

arg_boxplot <- ggplot(carbs_ARGs_df, aes(x = Category, y = Total_ARGs, fill = Category)) +
  geom_boxplot() +
  labs(x = "Blactam Category", y = "ARG copy numbers") +
  theme_minimal() 
arg_boxplot

## Check for statiscal difference 
wilcox_carb1 <- wilcox.test(Total_Transporters ~ Category, data = carbs_ARGs_df,
                             exact = TRUE, correct = TRUE, 
                             p.adjust.method = "bonferroni", 
                             conf.int = TRUE, conf.level = 0.95)
wilcox_carb1 # W = 19918, p-value = 0.0006478; 95% CI: 0.9999955 2.9999665

## Group observations into two levels (genomes with 0-14 ARG copies and with 15-30 ARG copies) and perform stat test
carbs_ARGs_df$Total_ARGs_group <- cut(carbs_ARGs_df$Total_ARGs, breaks = c(-Inf, 13, 30), labels = c("ARGs < median", "ARGs > median"))

wilcox_carb2 <- wilcox.test(carbs_ARGs_df$Total_Transporters ~ carbs_ARGs_df$Total_ARGs_group, 
            p.adjust.method = "bonferroni",
            exact = TRUE, 
            correct = TRUE,
            conf.int = TRUE,
            conf.level = 0.95) 
wilcox_carb2 # W = 19508, p-value < 2.2e-16; 95% CI: -8.276209e-05 -3.934423e-06
p_value_carb2 = 2.2e-16

ggplot(carbs_ARGs_df, aes(x = Total_ARGs_group, y = Total_Transporters, fill = Total_ARGs_group)) +
  geom_boxplot() +
  labs(x = "Category", y = "Carb. Transporter Copies") +
  theme_minimal() +
  # Add significance annotation
  geom_signif(comparisons = list(c("ARGs < median", "ARGs > median")), 
              textsize = 6, 
              map_signif_level = TRUE) +
  # Add p-value annotation
  annotate("text", x = 1.5, y = max(carbs_ARGs_df$Total_Transporters + carbs_ARGs_df$Total_ARGs), 
           label = paste("p-value:", round(p_value_carb4, 3)), 
           hjust = 0.5, vjust = -0.5, size = 4)

write.table(carbs_ARGs_df, file = "carbs_ARGs_df.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
























