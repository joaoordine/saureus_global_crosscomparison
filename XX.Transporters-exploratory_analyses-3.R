# Carbohydrate Transporters explatory analyses 3 - processing hmmer output files in R and plotting 

library(dplyr)

## Removing undesired things from the protein name and genome code 

carb_transporters_processed$Transporter <- gsub("_align", "", carb_transporters_processed$Transporter)
carb_transporters_processed$Transporter <- gsub("uniprotkb_", "", carb_transporters_processed$Transporter)
carb_transporters_processed$Transporter <- sub("_.*", "", carb_transporters_processed$Transporter) # remove everything after underscore


carb_transporters_processed$Genome <- gsub("glucose_", "", carb_transporters_processed$Genome)
carb_transporters_processed$Genome <- gsub("sia_", "", carb_transporters_processed$Genome)
carb_transporters_processed$Genome <- gsub("polysia_", "", carb_transporters_processed$Genome)
carb_transporters_processed$Genome <- gsub("arabinose_", "", carb_transporters_processed$Genome)
carb_transporters_processed$Genome <- gsub("galactose_", "", carb_transporters_processed$Genome)
carb_transporters_processed$Genome <- gsub("mannitol_", "", carb_transporters_processed$Genome)
carb_transporters_processed$Genome <- gsub("mannose_", "", carb_transporters_processed$Genome)
carb_transporters_processed$Genome <- gsub("ribose_", "", carb_transporters_processed$Genome)
carb_transporters_processed$Genome <- gsub("sucrose_", "", carb_transporters_processed$Genome)
carb_transporters_processed$Genome <- gsub("lactose_", "", carb_transporters_processed$Genome)

allCARB_transp_df <- carb_transporters_processed 
View(allCARB_transp_df)

## Calculate the coverage of each hit
allCARB_transp_df <- allCARB_transp_df %>%
  mutate(Coverage = ( Target_Length / Transporter_Length) * 100)

## Group rows by unique genomes (col 7)
allCARB_transp_grouped <- allCARB_transp_df %>%
  group_by(Genome)

## Filter out hits with low coverages and bad evalues 
allCARB_transp_grouped <- allCARB_transp_grouped %>%
  filter(E_Value <= 0.01) # filter for e-value below 0.01
allCARB_transp_grouped <- allCARB_transp_grouped %>%
  filter(Coverage > 70) # filter for coverage above 70%

## Removing duplicated hits (i.e. a hit for more than one annotated protein (prokka) in the same genome)	
unique_rows <- allCARB_transp_df %>%
  distinct(Target_Code, .keep_all = TRUE)

CARB_transp_filt <- unique_rows
View(CARB_transp_filt)

write.table(CARB_transp_filt, file = "CARB_transporters_filt.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

## Group the dataframe by Transporter and calculate summary statistics
transporter_summary <- CARB_transp_filt %>%
  group_by(Transporter) %>%
  summarise(
    Total_Transporters = n_distinct(Target_Code),  # Count unique transporters per genome
    Max_Score = max(Score),                         # Maximum score per genome
    Min_Score = min(Score),                         # Minimum score per genome
    Mean_Score = mean(Score),                       # Mean score per genome
    Median_Score = median(Score),                    # Median score per genome
    Max_E_Value = max(E_Value),                     # Maximum E_Value per genome
    Min_E_Value = min(E_Value),                     # Minimum E_Value per genome
    Mean_E_Value = mean(E_Value),                   # Mean E_Value per genome
    Median_E_Value = median(E_Value),               # Median E_Value per genome
    Max_Coverage = max(Coverage),                   # Maximum Coverage per genome
    Min_Coverage = min(Coverage),                   # Minimum Coverage per genome
    Mean_Coverage = mean(Coverage),                 # Mean Coverage per genome
    Median_Coverage = median(Coverage)              # Median Coverage per genome
  ) %>%
  ungroup()  # Remove grouping for further operations

View(transporter_summary)

write.table(transporter_summary, file = "carb_transporter_summary.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

## Group the dataframe by Genome and calculate summary statistics
genome_summary <- CARB_transp_filt %>%
  group_by(Genome) %>%
  summarise(
    Total_Transporters = n_distinct(Target_Code),  # Count unique transporters per genome
    Max_Score = max(Score),                         # Maximum score per genome
    Min_Score = min(Score),                         # Minimum score per genome
    Mean_Score = mean(Score),                       # Mean score per genome
    Median_Score = median(Score),                    # Median score per genome
    Max_E_Value = max(E_Value),                     # Maximum E_Value per genome
    Min_E_Value = min(E_Value),                     # Minimum E_Value per genome
    Mean_E_Value = mean(E_Value),                   # Mean E_Value per genome
    Median_E_Value = median(E_Value),               # Median E_Value per genome
    Max_Coverage = max(Coverage),                   # Maximum Coverage per genome
    Min_Coverage = min(Coverage),                   # Minimum Coverage per genome
    Mean_Coverage = mean(Coverage),                 # Mean Coverage per genome
    Median_Coverage = median(Coverage)              # Median Coverage per genome
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
wilcox_carb <- pairwise.wilcox.test(carbs_ARGs_df$Total_Transporters, carbs_ARGs_df$Total_ARGs, p.adjust.method = "fdr")
wilcox_carb # pairwise comparsion of ARG number x Transporter number

wilcox_carb2 <- wilcox.test(carbs_ARGs_df$Total_Transporters, 
            carbs_ARGs_df$Total_ARGs, 
            p.adjust.method = "bonferroni",
            paired = TRUE, 
            exact = TRUE, 
            correct = TRUE,
            conf.int = TRUE,
            conf.level = 0.95)
wilcox_carb2 # W = 75554, p-value < 2.2e-16

wilcox_carb3 <- wilcox.test(Total_Transporters ~ Category, data = carbs_ARGs_df,
                             exact = TRUE, correct = TRUE, 
                             p.adjust.method = "bonferroni", 
                             conf.int = TRUE, conf.level = 0.95)
wilcox_carb3 # W = 225188, p-value = 0.1469

## Group observations into two levels (genomes with 0-14 ARG copies and with 15-30 ARG copies) and perform stat test
carbs_ARGs_df$Total_ARGs_group <- cut(carbs_ARGs_df$Total_ARGs, breaks = c(-Inf, 11, 30), labels = c("ARGs < median", "ARGs > median"))

wilcox_carb4 <- wilcox.test(carbs_ARGs_df$Total_Transporters ~ carbs_ARGs_df$Total_ARGs_group, 
            p.adjust.method = "bonferroni",
            exact = TRUE, 
            correct = TRUE,
            conf.int = TRUE,
            conf.level = 0.95) 
wilcox_carb4 # W = 334730, p-value = 0.001049; 95 percent confidence interval: -8.276209e-05 -3.934423e-06
p_value_carb4 = 0.001049

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

























































