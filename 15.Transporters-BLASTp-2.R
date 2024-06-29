# Carbohydrate Transporters exploratatory analyses - processing BLASTp in R
library(data.table) 
setwd("/temporario2/11217468/projects/saureus_global/transporters_blast/blastp_output")
getwd()
carb_transporters_processed <- fread("carb_transporters_processed.tsv", sep = "\t")
print(head(carb_transporters_processed))

library(dplyr)

## Group rows by unique transporters (col 7)
allCARB_transp_grouped <- carb_transporters_processed %>%
  group_by(Annotated_Protein)

## Filter out hits with low identity and bad evalues 
CARB_transp_filt <- allCARB_transp_grouped %>%
  filter(E_Value <= 0.001 & Identity > 25 & Bitscore > 50) # Pearson, William R. "An introduction to sequence similarity (“homology”) searching." Current protocols in bioinformatics 42.1 (2013): 3-1. doi:10.1002/0471250953.bi0301s42.

## Removing duplicated hits (i.e. a hit for more than one annotated protein (prokka) in the same genome)	
CARB_transp_filt <- CARB_transp_filt %>%
  distinct(Genome, .keep_all = TRUE)

View(CARB_transp_filt)
write.table(CARB_transp_filt, file = "CARB_transporters_filt.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

## Make a correspondence between transporter codes and their corresponding carbohydrate 
codes_path <- "/temporario2/11217468/projects/saureus_global/transporters_blast/blastp_output/transporter_codes/transporters_codes.tsv"
carb_transporters_codes <- read.table(codes_path, sep = "\t")
transporters_codes <- as.data.frame(carb_transporters_codes)

library(dplyr)
library(stringr)
### Extract codes between bars "|" in Transporter_DB column
CARB_transp_filt$Transporter_DB <- str_extract(CARB_transp_filt$Transporter_DB, "(?<=\\|)[^|]+(?=\\|)") 

### Create a lookup table for transporter codes and their corresponding headers
lookup_table <- transporters_codes %>%
  pivot_longer(cols = everything(), names_to = "Header", values_to = "Code")

### Merge lookup table with CARB_transporters_filt to replace codes with headers
CARB_transporters_mapped <- CARB_transp_filt %>%
  left_join(lookup_table, by = c("Transporter_DB" = "Code"))

### Rename the new column 
CARB_transporters_mapped <- CARB_transporters_mapped %>%
  rename(Sugar = Header) %>%
  mutate(Sugar = if_else(is.na(Sugar), "polysia_transporters", Sugar))

## Group the dataframe by Transporter and calculate summary statistics
transporter_summary <- CARB_transporters_mapped %>%
  group_by(Sugar) %>%
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
genome_summary <- CARB_transporters_mapped %>%
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
carbs_ARGs_df <- merge(abricate_blactams, genome_summary, by = "Genome", all.x = FALSE)
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
  theme(axis.text.x = element_blank())  # Remove x-axis labels for better readability
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
wilcox_carb1 # W = 133360, p-value = 2.453e-08; 95% CI: 1.000017 2.000064
p_value_carb1 = 2.453e-08

ggplot(carbs_ARGs_df, aes(x = Category, y = Total_Transporters, fill = Category)) +
  geom_boxplot() +
  labs(x = "Category", y = "Carb. Transporter Copies") +
  theme_minimal() +
  # Add significance annotation
  geom_signif(comparisons = list(c("MRSA", "MSSA")), 
              textsize = 6, 
              map_signif_level = TRUE) +
  # Add p-value annotation
  annotate("text", x = 1.5, y = max(carbs_ARGs_df$Total_Transporters), 
           label = paste("p-value > 0.001"), 
           hjust = 0.5, vjust = -0.5, size = 4)

## Combining b_lactams table with the CARB_transporters_mapped df
full_carbs_ARGs_df <- merge(abricate_blactams, CARB_transporters_mapped, by = "Genome", all.x = TRUE)
View(full_carbs_ARGs_df)

## Group by Genome and Sugar, and count the occurrences
grouped_df <- full_carbs_ARGs_df %>%
  group_by(Genome, Sugar, Category) %>%
  summarise(Transporter_Count = n())

## Subset df and repeat this for each unique sugar transporter 
arabinose_ARGs <- grouped_df %>% filter(Sugar == "arabinose_transporters")
mannitol_ARGs <- grouped_df %>% filter(Sugar == "mannitol_transporters")
sia_ARGs <- grouped_df %>% filter(Sugar == "sia_transporters")
polysia_ARGs <- grouped_df %>% filter(Sugar == "polysia_transporters")
galactose_ARGs <- grouped_df %>% filter(Sugar == "galactose_transporters")
mannose_ARGs <- grouped_df %>% filter(Sugar == "mannose_transporters")
sucrose_ARGs <- grouped_df %>% filter(Sugar == "sucrose_transporters")
glucose_ARGs <- grouped_df %>% filter(Sugar == "glucose_transporters")
lactose_ARGs <- grouped_df %>% filter(Sugar == "lactose_transporters")
ribose_ARGs <- grouped_df %>% filter(Sugar == "ribose_transporters")

arabinose_wilcox <- wilcox.test(Transporter_Count ~ Category, data = arabinose_ARGs,
                             exact = TRUE, correct = TRUE, 
                             p.adjust.method = "bonferroni", 
                             conf.int = TRUE, conf.level = 0.95)
arabinose_wilcox # W = 117540, p-value = 0.02918; 95% CI: 6.533734e-05 4.505734e-05

arabinose_boxplot <- ggplot(arabinose_ARGs, aes(x = Category, y = Transporter_Count, fill = Category)) +
  geom_boxplot() +
  labs(x = "Category", y = "Arabinose Transporter Copies") +
  theme_minimal() +
  # Add significance annotation
  geom_signif(comparisons = list(c("MRSA", "MSSA")), 
              textsize = 4, 
              map_signif_level = TRUE) +
  # Add p-value annotation
  annotate("text", x = 1.5, y = max(arabinose_ARGs$Transporter_Count), 
           label = paste("p-value = 0.029"), 
           hjust = 0.5, vjust = -1.2, size = 4)
arabinose_boxplot

sia_wilcox <- wilcox.test(Transporter_Count ~ Category, data = sia_ARGs,
                             exact = TRUE, correct = TRUE, 
                             p.adjust.method = "bonferroni", 
                             conf.int = TRUE, conf.level = 0.95)
sia_wilcox # W = 129172, p-value = 3.589e-08; 95% CI: 1.602656e-05 1.518874e-05

sia_boxplot <- ggplot(sia_ARGs, aes(x = Category, y = Transporter_Count, fill = Category)) +
  geom_boxplot() +
  labs(x = "Category", y = "Sialic acid Transporter Copies") +
  theme_minimal() +
  # Add significance annotation
  geom_signif(comparisons = list(c("MRSA", "MSSA")), 
              textsize = 4, 
              map_signif_level = TRUE) +
  # Add p-value annotation
  annotate("text", x = 1.5, y = max(sia_ARGs$Transporter_Count), 
           label = paste("p-value > 0.001"), 
           hjust = 0.5, vjust = -1.2, size = 4)
sia_boxplot

mannitol_wilcox <- wilcox.test(Transporter_Count ~ Category, data = mannitol_ARGs,
                             exact = TRUE, correct = TRUE, 
                             p.adjust.method = "bonferroni", 
                             conf.int = TRUE, conf.level = 0.95)
mannitol_wilcox # W = 120170, p-value = 0.005939; 95% CI: 2.790275e-05 7.575369e-05

mannitol_boxplot <- ggplot(mannitol_ARGs, aes(x = Category, y = Transporter_Count, fill = Category)) +
  geom_boxplot() +
  labs(x = "Category", y = "Mannitol Transporter Copies") +
  theme_minimal() +
  # Add significance annotation
  geom_signif(comparisons = list(c("MRSA", "MSSA")), 
              textsize = 4, 
              map_signif_level = TRUE) +
  # Add p-value annotation
  annotate("text", x = 1.5, y = max(mannitol_ARGs$Transporter_Count), 
           label = paste("p-value > 0.001"), 
           hjust = 0.5, vjust = -1.2, size = 4)
mannitol_boxplot

polysia_wilcox <- wilcox.test(Transporter_Count ~ Category, data = polysia_ARGs,
                             exact = TRUE, correct = TRUE, 
                             p.adjust.method = "bonferroni", 
                             conf.int = TRUE, conf.level = 0.95)
polysia_wilcox # W = 110352, p-value = 0.1606; 95% CI: -7.522039e-05  4.373916e-06

polysia_boxplot <- ggplot(polysia_ARGs, aes(x = Category, y = Transporter_Count, fill = Category)) +
  geom_boxplot() +
  labs(x = "Category", y = "Polysialic acid Transporter Copies") +
  theme_minimal() +
  # Add significance annotation
  geom_signif(comparisons = list(c("MRSA", "MSSA")), 
              textsize = 4, 
              map_signif_level = TRUE) +
  # Add p-value annotation
  annotate("text", x = 1.5, y = max(polysia_ARGs$Transporter_Count), 
           label = paste("p-value = 0.161"), 
           hjust = 0.5, vjust = -1.2, size = 4)
polysia_boxplot

galactose_wilcox <- wilcox.test(Transporter_Count ~ Category, data = galactose_ARGs,
                             exact = TRUE, correct = TRUE, 
                             p.adjust.method = "bonferroni", 
                             conf.int = TRUE, conf.level = 0.95)
galactose_wilcox # W = 112826, p-value = 0.006195; 95% CI:  -7.444961e-05  5.228018e-05

galactose_boxplot <- ggplot(galactose_ARGs, aes(x = Category, y = Transporter_Count, fill = Category)) +
  geom_boxplot() +
  labs(x = "Category", y = "Galactose Transporter Copies") +
  theme_minimal() +
  # Add significance annotation
  geom_signif(comparisons = list(c("MRSA", "MSSA")), 
              textsize = 4, 
              map_signif_level = TRUE) +
  # Add p-value annotation
  annotate("text", x = 1.5, y = max(galactose_ARGs$Transporter_Count), 
           label = paste("p-value > 0.001"), 
           hjust = 0.5, vjust = -1.2, size = 4)
galactose_boxplot

mannose_wilcox <- wilcox.test(Transporter_Count ~ Category, data = mannose_ARGs,
                             exact = TRUE, correct = TRUE, 
                             p.adjust.method = "bonferroni", 
                             conf.int = TRUE, conf.level = 0.95)
mannose_wilcox # W = 107892, p-value = 0.9541; 95% CI:  -7.993303e-05  2.373064e-05

mannose_boxplot <- ggplot(mannose_ARGs, aes(x = Category, y = Transporter_Count, fill = Category)) +
  geom_boxplot() +
  labs(x = "Category", y = "Mannose Transporter Copies") +
  theme_minimal() +
  # Add significance annotation
  geom_signif(comparisons = list(c("MRSA", "MSSA")), 
              textsize = 4, 
              map_signif_level = TRUE) +
  # Add p-value annotation
  annotate("text", x = 1.5, y = max(mannose_ARGs$Transporter_Count), 
           label = paste("p-value = 0.945"), 
           hjust = 0.5, vjust = -1.2, size = 4)
mannose_boxplot

glucose_wilcox <- wilcox.test(Transporter_Count ~ Category, data = glucose_ARGs,
                             exact = TRUE, correct = TRUE, 
                             p.adjust.method = "bonferroni", 
                             conf.int = TRUE, conf.level = 0.95)
glucose_wilcox #W = 109246, p-value = 0.7497; 95% CI:  -2.777387e-05  3.076857e-07

glucose_boxplot <- ggplot(glucose_ARGs, aes(x = Category, y = Transporter_Count, fill = Category)) +
  geom_boxplot() +
  labs(x = "Category", y = "glucose Transporter Copies") +
  theme_minimal() +
  # Add significance annotation
  geom_signif(comparisons = list(c("MRSA", "MSSA")), 
              textsize = 4, 
              map_signif_level = TRUE) +
  # Add p-value annotation
  annotate("text", x = 1.5, y = max(glucose_ARGs$Transporter_Count), 
           label = paste("p-value = 0.749"), 
           hjust = 0.5, vjust = -1.2, size = 4)
glucose_boxplot

sucrose_wilcox <- wilcox.test(Transporter_Count ~ Category, data = sucrose_ARGs,
                             exact = TRUE, correct = TRUE, 
                             p.adjust.method = "bonferroni", 
                             conf.int = TRUE, conf.level = 0.95)
sucrose_wilcox # W = 108362, p-value = 0.7873; 95% CI:   -5.904733e-05  5.947456e-05

sucrose_boxplot <- ggplot(sucrose_ARGs, aes(x = Category, y = Transporter_Count, fill = Category)) +
  geom_boxplot() +
  labs(x = "Category", y = "sucrose Transporter Copies") +
  theme_minimal() +
  # Add significance annotation
  geom_signif(comparisons = list(c("MRSA", "MSSA")), 
              textsize = 4, 
              map_signif_level = TRUE) +
  # Add p-value annotation
  annotate("text", x = 1.5, y = max(sucrose_ARGs$Transporter_Count), 
           label = paste("p-value = 0.787"), 
           hjust = 0.5, vjust = -1.2, size = 4)
sucrose_boxplot

lactose_wilcox <- wilcox.test(Transporter_Count ~ Category, data = lactose_ARGs,
                             exact = TRUE, correct = TRUE, 
                             p.adjust.method = "bonferroni", 
                             conf.int = TRUE, conf.level = 0.95)
lactose_wilcox # W = 124464, p-value = 0.0001783; 95% CI:    2.766172e-05 2.888186e-05

lactose_boxplot <- ggplot(lactose_ARGs, aes(x = Category, y = Transporter_Count, fill = Category)) +
  geom_boxplot() +
  labs(x = "Category", y = "lactose Transporter Copies") +
  theme_minimal() +
  # Add significance annotation
  geom_signif(comparisons = list(c("MRSA", "MSSA")), 
              textsize = 4, 
              map_signif_level = TRUE) +
  # Add p-value annotation
  annotate("text", x = 1.5, y = max(lactose_ARGs$Transporter_Count), 
           label = paste("p-value > 0.001"), 
           hjust = 0.5, vjust = -1.2, size = 4)
lactose_boxplot

## Plotting all sugars with p < 0.05

library(ggpubr)
ggarrange(arabinose_boxplot, mannitol_boxplot, galactose_boxplot, lactose_boxplot, sia_boxplot,
          labels = c("A", "B", "C", "D", "E"),
          font.label = list(size = 18),
          ncol = 2, nrow = 3,
          legend = "right", 
          common.legend = TRUE)

write.table(carbs_ARGs_df, file = "carbs_ARGs_df.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(full_carbs_ARGs_df, file = "carbs_ARGs_df.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

## Spearman rank correlation test and visualization 
cor_ARG_transporters <- cor.test(carbs_ARGs_df$Total_Transporters, carbs_ARGs_df$Total_ARGs, method = "spearman") # S = 188171246, p-value < 2.2e-16; rho 0.2710382

ggplot(carbs_ARGs_df, aes(x = Total_Transporters, y = Total_ARGs, color = Category)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "auto", se = TRUE, color = "black", linetype = "solid") + # `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'
  labs(x = "Carbohydrate Transporter copies/genome", y = "ARG copies/genome", color = "Category", 
  subtitle = paste("Spearman's rho=", round(cor_ARG_transporters$estimate, 3),
  "\np-value>0.001")) +
  theme_minimal()   

## Calculate odds-ratio - method: median-unbiased estimate & mid-p exact CI  --- correct the code below 
library(epitools)

### Creating my matrix 
copy_number <- c('Genomes Above Median', 'Genomes Below Median')
resistance <- c('MRSA', 'MSSA')
OR_nanT <- matrix(c(594, 292, 127, 103), nrow=2, ncol=2, byrow=TRUE)
dimnames(OR_nanT) <- list('Copy Number'=copy_number, 'Resistance Category'=resistance)

oddsratio(OR_nanT) # odds ratio with 95% C.I. 1.649331 (1.226117 - 2.215623) / midp.exact 0.0009776512 

OR_kpsT <- matrix(c(363, 380, 29, 116), nrow=2, ncol=2, byrow=TRUE)
dimnames(OR_kpsT) <- list('Copy Number'=copy_number, 'Resistance Category'=resistance)

oddsratio(OR_kpsT) # odds ratio with 95% C.I. 3.80132 (2.499973 - 5.955532) / midp.exact 4.052003e-11

## Forest plot - odds ratio 
transporter <- c('kpsT', 'nanT')
col_names <- c('Sia Transporter', 'Odds Ratio', 'CI Low', 'CI High', 'Significance')
OR_plot <- matrix(c('kpsT', 3.80, 2.50, 5.95, 4.052003E-11, 'nanT', 1.65, 1.33, 2.21, 0.0009776512), nrow=2, ncol=5, byrow=TRUE)
dimnames(OR_plot) <- list('Sia Transporters'=transporter, 'Column Names'=col_names)

write.table(OR_plot, file = "OR_sia_transporters.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

OR_data <- read.table("OR_sia_transporters.tsv", header = TRUE, sep = "\t")

ggplot(OR_data, aes(x = Odds_Ratio, y = Sia_Transporter)) +
  geom_errorbarh(aes(xmax = CI_High, xmin = CI_Low), size = 0.5, height = 0.2, color = 'gray50') +
  geom_point(size = 3.5) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = c("top")) +  # Set legend position
  ylab('') +
  xlab('Odds ratio') +
  #ggtitle('Odds Ratio for Sia Transporters') +
  geom_text(aes(label = Odds_Ratio), nudge_y = 0.25) +  # Add odds ratio values
  geom_vline(xintercept = 1, linetype = "dashed", color = "red")  # Add dashed line at x = 1.0




















