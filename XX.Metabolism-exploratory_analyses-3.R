# Metabolism explatory analyses 3 - processing hmmer output files in R and plotting 

library(dplyr)

## Removing undesired things from the protein name and genome code 
sym_processed$Transporter <- gsub("_align", "", sym_processed$Transporter)
sym_processed$Genome <- gsub("sym_", "", sym_processed$Genome)

## Calculate the coverage of each hit
sym_processed <- sym_processed %>%
  mutate(Coverage = ( Target_Length / Transporter_Length) * 100)
  
## Removing duplicated hits (i.e. more than one prokka annotated protein)
sorted_processed <- sym_processed %>% 
  arrange(Target_Code, desc(Coverage)) # Sort the dataframe by column 1 (Target_Code) and then by column 6 (Score) in descending order

processed_max_score <- sorted_processed %>% 
  group_by(Target_Code) %>% 
  slice_max(order_by = Coverage) # Keep only the rows with the highest coverage value for each unique value in column 1 (Target_Code)

processed_max_score <- processed_max_score %>%
  filter(E_Value <= 0.01) # filter for e-value below 0.01
  
sym_filt <- processed_max_score %>% 
  distinct(Target_Code, .keep_all = TRUE) # removing remaining duplicates with identical scores 
View(sym_filt)

write.table(sym_filt, file = "sym_filt_unique.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

## Group the dataframe by Genome and calculate summary statistics
genome_summary <- sym_filt %>%
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

sym_stats_per_genome <- genome_summary
write.table(sym_stats_per_genome, file = "sym_stats_per_genome.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

## Combining b_lactams table with the summary of each transporter

### Remove suffix and rename columns
abricate_blactams$Genome <- gsub("_genomic-ARGs.txt", "", abricate_blactams$FILE)
names(abricate_blactams)[1] <- "Genome"
names(abricate_blactams)[2] <- "Total_ARGs"

### Merge dfs 
sym_ARGs_df <- merge(abricate_blactams, sym_stats_per_genome, by = "Genome", all.x = TRUE)
sym_ARGs_df <- sym_ARGs_df[, c("Genome", "Total_Transporters", "Total_ARGs", "Category")] # filtering out for specific columns
sym_ARGs_df$Total_Transporters[is.na(sym_ARGs_df$Total_Transporters)] <- 0 # replace NAs with 0


## Make a scatter plot to check the distribution of your data 
library(ggplot2)
library(ggsignif)

ggplot(sym_ARGs_df, aes(x = Total_ARGs, y = Total_Transporters, color = Category)) +
  geom_point() +
  labs(x = "ARG copies/genome", y = "sym copies/genome", color = "Category") +
  theme_minimal()
    
## Make a box plot by b-lactam resistance category
boxplot(Total_Transporters ~ Category, data = sym_ARGs_df, col = c("blue", "red"), main = "Sym copy numbers by Category", xlab = "Category", ylab = "Sym copy numbers")
sym_outliers <- boxplot.stats(sym_ARGs_df$Total_Transporters)$out

wilcox_sym <- wilcox.test(Total_Transporters ~ Category, data = sym_ARGs_df)
p_value_sym <- wilcox_sym$p.value

ggplot(sym_ARGs_df, aes(x = Category, y = Total_Transporters, fill = Category)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("MRSA", "MSSA")), 
              textsize = 6, 
              map_signif_level = TRUE) +
  labs(title = "Sym Transporters by Category",
       x = "Category",
       y = "Sym Copy Number",
       fill = "Category") +
  theme_minimal() +
  annotate("text", x = 1.5, y = max(sym_ARGs_df$Total_Transporters), 
           label = paste("p-value:", round(p_value_sym, 4)), 
           hjust = 0.5, vjust = -0.5, size = 4)

write.table(sym_ARGs_df, file = "sym_ARGs_df.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

## Repeated the whole workflow for each of the transporters 

##############################


## Combining dataframes 
### Add Transporter_Name column to each dataframe
sym_stats_per_genome$Transporter_Name <- "sym"
nanT_stats_per_genome$Transporter_Name <- "nanT"
kpsT_stats_per_genome$Transporter_Name <- "kpsT"
kpsM_stats_per_genome$Transporter_Name <- "kpsM"

### Concatenate the dataframes
tranposters_stats_per_genome <- rbind(sym_stats_per_genome, nanT_stats_per_genome, kpsT_stats_per_genome, kpsM_stats_per_genome)

View(tranposters_stats_per_genome)

## Creating a transporter count dataframe 
library(dplyr)
library(tidyr)

### Group by Genome and Transporter_Name, then summarize to get the sum of Total_Transporters for each combination
transporters_sum <- tranposters_stats_per_genome %>%
  group_by(Genome, Transporter_Name) %>%
  summarize(Total_Transporters_Sum = sum(Total_Transporters, na.rm = TRUE))

### Pivot the table to get Transporter_Name as columns
transporters_sum_pivot <- pivot_wider(transporters_sum, names_from = Transporter_Name, values_from = Total_Transporters_Sum)

### Fill NA values with 0
transporters_sum_pivot[is.na(transporters_sum_pivot)] <- 0

### Create the whole_transporters_count dataframe
whole_transporters_count <- transporters_sum_pivot %>%
  mutate(Total_Transporters_kpsM = kpsM,
         Total_Transporters_kpsT = kpsT,
         Total_Transporters_nanT = nanT,
         Total_Transporters_sym = sym) %>%
  select(Genome, Total_Transporters_kpsM, Total_Transporters_kpsT, Total_Transporters_nanT, Total_Transporters_sym)

View(whole_transporters_count)

## Mege transporter count dataframe with AMR profile against b-lactams
whole_transporters_count <- merge(abricate_blactams, whole_transporters_count, by = "Genome", all.x = TRUE)
whole_transporters_count <- whole_transporters_count[, c("Genome", "Total_Transporters_kpsM", "Total_Transporters_kpsT", "Total_Transporters_nanT", "Total_Transporters_sym", "Total_ARGs", "Category")] # filtering out for specific columns
whole_transporters_count[is.na(whole_transporters_count)] <- 0 # replace NAs with 0

write.table(whole_transporters_count, file = "whole_transporters_count.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

## Visualize data comparing individual transporters number x ARGs number
ggplot(whole_transporters_count, aes(x = Total_ARGs, y = Total_Transporters_kpsM, color = Category)) +
  geom_point() +
  labs(x = "ARG copies/genome", y = "kpsM copies/genome", color = "Category") +
  theme_minimal()
  
ggplot(whole_transporters_count, aes(x = Total_ARGs, y = Total_Transporters_kpsT, color = Category)) +
  geom_point() +
  labs(x = "ARG copies/genome", y = "kpsT copies/genome", color = "Category") +
  theme_minimal()
  
ggplot(whole_transporters_count, aes(x = Total_ARGs, y = Total_Transporters_nanT, color = Category)) +
  geom_point() +
  labs(x = "ARG copies/genome", y = "nanT copies/genome", color = "Category") +
  theme_minimal()  
  
## Correlation analysis with XICOR - Computes robust association measures that do not presuppose linearity. 
###The xi correlation (xicor) is based on cross correlation between ranked increments. The reference for the methods implemented here is Chatterjee, Sourav (2020) 
install.packages("XICOR")
library(XICOR)

xicor_result_nanT <- xicor(whole_transporters_count$Total_Transporters_nanT, whole_transporters_count$Total_ARGs) # Compute xi correlation (xicor)
print(xicor_result_nanT)

## Calculate the median of each variable 
### Splitting the dataframe based on the "Category" column
MRSA_transporters_count <- subset(whole_transporters_count, Category == "MRSA")
MSSA_transporters_count <- subset(whole_transporters_count, Category == "MSSA")

exclude_cols <- c(1, 7) # define columns to exclude during calculating medians 
medians_MRSA <- sapply(MRSA_transporters_count[, -exclude_cols], median)
medians_MSSA <- sapply(MSSA_transporters_count[, -exclude_cols], median) # calculate the meadian 

print(medians_MRSA)
print(medians_MSSA)

write.table(MRSA_transporters_count, file = "MRSA_transporters_count.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(MSSA_transporters_count, file = "MSSA_transporters_count.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

## Calculate odds-ratio - only for kpsT and nanT - method: median-unbiased estimate & mid-p exact CI
install.packages('epitools')
library(epitools)

### Creating my matrix 
copy_number <- c('Genomes Above Median', 'Genomes Below Median')
resistance <- c('MRSA', 'MSSA')
OR_nanT <- matrix(c(456, 153, 603, 118), nrow=2, ncol=2, byrow=TRUE)
dimnames(OR_nanT) <- list('Copy Number'=copy_number, 'Resistance Category'=resistance)

oddsratio(OR_nanT) # odds ratio with 95% C.I. 0.5836476 (0.4452102 - 0.7634528) / midp.exact 8.283424e-05 

OR_kpsT <- matrix(c(554, 129, 552, 107), nrow=2, ncol=2, byrow=TRUE)
dimnames(OR_kpsT) <- list('Copy Number'=copy_number, 'Resistance Category'=resistance)

oddsratio(OR_kpsT) # odds ratio with 95% C.I. 0.8327422 (0.6271817 - 1.103854) / midp.exact 0.2033205

## Forest plot - odds ratio 
transporter <- c('kpsT', 'nanT')
col_names <- c('Sia Transporter', 'Odds Ratio', 'CI Low', 'CI High', 'Significance')
OR_plot <- matrix(c('kpsT', 0.83, 0.62, 1.10, 0.2033205, 'nanT', 0.58, 0.44, 0.76, 0.00008283424), nrow=2, ncol=5, byrow=TRUE)
dimnames(OR_plot) <- list('Sia Transporters'=transporter, 'Column Names'=col_names)

write.table(OR_plot, file = "OR_sia_transporters.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

OR_data <- read.table("OR_sia_transporters.tsv", header = TRUE, sep = "\t")

#OR_data <- OR_plot

ggplot(OR_data, aes(x = `Odds Ratio`, y = `Sia Transporter`, color = `Significance`)) +
  geom_errorbarh(aes(xmax = `CI High`, xmin = `CI Low`), size = 0.5, height = 0.2, color = 'gray50') +
  geom_point(size = 3.5) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = c("top")) +  # Set legend position
  ylab('') +
  xlab('Odds ratio') +
  #ggtitle('Odds Ratio for Sia Transporters') +
  geom_text(aes(label = `Odds Ratio`), nudge_y = 0.25) +  # Add odds ratio values
  geom_vline(xintercept = 5.9, linetype = "dashed", color = "red")  # Add dashed line at x = 1.0
  #coord_flip()  # Flip the plot












