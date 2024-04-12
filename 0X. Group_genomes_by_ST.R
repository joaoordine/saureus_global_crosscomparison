# Grouping genomes based on ST

## Imported genome_STs.tsv into R; I added a header line by hand "Genome ST" before importing it 
library(dplyr)

## Grouping genomes based on ST and count how many genomes there are in each ST category
grouped_STs <- genome_STs %>%
  group_by(ST) %>%
  summarize(
  Genomes = paste(Genome, collapse = ","),
  Count = n()  
  ) %>%
  ungroup()
View(grouped_STs)

## Save this table 
write.table(grouped_STs, file = "/home/strawberry/Documents/Saureus_genomes/grouped_genomes-ST.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
