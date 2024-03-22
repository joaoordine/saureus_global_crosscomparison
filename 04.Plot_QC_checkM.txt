# Plotting genomes QC in R - scatter plots: completeness x contamination

## import dataset into Rstudio
library(ggplot2)
genomes <- genome_quality_filt_Whead

## Set cutoff values
completeness_cutoff <- 97
contamination_cutoff <- 3

## Create a new column for color based on conditions
genomes$Color <- ifelse(genomes$Completeness > completeness_cutoff & genomes$Contamination < contamination_cutoff, "Above Cutoff", "Below Cutoff")

## Figure 1A
### Create scatter plot
A <- ggplot(genomes, aes(x = Completeness, y = Contamination, color = Color)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Above Cutoff" = "#8A9A5B", "Below Cutoff" = "#FF00FF")) +
  labs(title = "Completeness vs. Contamination Scatter Plot", x = "Completeness", y = "Contamination") +
  theme_minimal() +
  expand_limits(x = 0, y = 0)  # Set the axes to start at zero

## Figure 1B
### Create scatter plot
B <- ggplot(genomes, aes(x = Completeness, y = Contamination, color = Color)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Above Cutoff" = "#8A9A5B", "Below Cutoff" = "#FF00FF")) +
  labs(title = "Completeness vs. Contamination Scatter Plot", x = "Completeness", y = "Contamination") +
  theme_minimal()
  
## Combine figures 
options(repr.plot.width=12, repr.plot.height=12)

library(ggpubr)

ggarrange(A, B, 
          labels = c("A", "B"),
          font.label = list(size = 18),
          ncol = 2, nrow = 1)
