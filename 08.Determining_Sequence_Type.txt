# Multi-locus sequence typing - mlst v2.16.1

## Installing package 
source ~/.bashrc
conda init --all
conda activate bioinfo
 conda install bioconda::mlst
 
## Getting S. aureus pubMLST scheme 
mlst --list | grep saureus 

## Running the tool
for file in /temporario2/11217468/projects/saureus_global/QC_filtered_genome_files/*.fna; do
    mlst --scheme saureus "$file" >> /temporario2/11217468/projects/saureus_global/mlst_results.txt
done

## Processing mlst output - splitting it into two separate files (one with genome ID and its corresponding ST and another with the genome ID and a presence/absence matrix of all genes detected during MLST)

### Genome ID and ST
awk -F'\t' '/temporario/ { print $1, $3 }' mlst_results.txt > genome_STs.txt
sed -i 's|.*/\([^/]*\)_genomic\.fna|\1|' genome_STs.txt # searches for lines in your file that match the pattern */<identifier>_genomic.fna and replaces the entire line with just the <identifier>

### Genome ID and allelle ID
awk -F'\t' '/temporario/ { print $1, $4, $5, $6, $7, $8, $9, $10 }' mlst_results.txt > genome_allele_ST.txt
sed -i 's|.*/\([^/]*\)_genomic\.fna|\1|' genome_allele_ST.txt

### Frequency table/plot of each ST identified
sort -k2n genome_STs.txt -o genome_STs_sorted.txt
awk '{print $2}' genome_STs_sorted.txt | uniq -c > just_STs.txt
vim header # Count Sequence_Type
cat header.txt just_STs.txt > just_STs_W_header.txt
sed -i 's/-/Unclassified/g' just_STs_W_header.txt # replacing '-' character for something else so there's isn't any problems importing the dataset into R
awk 'NR <= 2 {print; next} {print $1, "ST_"$2}' just_STs_W_header.txt > just_STs_W_header_R.txt # skip the first two rows, print them as they are, and for the rest of the rows, it will add 'ST_' in front of the numbers in the second column


