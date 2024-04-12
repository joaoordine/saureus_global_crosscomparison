# Retrieving S. aureus complete genomes from NCBI - March 18th 2024

## Getting the accessions list and processing it
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt # Downloading all bacterial assembled genomes from NCBI
wget  ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt	# Checking which ones are my columns of interest: 1, 2, 3, 6, 7, 8, 9, 11, 12, 20, 26, 27, 28, 35, 36
cut -f1,2,3,6,7,8,9,11,12,20,26,27,28,35,36 assembly_summary.txt > filtered_genomes.txt ### selected only the columns with useful information and created a separated file 

## Selecting only the complete genomes and reasuring that all genomes there are actually only the complete ones
head ./filtered_genomes.txt
grep -c 'Complete' filtered_genomes.txt # 49196 complete bacterial genomes
grep 'Complete' ./filtered_genomes.txt >> ./complete_genomes.txt      
grep -c 'Complete' complete_genomes.txt # 49196
grep -c 'Scafold' complete_genomes.txt # 0
grep -c 'Contig' complete_genomes.txt # 0

## Selecting only S. aureus genomes
grep -c 'Staphylococcus aureus' complete_genomes.txt ### checking how many S. aureus genomes there are on the file -> # 1815
grep 'Staphylococcus aureus' ./complete_genomes.txt >> ./Saureus_complete_genomes.txt
grep -c 'Staphylococcus aureus' Saureus_complete_genomes.txt # 1815 
grep -c 'Scafold' Saureus_complete_genomes.txt # 0
grep -c 'Contig' Saureus_complete_genomes.txt # 0
grep -c 'Complete' Saureus_complete_genomes.txt # 1815

## Adding header  
head -2 filtered_genomes.txt > header ### creating a file that contains only the header of my columns of interest 
cat header Saureus_complete_genomes.txt > Saureus_complete_genomesWhead.txt 

## Creating a loop to read the FTP paths from the file and download all corresponding genomes 

### Moving all ftp paths to a separated file  
while IFS=$'\t' read -r _ _ _ _ _ _ _ _ _ ftp_path _ _ _ _ _; do     
dir_path="${ftp_path%/*}";     
filename="${ftp_path##*/}";     
modified_path="$dir_path/$filename/${filename}_genomic.fna.gz";
echo "$modified_path" >> modified_paths.txt; done < Saureus_complete_genomesWhead.txt 
sed '1,2d' modified_paths.txt > new_modified_paths.txt # removing the first two lines (only headers, don't contain any ftp paths)
wc -l new_modified_paths.txt # 1815

### Read the modified FTP paths file and download genomes
while IFS= read -r ftp_path; do
    wget "$ftp_path"; done < new_modified_paths.txt

### Extract gzip files and move them to another directory 
mkdir genome_files
for file in ./*_genomic.fna.gz; do
    gzip -d "$file"    
    filename="${file%.*}"    
    mv "$filename" genome_files/; done

### Checking how many files were downloaded and moved to genome_files
ls | grep -c GCA # 1813












