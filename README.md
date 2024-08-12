# Global cross-genome comparison of carbohydrate transporters in Staphylococcus aureus with distinct antimicrobial susceptibility profiles
We aimed to cross-compare genetic elements related to AMR and carbohydrate transport in the global S. aureus genomic population. Here, you'll find all scripts used during the development of the project and its analyses. 
Manuscript pre-print: (replace later with article doi) 

# Genome selection and curation

## 1) Genome selection 

### Getting the accessions list and processing it
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt # Downloading all bacterial assembled genomes from NCBI 
wget  ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt	# Checking which ones are my columns of interest: 1, 2, 3, 6, 7, 8, 9, 11, 12, 20, 26, 27, 28, 35, 36
cut -f1,2,3,6,7,8,9,11,12,20,26,27,28,35,36 assembly_summary.txt > filtered_genomes.txt ### selected only the columns with useful information and created a separated file
```

### Selecting only the complete genomes and reassuring that all genomes there are actually only the complete ones
```
head ./filtered_genomes.txt 
grep -c 'Complete' filtered_genomes.txt # 49196 complete bacterial genomes  
grep 'Complete' ./filtered_genomes.txt >> ./complete_genomes.txt    
grep -c 'Complete' complete_genomes.txt # 49196  
grep -c 'Scafold' complete_genomes.txt # 0  
grep -c 'Contig' complete_genomes.txt # 0
```

### Selecting only S. aureus genomes
```
grep -c 'Staphylococcus aureus' complete_genomes.txt ### checking how many S. aureus genomes there are on the file -> # 1815  
grep 'Staphylococcus aureus' ./complete_genomes.txt >> ./Saureus_complete_genomes.txt  
grep -c 'Staphylococcus aureus' Saureus_complete_genomes.txt # 1815  
grep -c 'Scafold' Saureus_complete_genomes.txt # 0  
grep -c 'Contig' Saureus_complete_genomes.txt # 0  
grep -c 'Complete' Saureus_complete_genomes.txt # 1815
```

### Adding header  
```
head -2 filtered_genomes.txt > header ### creating a file that contains only the header of my columns of interest  
cat header Saureus_complete_genomes.txt > Saureus_complete_genomesWhead.txt
```

### Creating a loop to read the FTP paths from the file and download all corresponding genomes 

#### Moving all ftp paths to a separated file  
```
while IFS=$'\t' read -r _ _ _ _ _ _ _ _ _ ftp_path _ _ _ _ _; do      
 dir_path="${ftp_path%/*}";      
 filename="${ftp_path##*/}";      
 modified_path="$dir_path/$filename/${filename}_genomic.fna.gz";
 echo "$modified_path" >> modified_paths.txt; done < Saureus_complete_genomesWhead.txt 
 sed '1,2d' modified_paths.txt > new_modified_paths.txt # removing the first two lines (only headers, don't contain any ftp paths) 
 wc -l new_modified_paths.txt # 1815
``` 

#### Read the modified FTP paths file and download genomes
```
while IFS= read -r ftp_path; do 
     wget "$ftp_path"; done < new_modified_paths.txt
``` 

#### Extract gzip files and move them to another directory 
```
mkdir genome_files 
 for file in ./*_genomic.fna.gz; do 
     gzip -d "$file"     
     filename="${file%.*}"     
     mv "$filename" genome_files/; done
``` 

#### Checking how many files were downloaded and moved to genome_files
``` 
ls | grep -c GCA # 1813
``` 

### Setting up my conda environments

#### Initiating miniconda 
```
export PATH="/temporario/11217468/miniconda3/bin:$PATH" 
 conda init bash 
 conda create --name bioinfo 
 conda activate bioinfo
``` 

#### Downloading required packages available through bioconda (all have python version 3.8)
```
conda install bioconda::blast 
 conda install bioconda::prokka 
 conda install bioconda::mlst 
 conda deactivate
``` 

#### Creating separated environments for specific packages I know will conflict with those aforementioned
```
conda create --name checkm
conda activate checkm
conda install -c bioconda checkm-genome
conda deactivate
``` 

## 2) Quality control of downloaded genomes -checkM v1.2.2
``` 
source </temporario2/11217468/miniconda3/bin/activate checkm> # replace path accordingly
checkm lineage_wf genome_files checkm_output -t 5
``` 

### Output visualization*
```
awk -F',' '/^GCA_/ { print $1, $11, $12 }' bin_stats_ext.tsv > genome_quality.txt # extracting only those columns I'm interested. Note that I used comma as a separator, which will give me as the first column not only the genome GCA identification, but also the taxonomic annotation, which we'll remove with the command below
awk -F' ' '/^GCA_/ { print $1, $6, $8 }' genome_quality.txt > genome_quality_filtered.txt
echo "Genome_ID Completeness Contamination" > genome_quality_filt_Whead.txt # creating a file with the headers
cat genome_quality_filtered.txt >>  genome_quality_filt_Whead.txt
``` 
* more easily visualized in the end of the slurm output file (slurm-JOBID) or also in the 'bin_stats_ext.tsv' (contains along with completeness and contamination parameters, other metrics related to the quality of the assembled genome, such as N50, mean contig, etc)


## 3) Selecting high-quality genomes
``` 
awk '{print $1".fna", $2, $3}' genome_quality_filt_Whead.txt > fna_quality_filtered_genomes.txt # Discarding genomes that didn't meet the quality criteria ; adding .fna in the end of each Genome_ID so that it`ll have the same name as its corresponding file
```

### Set the criteria values
```
completeness_cutoff=97
contamination_cutoff=3

while read -r genome completeness contamination; do
    # Check if completeness is greater than 97 and contamination is less than 3
    if (( $(echo "$completeness > $completeness_cutoff" | bc -l) )) && (( $(echo "$contamination < $contamination_cutoff" | bc -l) )); then
        mv "genome_files/$genome" "QC_filtered_genome_files" # Copy the genome file from genome_files to the new directory
        echo "Copied $genome to QC_filtered_genome_files"
    fi
done < "fna_quality_filtered_genomes.txt" # Loop through each line in the quality filtered genomes file
```

### Count how many high quality genomes there are
```
cd QC_filtered_genome_files
ls | grep -c GCA_
```

## 4) QC plot
Check script 01.Plot_QC_checkM.ipynb

## 5) Genome annotation - Prokka v1.13
```
for genome_file in /temporario2/11217468/projects/saureus_global/QC_filtered_genome_files/*.fna; do # Iterate over each genome file in the directory
    genome_name=$(basename "$genome_file" .fna)     # Get the genome name from the file name

    index=1 # Initialize an index to make the output directory name unique
    output_dir="${genome_name}_prokka_out" # Construct the initial output directory name
    
    ### Check if the output directory already exists, and if so, increment the index. Then, run prokka 
    while [ -e "$output_dir" ]; do
        index=$((index + 1))
        output_dir="${genome_name}_prokka_out_$index"
    done

    prokka --outdir prokka_output/"$output_dir" --prefix "$genome_name" "$genome_file" --cpus 0 && echo "Annotated $genome_name"
done
```

## 6) Multi-locus sequence typing - mlst v2.16.1
 
```
mlst --list | grep saureus # Getting S. aureus pubMLST scheme 

for file in /temporario2/11217468/projects/saureus_global/QC_filtered_genome_files/*.fna; do
    mlst --scheme saureus "$file" >> /temporario2/11217468/projects/saureus_global/mlst_results.txt
done
```

### Processing mlst output 
Splitting output into two separate files (one with genome ID and its corresponding ST and another with the genome ID and a presence/absence matrix of all genes detected during MLST)

```
awk -F'\t' '/temporario/ { print $1, $3 }' mlst_results.txt > genome_STs.txt
sed -i 's|.*/\([^/]*\)_genomic\.fna|\1|' genome_STs.txt # searches for lines in your file that match the pattern */<identifier>_genomic.fna and replaces the entire line with just the <identifier>
```
``` 
awk -F'\t' '/temporario/ { print $1, $4, $5, $6, $7, $8, $9, $10 }' mlst_results.txt > genome_allele_ST.txt
sed -i 's|.*/\([^/]*\)_genomic\.fna|\1|' genome_allele_ST.txt
```

####  Frequency table/plot of each ST identified
```
sort -k2n genome_STs.txt -o genome_STs_sorted.txt
awk '{print $2}' genome_STs_sorted.txt | uniq -c > just_STs.txt
vim header # Count Sequence_Type
cat header.txt just_STs.txt > just_STs_W_header.txt
sed -i 's/-/Unclassified/g' just_STs_W_header.txt # replacing '-' character for something else so there's isn't any problems importing the dataset into R
awk 'NR <= 2 {print; next} {print $1, "ST_"$2}' just_STs_W_header.txt > just_STs_W_header_R.txt # skip the first two rows, print them as they are, and for the rest of the rows, it will add 'ST_' in front of the numbers in the second column
```

## 7) MLST plot 
Check script 02.Visualizing_MLST_output.ipynb

### Grouping genomes based on ST
Check script 03. Group_genomes_by_ST.ipynb

# Genomic analyses of antibiotic resistance

## 8) SCCmec typing - in silico multiplex PCR - BLAST-detection of specific SCCmec typing primers, according to IWG-SCC - www.sccmec.org
```
vim SCCmec_primers.fasta # Create primers file 

genome_dir="/temporario2/11217468/projects/saureus_global/QC_filtered_genome_files" # Set the directory containing genome files

makeblastdb -in /temporario2/11217468/projects/saureus_global/sccmec/SCCmec_primers.fasta -dbtype nucl -out db/SCCmecprimer_db # Create a BLAST database

SCCmecprimer_db="/temporario2/11217468/projects/saureus_global/sccmec/db/SCCmecprimer_db" # Set the directory containing the primer database

mkdir BLASTn_output # to store the results 
for genome_file in "${genome_dir}"/*.fna; do
    genome_name=$(basename "$genome_file" _genomic.fna)
    output_file="/temporario2/11217468/projects/saureus_global/sccmec/BLASTn_output/${genome_name}_primer_hits.txt"
    blastn -query "$genome_file" -db "$SCCmecprimer_db" -word_size 18 -out "$output_file" -outfmt "6 qseqid sseqid qstart qend length pident evalue bitscore"
done
```

### Checking how many empty files I have
```
find /temporario2/11217468/projects/saureus_global/sccmec/BLASTn_output/ -type f -empty | wc -l # 783
```

## 9) Processing output 
Make tabular (True/False) summary of BLASTn output files all genomes

```
echo -e "Sample\tTYPE_I\tTYPE_II\tTYPE_III\tTYPE_IVa\tTYPE_IVb\tTYPE_IVc\tTYPE_IVd\tTYPE_IVe\tTYPE_V\tccr4\tmecA147" > //temporario2/11217468/projects/saureus_global/sccmec//summary_SCCmectable.txt # Create the header 

for blastn_result in /temporario2/11217468/projects/saureus_global/sccmec/BLASTn_output/*.txt; do
    ## Extract the sample name from the filename
    sample=$(basename "$blastn_result" | sed 's/_primer_hits\.txt//')

    ### Initialize variables for each column to False
    TYPE_I="False"
    TYPE_II="False"
    TYPE_III="False"
    TYPE_IVa="False"
    TYPE_IVb="False"
    TYPE_IVc="False"
    TYPE_IVd="False"
    TYPE_IVe="False"
    TYPE_V="False"
    ccr4="False"
    mecA147="False"

    ### Check if hits are found in the second column of the blastn result
    if grep -q -e "TYPE_I" "$blastn_result"; then
        TYPE_I="True"
    fi

    if grep -q -e "TYPE_II" "$blastn_result"; then
        TYPE_II="True"
    fi

    if grep -q -e "TYPE_III" "$blastn_result"; then
        TYPE_III="True"
    fi

    if grep -q -e "TYPE_IVa" "$blastn_result"; then
        TYPE_IVa="True"
    fi

    if grep -q -e "TYPE_IVb" "$blastn_result"; then
        TYPE_IVb="True"
    fi

    if grep -q -e "TYPE_IVc" "$blastn_result"; then
        TYPE_IVc="True"
    fi

    if grep -q -e "TYPE_IVd" "$blastn_result"; then
        TYPE_IVd="True"
    fi

    if grep -q -e "TYPE_IVe" "$blastn_result"; then
        TYPE_IVe="True"
    fi

    if grep -q -e "TYPE_V" "$blastn_result"; then
        TYPE_V="True"
    fi

    if grep -q -e "ccr4" "$blastn_result"; then
        ccr4="True"
    fi

    if grep -q -e "mecA147" "$blastn_result"; then
        mecA147="True"
    fi

    ### Output the tabular summary for this sample
    echo -e "$sample\t$TYPE_I\t$TYPE_II\t$TYPE_III\t$TYPE_IVa\t$TYPE_IVb\t$TYPE_IVc\t$TYPE_IVd\t$TYPE_IVe\t$TYPE_V\t$ccr4\t$mecA147" >> /temporario2/11217468/projects/saureus_global/sccmec/summary_SCCmectable.txt
done

sed 's/\t/,/g' summary_SCCmectable.txt > summary_SCCmectable.csv # to make it more easily to import in R/excel 
```

## 9) Visualizing SCCmec types across S. aureus genomes in R
Check script 04.Visualizing_SCCmec_groups.R

## 10) ABRICATE analysis 
```
mkdir abricate_output
```

### Specify the directory containing your genome files
```
genome_directory="/temporario2/11217468/projects/saureus_global/QC_filtered_genome_files"
```

### Get the genome name from the file name and Run ABRICATE on the genome file
```
output_directory="/temporario2/11217468/projects/saureus_global/abricate_output"
```

### Retrieving ARGs aligning against CARD database
```
for genome_file in "$genome_directory"/*.fna; do
    genome_name=$(basename "$genome_file" .fna)
    abricate --db card "$genome_file" > "$output_directory/$genome_name-ARGs.txt"
    && echo "Analyzed $genome_name"
done
abricate --summary *ARGs.txt > abricate-ARGs_summary.txt # sum up ABRICATE results in a boolean table with coverage of each ARG in that genome
```

### Retrieving VGs aligning against VFDB database
```
for genome_file in "$genome_directory"/*.fna; do
    genome_name=$(basename "$genome_file" .fna)
    abricate --db vfdb "$genome_file" > "$output_directory/$genome_name-VGs.txt"
    && echo "Analyzed $genome_name"
done
abricate --summary *VGs.txt > abricate-VGs_summary.txt 
```

### Retrieving plasmids aligning against PlasmidFinder database
```
for genome_file in "$genome_directory"/*.fna; do
    genome_name=$(basename "$genome_file" .fna)
    abricate --db plasmidfinder "$genome_file" > "$output_directory/$genome_name-plasmids.txt"
    && echo "Analyzed $genome_name"
done
abricate --summary *plasmids.txt > abricate-plasmids_summary.txt 
```

## 11) Visualizing ABRICATE summary output in R
Check script 05.Visualizing_ABRICATE_output.ipynb

## Genomic analyses of carbohydrate transporters 

## 12) Carbohydrate Transporters exploratatory analyses - BLASTp 

### Concatenate all sugar transporter files into a single file
```
cat *.fasta > concat_transporters.fasta
scp concat_transporters.fasta ../../transporters_blast/ # did the same with annotated files from prokka for all genomes
```

### Create a BLAST database from the combined file
```
mkdir db_sugarT
makeblastdb -in concat_transporters.fasta -dbtype prot -out sugarT_db
mv sugarT_db* db_sugarT/ 
```

### Specifying the path to genomes database
```
export BLASTDB="/temporario2/11217468/projects/saureus_global/transporters_blast/db_sugarT"
```

### Running blastn on multiple queries
```
mkdir -p blastp_output
```

### Setting up the path to our protein sequences and to our genomic database
``` genome_dir="/temporario2/11217468/projects/saureus_global/transporters_blast/all_faa_files"
transporter_db="/temporario2/11217468/projects/saureus_global/transporters_blast/db_sugarT/sugarT_db"
```

``` cd "$genome_dir"

for genome_file in "${genome_dir}"/*.faa; do
genome_name=$(basename "$genome_file" _genomic.faa)
genome_no_ext="${genome_name%_genomic.faa}" # Remove the file extension
output_file="/temporario2/11217468/projects/saureus_global/transporters_blast/blastp_output/${genome_no_ext}_hits.tsv"
blastp -query "$genome_file" -db "$transporter_db" -out "$output_file" -evalue 0.01 -num_threads 10 -outfmt "6 qseqid sseqid sstart send length pident evalue bitscore sseq"
done
```

## 13) Processing BLASTp outputs

### Adding the genome code as last column for each rowbox
``` for file in *_hits.tsv; do
    genome_code="${file%%_hits.tsv}" # Extract genome code from the file name
    out_file="${file%%_hits.tsv}_filt.tsv"   
    awk -v genome_code="$genome_code" -F'\t' 'BEGIN {OFS="\t"} {print $0, genome_code}' "$file" > "$out_file" # Append file name as a new column to each row
done
``` 
### Concatenate files
``` cat *_filt.tsv > carb_transporters_concat.tsv
```

### Make a header and append it to the concatenated file
```
echo -e "Annotated_Protein\tTransporter_DB\tStart\tEnd\tLength\\tIdentity\tE_Value\tBitscore\tSequence\tGenome" > blast_header.tsv
cat blast_header.tsv carb_transporters_concat.tsv > carb_transporters_processed.tsv
```

### Import file into R and keep processing 
Check script 06.Transporters-BLASTp-2.ipynb

# Phylogenomic and pangenome analyses

## 14) Phylogenomics analysis - PhyloPhlan v3.1.1   
Reconstructs strain-level phylogenies from among the closest species using clade specific maximally informative markers

```
conda install -c bioconda phylophlan
```

### Setting up database with specific S. aureus marker genes
```
phylophlan_setup_database \
    -g s__Staphylococcus_aureus \
    -o Saureus_markers_out \
    --verbose 2>&1 | tee logs/phylophlan_setup_database.log
#### Counting the number of marker genes identified 
grep -c ">s--Staphylococcus-aureus" Saureus_markers_out/s__Staphylococcus_aureus.faa # 1405
```

### Generating the configuration file
```
phylophlan_write_config_file \
    -o Saureus_phylog_config.cfg \
    -d a \
    --force_nucleotides \
    --db_aa diamond \
    --map_aa diamond \
    --map_dna diamond \
    --msa mafft \
    --trim trimal \
    --tree1 fasttree \
    --tree2 raxml
```

### Creating phylogenetic tree
```
phylophlan \
    -i QC_filtered_genome_files \
    -o PhyloPhlan_Tree_Output \
    -d Saureus_markers_out \
    --trim greedy \
    --not_variant_threshold 0.99 \
    --remove_fragmentary_entries \
    --fragmentary_threshold 0.67 \
    -t a \
    -f Saureus_phylog_config.cfg \
    --diversity low \
    --force_nucleotides \
    --nproc 20 \
    --verbose 2>&1 | tee logs/phylophlan__output_phylophlan.log
```

## 15) Reinferring ML tree from reduced alignment and computing bootstrap - improve the accuracy and phylogenetic representation of the tree - IQ-TREE multicore v2.1.4

```
iqtree -s QC_filtered_genome_files_concatenated.aln.reduced -m TEST -nt AUTO -seed 2585 -B 1000 # Standard model selection by ModelFinderPlus followed by tree inference, automatic detection of best number of threads to be used and 1,000 Replicates for ultrafast bootstrap
```

- Best-fit model: GTR+F+G4 chosen according to BIC (Bayesian Information Criterion) score

## 16) Tree Annotation: IToL
Scripts used for tree annotation in IToL are in the "data/itol_annotations-scripts" folder

## 17) Pan-resistome inferral - PanViTa v1.0 
Clustermap with data grouped with Euclidean distance measure  

### Installing tool and its dependencies 
```
git clone https://github.com/dlnrodrigues/panvita.git
cd panvita
python3 panvita.py -u
python3 panvita.py -h
```

### Moving all gbk files from Prokka into a single directory 
```
cd prokka_output
mkdir all_gbk_files
for dir in GCA_*; do
    scp /temporario2/11217468/projects/saureus_global/prokka_output/"$dir"/*.gbk all_gbk_files && 
    echo "Copy of $dir was sent"
done 
```

### If any dependencie isn't installed automatically, run the following
```
####conda install -c anaconda wget
####conda install seaborn
####conda install pandas
####conda install -c conda-forge matplotlib
####conda install -c anaconda basemap
```

### Run the tool 
```
input_dir="/temporario2/11217468/projects/saureus_global/prokka_output/all_gbk_files"

python3 panvita.py -card -vfdb -bacmet -i 80 -c 80 "$input_dir"/*.gbk 
```

## 18) Pan-genome inferral ------- ROARY v3.13.0
```
conda install bioconda::roary
```

### Moving all gff files from Prokka into a single directory 
```
mkdir all_gff_files
for dir in GCA_*; do
    scp /temporario2/11217468/projects/saureus_global/prokka_output/"$dir"/*.gff all_gff_files && 
    echo "Copy of $dir was sent"
done 
```

### Inferring the pangenome
```
mkdir roary_output  
input_dir="/temporario2/11217468/projects/saureus_global/prokka_output/all_gff_files"
roary -e -n -v -f roary_output -p 5 "$input_dir"/*.gff 
```

### Creating plots 
```
python roary_plots.py core_gene_alignment.nwk gene_presence_absence.csv
```

## Analyzing specific genes annotated in the pangenome - check Pangenome sequence analysis in https://github.com/microgenomics/tutorials/blob/master/pangenome.md

## 19) Replotting PanViTa output 
Check script 07.Replotting_Panvita-output.ipynb

## 20) Plot Roary outputs in R
Check script 08.Plot_roary_output.ipynb

## 21) Process tables to annotate phylogenomic tree
Check script 10.Process_table_tree_annot.ipynb
















 
