# Discarding genomes that didn't meet the quality criteria 
awk '{print $1".fna", $2, $3}' genome_quality_filt_Whead.txt > fna_quality_filtered_genomes.txt # adding .fna in the end of each Genome_ID so that it`ll have the same name as its corresponding file

## Set the criteria values
completeness_cutoff=97
contamination_cutoff=3
# Loop through each line in the quality filtered genomes file
while read -r genome completeness contamination; do
    # Check if completeness is greater than 97 and contamination is less than 3
    if (( $(echo "$completeness > $completeness_cutoff" | bc -l) )) && (( $(echo "$contamination < $contamination_cutoff" | bc -l) )); then
        mv "genome_files/$genome" "QC_filtered_genome_files" # Copy the genome file from genome_files to the new directory
        echo "Copied $genome to QC_filtered_genome_files"
    fi
done < "fna_quality_filtered_genomes.txt"

## Count how many high quality genomes there are
cd QC_filtered_genome_files
ls | grep -c GCA_


