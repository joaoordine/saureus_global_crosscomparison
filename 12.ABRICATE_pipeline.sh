# ABRICATE analysis 
mkdir abricate_output

## Specify the directory containing your genome files
genome_directory="/temporario2/11217468/projects/saureus_global/QC_filtered_genome_files"

## Get the genome name from the file name and Run ABRICATE on the genome file
output_directory="/temporario2/11217468/projects/saureus_global/abricate_output"

### Retrieving ARGs aligning against CARD database
for genome_file in "$genome_directory"/*.fna; do
    genome_name=$(basename "$genome_file" .fna)
    abricate --db card "$genome_file" > "$output_directory/$genome_name-ARGs.txt"
    && echo "Analyzed $genome_name"
done
abricate --summary *ARGs.txt > abricate-ARGs_summary.txt # sum up ABRICATE results in a boolean table with coverage of each ARG in that genome

### Retrieving VGs aligning against VFDB database
for genome_file in "$genome_directory"/*.fna; do
    genome_name=$(basename "$genome_file" .fna)
    abricate --db vfdb "$genome_file" > "$output_directory/$genome_name-VGs.txt"
    && echo "Analyzed $genome_name"
done
abricate --summary *VGs.txt > abricate-VGs_summary.txt 

### Retrieving plasmids aligning against PlasmidFinder database
for genome_file in "$genome_directory"/*.fna; do
    genome_name=$(basename "$genome_file" .fna)
    abricate --db plasmidfinder "$genome_file" > "$output_directory/$genome_name-plasmids.txt"
    && echo "Analyzed $genome_name"
done
abricate --summary *plasmids.txt > abricate-plasmids_summary.txt 





