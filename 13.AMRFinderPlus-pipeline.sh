# AMRFinderPlus v3.10.24 pipeline

# Specify the directory containing your genome files and the output directory
genome_directory="/temporario2/11217468/projects/saureus_global/QC_filtered_genome_files"
output_directory="/temporario2/11217468/projects/saureus_global/amrfinderplus_output"

mkdir amrfinderplus_output

# Update database 
amrfinder -u

# Loop through each genome file in the directory
for genome_file in "$genome_directory"/*.fna; do
    # Get the genome name from the file name
    genome_name=$(basename "$genome_file" .fna)
    
    # Run amrfinder on the genome file and direct the output to the output directory
    amrfinder -n "$genome_file" --organism Staphylococcus_aureus --plus -o "$output_directory/$genome_name.fna"
    
    echo "Processed $genome_name"
done

# Filtering the output files by Element type --> selecting only environmental stress resistance genes (AMR/VF genes will be obtained from ABRICATE)
for file in "$output_directory"/*.fna; do
    output_file="${file%.fna}_filtered.tsv" # each output file will be named according to the file (with .fna extension removed) and added by filtered
    awk -F'\t' 'NR==1 || $9 == "STRESS"' "$file" > "$output_file"
done

mkdir filtered_AMRFind_output
mv *_filtered.tsv filtered_AMRFind_output/


