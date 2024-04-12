# Genome annotation - Prokka v1.13

## Iterate over each genome file in the directory
for genome_file in /temporario2/11217468/projects/saureus_global/QC_filtered_genome_files/*.fna; do
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
