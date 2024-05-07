# Carbohydrate Transporters Exploratory Analysis
source ~/.bashrc
conda init --all
conda activate bioinfo

## Download UNIPROT sequences
### Searched terms (all filtered for "Bacteria")
#wget -O carbohydrate_transporters.fasta.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%22carbohydrate+transporter%22+AND+%28taxonomy_id%3A2%29%29"
#wget -O sugar_transporters.fasta.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%22sugar+transporter%22+AND+%28taxonomy_id%3A2%29%29"
#wget -O sia_transporters.fasta.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%22sialic+acid+transporter%22+AND+%28taxonomy_id%3A2%29%29"
#wget -O polysia_transporters.fasta.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%22polysialic+acid+transporter%22+AND+%28taxonomy_id%3A2%29%29"

## Extracting gunzziped files 
for file in *.gz; do
    gunzip "$file"
done

## Removing duplicated sequences with CD-HIT v4.8.1

mkdir -p ../filtered_sequences
cd uniprot_seqs

for file in *.fasta; do
    output_file="/temporario2/11217468/projects/saureus_global/transporters/filtered_sequences/${file%.*}_filt.fasta"
    cd-hit -i "$file" -o "$output_file" -c 1.0 # using 100% identity as threashold 
done

cd ../filtered_sequences
mkdir -p clusters
mv *.fasta.clstr clusters/

## Multiple sequence alignment of filtered sequences with MAFFT
mkdir -p ../aligned_proteins

for file in *.fasta; do
    filename=$(basename "$file")  # Extract the filename
    filename_no_ext="${filename%_filt.fasta}"  # Remove the file extension
    output_file="/temporario2/11217468/projects/saureus_global/transporters/aligned_proteins/${filename_no_ext}_align.fasta"
    mafft --thread 5 "${file}" > "${output_file}" 
done
cd ../aligned_proteins

## Create the HMM models with HMMER 
mkdir -p ../hmm_profiles

for file in *.fasta; do
    filename=$(basename "$file")  # Extract the filename
    filename_no_ext="${filename%_align.fasta}"  # Remove the file extension
    output_file="/temporario2/11217468/projects/saureus_global/transporters/hmm_profiles/${filename_no_ext}.hmm"
    
    hmmbuild "${output_file}" "${file}"
done

## Search for homologous proteins and processing alignmnet output
mkdir -p ../hmmersearch_output

### Define the directories containing genomes and hmm profiles 
hmm_dir="/temporario2/11217468/projects/saureus_global/transporters/hmm_profiles"
genome_dir="/temporario2/11217468/projects/saureus_global/prokka_output/all_faa_files"

### Loop through hmm files in their directory
for hmm_file in "$hmm_dir"/*.hmm; do
    hmm_name=$(basename "$hmm_file" .hmm)    # Extract protein name (without .hmm extension)
    for dir_file in "$genome_dir"/*.faa; do   # Loop through genome files 
        seq_name=$(basename "$dir_file" .faa) # Extract genome name (without .faa extension)
        # Debug: Print sequence name
        echo "Processing sequence: $seq_name"
        # Define output file path
        output_file="/temporario2/11217468/projects/saureus_global/transporters/hmmersearch_output/${hmm_name}_${seq_name}_output.tsv"
        # Perform hmmsearch and redirect the output
        hmmsearch --cpu 5 --domtblout "$output_file" "$hmm_file" "$dir_file"  
        echo "Executed hmmsearch for $hmm_name against $seq_name. Output saved to $output_file."
    done
done

### Checking out how many hits I got 
ls | grep -c .tsv # 17140







































