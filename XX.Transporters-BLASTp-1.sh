# Carbohydrate Transporters exploratatory analyses - BLASTp 

## Concatenate all sugar transporter files into a single file
cat *.fasta > concat_transporters.fasta
scp concat_transporters.fasta ../../transporters_blast/ # did the same with annotated files from prokka for all genomes

## Create a BLAST database from the combined file
mkdir db_sugarT
makeblastdb -in concat_transporters.fasta -dbtype prot -out sugarT_db
mv sugarT_db* db_sugarT/ 

## Specifying the path to genomes database
export BLASTDB="/temporario2/11217468/projects/saureus_global/transporters_blast/db_sugarT"

## Running blastn on multiple queries
mkdir -p blastp_output

## Setting up the path to our protein sequences and to our genomic database
genome_dir="/temporario2/11217468/projects/saureus_global/transporters_blast/all_faa_files"
transporter_db="/temporario2/11217468/projects/saureus_global/transporters_blast/db_sugarT/sugarT_db"

cd "$genome_dir"

for genome_file in "${genome_dir}"/*.faa; do
genome_name=$(basename "$genome_file" _genomic.faa)
genome_no_ext="${genome_name%_genomic.faa}" # Remove the file extension
output_file="/temporario2/11217468/projects/saureus_global/transporters_blast/blastp_output/${genome_no_ext}_hits.tsv"
blastp -query "$genome_file" -db "$transporter_db" -out "$output_file" -evalue 0.01 -num_threads 10 -outfmt "6 qseqid sseqid sstart send length pident evalue bitscore sseq"
done

# Processing BLASTp outputs

## Adding the genome code as last column for each row
for file in *_hits.tsv; do
    genome_code="${file%%_hits.tsv}" # Extract genome code from the file name
    out_file="${file%%_hits.tsv}_filt.tsv"   
    awk -v genome_code="$genome_code" -F'\t' 'BEGIN {OFS="\t"} {print $0, genome_code}' "$file" > "$out_file" # Append file name as a new column to each row
done

## Concatenate files
cat *_filt.tsv > carb_transporters_concat.tsv

## Make a header and append it to the concatenated file
echo -e "Annotated_Protein\tTransporter_DB\tStart\tEnd\tLength\\tIdentity\tE_Value\tBitscore\tSequence\tGenome" > blast_header.tsv
cat blast_header.tsv carb_transporters_concat.tsv > carb_transporters_processed.tsv

## Import file into R and keep processing 



















