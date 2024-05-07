# Carbohydrate Transporters exploratatory analyses - processing HMMsearch outputs

#!/bin/bash

mkdir -p raw_files
mv *_output.tsv raw_files

## Remove rows starting with # 
for file in raw_files/sia*.tsv; do
    grep -v '^#' "$file" > "${file%.tsv}_filtered.tsv"    
done
cd raw_files
mv -f *_filtered.tsv ../
cd ..

## Making all elements in the output separed by tabs 
for file in *.tsv; do
    # Replace spaces with tabs in each file
    sed -i 's/ \+/\t/g' "$file"
done

## Adding the file name as last column for each row
for file in *_filtered.tsv; do
    genome_code="${file%%_genomic_output_filtered.tsv}" # Extract genome code from the file name
    out_file="${file%_genomic_output_filtered.tsv}_filt2.tsv"   
    awk -v genome_code="$genome_code" -F'\t' 'BEGIN {OFS="\t"} {print $0, genome_code}' "$file" > "$out_file" # Append file name as a new column to each row
done

## Filter out for columns of interest - cols 1,3,4,6,7,8 
for file in ./*_filt2.tsv; do
    filtered_file="${file%_filt2.tsv}_filt3.tsv"
    awk -v OFS='\t' '{print $1, $3, $4, $6, $7, $8, $NF}' "$file" > "$filtered_file" # Filter for columns of interest and the last column, then save to filtered_file filtered_file
done

## Concatenate files
cat *_filt3.tsv > carb_transporters_concat.tsv

## Make a header and append it to the concatenated file
echo -e "Target_Code\tTarget_Length\tTransporter\tTransporter_Length\tE_Value\tScore\tGenome" > hmm_header.tsv
cat hmm_header.tsv carb_transporters_concat.tsv > carb_transporters_processed.tsv

