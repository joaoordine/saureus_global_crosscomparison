# Processing the output files in `filteredAMRFind_output´ - only the columns of interest and create a matrix 

### Create a directory to store the double-filtered output
mkdir -p double_filt_output

### Directory containing my previously filtered AMRFinderPlus output files
filt_output_dir="/temporario2/11217468/projects/saureus_global/amrfinderplus_output/filtered_AMRFind_output"

### Loop
for genome_file in "${filt_output_dir}"/*.tsv; do
    file_name=$(basename "$genome_file" .tsv)
output_file="/temporario2/11217468/projects/saureus_global/amrfinderplus_output/filtered_AMRFind_output/double_filt_output/${file_name}_filt2.tsv"
cut -f2,6,7,9,10,16,17 "$genome_file" > "$output_file"
done

## Combine column 2 of all double_filt.tsv files into a single file
cut -f2 /temporario2/11217468/projects/saureus_global/amrfinderplus_output/filtered_AMRFind_output/double_filt_output/*.tsv > all_genes_amrfinder.tsv
wc -l all_genes_amrfinder.tsv # checking out - 577 environmental stress resistance genes 

## Remove duplicate gene symbols
sort all_genes_amrfinder.tsv | uniq > unique_genes_amrfinder.tsv
wc -l unique_genes_amrfinder.tsv # checking out - 15 unique_genes_amrfinder.tsv

## Extract file names and create the header
header="File $(cat unique_genes_amrfinder.tsv | tr '\n' '\t')"
echo -e "$header" > presence_absence_matrix.tsv
awk -F'\t' '{print NF; exit}' presence_absence_matrix.tsv # checking out how many columns there are to see if there`s no missing gene

## Iterate through double_filt.tsv files
for file in /temporario2/11217468/projects/saureus_global/amrfinderplus_output/filtered_AMRFind_output/double_filt_output/*.tsv; do
    # Get the file name (without extension)
    file_name=$(basename "$file" .tsv)
   
    ### Initialize the matrix line with the file name
    matrix_line="$file_name"

    ### Read each gene symbol from unique_genes_amrfinder.tsv
    while IFS= read -r gene; do
        # Check if the gene is present in the current file (all_genes_amrfinder.tsv)
        if grep -q "$gene" "$file"; then
            matrix_line+="\t1"
        else
            matrix_line+="\t0"
        fi
    done < unique_genes_amrfinder.tsv
   
    echo -e "$matrix_line" >> presence_absence_matrix.tsv
done

grep -c GCA presence_absence_matrix.tsv # checking out how many files I have - 162 genomes 

mv presence_absence_matrix.tsv eSRGs_bool-table.tsv

## Imported table into excel and mannually correct columns and rows with problems before importing it into R - file name "amrfinder-eSRGs_corrected.csv"

