# Pan-genome inferral ------- ROARY v3.13.0
conda install bioconda::roary

## Moving all gff files from Prokka into a single directory 
mkdir all_gff_files
for dir in GCA_*; do
    scp /temporario2/11217468/projects/saureus_global/prokka_output/"$dir"/*.gff all_gff_files && 
    echo "Copy of $dir was sent"
done 

## Inferring the pangenome
mkdir roary_output  
input_dir="/temporario2/11217468/projects/saureus_global/prokka_output/all_gff_files"
roary -e -n -v -f roary_output -p 5 "$input_dir"/*.gff 

## Creating plots 
python roary_plots.py core_gene_alignment.nwk gene_presence_absence.csv

## Analyzing specific genes annotated in the pangenome - check Pangenome sequence analysis in https://github.com/microgenomics/tutorials/blob/master/pangenome.md

