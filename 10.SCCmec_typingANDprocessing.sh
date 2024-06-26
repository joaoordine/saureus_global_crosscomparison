# SCCmec typing - in silico multiplex PCR - BLAST-detection of specific SCCmec typing primers, according to IWG-SCC - www.sccmec.org

## Create primers file 
vim SCCmec_primers.fasta

## Set the directory containing genome files
genome_dir="/temporario2/11217468/projects/saureus_global/QC_filtered_genome_files"

## Create a BLAST database
makeblastdb -in /temporario2/11217468/projects/saureus_global/sccmec/SCCmec_primers.fasta -dbtype nucl -out db/SCCmecprimer_db

## Set the directory containing the primer database
SCCmecprimer_db="/temporario2/11217468/projects/saureus_global/sccmec/db/SCCmecprimer_db" 

## Iterate through each genome file and run blastn
mkdir BLASTn_output # to store the results 
for genome_file in "${genome_dir}"/*.fna; do
    genome_name=$(basename "$genome_file" _genomic.fna)
    output_file="/temporario2/11217468/projects/saureus_global/sccmec/BLASTn_output/${genome_name}_primer_hits.txt"
    blastn -query "$genome_file" -db "$SCCmecprimer_db" -word_size 18 -out "$output_file" -outfmt "6 qseqid sseqid qstart qend length pident evalue bitscore"
done

## Checking how many empty files I have
find /temporario2/11217468/projects/saureus_global/sccmec/BLASTn_output/ -type f -empty | wc -l # 783

## Processing output - make tabular (True/False) summary of BLASTn output files all genomes
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
