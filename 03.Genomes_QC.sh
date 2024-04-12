# Quality control of downloaded genomes -checkM v1.2.2
source /temporario2/11217468/miniconda3/bin/activate checkm
checkm lineage_wf genome_files checkm_output -t 5

## Output visualization: more easily visualized in the end of the slurm output file (slurm-JOBID) or also in the 'bin_stats_ext.tsv' (contains along with completeness and contamination parameters, other metrics related to the quality of the assembled genome, such as N50, mean contig, etc)
awk -F',' '/^GCA_/ { print $1, $11, $12 }' bin_stats_ext.tsv > genome_quality.txt # extracting only those columns I'm interested. Note that I used comma as a separator, which will give me as the first column not only the genome GCA identification, but also the taxonomic annotation, which we'll remove with the command below
awk -F' ' '/^GCA_/ { print $1, $6, $8 }' genome_quality.txt > genome_quality_filtered.txt
echo "Genome_ID Completeness Contamination" > genome_quality_filt_Whead.txt # creating a file with the headers
cat genome_quality_filtered.txt >>  genome_quality_filt_Whead.txt

# Quality control of downloaded genomes -QUAST v5.2.0
conda create --name quast
conda activate quast
conda install -c bioconda quast
/temporario2/11217468/miniconda3/envs/quast/bin/quast.py *_genomic.fna



