# Phylogenomics analysis --------------------------- PhyloPhlan v3.1.1   - reconstructs strain-level phylogenies from among the closest species using clade specific maximally informative markers
conda install -c bioconda phylophlan

## Setting up database with specific S. aureus marker genes
phylophlan_setup_database \
    -g s__Staphylococcus_aureus \
    -o Saureus_markers_out \
    --verbose 2>&1 | tee logs/phylophlan_setup_database.log
### Counting the number of marker genes identified 
grep -c ">s--Staphylococcus-aureus" Saureus_markers_out/s__Staphylococcus_aureus.faa # 1405

## Generating the configuration file
phylophlan_write_config_file \
    -o Saureus_phylog_config.cfg \
    -d a \
    --force_nucleotides \
    --db_aa diamond \
    --map_aa diamond \
    --map_dna diamond \
    --msa mafft \
    --trim trimal \
    --tree1 fasttree \
    --tree2 raxml

## Creating phylogenetic tree
phylophlan \
    -i QC_filtered_genome_files \
    -o PhyloPhlan_Tree_Output \
    -d Saureus_markers_out \
    --trim greedy \
    --not_variant_threshold 0.99 \
    --remove_fragmentary_entries \
    --fragmentary_threshold 0.67 \
    -t a \
    -f Saureus_phylog_config.cfg \
    --diversity low \
    --force_nucleotides \
    --nproc 20 \
    --verbose 2>&1 | tee logs/phylophlan__output_phylophlan.log

# Obtaining bootstrap values for the phylogenetic tree - since PhyloPhlan does not generate bootstrap values by itself

## Generate Bootstrap Replicates
raxmlHPC-PTHREADS -m GTRCAT -p 51234 -N 100 -s QC_filtered_genome_files_concatenated.aln.reduced -n bootstrap -T 5 -w /temporario2/11217468/projects/saureus_global/phylogeny/PhyloPhlan_Tree_Output/bootstrap-trees

## Perform Bootstrapping with RAxML
raxmlHPC-PTHREADS -m GTRCAT -p 51234 -f b -t RAxML_bestTree.QC_filtered_genome_files_refined.tre -z RAxML_bootstrap.boot -n output_bootstrap.tre -T 5 -w /temporario2/11217468/projects/saureus_global/phylogeny/PhyloPhlan_Tree_Output/bootstrap-trees


# Tree Annotation: IToL




