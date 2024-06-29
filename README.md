# Global cross-genome comparison of carbohydrate transporters in Staphylococcus aureus with distinct antimicrobial susceptibility profiles
We aimed to cross-compare genetic elements related to AMR and carbohydrate transport in the global S. aureus genomic population. Here, you'll find all scripts used during the development of the project and its analyses.  

## Genome selection and curation
Scripts: 
- 01.DownloadingANDprocessing_genomes.sh
- 02.Settingup_conda_envs.sh
- 03.Genomes_QC.sh
- 04.Plot_QC_checkM.R
- 05.Selecting_highquality_genomes.sh
- 06.Genome_annotation.sh
- 07.Determining_Sequence_Type.sh
- 08.Visualizing_MLST_output.R
- 09.Group_genomes_by_ST.R

## Genomic analyses of antibiotic resistance
Scripts:
- 10.SCCmec_typingANDprocessing.sh
- 11.Visualizing_SCCmec_groups.R
- 12.ABRICATE_pipeline.sh
- 13.Visualizing_ABRICATE_output.txt

## Genomic analyses of carbohydrate transporters 
Scripts:
- 14.Transporters-BLASTp-1.sh
- 15.Transporters-BLASTp-2.R

## Phylogenomic and pangenome analyses
Scripts:
- 16.Phylogenomics-analysis.sh
- 17.Panresistome_infer.sh
- 18.Pangenome_infer.sh
- 19.Replotting_Panvita-output.R
- 20.Plot_Roary_outputs.R

