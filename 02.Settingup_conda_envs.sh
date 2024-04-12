# Setting up my conda environments

## Initiating miniconda 
export PATH="/temporario/11217468/miniconda3/bin:$PATH"
conda init bash
conda create --name bioinfo
conda activate bioinfo

## Downloading required packages available through bioconda (all have python version 3.8)
conda install bioconda::blast
conda install bioconda::prokka
conda install bioconda::bowtie2  
conda install bioconda::diamond  
conda install bioconda::cd-hit
conda install bioconda::hmmer
conda install bioconda::quast
conda deactivate 

## Creating separated environments for specific packages I know will conflict with those aforementioned
conda create --name checkm
conda activate checkm
conda install -c bioconda checkm-genome
conda deactivate 

 



