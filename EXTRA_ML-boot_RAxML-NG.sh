## Reinferring ML tree from reduced alignment and computing bootstrap - improve the accuracy and phylogenetic representation of the tree - RAxML-NG v. 1.2.2

### Preparing the alignment and checking it  
raxml-ng --parse --msa QC_filtered_genome_files_concatenated.aln.reduced --model GTR+G --prefix correct ## Info- Model: GTR+FO+G4m; Alignment sites / patterns: 119718 / 53739; Gaps: 1.15 %; Invariant sites: 0.02 %; * Estimated memory requirements: 19943 MB; * Recommended number of threads / MPI processes: 14

### Perform all-in-one (ML search + bootstrapping) 
### --------------- The following lines as a whole should be submitted as a job in a HPC ---------------
#! /bin/bash 
#SBATCH --job-name=raxml-ng-mpi
#SBATCH --partition=SP2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=192:00:00

echo $SLURM_JOB_ID 
echo $SLURM_SUBMIT_DIR
echo $SLURM_JOB_NODELIST
echo $SLURM_NTASKS  

#export OMP_NUM_THREADS=20

module load gnu8
module load openmpi3
module load raxml-ng-mpi

# Before running the script, unload the gnu module with 'module unload gnu/5.4.0'. Then, you can simply run the script. 

# Execute RAxML-NG with MPI parallelization (starting trees: 20 random and 20 parsimony/ model 
mpirun -np $SLURM_NTASKS raxml-ng-mpi --all --tree rand{20},pars{20} --msa correct.raxml.rba --model GTR+FO+G4m --prefix correct --seed 88585 --threads 10 --workers 2

### ------------------------------------------------------
# OBS: by default in --all parameter: up to 1000 bootstrap replicates will be generated, with the sufficient number of replicates determined automatically using the autoMRE bootstrap convergence test (so-called bootstopping, Pattengale et al. 2009)
