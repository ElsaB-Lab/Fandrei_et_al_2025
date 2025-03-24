#!/bin/bash

#SBATCH -p mediumq
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -o sc_dna-%u-%N-%j.out
#SBATCH -e sc_dna-%u-%N-%j.err

module purge
module load anaconda3/2020-11
source /mnt/beegfs/software/anaconda3/2020-11/etc/profile.d/conda.sh
conda activate snakemake
module load singularity

python --version
snakemake --version
singularity --version
 
path_to_configfile="/mnt/beegfs/scratch/d_fandrei/config_gr.yaml"
path_to_pipeline="/mnt/beegfs/pipelines/single-cell_dna/1.1"

snakemake --profile ${path_to_pipeline}/profiles/slurm -s ${path_to_pipeline}/Snakefile --configfile ${path_to_configfile}

conda deactivate
