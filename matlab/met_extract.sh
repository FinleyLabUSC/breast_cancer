#!/bin/bash
#SBATCH --account=sfinley_500
#SBATCH --job-name=matlab_met_extraction
#SBATCH --partition=largemem
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --output=./eao/output-%A.out
#SBATCH --error=./eao/error-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="kingtyle@usc.edu"

module purge
module load matlab/2024a

matlab -nodisplay -batch 'extract_clusters_complete'