#!/bin/bash
# tells the system that this is a bash file
#SBATCH --partition=courses    # choose from debug, express, or short
#SBATCH --job-name=mono        # change this name to be informative for what you are running (eg. name of key script)
#SBATCH --time=08:00:00       # max time to run in hh:mm:ss, must be at or below max for partition
#SBATCH -N 6                   # nodes requested
#SBATCH -n 6
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G
#SBATCH --output=logs/%j.output     # where to direct standard output; will be jobnumber.output
#SBATCH --error=logs/%j.error       # where to direct standard error; will be jobnumber.error

# Automated setup script for Seurat analysis environment

# Step 1: Create conda environment from YAML specification
conda env create --name monocle3 -f ./config/monocle3_env.yml

# Step 2: Activate the seurat conda environment
eval "$(conda shell.bash hook)"
conda activate monocle3

# Step 3: Restore R packages from renv lockfile
#Rscript -e "renv::restore(lockfile = 'config/renv.lock', prompt = FALSE)"

Rscript /home/forbes.ai/capstone_dir/Monocle3/mono_full.R
