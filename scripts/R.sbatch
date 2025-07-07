#!/bin/bash
#SBATCH --job-name=run_R_script      # Job name
#SBATCH --output=R_script_%j.out     # Output file (%j will be replaced with job ID)
#SBATCH --error=R_script_%j.err      # Error file (%j will be replaced with job ID)
#SBATCH --time=01:00:00              # Time limit (hh:mm:ss)
#SBATCH --partition=standard         # Partition name
#SBATCH --ntasks=1                   # Number of tasks
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=4G                     # Memory per node

# Load the R module
module load build-env/f2022
module load r-bundle-bioconductor/3.18-foss-2023a-r-4.3.2

# Run the R script
Rscript scripts/R.R