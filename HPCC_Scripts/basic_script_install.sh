#!/bin/bash
## Example of R job
#SBATCH --job-name=test_R
#SBATCH --partition=haswell
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=katelynn.lankowicz@umaine.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1gb
#SBATCH --time=4:00:00
#SBATCH --output=test_%j.log

## Load application environment
module load R/4.2.2

## Navigate to folder
cd VAST_install

## Run application commands
Rscript --vanilla installing_packages.R
