#!/bin/bash
#SBATCH --job-name=coheritability
#SBATCH --time=180:00:00
#SBATCH --mail-type=end,fail
#SBATCH --mem=4g
#SBATCH --cpus-per-task=1
#SBATCH --array=1-8410

R CMD BATCH --no-save --no-restore 5_coheritability_alldat.R heri.Rout