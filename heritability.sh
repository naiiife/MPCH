#!/bin/bash
#SBATCH --job-name=heritability
#SBATCH --time=240:00:00
#SBATCH --mail-type=end,fail
#SBATCH --mem=4g
#SBATCH --cpus-per-task=1
#SBATCH --array=1-290

R CMD BATCH --no-save --no-restore 2_heritability_alldat.R heri.Rout