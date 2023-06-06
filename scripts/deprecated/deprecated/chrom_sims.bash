#!/bin/sh
#SBATCH -J chrom_sims
#SBATCH --time=2-00:00:00
#SBATCH --partition=xlong
#SBATCH -n1
#SBATCH --mem=32000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maxchin@tamu.edu
#SBATCH -o chrom_sims-%j.out
#SBATCH -e chrom_sims-%j.err

#load R modules
module load GCC/10.2.0
module load OpenMPI/4.0.5
module load R_tamu/4.0.3

Rscript simChrom_tests.R

