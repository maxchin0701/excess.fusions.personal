#!/bin/sh
#SBATCH -J simmap_test
#SBATCH --time=4-00:00:00
#SBATCH --partition=xlong
#SBATCH -n1
#SBATCH --mem=32000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maxchin@tamu.edu
#SBATCH -o simmap_test-%j.out
#SBATCH -e simmap_test-%j.err

#load R modules
module load GCC/10.2.0
module load OpenMPI/4.0.5
module load R_tamu/4.0.3

#Set temporary output library
export R_LIBS=$SCRATCH/R_LIBS_USER/OpenMPI-4.0.5-GCC-10.2.0

#Set LC
unset LC_ALL
export LC_CTYPE="en_US.UTF-8"
export LC_COLLATE="en_US.UTF-8"
export LC_TIME="en_US.UTF-8"
export LC_MESSAGES="en_US.UTF-8"
export LC_MONETARY="en_US.UTF-8"
export LC_PAPER="en_US.UTF-8"
export LC_MEASUREMENT="en_US.UTF-8"

#Run script
Rscript simmap_test.R




