#!/bin/sh
#SBATCH -J subtree_simmap
#SBATCH --time=10-00:00:00
#SBATCH --partition=xlong
#SBATCH --ntasks=6
#SBATCH --ntasks-per-node=6
#SBATCH --mem=75G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maxchin@tamu.edu
#SBATCH -o subtree_simmap-%j.out
#SBATCH -e subtree_simmap-%j.err

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
Rscript subtree_simmap.R
