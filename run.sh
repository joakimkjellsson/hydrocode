#!/bin/bash

#SBATCH -N1-1
#SBATCH -n 8
#SBATCH --mem=8000
#SBATCH -t 8:00:00
#SBATCH -J T2
#SBATCH -o hydro.out
#SBATCH -e hydro.out

ulimit -s unlimited
export OMP_NUM_THREADS=8
export OMP_STACKSIZE=1g

##
## Run it
##
make
time ./psi.x

