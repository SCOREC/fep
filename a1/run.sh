#!/bin/bash -e
#SBATCH --job-name=fepsA1
#SBATCH --partition=debug
#SBATCH -t 00:10:00
#SBATCH -n 1

exe=$1; shift
args=("$@")
srun $exe ${args[@]}
