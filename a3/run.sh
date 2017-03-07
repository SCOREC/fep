#!/bin/bash -e
#SBATCH --job-name=fepsA3
#SBATCH --partition=debug
#SBATCH -t 00:10:00
#SBATCH -n 2

exe=$1; shift
args=("$@")
srun $exe ${args[@]}
