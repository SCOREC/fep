#!/bin/bash -e
#SBATCH --job-name=fepsA4
#SBATCH --partition=debug
#SBATCH -t 00:10:00
#SBATCH -n 1 

source erp_env_setup.sh

exe=$1; shift
args=("$@")
srun --mpi=pmi2 --cpu_bind=cores $exe ${args[@]}
