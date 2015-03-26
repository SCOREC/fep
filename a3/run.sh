#!/bin/bash -e
#SBATCH --job-name=fepsA3
#SBATCH --partition=debug
#SBATCH -t 00:10:00
#SBATCH -n 2

MY_RUNJOB_OPTS="PAMI_MAX_COMMTHREADS=0 PAMID_ASYNC_PROGRESS=0 PAMID_CONTEX_POST=1 BG_THREADMODEL=1"
exe=$1; shift
args=("$@")
srun --runjob-opts="--envs $MY_RUNJOB_OPTS" $exe ${args[@]}
