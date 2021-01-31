#!/bin/bash -e
#SBATCH --job-name=fepsA2
#SBATCH --partition=debug
#SBATCH -t 00:10:00
#SBATCH -n 1

gcc9libs=/gpfs/u/software/erp-rhel7/modulefiles/.gcc/9.1.0/1/lib64
export LD_LIBRARY_PATH=$gcc9libs:$LD_LIBRARY_PATH
module use /gpfs/u/software/erp-spack-install/vDev7a99c49_gcc910_1/lmod/linux-centos7-x86_64/gcc/9.1.0/
module load openmpi/4.0.1
module load pumi/master-int32-sim

exe=$1; shift
args=("$@")
srun --mpi=pmi2 --cpu_bind=cores $exe ${args[@]}
