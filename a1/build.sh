#!/bin/bash -e
[ $# -lt 1 ] && echo "Usage $0 <name of binary> [verbose=on|off]" && exit 1
module use /gpfs/u/software/erp-spack-install/vDev7a99c49_gcc910_1/lmod/linux-centos7-x86_64/gcc/9.1.0/
module load gcc/9.1.0/1
module load openmpi/4.0.1
module load cmake/3.19.2
module load pumi/master-int32-sim

[ ! -e build ] && echo "ERROR: build directory does not exist. Run setup.sh" && exit 1
set -x
cd build
verbose=""
[ "$2" = "on" ] && verbose="VERBOSE=1"
make $verbose $1
cd -
set +x
