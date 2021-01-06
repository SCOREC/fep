#!/bin/bash -e
module use /gpfs/u/software/erp-spack-install/v0160_gcc910_1/lmod/linux-centos7-x86_64/gcc/9.1.0/
module load gcc/9.1.0/1
module load openmpi/4.0.1
module load cmake/3.19.2
module load pumi/master-int32-sim

set -x
[ ! -e build ] && mkdir build
cd build
flags="-g -O2"
cmake -DCMAKE_CXX_COMPILER=`which mpicxx` -DCMAKE_CXX_FLAGS="$flags" ../
cd -
set +x
