#!/bin/bash -e
module use /gpfs/u/software/erp-spack-install/v0190_0/lmod/linux-centos7-x86_64/gcc/9.1.0
module unuse /gpfs/u/software/erp-spack-install/lmod/linux-centos7-x86_64/Core
module load gcc/9.1.0/1
module load openmpi/3.1.6
module load cmake/3.24.3
module load pumi/master-int32-sim

set -x
[ ! -e build ] && mkdir build
cd build
flags="-g -O2"
cmake -DCMAKE_CXX_COMPILER=`which mpicxx` -DCMAKE_CXX_FLAGS="$flags" ../
cd -
set +x
