#!/bin/bash -e
[ $# -lt 1 ] && echo "Usage $0 <name of binary> [verbose=on|off]" && exit 1
module use /gpfs/u/software/erp-spack-install/v0190_0/lmod/linux-centos7-x86_64/gcc/9.1.0
module unuse /gpfs/u/software/erp-spack-install/lmod/linux-centos7-x86_64/Core
module load gcc/9.1.0/1
module load openmpi/3.1.6
module load cmake/3.24.3
module load pumi/master-int32-sim


[ ! -e build ] && echo "ERROR: build directory does not exist. Run setup.sh" && exit 1
set -x
cd build
verbose=""
[ "$2" = "on" ] && verbose="VERBOSE=1"
make $verbose $1
cd -
set +x
