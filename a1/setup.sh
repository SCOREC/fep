#!/bin/bash -e
module load cmake xl proprietary/simmetrix/simModSuite proprietary/core-sim/xl
set -x
[ ! -e build ] && mkdir build
cd build
flags='-g -O2'
cmake -DCMAKE_CXX_COMPILER=`which mpicxx` -DCMAKE_CXX_FLAGS=$flags ../
cd -
set +x
