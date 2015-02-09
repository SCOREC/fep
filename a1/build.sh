#!/bin/bash -e
x=$(module load xl 2>&1)
if [ ! -z "$x" ]; then
 echo $x
 echo "Module conflict.... exiting."
 exit 1
fi
module load xl
pumi=/gpfs/u/barn/FEP2/shared/pumi/dbg/lib/pkgconfig
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$pumi
set -x
CXXFLAGS="$CXXFLAGS -g "
mpicxx $1.cc $CXXFLAGS `pkg-config --cflags --libs libmds` -o $1
set +x
