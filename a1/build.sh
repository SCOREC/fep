#!/bin/bash -e
[ $# -lt 1 ] && echo "Usage $0 <name of binary> [verbose=on|off]" && exit 1
module load xl proprietary/simmetrix/simModSuite proprietary/core-sim/xl
[ ! -e build ] && echo "ERROR: build directory does not exist. Run setup.sh" && exit 1
set -x
cd build
verbose=""
[ "$2" = "on" ] && verbose="VERBOSE=1"
make $verbose $1
cd -
set +x
