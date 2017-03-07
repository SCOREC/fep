#!/bin/bash -e
[ $# -gt 1 ] && echo "Usage $0 [verbose=on|off]" && exit 1
module load xl proprietary/simmetrix/simModSuite proprietary/core-sim/xl
[ ! -e build ] && echo "ERROR: build directory does not exist. Run setup.sh" && exit 1
set -x
cd build
verbose=""
[ "$1" = "on" ] && verbose="VERBOSE=1"
make $verbose
cd -
set +x
