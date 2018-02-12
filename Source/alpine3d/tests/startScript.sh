#!/bin/bash
#warning: the order IS important (they depend on each other)

log_echo() {
        msg=$1
        me=${USER}
        datum=$(date "+%Y-%m-%dT%H:%M:%S")
        printf "[${datum}] [${me}] ${msg}\n"
}

#get expanded path to ALPINE3D's root, ie one level above this script
tmp=`dirname $0`
tmp="${tmp}/../"
cd ${tmp}
A3D_ROOT=`pwd`
cd -

#run the tests
log_echo "Starting ALPINE3D testing with A3D_ROOT=${A3D_ROOT}"
cd ${A3D_ROOT}

make distclean
/usr/bin/ctest -S ${A3D_ROOT}/tests/startScriptCoverage.cmake -V > ${A3D_ROOT}/tests/startScriptCoverage.log 2>&1

make distclean
/usr/bin/ctest -S ${A3D_ROOT}/tests/startScriptValgrind.cmake -V > ${A3D_ROOT}/tests/startScriptValgrind.log 2>&1

make distclean
/usr/bin/ctest -S ${A3D_ROOT}/tests/startScriptNightly.cmake -V > ${A3D_ROOT}/tests/startScriptNightly.log 2>&1

#so the up-to-date documentation can be available on the server
make doc

log_echo "ALPINE3D testing done!"
