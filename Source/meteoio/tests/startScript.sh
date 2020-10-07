#!/bin/bash
#warning: the order IS important (they depend on each other)

log_echo() {
        msg=$1
        me=${USER}
        datum=$(date "+%Y-%m-%dT%H:%M:%S")
        printf "[${datum}] [${me}] ${msg}\n"
}

#get expanded path to MeteoIO's root, ie one level above this script
tmp=`dirname $0`
tmp="${tmp}/../"
cd ${tmp}
MIO_ROOT=`pwd`

#run the tests
log_echo "Starting MeteoIO testing"

make distclean
/usr/bin/ctest -S ${MIO_ROOT}/tests/startScriptCoverage.cmake -V > ${MIO_ROOT}/tests/startScriptCoverage.log 2>&1

make distclean
/usr/bin/ctest -S ${MIO_ROOT}/tests/startScriptValgrind.cmake -V > ${MIO_ROOT}/tests/startScriptValgrind.log 2>&1

make distclean
/usr/bin/ctest -S ${MIO_ROOT}/tests/startScriptNightly.cmake -V > ${MIO_ROOT}/tests/startScriptNightly.log 2>&1

#so the up-to-date documentation can be available on the server
make doc

log_echo "MeteoIO testing done!"
