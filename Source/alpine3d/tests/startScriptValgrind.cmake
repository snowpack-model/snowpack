# simple Script that makes Nigthly build with cTest and publie it on CDash with valgrind informations.
# This script need to be started be cron. The command for the cron is [see next section]:
# "ctest -S startScriptValgrind.cmake"

# Set Information for build
SET(CTEST_SITE "srv7")
SET(CTEST_BUILD_NAME "Linux_Standart_Valgrind")
SET(CTEST_BUILD_CONFIGURATION "release")

# set for the thest source and binary directories
SET(CTEST_SOURCE_DIRECTORY .)
SET(CTEST_BINARY_DIRECTORY .)

# set SVN command
SET(CTEST_SVN_COMMAND /usr/bin/svn)

# set Generator
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

#Configure valgrind
SET(CTEST_MEMORYCHECK_COMMAND "/usr/bin/valgrind")
#SET(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all --log-file=suppressions.log")
SET(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-reachable=yes --error-limit=no")
SET(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "${CMAKE_A3D_ROOT}/tests/bash_suppressions.txt")

#set cMake command to be used
SET(CTEST_CMAKE_COMMAND "\"${CMAKE_EXECUTABLE_NAME}\"")

# run tests
ctest_start(Nightly)
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}")
ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" OPTIONS "-DBUILD_TESTING=ON BUILD_TESTING_WITH_COVERAGE=OFF")
ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" )
ctest_memcheck(BUILD "${CTEST_BINARY_DIRECTORY}" INCLUDE "simple")
ctest_submit()
