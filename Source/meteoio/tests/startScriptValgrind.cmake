# simple Script that makes Nigthly build with cTest and publie it on CDash with valgrind informations.
# This script need to be started be cron. The command for the cron is [see next section]:
# "ctest -S startScriptValgrind.cmake"

# Set Information for build
SET(CTEST_SITE "srv7")
SET(CTEST_BUILD_NAME "Linux_Standart_Valgrind")
SET(CTEST_BUILD_CONFIGURATION "debug")

# set for the thest source and binary directories
SET(CTEST_SOURCE_DIRECTORY .)
SET(CTEST_BINARY_DIRECTORY .)

# set SVN command
SET(CTEST_SVN_COMMAND /usr/bin/svn)

# set Generator
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

#Configure valgrind
SET(CTEST_MEMORYCHECK_COMMAND "/usr/bin/valgrind")
SET(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--log-file=vlogs --leak-check=full")
#SET(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE) # to say, what valgrind should not "test" ;-), filter out the false positive errors of valgrind

#set cMake command to be used
SET(CTEST_CMAKE_COMMAND "\"${CMAKE_EXECUTABLE_NAME}\"")
SET(CTEST_BUILD_COMMAND "/usr/bin/make")

ctest_start(Nightly)
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}")
ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" OPTIONS "-DBUILD_TESTING=ON BUILD_TESTING_WITH_COVERAGE=OFF")
ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}")
ctest_memcheck(BUILD "${CTEST_BINARY_DIRECTORY}")
ctest_submit()
