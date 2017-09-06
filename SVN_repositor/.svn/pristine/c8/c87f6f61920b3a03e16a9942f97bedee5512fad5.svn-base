# simple Script that makes Experimental  build with cTest and publie it on CDash with coverage informations.
# This script need to be started be cron. The command for the cron is :
# "ctest -S startScriptCoverage.cmake"

# set for the thest source and binary directories
SET(CTEST_SOURCE_DIRECTORY .)
SET(CTEST_BINARY_DIRECTORY .)

# set cTest commands to be used
SET(CTEST_COMMAND "\"${CTEST_EXECUTABLE_NAME}\" -D Nightly")

#set cMake command to be used
SET(CTEST_CMAKE_COMMAND "\"${CMAKE_EXECUTABLE_NAME}\"")

#also possible to set initial cache values for config to set that the test are build
# BUILD NAME SET HERE DROUG CACHE.. OLD WAY TO OD BUT ONLY WORKING WAY
SET(CTEST_INITIAL_CACHE "
	BUILD_TESTING:BOOL=ON
	BUILD_TESTING_WITH_COVERAGE:BOOL=ON
	BUILDNAME:STRING=Linux_Coverage_noOptim
")
