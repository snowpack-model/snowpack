## REALY ??? somthing about old/new format of command error, not get it why....
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#############################
## Global Configuration    ##
#############################

# set Memory check command, more flags possible to be set
set(CTEST_MEMORYCHECK_COMMAND valgrind)

#get all include flags needed to make executables
SET(SHAREDNAME ${PROJECT_NAME})
SET(LIBRARIES ${SHAREDNAME} ${CMAKE_DL_LIBS})

# go back to the source dir to have all .h files, same tree/structure as would be installed
INCLUDE_DIRECTORIES(..)

###################
## Add Tests     ##
###################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To group test (for local run test or script) there are two possibelities
# 1. Use "ctest -R RecEx" where RecEx is a regular expretion matching the name
# 2. Group Test on Labels wiht Set Properties and run "ctest -L LabelName"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ADD_SUBDIRECTORY(meteo_reading_interpol)
ADD_SUBDIRECTORY(meteo_reading_no_interpol)
ADD_SUBDIRECTORY(sun)
ADD_SUBDIRECTORY(dem_reading)
ADD_SUBDIRECTORY(2D_interpolations)
ADD_SUBDIRECTORY(arrays)
ADD_SUBDIRECTORY(coords)
ADD_SUBDIRECTORY(stats)
