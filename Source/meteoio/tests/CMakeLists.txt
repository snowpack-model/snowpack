#############################
## Global Configuration    ##
#############################

# set Memory check command
set(CTEST_MEMORYCHECK_COMMAND valgrind)

#get all include flags needed to make executables
SET(SHAREDNAME ${PROJECT_NAME})
SET(LIBRARIES ${SHAREDNAME} ${CMAKE_DL_LIBS})

# go back to the source dir to have all .h files as would be installed
INCLUDE_DIRECTORIES(..)

###################
## Add Tests     ##
###################
ADD_SUBDIRECTORY(meteo_reading_interpol)
ADD_SUBDIRECTORY(meteo_reading_no_interpol)
ADD_SUBDIRECTORY(sun)
ADD_SUBDIRECTORY(dem_reading)
ADD_SUBDIRECTORY(2D_interpolations)
ADD_SUBDIRECTORY(arrays)
ADD_SUBDIRECTORY(coords)
ADD_SUBDIRECTORY(stats)
