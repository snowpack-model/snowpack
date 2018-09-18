#############################
## Global Configuration    ##
#############################

# set Memory check command
set(CTEST_MEMORYCHECK_COMMAND valgrind)

#get the proper MeteoIO library
IF(BUILD_SHARED_LIBS)
	SET(METEOIO_LIBRARIES ${PROJECT_NAME})
ELSE(BUILD_SHARED_LIBS)
	IF(BUILD_STATIC_LIBS)
		SET(METEOIO_LIBRARIES "${PROJECT_NAME}_STATIC")
	ELSE(BUILD_STATIC_LIBS)
		MESSAGE(SEND_ERROR "Not building MeteoIO, the test won't be able to build")
	ENDIF(BUILD_STATIC_LIBS)
ENDIF(BUILD_SHARED_LIBS)

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
