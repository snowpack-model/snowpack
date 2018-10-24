#############################
## Global Configuration    ##
#############################

# set Memory check command
set(CTEST_MEMORYCHECK_COMMAND valgrind)

#get the proper Snowpack library
IF(BUILD_SHARED_LIBS)
	SET(LIBALPINE3D_LIBRARY ${PROJECT_NAME})
ELSE(BUILD_SHARED_LIBS)
	IF(BUILD_STATIC_LIBS)
		SET(LIBALPINE3D_LIBRARY "${PROJECT_NAME}_STATIC")
	ELSE(BUILD_STATIC_LIBS)
		MESSAGE(SEND_ERROR "Not building Snowpack, the test won't be able to build")
	ENDIF(BUILD_STATIC_LIBS)
ENDIF(BUILD_SHARED_LIBS)
INCLUDE_DIRECTORIES(../)

FIND_PACKAGE(MeteoIO REQUIRED)
INCLUDE_DIRECTORIES(${METEOIO_INCLUDE_DIR})
FIND_PACKAGE(Libsnowpack REQUIRED)
INCLUDE_DIRECTORIES(${LIBSNOWPACK_INCLUDE_DIR})

#in order to link a new executable: 
#TARGET_LINK_LIBRARIES(${BINARY} ${LIBALPINE3D_LIBRARY} ${LIBSNOWPACK_LIBRARY} ${METEOIO_LIBRARY} ${EXTRA_LINKS})

###################
## Add Tests     ##
###################
ADD_SUBDIRECTORY(simple)
ADD_SUBDIRECTORY(basics)
