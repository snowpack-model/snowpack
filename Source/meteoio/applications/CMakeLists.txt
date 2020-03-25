INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")
INCLUDE("${CMAKE_MODULE_PATH}/BuildVersion.cmake")
BuildVersion()

# IF(CMAKE_COMPILER_IS_GNUCXX)
# 	SET(PROFILING "-pg -fprofile-arcs") #add ${PROFILING} to the CFLAGS when necessary
# 	EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
# 	IF(GCC_VERSION VERSION_GREATER 4.8 OR GCC_VERSION VERSION_EQUAL 4.8)
# 		# if set to ON, all binaries depending on the library have to be compiled the same way. 
# 		#Then, do an "export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.4" and run with "ASAN_OPTIONS=symbolize=1"
# 		SET(LEAKS_CHECK OFF CACHE BOOL "Set to ON to dynamically check for memory corruption (and do the same for applications linked with MeteoIO)")
# 		IF (LEAKS_CHECK)
# 			SET(EXTRA "-fsanitize=address -fno-omit-frame-pointer")
# 		ENDIF(LEAKS_CHECK)
# 	ENDIF()
# 	IF(GCC_VERSION VERSION_GREATER 5.0 OR GCC_VERSION VERSION_EQUAL 5.0)
# 	OPTION(USE_OCCI "Set to ON when using Oracle's OCCI" OFF)
# 		IF (USE_OCCI) #HACK: current OCCI does not support the short strings optimizations of gcc>=5
# 			SET(EXTRA "-D_GLIBCXX_USE_CXX11_ABI=0 ${EXTRA}")
# 		ENDIF(USE_OCCI)
# 	ENDIF()
# 	SET(EXTRA "-g ${EXTRA}") #add debug symbols
# ENDIF(CMAKE_COMPILER_IS_GNUCXX)
# SET(CMAKE_CXX_FLAGS "${EXTRA}" CACHE STRING "" FORCE)

#Handle the missing getopt on Windows
IF(MSVC)
	SET(getopt_src getopt.c getopt_long.c)
ENDIF(MSVC)

#get the proper MeteoIO library
SET(METEOIO_LIBRARIES ${PROJECT_NAME})
#link_directories(${PROJECT_SOURCE_DIR}/lib)


#Handle extra libraries needed for GUI exceptions
IF(APPLE)
	SET(EXTRA_LINKS "-framework CoreServices")
ENDIF(APPLE)
IF(UNIX AND NOT HAIKU)
	SET(EXTRA_LINKS "dl;pthread")
ENDIF(UNIX AND NOT HAIKU)

#Prepare executable
ADD_EXECUTABLE(meteoio_timeseries meteoio_timeseries.cc ${getopt_src})
TARGET_LINK_LIBRARIES(meteoio_timeseries ${METEOIO_LIBRARIES})
SET_TARGET_PROPERTIES(meteoio_timeseries PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
	CLEAN_DIRECT_OUTPUT 1
	OUTPUT_NAME "meteoio_timeseries")
INSTALL(TARGETS meteoio_timeseries RUNTIME DESTINATION bin COMPONENT exe)
