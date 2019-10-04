CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
#in case the user has copied the examples directory, we fallback on to the local copy of FindMeteoIO
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../../tools/cmake/" "${CMAKE_SOURCE_DIR}/cmake")

PROJECT(meteoio_examples)

FIND_PACKAGE(MeteoIO)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${METEOIO_INCLUDE_DIR})

IF(CMAKE_COMPILER_IS_GNUCXX)
	SET(PROFILING "-pg -fprofile-arcs") #add ${PROFILING} to the CFLAGS when necessary
	EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
	IF(GCC_VERSION VERSION_GREATER 4.8 OR GCC_VERSION VERSION_EQUAL 4.8)
		# if set to ON, all binaries depending on the library have to be compiled the same way. 
		#Then, do an "export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.4" and run with "ASAN_OPTIONS=symbolize=1"
		SET(LEAKS_CHECK OFF CACHE BOOL "Set to ON to dynamically check for memory corruption (and do the same for applications linked with MeteoIO)")
		IF (LEAKS_CHECK)
			SET(EXTRA "-fsanitize=address -fno-omit-frame-pointer")
		ENDIF(LEAKS_CHECK)
	ENDIF()
	IF(GCC_VERSION VERSION_GREATER 5.0 OR GCC_VERSION VERSION_EQUAL 5.0)
	OPTION(USE_OCCI "Set to ON when using Oracle's OCCI" OFF)
		IF (USE_OCCI) #HACK: current OCCI does not support the short strings optimizations of gcc>=5
			SET(EXTRA "-D_GLIBCXX_USE_CXX11_ABI=0 ${EXTRA}")
		ENDIF(USE_OCCI)
	ENDIF()
	SET(EXTRA "-g ${EXTRA}") #add debug symbols
ENDIF(CMAKE_COMPILER_IS_GNUCXX)
SET(CMAKE_CXX_FLAGS "${EXTRA}" CACHE STRING "" FORCE)

#####################################################################
# add all the executables
ADD_EXECUTABLE(2D_interpolations 2D_interpolations.cc)
TARGET_LINK_LIBRARIES(2D_interpolations ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(coordinates coordinates.cc)
TARGET_LINK_LIBRARIES(coordinates ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(data_converter data_converter.cc)
TARGET_LINK_LIBRARIES(data_converter ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(dem_reading dem_reading.cc)
TARGET_LINK_LIBRARIES(dem_reading ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(grid2d_reading grid2d_reading.cc)
TARGET_LINK_LIBRARIES(grid2d_reading ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(matrix matrix.cc)
TARGET_LINK_LIBRARIES(matrix ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(meteo_reading meteo_reading.cc)
TARGET_LINK_LIBRARIES(meteo_reading ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(sun sun.cc)
TARGET_LINK_LIBRARIES(sun ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(time time.cc)
TARGET_LINK_LIBRARIES(time ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(random_numbers random_numbers.cc)
TARGET_LINK_LIBRARIES(random_numbers ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(statistical_filters statistical_filters.cc)
TARGET_LINK_LIBRARIES(statistical_filters ${METEOIO_LIBRARIES})

#####################################################################
# add a distclean target
ADD_CUSTOM_TARGET(distclean make clean
	COMMAND cmake -E remove *~ *.o
	COMMAND cmake -E remove CMakeCache.txt cmake_install.cmake
	COMMAND cmake -E remove_directory CMakeFiles
	)
