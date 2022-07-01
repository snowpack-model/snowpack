#SPDX-License-Identifier: LGPL-3.0-or-later
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")
INCLUDE("${CMAKE_MODULE_PATH}/BuildVersion.cmake")
BuildVersionGIT()
SET(BINARY "meteoio_timeseries")

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
	SET(getopt_src getopt.c)
ENDIF(MSVC)

#get the proper MeteoIO library
IF(BUILD_SHARED_LIBS)
	SET(METEOIO_LIBRARIES ${PROJECT_NAME})
ELSE(BUILD_SHARED_LIBS)
	IF(BUILD_STATIC_LIBS)
		SET(METEOIO_LIBRARIES "${PROJECT_NAME}_STATIC")
	ELSE(BUILD_STATIC_LIBS)
		MESSAGE(SEND_ERROR "Not building MeteoIO, the standalone application won't be able to build")
	ENDIF(BUILD_STATIC_LIBS)
ENDIF(BUILD_SHARED_LIBS)
INCLUDE_DIRECTORIES(../)

#Handle extra libraries needed for GUI exceptions
IF(APPLE)
	SET(EXTRA_LINKS "-framework CoreServices")
ENDIF(APPLE)
IF(UNIX AND NOT HAIKU)
	SET(EXTRA_LINKS "dl;pthread")
ENDIF(UNIX AND NOT HAIKU)
IF (LEAKS_CHECK)
	SET(EXTRA_LINKS "${EXTRA_LINKS} -fsanitize=address -fno-omit-frame-pointer -static-libasan")
ENDIF(LEAKS_CHECK)

#Prepare executable
ADD_EXECUTABLE(${BINARY} meteoio_timeseries.cc ${getopt_src})
TARGET_LINK_LIBRARIES(${BINARY} ${METEOIO_LIBRARIES} ${EXTRA_LINKS})
IF(BUILD_STATIC_LIBS)
	TARGET_LINK_OPTIONS(${BINARY} PRIVATE -static-libgcc -static-libstdc++ -static)
ENDIF(BUILD_STATIC_LIBS)
SET_TARGET_PROPERTIES(${BINARY} PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
	CLEAN_DIRECT_OUTPUT 1
	OUTPUT_NAME "${BINARY}")
INSTALL(TARGETS ${BINARY} RUNTIME DESTINATION bin COMPONENT exe)
