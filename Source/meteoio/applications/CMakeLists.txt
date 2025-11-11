#SPDX-License-Identifier: LGPL-3.0-or-later
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")
INCLUDE("${CMAKE_MODULE_PATH}/BuildVersion.cmake")
BuildVersionGIT()
SET(MIO_BINARY "meteoio_timeseries")

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
IF (LEAKS_CHECK)
	SET(EXTRA_LINKS "${EXTRA_LINKS} -fsanitize=address -fno-omit-frame-pointer -static-libasan")
ENDIF(LEAKS_CHECK)

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

#Plateform specific options
IF(APPLE)
	SET(TARGET_ARGS ${APP_ICON_MACOSX})
	#Handle extra libraries needed for GUI exceptions
	SET(EXTRA_LINKS "-framework CoreServices")
ENDIF(APPLE)
IF(WIN32)
	#Handle the missing getopt on Windows
	IF(MSVC)
		SET(getopt_src getopt.c)
	ENDIF(MSVC)
	SET(TARGET_ARGS ${getopt_src} ${APP_ICON_RESOURCE_WINDOWS})
ENDIF(WIN32)

IF(UNIX AND NOT HAIKU)
	SET(EXTRA_LINKS "dl;pthread")
ENDIF(UNIX AND NOT HAIKU)

SET(meteoio-app_sources
	meteoio_timeseries.cc
	${getopt_src}
)

#Prepare executable
ADD_EXECUTABLE(${MIO_BINARY} ${meteoio-app_sources} ${TARGET_ARGS})
TARGET_LINK_LIBRARIES(${MIO_BINARY} ${METEOIO_LIBRARIES} ${EXTRA_LINKS})
IF(BUILD_STATIC_LIBS)
	TARGET_LINK_OPTIONS(${MIO_BINARY} PRIVATE -static-libgcc -static-libstdc++ -static)
ENDIF(BUILD_STATIC_LIBS)

SET_TARGET_PROPERTIES(${MIO_BINARY} PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
	CLEAN_DIRECT_OUTPUT 1
	OUTPUT_NAME "${MIO_BINARY}")

INSTALL(TARGETS ${MIO_BINARY} RUNTIME DESTINATION ${EXE_DEST} ${INSTALL_WIN32} COMPONENT exe)
