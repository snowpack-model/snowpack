#SPDX-License-Identifier: LGPL-3.0-or-later
CMAKE_MINIMUM_REQUIRED(VERSION 3.1...3.27.2)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/tools/cmake/")

PROJECT(meteoio)
# Version of library
IF(DEFINED ENV{VERSION_MAJOR} AND DEFINED ENV{VERSION_MINOR} AND DEFINED ENV{VERSION_PATCH})
	SET(VERSION_MAJOR $ENV{VERSION_MAJOR} CACHE INTERNAL "Major version")
	SET(VERSION_MINOR $ENV{VERSION_MINOR} CACHE INTERNAL "Minor version")
	SET(VERSION_PATCH $ENV{VERSION_PATCH} CACHE INTERNAL "Patch version")
ELSE()
	SET(VERSION_MAJOR 2 CACHE INTERNAL "Major version")
	SET(VERSION_MINOR 11 CACHE INTERNAL "Minor version")
	SET(VERSION_PATCH 0 CACHE INTERNAL "Patch version")
ENDIF()


SET(LIBPREFIX "lib") #make sure ALL plugins are libXXX for ALL plateforms
SET(CMAKE_IMPORT_LIBRARY_PREFIX "lib") #for Visual C++ import libraries

###########################################################
#platform specific. First, compiler specific flags, then system specific flags
OPTION(BUILD_SHARED_LIBS "Set to ON to build dynamic libraries" ON)
IF(WIN32)
	SET(BUILD_STATIC_LIBS ON) #otherwise the lib files are not generated
ELSE(WIN32)
	OPTION(BUILD_STATIC_LIBS "Set to ON to build static libraries" OFF)
ENDIF(WIN32)

###########################################################
#set compiler flags according to the detected compiler and user options
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED True)
INCLUDE("${CMAKE_SOURCE_DIR}/tools/cmake/SetCompilerOptions.cmake")
SET_COMPILER_OPTIONS()

#The following pre-processor variables are automatically defined:
# __linux, __GNUC__, __CYGWIN__, _WIN32, __MINGW32__, __APPLE__
IF(UNIX)
	INCLUDE (CheckStructHasMember)
	CHECK_STRUCT_HAS_MEMBER("struct stat" st_flags sys/stat.h HAVE_STRUCT_STAT_ST_FLAGS)
	IF(HAVE_STRUCT_STAT_ST_FLAGS) 
		SET(EXTRA "-DHAVE_STRUCT_STAT_ST_FLAGS ${EXTRA}")
	ENDIF()
	
	OPTION(USE_RPATH "Use rpath (not recommended for Fedora)" ON)
	MARK_AS_ADVANCED(FORCE USE_RPATH)
	if(USE_RPATH)
		#set a full rpath for both build & install tree
		SET(CMAKE_SKIP_BUILD_RPATH FALSE)
		SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
		# add the directory of the binary to the rpath
		#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'$ORIGIN'")
		IF(APPLE)
			SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'@executable_path',-rpath,'@executable_path/../lib'")
			SET(CMAKE_OSX_DEPLOYMENT_TARGET "10.13" CACHE STRING "Minimum OS X deployment version")
		ELSE(APPLE)
			SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'$ORIGIN',-rpath,'$ORIGIN'/../lib")
		ENDIF(APPLE)

		# add the automatically determined parts of the RPATH
		SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
		#SET(CMAKE_INSTALL_RPATH $ORIGIN)
		# the RPATH to be used when installing, but only if it's not a system directory
		LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
		IF("${isSystemDir}" STREQUAL "-1")
			SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
		ENDIF("${isSystemDir}" STREQUAL "-1")
	ELSE(USE_RPATH)
		SET(CMAKE_SKIP_BUILD_RPATH TRUE)
	ENDIF(USE_RPATH)

	IF(APPLE)
		SET(DLL_EXT "dylib*")
		SET(STAT_EXT "a")
		SET(EXTRA_WARNINGS "${EXTRA_WARNINGS} -Wno-sign-conversion") #for now, keep clang quiet
		SET(CMAKE_SHARED_LINKER_FLAGS "-framework CoreFoundation")
		#set(CMAKE_OSX_SYSROOT iphoneos2.2.1) #for IOS
		#set(CMAKE_OSX_ARCHITECTURES $(ARCHS_STANDARD_32_BIT)) #for IOS
	ELSE(APPLE)
		SET(DLL_EXT "so*")
		SET(STAT_EXT "a")
	ENDIF(APPLE)
ELSE(UNIX)
	IF(WIN32)
		SET(DLL_EXT "dll")
		SET(STAT_EXT "a")
	ENDIF(WIN32)
ENDIF(UNIX)

IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE "release" CACHE STRING "Choose the type of build, options are: debug release coverage cxxflags" FORCE)
	SET(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE})
ENDIF(NOT CMAKE_BUILD_TYPE)

###########################################################
#plugins and applications
SET(METEOIO_TIMESERIES ON CACHE BOOL "Compilation METEOIO_TIMESERIES standalone ON or OFF")

SET(PLUGIN_A3DIO ON CACHE BOOL "Compilation A3DIO ON or OFF")
SET(PLUGIN_ALPUG ON CACHE BOOL "Compilation ALPUG ON or OFF")
SET(PLUGIN_ARCIO ON CACHE BOOL "Compilation ARCIO ON or OFF")
SET(PLUGIN_ARGOSIO OFF CACHE BOOL "Compilation ARGOSIO ON or OFF")
SET(PLUGIN_ARPSIO ON CACHE BOOL "Compilation ARPSIO ON or OFF")
SET(PLUGIN_COSMOXMLIO OFF CACHE BOOL "Compilation COSMOXMLIO ON or OFF")
SET(PLUGIN_CSVIO ON CACHE BOOL "Compilation CsvIO ON or OFF")
SET(PLUGIN_DBO OFF CACHE BOOL "Compilation DBO ON or OFF")
SET(PLUGIN_GEOTOPIO ON CACHE BOOL "Compilation GEOTOPIO ON or OFF")
SET(PLUGIN_GOESIO ON CACHE BOOL "Compilation GOESIO ON or OFF")
SET(PLUGIN_GRASSIO ON CACHE BOOL "Compilation GRASSIO ON or OFF")
SET(PLUGIN_GRIBIO OFF CACHE BOOL "Compilation GRIBIO ON or OFF")
SET(PLUGIN_iCSVIO ON CACHE BOOL "Compilation iCSVIO ON or OFF")
SET(PLUGIN_IMISIO OFF CACHE BOOL "Compilation IMISIO ON or OFF")
SET(PLUGIN_METEOBLUE OFF CACHE BOOL "Compilation METEOBLUE ON or OFF")
SET(PLUGIN_MYSQLIO OFF CACHE BOOL "Compilation MYSQLIO ON or OFF")
SET(PLUGIN_NETCDFIO OFF CACHE BOOL "Compilation NETCDFIO ON or OFF")
SET(PLUGIN_OSHDIO OFF CACHE BOOL "Compilation OshdIO ON or OFF")
SET(PLUGIN_PGMIO ON CACHE BOOL "Compilation PGMIO ON or OFF")
SET(PLUGIN_PMODIO OFF CACHE BOOL "Compilation PMODIO ON or OFF")
SET(PLUGIN_PNGIO OFF CACHE BOOL "Compilation PNGIO ON or OFF")
SET(PLUGIN_PSQLIO OFF CACHE BOOL "Compilation PSQLIO ON or OFF")
SET(PLUGIN_SASEIO OFF CACHE BOOL "Compilation SASEIO ON or OFF")
SET(PLUGIN_SMETIO ON CACHE BOOL "Compilation SMETIO ON or OFF")
SET(PLUGIN_SNIO ON CACHE BOOL "Compilation SNIO ON or OFF")
SET(PLUGIN_SYNTHIO ON CACHE BOOL "Compilation SYNTHIO ON or OFF")
SET(PLUGIN_ZRXPIO OFF CACHE BOOL "Compilation ZRXPIO ON or OFF")
SET(PROJ OFF CACHE BOOL "Use PROJ for the class MapProj ON or OFF")
###########################################################
## Additional Options
###########################################################
SET(LIMIT_OUT_ACCESS OFF CACHE BOOL "Limit write access of the library to CWD")
MARK_AS_ADVANCED(FORCE LIMIT_OUT_ACCESS)
# SET(LIMIT_BASE_DIR ".")
IF(LIMIT_OUT_ACCESS)
	IF(MSVC)	#TODO: replace by add_compile_definitions(LIMIT_WRITE_ACCESS) when moving to cmake >=3.12 released 11.2018
		ADD_DEFINITIONS(/LIMIT_WRITE_ACCESS)
	ELSE(MSVC)
		ADD_DEFINITIONS(-DLIMIT_WRITE_ACCESS)
	ENDIF(MSVC)
	SET(LIMIT_BASE_DIR "." CACHE STRING "Limit write access of the library to the given directory")
	IF(MSVC)	#TODO: replace by add_compile_definitions("LIMIT_BASE_DIR=\"${LIMIT_BASE_DIR}\"") when moving to cmake >=3.12 released 11.2018
		ADD_DEFINITIONS(/LIMIT_BASE_DIR=\"${LIMIT_BASE_DIR}\")
	ELSE(MSVC)
		ADD_DEFINITIONS(-DLIMIT_BASE_DIR=\"${LIMIT_BASE_DIR}\")
	ENDIF(MSVC)
ENDIF()

IF(CMAKE_CXX_COMPILER_ID MATCHES "^GNU$")
	IF(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0 OR CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 5.0)
		IF(PLUGIN_IMISIO) #HACK: current OCCI does not support the short strings optimizations of gcc>=5
			SET(EXTRA "-D_GLIBCXX_USE_CXX11_ABI=0 ${EXTRA}")
		ENDIF(PLUGIN_IMISIO)
	ENDIF()
ENDIF()

#show exception messages in a graphical message box
SET(GUI_EXCEPTIONS OFF CACHE BOOL "Show a message box with exceptions texts ON or OFF")

###########################################################
#finally, SET compile flags
SET(CMAKE_CXX_FLAGS "${_VERSION} ${ARCH} ${EXTRA}" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_RELEASE "${OPTIM}" CACHE STRING "" FORCE)

SET(CMAKE_CXX_FLAGS_DEBUG "${DEBUG} ${WARNINGS} ${EXTRA_WARNINGS}" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_CXXFLAGS "$ENV{CXXFLAGS}" CACHE STRING "" FORCE)
MARK_AS_ADVANCED(FORCE CMAKE_CXX_FLAGS_CXXFLAGS)

#Add subdirectory, so they can inherit all variables defined above
ADD_SUBDIRECTORY(meteoio/dataClasses)
ADD_SUBDIRECTORY(meteoio/dataGenerators)
ADD_SUBDIRECTORY(meteoio/plugins)
ADD_SUBDIRECTORY(meteoio/meteoLaws)
ADD_SUBDIRECTORY(meteoio/meteoFilters)
ADD_SUBDIRECTORY(meteoio/meteoResampling)
ADD_SUBDIRECTORY(meteoio/gridResampling)
ADD_SUBDIRECTORY(meteoio/spatialInterpolations)
ADD_SUBDIRECTORY(meteoio)

IF(METEOIO_TIMESERIES)
	#it seems that the only reliable way is to test-compile...
	INCLUDE (CheckCXXSourceCompiles)
	CHECK_CXX_SOURCE_COMPILES("
	#include <fenv.h>
	int main() {
	feenableexcept(FE_DIVBYZERO);
	return 0;
	}" HAVE_FEENABLE)

	IF(HAVE_FEENABLE)
		SET(DEBUG_ARITHM ON CACHE BOOL "Force-crash the application if doing an arithmetic exception")
	ENDIF(HAVE_FEENABLE)
	ADD_SUBDIRECTORY(applications)
ENDIF(METEOIO_TIMESERIES)

###########################################################
#for the install target
FILE(GLOB hfiles "meteoio/*.h")
FILE(GLOB hdatafiles "meteoio/dataClasses/*.h")
FILE(GLOB hgeneratorsfiles "meteoio/dataGenerators/*.h")
FILE(GLOB hresamplingfiles "meteoio/meteoResampling/*.h")
FILE(GLOB hgridresamplingfiles "meteoio/gridResampling/*.h")
FILE(GLOB hspatialfiles "meteoio/spatialInterpolations/*.h")
FILE(GLOB hpluginfiles "meteoio/plugins/*.h")
FILE(GLOB hlawsfiles "meteoio/meteoLaws/*.h")
FILE(GLOB hfilterfiles "meteoio/meteoFilters/*.h")
FILE(GLOB hstatsfiles "meteoio/meteoStats/*.h")

INSTALL(FILES ${hfiles} DESTINATION include/meteoio COMPONENT headers)
INSTALL(FILES ${hdatafiles} DESTINATION include/meteoio/dataClasses COMPONENT headers)
INSTALL(FILES ${hgeneratorsfiles} DESTINATION include/meteoio/dataGenerators COMPONENT headers)
INSTALL(FILES ${hresamplingfiles} DESTINATION include/meteoio/meteoResampling COMPONENT headers)
INSTALL(FILES ${hgridresamplingfiles} DESTINATION include/meteoio/gridResampling COMPONENT headers)
INSTALL(FILES ${hspatialfiles} DESTINATION include/meteoio/spatialInterpolations COMPONENT headers)
INSTALL(FILES ${hpluginfiles} DESTINATION include/meteoio/plugins COMPONENT headers)
INSTALL(FILES ${hlawsfiles} DESTINATION include/meteoio/meteoLaws COMPONENT headers)
INSTALL(FILES ${hfilterfiles} DESTINATION include/meteoio/meteoFilters COMPONENT headers)
INSTALL(FILES ${hstatsfiles} DESTINATION include/meteoio/meteoStats COMPONENT headers)

ADD_CUSTOM_TARGET(distclean make clean
	COMMAND cmake -E remove CMakeCache.txt CPackSourceConfig.cmake CPackConfig.cmake CTestTestfile.cmake DartConfiguration.tcl cmake_install.cmake cmake_uninstall.cmake install_manifest.txt Makefile
	COMMAND cmake -E remove "lib/libmeteo*.${DLL_EXT}" "lib/libmeteo*.${STAT_EXT}" "lib/libmeteo*.lib"
	COMMAND cmake -E remove */CMakeCache.txt */CTestTestfile.cmake */cmake_install.cmake
	COMMAND cmake -E remove */*/CMakeCache.txt */*/CTestTestfile.cmake */*/cmake_install.cmake
	COMMAND cmake -E remove_directory CMakeFiles
	COMMAND cmake -E remove_directory Testing
	COMMAND cmake -E remove_directory meteoio/CMakeFiles
	COMMAND cmake -E remove_directory meteoio/dataClasses/CMakeFiles
	COMMAND cmake -E remove_directory meteoio/dataGenerators/CMakeFiles
	COMMAND cmake -E remove_directory meteoio/meteoFilters/CMakeFiles
	COMMAND cmake -E remove_directory meteoio/meteoLaws/CMakeFiles
	COMMAND cmake -E remove_directory meteoio/meteoResampling/CMakeFiles
	COMMAND cmake -E remove_directory meteoio/gridResampling/CMakeFiles
	COMMAND cmake -E remove_directory meteoio/spatialInterpolations/CMakeFiles
	COMMAND cmake -E remove_directory meteoio/plugins/CMakeFiles
	COMMAND cmake -E remove_directory tests/CMakeFiles
	COMMAND cmake -E remove_directory tests/Testing
	COMMAND cmake -E remove_directory tests/meteo_reading_interpol/CMakeFiles
	COMMAND cmake -E remove_directory tests/meteo_reading_interpol/Testing
	COMMAND cmake -E remove_directory tests/meteo_reading_no_interpol/CMakeFiles
	COMMAND cmake -E remove_directory tests/meteo_reading_no_interpol/Testing
	COMMAND cmake -E remove_directory tests/sun/CMakeFiles
	COMMAND cmake -E remove_directory tests/sun/Testing
	COMMAND cmake -E remove tests/sun/generated_results.txt
	COMMAND cmake -E remove_directory tests/dem_reading/CMakeFiles
	COMMAND cmake -E remove_directory tests/dem_reading/Testing
	COMMAND cmake -E remove tests/dem_reading/AZI.asc tests/dem_reading/DEM.asc tests/dem_reading/SLOPE.asc
	COMMAND cmake -E remove_directory tests/2D_interpolations/CMakeFiles
	COMMAND cmake -E remove_directory tests/2D_interpolations/Testing
	COMMAND cmake -E remove tests/2D_interpolations/2009-01-19T12.00_HNW.asc tests/2D_interpolations/2009-01-19T12.00_RH.asc tests/2D_interpolations/2009-01-19T12.00_RSWR.asc tests/2D_interpolations/2009-01-19T12.00_TA.asc
)

###########################################################
#for the uninstall target
CONFIGURE_FILE(
	"${PROJECT_SOURCE_DIR}/tools/cmake/cmake_uninstall.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
	IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

#Build the documentation if doxygen is found
INCLUDE(FindDoxygen)
FIND_PACKAGE(Doxygen)
if (DOXYGEN_FOUND)
	ADD_CUSTOM_TARGET(doc COMMAND ${CMAKE_COMMAND} -E env MIO_VERSION="${_versionString}" "doxygen" ${PROJECT_SOURCE_DIR}/meteoio/config.dox)
ENDIF (DOXYGEN_FOUND)
#Install the documentation if the user desires it
SET(INSTALL_DOC OFF CACHE BOOL "Install documentation ON or OFF")
IF(INSTALL_DOC)
	IF(WIN32)
		SET(DOC_PATH "doc")
	ELSE(WIN32)
		SET(DOC_PATH "share/doc/${CMAKE_PROJECT_NAME}")
	ENDIF(WIN32)
	SET(DOC_FILES License.txt meteoio/icons/icon.ico doc/SMET_specifications.pdf)
	IF (DOXYGEN_FOUND)
		INSTALL(DIRECTORY doc/html DESTINATION ${DOC_PATH} COMPONENT documentation)
		LIST(APPEND DOC_FILES doc/MeteoIO_doc.html)
	ELSE (DOXYGEN_FOUND)
		MESSAGE("Doxygen has not been found and is needed to build online documentation. Please either install doxygen or download a documentation-only package!")
	ENDIF (DOXYGEN_FOUND)

	INSTALL(FILES ${DOC_FILES} DESTINATION ${DOC_PATH} COMPONENT documentation)
	INSTALL(DIRECTORY "doc/examples" USE_SOURCE_PERMISSIONS DESTINATION ${DOC_PATH} COMPONENT examples PATTERN ".svn" EXCLUDE)
ENDIF(INSTALL_DOC)

###########################################################
## CPACK section
###########################################################
SET(CPACK_BUNDLE_ALL OFF CACHE BOOL "Pack all necessary libs and tools for a standalone release")
IF(CPACK_BUNDLE_ALL)
	SET(ALL_RX_PERM OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
	SET(ALL_RO_PERM OWNER_READ GROUP_READ WORLD_READ)
	
	#standard libraries must be installed for gcc on Windows
	IF(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "^GNU$")
		SET(CPACK_BUNDLE_STDLIBS ON CACHE BOOL "Pack all necessary standard libs for a standalone release")
		IF(CPACK_BUNDLE_STDLIBS)
			FILE(GLOB winpthreads "${CMAKE_SOURCE_DIR}/bin/libwinpthread*.dll")
			FILE(GLOB stdlibs "${CMAKE_SOURCE_DIR}/bin/libstdc++*.dll")
			FILE(GLOB libgcc "${CMAKE_SOURCE_DIR}/bin/libgcc_s_seh*.dll")
			IF("${stdlibs}" STREQUAL "")
				MESSAGE("STDLIBS not found! Please make sure the dlls from your compiler for libwinpthread*, libstdc++* and libgcc_s_seh* are in ${CMAKE_SOURCE_DIR}/bin")
			ENDIF("${stdlibs}" STREQUAL "")
			INSTALL(FILES ${winpthreads} ${stdlibs} ${libgcc} DESTINATION bin PERMISSIONS ${ALL_RX_PERM} COMPONENT exe)
		ENDIF(CPACK_BUNDLE_STDLIBS)
	ENDIF(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "^GNU$")
ENDIF(CPACK_BUNDLE_ALL)

SET(CPACK_PACKAGE_DESCRIPTION "Meteorological preprocessing library")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Meteorological Input/Output preprocessing library: it covers data retrieval, data filtering, data transformations and spatial interpolations. Please refer to the included documentation and the included code examples!")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/Welcome.txt")
SET(CPACK_PACKAGE_NAME "MeteoIO")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}")
SET(CPACK_PACKAGE_CONTACT "Mathias Bavay <Mathias.Bavay@wsl.ch>")
SET(CPACK_PACKAGE_VENDOR "WSL Institute for Snow and Avalanche Research")
SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/meteoio/icons\\\\icon.ico")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/License.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/Welcome.txt")
SET(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/Readme.txt")

SET(CPACK_PACKAGE_VERSION_MAJOR "${VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${VERSION_PATCH}")
SET(VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_ARCH "${CMAKE_SYSTEM_PROCESSOR}")
SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${VERSION}-${CPACK_ARCH}")

#installation by components for formats that support it
#SET(CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE 1)
SET(CPACK_COMPONENTS_ALL libraries headers documentation examples exe resources)
SET(CPACK_COMPONENT_RESOURCES_DISPLAY_NAME "Inishell resources")
SET(CPACK_COMPONENT_RESOURCES_DESCRIPTION "Additional files to use Inishell")
SET(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "MeteoIO library")
SET(CPACK_COMPONENT_LIBRARIES_DESCRIPTION "The MeteoIO library itself")
SET(CPACK_COMPONENT_LIBRARIES_GROUP Runtime)
SET(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "Header files")
SET(CPACK_COMPONENT_HEADERS_DESCRIPTION "Necessary C++ header files for developing applications using MeteoIO")
SET(CPACK_COMPONENT_HEADERS_GROUP Development)
SET(CPACK_COMPONENT_HEADERS_DEPENDS libraries)
SET(CPACK_COMPONENT_EXE_DISPLAY_NAME "MeteoIO's applications")
SET(CPACK_COMPONENT_EXE_DESCRIPTION "MeteoIO's standalone applications")
SET(CPACK_COMPONENT_EXE_GROUP Runtime)
SET(CPACK_COMPONENT_RUNTIME_GROUP "Runtime")
SET(CPACK_COMPONENT_GROUP_RUNTIME_DESCRIPTION "What you need to run an application based on MeteoIO")
IF(INSTALL_DOC)
	SET(CPACK_COMPONENT_DOCUMENTATION_DISPLAY_NAME "Documentation")
	SET(CPACK_COMPONENT_DOCUMENTATION_DESCRIPTION "Documentation of the MeteoIO usage and API both for the end user and the developer using MeteoIO in his/her application")
	SET(CPACK_COMPONENT_DOCUMENTATION_GROUP Docs)
	SET(CPACK_COMPONENT_EXAMPLES_DISPLAY_NAME "Examples")
	SET(CPACK_COMPONENT_EXAMPLES_DESCRIPTION "C++ Source code of various simple examples of applications with their example data")
	SET(CPACK_COMPONENT_EXAMPLES_DEPENDS libraries headers)
	SET(CPACK_COMPONENT_EXAMPLES_GROUP Docs)
	SET(CPACK_COMPONENT_DOCS_GROUP "Documentation")
	SET(CPACK_COMPONENT_GROUP_DOCS_DESCRIPTION "How to use MeteoIO as an end user and as a developer")
ENDIF(INSTALL_DOC)

SET(CPACK_COMPONENT_RUNTIME_GROUP "Runtime")
SET(CPACK_COMPONENT_GROUP_RUNTIME_DESCRIPTION "What you need to run an application based on MeteoIO")
SET(CPACK_COMPONENT_DEVELOPMENT_GROUP "Development")
SET(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION "What you need to develop an application based on MeteoIO")

IF(UNIX)
	IF(APPLE)	#OsX packages: IMG
		SET(CPACK_GENERATOR "DragNDrop")
		SET(CPACK_PACKAGING_INSTALL_PREFIX "/${CPACK_PACKAGE_NAME}")
		SET(CPACK_DMG_FORMAT "UDBZ")
		#SET(CPACK_GENERATOR "Bundle")
		SET(CPACK_BUNDLE_PLIST "${CMAKE_SOURCE_DIR}/tools/meteoio.plist")
		SET(MACOSX_BUNDLE_INFO_STRING "${CPACK_PACKAGE_DESCRIPTION_SUMMARY}")
		SET(CPACK_BUNDLE_ICON "${CMAKE_SOURCE_DIR}/meteoio/icons/icon.icns")
		SET(MACOSX_BUNDLE_ICON_FILE "${CPACK_BUNDLE_ICON}")
		SET(MACOSX_BUNDLE_GUI_IDENTIFIER "${CPACK_PACKAGE_DESCRIPTION}")
		SET(MACOSX_BUNDLE_LONG_VERSION_STRING "${VERSION}")
		SET(CPACK_BUNDLE_NAME "${CPACK_PACKAGE_NAME}")
		SET(MACOSX_BUNDLE_BUNDLE_NAME "${CPACK_PACKAGE_NAME}")
		SET(MACOSX_BUNDLE_SHORT_VERSION_STRING "${VERSION}")
		SET(MACOSX_BUNDLE_BUNDLE_VERSION "0")
		SET(MACOSX_BUNDLE_COPYRIGHT "${CPACK_PACKAGE_VENDOR}")
		#SET(CPACK_PACKAGING_INSTALL_PREFIX "/Applications/${CPACK_PACKAGE_NAME}")
		SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/meteoio/icons/icon.icns")
		SET_SOURCE_FILES_PROPERTIES(${CPACK_BUNDLE_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
	ELSE(APPLE)	#Linux packages: rmp, deb and tgz
		SET(CPACK_GENERATOR "DEB;RPM;TGZ" CACHE STRING "Which packages to build")
		MARK_AS_ADVANCED(FORCE CPACK_GENERATOR)
		SET(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
		SET(CPACK_DEBIAN_PACKAGE_NAME "meteoio")
		SET(CPACK_DEBIAN_PACKAGE_SECTION "devel")
		SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6") #objdump -p lib/libmeteoio.so | grep NEEDED | tr -s ' ' | cut -d' ' -f3 | tr '\n' ','
		#SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386") #dpkg --print-architecture
		IF("${CPACK_ARCH}" STREQUAL "x86_64")
			SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
		ELSE("${CPACK_ARCH}" STREQUAL "x86_64")
			SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "${CPACK_ARCH}")
		ENDIF("${CPACK_ARCH}" STREQUAL "x86_64")
		SET(CPACK_RPM_PACKAGE_NAME "meteoio")
		SET(CPACK_RPM_PACKAGE_LICENSE "LGPLv3")
		SET(CPACK_RPM_PACKAGE_REQUIRES "libstdc++") #simply libproj for OpenSuse
		SET(CPACK_RPM_PACKAGE_ARCHITECTURE "${CPACK_ARCH}")
		SET(CPACK_RPM_PACKAGE_GROUP "Development/Libraries")
		SET(CPACK_RPM_SPEC_INSTALL_POST "/bin/true") #disable running strip on the binaries (fails on some distributions and makes debugging much harder)
	ENDIF(APPLE)
ELSE(UNIX)
	IF(WIN32)	#Windows packages
		SET(CPACK_GENERATOR "NSIS")
		SET(CPACK_NSIS_MUI_ICON "${CPACK_PACKAGE_ICON}")
		SET(CPACK_PACKAGE_ICON "${CPACK_PACKAGE_ICON}")
		SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME}")
		SET(CPACK_NSIS_MODIFY_PATH ON)
		
		#Main links (they can get overwritten later)
		SET(MENU_LINKS_NIVIZ "https://run.niviz.org" "NIVIZ online visualization" )
		SET(MENU_LINKS_METEOIO "https://meteoio.slf.ch" "MeteoIO Home page")
		IF(INSTALL_DOC)
			LIST(APPEND MENU_LINKS_METEOIO
			"${DOC_PATH}/MeteoIO_doc.html" "MeteoIO documentation"
			"${DOC_PATH}/examples" "MeteoIO examples"
			)
		ENDIF(INSTALL_DOC)
		#SET(CPACK_NSIS_EXTRA_INSTALL_COMMANDS "
		#	WriteRegStr HKLM 'SOFTWARE\\\\Microsoft\\\\Windows\\\\CurrentVersion\\\\App Paths\\\\meteoio_timeseries' 'c:\\\\Program Files\\\\MeteoIO\\\\bin\\\\meteoio_timeseries.exe'
		#	WriteRegStr HKLM 'SOFTWARE\\\\Microsoft\\\\Windows\\\\CurrentVersion\\\\App Paths\\\\meteoio_timeseries\\\\path' 'c:\\\\Program Files\\\\MeteoIO\\\\lib'
		#	")
		
		#all other links: inishell
		IF(CPACK_BUNDLE_ALL)
# 			SET(MENU_LINKS_INISHELL "bin/inishell.exe" "MeteoIO's Inishell")
		ENDIF(CPACK_BUNDLE_ALL)
		
		SET(CPACK_NSIS_MENU_LINKS ${MENU_LINKS_METEOIO} ${MENU_LINKS_INISHELL} ${MENU_LINKS_NIVIZ})
		
		#info for the add/uninstall software administration
		SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_NAME}-${VERSION}")
		SET(CPACK_NSIS_INSTALLED_ICON_NAME "doc\\\\icon.ico")
		SET(CPACK_NSIS_URL_INFO_ABOUT "http://www.slf.ch")
		SET(CPACK_NSIS_HELP_LINK "https://meteoio.slf.ch")
		SET(CPACK_NSIS_CONTACT "${CPACK_PACKAGE_CONTACT}")
		SET(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL "ON")
	ENDIF(WIN32)
ENDIF(UNIX)

#Source packages
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${VERSION}" CACHE INTERNAL "tarball basename")
SET(CPACK_SOURCE_GENERATOR TGZ)
# The following components are regex's to match anywhere (unless anchored)
# in absolute path + filename to find files or directories to be excluded
# from source tarball.
SET(CPACK_SOURCE_IGNORE_FILES
	"~$"
	"\\\\.svn/"
	"/CMakeFiles/"
	"/CMakeCache*"
	"/cmake_install\\\\.cmake$"
	"/cmake_uninstall\\\\.cmake$"
	"^_CPack_.*/"
	"/CPack*"
	"\\\\.o$"
	"/m4/"
)

INCLUDE(CPack)

###########################################################
## CTest section
###########################################################
OPTION(BUILD_TESTING "Build the tests" OFF)
IF(BUILD_TESTING)
	# Special build type configuration for Coverage Testing.
	#Coverage testing only supported by GCC and a few commercial tools
	IF(CMAKE_COMPILER_IS_GNUCXX)
		OPTION(BUILD_TESTING_WITH_COVERAGE "Use gcc flags to make coverage tests" OFF)
		IF(BUILD_TESTING_WITH_COVERAGE)
			# force the correct flags to make coverage trought a build type
			SET(CMAKE_BUILD_TYPE coverage)
			SET(CMAKE_CXX_FLAGS_COVERAGE "${OPTIM} -ftest-coverage -fprofile-arcs" CACHE STRING "" FORCE)
			MARK_AS_ADVANCED(FORCE CMAKE_CXX_FLAGS_COVERAGE)
			# Force debug mode, so that no optimisation is made to get "Correct" results of coverage testing
			SET(CMAKE_BUILD_TYPE "coverage" CACHE STRING "Choose the type of build, options are: debug release coverage " FORCE)
			SET(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE})
		ENDIF(BUILD_TESTING_WITH_COVERAGE)
	ENDIF(CMAKE_COMPILER_IS_GNUCXX)

	INCLUDE(CTest) # This makes ENABLE_TESTING() and gives support for Dashboard
	ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTING)
