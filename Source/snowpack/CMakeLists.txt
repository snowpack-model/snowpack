#SPDX-License-Identifier: LGPL-3.0-or-later
CMAKE_MINIMUM_REQUIRED(VERSION 3.1...3.31.6)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/tools/cmake/")

PROJECT(snowpack)
# Version of library
IF(DEFINED ENV{VERSION_MAJOR} AND DEFINED ENV{VERSION_MINOR} AND DEFINED ENV{VERSION_PATCH})
	SET(VERSION_MAJOR $ENV{VERSION_MAJOR} CACHE INTERNAL "Major version")
	SET(VERSION_MINOR $ENV{VERSION_MINOR} CACHE INTERNAL "Minor version")
	SET(VERSION_PATCH $ENV{VERSION_PATCH} CACHE INTERNAL "Patch version")
ELSE()
	SET(VERSION_MAJOR 3 CACHE INTERNAL "Major version")
	SET(VERSION_MINOR 7 CACHE INTERNAL "Minor version")
	SET(VERSION_PATCH 0 CACHE INTERNAL "Patch version")
ENDIF()

SET(LIBPREFIX "lib") #make sure ALL plugins are libXXX for ALL plateforms
SET(CMAKE_IMPORT_LIBRARY_PREFIX "lib") #for Visual C++ import libraries
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
	SET(LIB_DEST "lib")

	OPTION(USE_RPATH "Use rpath (not recommended for Fedora)" ON)
	MARK_AS_ADVANCED(FORCE USE_RPATH)
	if(USE_RPATH)
		#set a full rpath for both build & install tree
		SET(CMAKE_SKIP_BUILD_RPATH FALSE)
		SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
		# add the directory of the binary to the rpath
		IF(APPLE)
			SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'@executable_path',-rpath,'@executable_path/../lib'")
			SET(CMAKE_OSX_DEPLOYMENT_TARGET "10.13" CACHE STRING "Minimum OS X deployment version")
		ELSE(APPLE)
			SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'$ORIGIN',-rpath,'$ORIGIN'/../lib")
		ENDIF(APPLE)
		
		# add the automatically determined parts of the RPATH
		SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
		# the RPATH to be used when installing, but only if it's not a system directory
		LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
		IF("${isSystemDir}" STREQUAL "-1")
			SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
		ENDIF("${isSystemDir}" STREQUAL "-1")
	ELSE(USE_RPATH)
		SET(CMAKE_SKIP_BUILD_RPATH TRUE)
	ENDIF(USE_RPATH)

	IF(APPLE)
		SET(DLL_EXT "dylib")
		SET(STAT_EXT "a")
		#set(CMAKE_OSX_SYSROOT iphoneos2.2.1) #for IOS
		#set(CMAKE_OSX_ARCHITECTURES $(ARCHS_STANDARD_32_BIT)) #for IOS
	ELSE(APPLE)
		SET(DLL_EXT "so*")
		SET(STAT_EXT "a")
	ENDIF(APPLE)
ELSE(UNIX)
	IF(WIN32)
		SET(LIB_DEST "bin") #since for Windows, dlls have to be alongside the exe
		#SET(PLATFORM "-D_WIN32 -U__STRICT_ANSI__")
		ADD_DEFINITIONS(-U__STRICT_ANSI__)
		SET(DLL_EXT "dll")
		SET(STAT_EXT "a")
	ENDIF(WIN32)
ENDIF(UNIX)

IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE "release" CACHE STRING "Choose the type of build, options are: debug release coverage cxxflags" FORCE)
	SET(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE})
ENDIF(NOT CMAKE_BUILD_TYPE)

FIND_PACKAGE(MeteoIO REQUIRED)

###########################################################
#plugins and applications
SET(SNOWPACK_STANDALONE ON CACHE BOOL "Compilation SNOWPACK standalone ON or OFF")
SET(PROFEVAL OFF CACHE BOOL "Compilation PROFEVAL ON or OFF")

SET(ENABLE_LAPACK OFF CACHE BOOL "Compile with the CLAPACK library?")
SET(PLUGIN_IMISIO OFF CACHE BOOL "Compilation IMISDBIO ON or OFF - only relevant for SLF")
SET(PLUGIN_CAAMLIO ON CACHE BOOL "Compilation CAAMLIO ON or OFF to read CAAML profiles")
###########################################################
## Additional Options
###########################################################
IF(GCC_VERSION VERSION_GREATER 5.0 OR GCC_VERSION VERSION_EQUAL 5.0)
	IF(PLUGIN_IMISIO) #HACK: current OCCI does not support the short strings optimizations of gcc>=5
		SET(EXTRA "-D_GLIBCXX_USE_CXX11_ABI=0 ${EXTRA}")
	ELSE(PLUGIN_IMISIO)
		SET(FORCE_C03_ABI OFF CACHE BOOL "Compatibility with Oracle's OCCI")
		IF (FORCE_C03_ABI)
			SET(EXTRA "-D_GLIBCXX_USE_CXX11_ABI=0 ${EXTRA}")
		ENDIF(FORCE_C03_ABI)
	ENDIF(PLUGIN_IMISIO)
ENDIF()

###########################################################
#finally, SET compile flags
SET(CMAKE_CXX_FLAGS "${_VERSION} ${ARCH} ${EXTRA}" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_RELEASE "${OPTIM}" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_DEBUG "${DEBUG} ${WARNINGS} ${EXTRA_WARNINGS}" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_CXXFLAGS "$ENV{CXXFLAGS}" CACHE STRING "" FORCE)
MARK_AS_ADVANCED(FORCE CMAKE_CXX_FLAGS_CXXFLAGS)

#Add subdirectory, so they can inherit all variables defined above
ADD_SUBDIRECTORY(snowpack)

IF(SNOWPACK_STANDALONE)
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
	ADD_SUBDIRECTORY(applications/snowpack)
ENDIF(SNOWPACK_STANDALONE)

IF(PROFEVAL)
	ADD_SUBDIRECTORY(applications/profeval)
ENDIF(PROFEVAL)

###########################################################
#for the install target
FILE(GLOB hfiles "snowpack/*.h")
FILE(GLOB hpluginfiles "snowpack/plugins/*.h")
FILE(GLOB pugifiles "snowpack/plugins/pugixml/pugi*.?pp")
FILE(GLOB hcorefiles "snowpack/snowpackCore/*.h")
INSTALL(FILES ${hfiles} DESTINATION include/snowpack COMPONENT headers)
INSTALL(FILES ${hpluginfiles} DESTINATION include/snowpack/plugins COMPONENT headers)
INSTALL(FILES ${pugifiles} DESTINATION include/snowpack/plugins/pugixml COMPONENT headers)
INSTALL(FILES ${hcorefiles} DESTINATION include/snowpack/snowpackCore COMPONENT headers)

ADD_CUSTOM_TARGET(distclean make clean
	COMMAND cmake -E remove CMakeCache.txt CTestTestfile.cmake DartConfiguration.tcl CPack*.cmake cmake_*.cmake install_manifest.txt Makefile CMakeDoxy*
	COMMAND cmake -E remove "lib/libsnow*.${DLL_EXT}" "lib/libsnow*.${STAT_EXT}" "lib/libsnow*.lib"
	COMMAND cmake -E remove */CMakeCache.txt */CTestTestfile.cmake */cmake_install.cmake
	COMMAND cmake -E remove *~ snowpack/*~ snowpack/*/*~ applications/snowpack/*~
	COMMAND cmake -E remove applications/snowpack/cmake_install.cmake applications/snowpack/Makefile  applications/snowpack/cmake_*.cmake
	COMMAND cmake -E remove */*/CMakeCache.txt */*/CTestTestfile.cmake */*/cmake_install.cmake
	COMMAND cmake -E remove_directory CMakeFiles snowpack/CMakeFiles snowpack/*/CMakeFiles applications/snowpack/CMakeFiles
	COMMAND cmake -E remove compile_commands.json
	COMMAND cmake -E remove_directory Testing
	COMMAND cmake -E remove_directory tests/CMakeFiles tests/*/CMakeFiles
	COMMAND cmake -E remove tests/Makefile tests/*/Makefile tests/CTestTestfile.cmake tests/*/CTestTestfile.cmake tests/cmake_install.cmake tests/*/cmake_install.cmake
	COMMAND cmake -E remove tests/*~ tests/*/*~
	COMMAND cmake -E remove tests/res1exp/output/* tests/res5exp/output/*
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
IF(CMAKE_VERSION VERSION_GREATER 3.9)
	FIND_PACKAGE(Doxygen OPTIONAL_COMPONENTS dot)
ELSE(CMAKE_VERSION VERSION_GREATER 3.9)
	FIND_PACKAGE(Doxygen)
ENDIF(CMAKE_VERSION VERSION_GREATER 3.9)
IF(DOXYGEN_FOUND)
	ADD_CUSTOM_TARGET(doc COMMAND ${CMAKE_COMMAND} -E env SN_VERSION="${_versionString}" "doxygen" ${PROJECT_SOURCE_DIR}/snowpack/config.dox)
ENDIF(DOXYGEN_FOUND)
#Install the documentation if the user desires it
SET(INSTALL_DOC OFF CACHE BOOL "Install documentation ON or OFF")
IF(INSTALL_DOC)
	IF(WIN32)
		SET(DOC_PATH "doc")
	ELSE(WIN32)
		SET(DOC_PATH "share/doc/${CMAKE_PROJECT_NAME}")
	ENDIF(WIN32)
	SET(DOC_FILES License.txt Readme.md snowpack/icons/icon.ico)
	IF (DOXYGEN_FOUND)
		INSTALL(DIRECTORY doc/html DESTINATION ${DOC_PATH} COMPONENT documentation)
		LIST(APPEND DOC_FILES doc/Snowpack_doc.html)
	ELSE (DOXYGEN_FOUND)
		MESSAGE("Doxygen has not been found and is needed to build online documentation. Please either install doxygen or download a documentation-only package!")
	ENDIF (DOXYGEN_FOUND)

	INSTALL(FILES ${DOC_FILES} DESTINATION ${DOC_PATH} COMPONENT documentation)
	INSTALL(DIRECTORY doc/examples USE_SOURCE_PERMISSIONS DESTINATION ${DOC_PATH} COMPONENT examples PATTERN ".svn" EXCLUDE)
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
 		        MESSAGE("STDLIBS not found! Please make sure they are in ${CMAKE_SOURCE_DIR}/bin")
			ENDIF("${stdlibs}" STREQUAL "")
			INSTALL(FILES ${winpthreads} ${stdlibs} ${libgcc} DESTINATION bin PERMISSIONS ${ALL_RX_PERM} COMPONENT exe)
		ENDIF(CPACK_BUNDLE_STDLIBS)
	ENDIF(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "^GNU$")
	
	#bundle MeteoIO
	FIND_PACKAGE(MeteoIO REQUIRED)
	IF(BUILD_SHARED_LIBS) #consider that if STATIC, we don't need the libs
		FILE(GLOB meteoio_libs "${METEOIO_ROOT}/lib/libmeteo*.${DLL_EXT}")
		IF(WIN32)
			INSTALL(FILES ${meteoio_libs} DESTINATION bin COMPONENT exe)
		ELSE(WIN32)
			INSTALL(FILES ${meteoio_libs} DESTINATION lib COMPONENT exe)
		ENDIF(WIN32)
		FILE(GLOB meteoio_bins "${METEOIO_ROOT}/bin/meteoio_timeserie*")
		INSTALL(FILES ${meteoio_bins} DESTINATION bin PERMISSIONS ${ALL_RX_PERM} COMPONENT exe)
		INSTALL(DIRECTORY ${METEOIO_INCLUDE_DIR}/meteoio/ DESTINATION include/meteoio COMPONENT headers 
		        FILES_MATCHING PATTERN "*.h" 
		        PATTERN "CMakeFiles" EXCLUDE 
		        PATTERN "icons" EXCLUDE 
		        PATTERN "jnative" EXCLUDE)
	ENDIF(BUILD_SHARED_LIBS)
	IF(INSTALL_DOC)
		#meteoio's documentation
		IF(WIN32)
			SET(METEOIO_INSTALLED_DOC "${METEOIO_ROOT}/doc")
			SET(METEOIO_DOC ${DOC_PATH}/meteoio)
		ELSE(WIN32)
			SET(METEOIO_INSTALLED_DOC "${METEOIO_ROOT}/share/doc/meteoio")
			SET(METEOIO_DOC ${DOC_PATH}/../meteoio)
		ENDIF(WIN32)
		IF(NOT IS_DIRECTORY "${METEOIO_INSTALLED_DOC}") #if we detected neighboring checked out version of MeteoIO
			SET(METEOIO_INSTALLED_DOC "${METEOIO_ROOT}/doc")
		ENDIF()
		INSTALL(DIRECTORY "${METEOIO_INSTALLED_DOC}/html" DESTINATION ${METEOIO_DOC} COMPONENT documentation)
		FILE(GLOB meteoio_docs "${METEOIO_INSTALLED_DOC}/*.*")
		INSTALL(FILES ${meteoio_docs} DESTINATION ${METEOIO_DOC} COMPONENT documentation)
	ENDIF(INSTALL_DOC)
ENDIF(CPACK_BUNDLE_ALL)

SET(CPACK_PACKAGE_DESCRIPTION "Snowpack Model")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "This is the SNOWPACK model ([Bartelt 2002], [Lehning 2002a], [Lehning 2002b]) for physical modeling of snow pack. This consists of a library and a program calling it so it can easily be integrated into other models.")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/Welcome.txt")
SET(CPACK_PACKAGE_NAME "Snowpack")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}")
SET(CPACK_PACKAGE_CONTACT "Mathias Bavay <Mathias.Bavay@wsl.ch>")
SET(CPACK_PACKAGE_VENDOR "WSL Institute for Snow and Avalanche Research")
SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/snowpack/icons\\\\icon.ico")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/License.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/Welcome.txt")
SET(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/Readme.md")

SET(CPACK_PACKAGE_VERSION_MAJOR "${VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${VERSION_PATCH}")
SET(VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_ARCH "${CMAKE_SYSTEM_PROCESSOR}")
SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}${PACKAGE_EXTRA}-${VERSION}-${CPACK_ARCH}")

#installation by components for formats that support it
#SET(CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE 1)
SET(CPACK_ARCHIVE_COMPONENT_INSTALL ON)
SET(CPACK_COMPONENTS_ALL documentation examples exe resources)
SET(CPACK_COMPONENT_RESOURCES_DISPLAY_NAME "Inishell resources")
SET(CPACK_COMPONENT_RESOURCES_DESCRIPTION "Additional files to use Inishell")
SET(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "Snowpack library")
SET(CPACK_COMPONENT_LIBRARIES_DESCRIPTION "The SNOWPACK library")
SET(CPACK_COMPONENT_LIBRARIES_GROUP Runtime)
SET(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "Header files")
SET(CPACK_COMPONENT_HEADERS_DESCRIPTION "Necessary C++ header files for developing applications using SNOWPACK")
SET(CPACK_COMPONENT_HEADERS_GROUP Development)
SET(CPACK_COMPONENT_HEADERS_DEPENDS libraries)
SET(CPACK_COMPONENT_DEVELOPMENT_GROUP "Development")
SET(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION "What you need to develop an application based on SNOWPACK")
SET(CPACK_COMPONENT_EXE_DISPLAY_NAME "Snowpack application")
SET(CPACK_COMPONENT_EXE_DESCRIPTION "The SNOWPACK standalone model")
SET(CPACK_COMPONENT_EXE_GROUP Runtime)
SET(CPACK_COMPONENT_RUNTIME_GROUP "Runtime")
SET(CPACK_COMPONENT_GROUP_RUNTIME_DESCRIPTION "What you need to run an application based on SNOWPACK")
IF(INSTALL_DOC)
	SET(CPACK_COMPONENT_DOCUMENTATION_DISPLAY_NAME "Documentation")
	SET(CPACK_COMPONENT_DOCUMENTATION_DESCRIPTION "Documentation of the SNOWPACK usage and API both for the end user and the developer using SNOWPACK in his/her application")
	SET(CPACK_COMPONENT_DOCUMENTATION_GROUP Docs)
	SET(CPACK_COMPONENT_EXAMPLES_DISPLAY_NAME "Examples")
	SET(CPACK_COMPONENT_EXAMPLES_DESCRIPTION "Various simulations with their example data")
	SET(CPACK_COMPONENT_EXAMPLES_GROUP Docs)
	#SET(CPACK_COMPONENT_EXAMPLES_DEPENDS libraries headers runtime)
	SET(CPACK_COMPONENT_DOCS_GROUP "Documentation")
	SET(CPACK_COMPONENT_GROUP_DOCS_DESCRIPTION "How to use Snowpack as a standalone model or as a library called by your own model")
ENDIF(INSTALL_DOC)

IF(UNIX)
	IF(APPLE)
		SET(CPACK_GENERATOR "DragNDrop")
		SET(CPACK_PACKAGING_INSTALL_PREFIX "/${CPACK_PACKAGE_NAME}")
		SET(CPACK_DMG_FORMAT "UDBZ")
		#SET(CPACK_GENERATOR "Bundle")
		SET(CPACK_BUNDLE_PLIST "${CMAKE_SOURCE_DIR}/tools/snowpack.plist")
		SET(MACOSX_BUNDLE_INFO_STRING "${CPACK_PACKAGE_DESCRIPTION_SUMMARY}")
		SET(CPACK_BUNDLE_ICON "${CMAKE_SOURCE_DIR}/snowpack/icons/icon.icns")
		SET(MACOSX_BUNDLE_ICON_FILE "${CPACK_BUNDLE_ICON}")
		SET(MACOSX_BUNDLE_GUI_IDENTIFIER "${CPACK_PACKAGE_DESCRIPTION}")
		SET(MACOSX_BUNDLE_LONG_VERSION_STRING "${VERSION}")
		SET(CPACK_BUNDLE_NAME "${CPACK_PACKAGE_NAME}")
		SET(MACOSX_BUNDLE_BUNDLE_NAME "${CPACK_PACKAGE_NAME}")
		SET(MACOSX_BUNDLE_SHORT_VERSION_STRING "${VERSION}")
		SET(MACOSX_BUNDLE_BUNDLE_VERSION "0")
		SET(MACOSX_BUNDLE_COPYRIGHT "${CPACK_PACKAGE_VENDOR}")
		#SET(CPACK_PACKAGING_INSTALL_PREFIX "/Applications/${CPACK_PACKAGE_NAME}")
		SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/snowpack/icons/icon.icns")
		SET_SOURCE_FILES_PROPERTIES(${CPACK_BUNDLE_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
		INSTALL(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink bin/snowpack /usr/bin/snowpack)")
	ELSE(APPLE)
		SET(CPACK_ARCHIVE_COMPONENT_INSTALL OFF)
		SET(CPACK_GENERATOR "DEB;RPM;TGZ" CACHE STRING "Which packages to build")
		MARK_AS_ADVANCED(FORCE CPACK_GENERATOR)
		SET(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
		SET(CPACK_DEBIAN_PACKAGE_NAME "snowpack")
		SET(CPACK_DEBIAN_PACKAGE_SECTION "science")
		#objdump -p lib/libmeteoio.so | grep NEEDED | tr -s ' ' | cut -d' ' -f3 | tr '\n' ','
		IF(CPACK_BUNDLE_ALL OR BUILD_STATIC_LIBS)
			SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6, libc6")
		ELSE(CPACK_BUNDLE_ALL OR BUILD_STATIC_LIBS)
			SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++, meteoio")
		ENDIF(CPACK_BUNDLE_ALL OR BUILD_STATIC_LIBS)
		#SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386") #dpkg --print-architecture
		IF("${CPACK_ARCH}" STREQUAL "x86_64")
			SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
		ELSE("${CPACK_ARCH}" STREQUAL "x86_64")
			SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "${CPACK_ARCH}")
		ENDIF("${CPACK_ARCH}" STREQUAL "x86_64")
		SET(CPACK_RPM_PACKAGE_NAME "snowpack")
		SET(CPACK_RPM_PACKAGE_LICENSE "LPGLv3")
		SET(CPACK_RPM_PACKAGE_REQUIRES "libstdc++, meteoio")
		SET(CPACK_RPM_PACKAGE_ARCHITECTURE "${CPACK_ARCH}")
		SET(CPACK_RPM_PACKAGE_GROUP "Applications/Science")
	ENDIF(APPLE)
ELSE(UNIX)
	IF(WIN32)
		SET(CPACK_GENERATOR "NSIS")
		SET(CPACK_NSIS_MUI_ICON "${CPACK_PACKAGE_ICON}")
		SET(CPACK_PACKAGE_ICON "${CPACK_PACKAGE_ICON}")
		SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME}")
		SET(CPACK_NSIS_MODIFY_PATH ON)
		
		#Main links (they can get overwritten later)
		SET(MENU_LINKS_NIVIZ "https://run.niviz.org" "NIVIZ online visualization" )
		SET(MENU_LINKS_SNOWPACK "https://snowpack.slf.ch" "SNOWPACK Home page")
		IF(INSTALL_DOC)
			LIST(APPEND MENU_LINKS_SNOWPACK
			"${DOC_PATH}/Snowpack_doc.html" "SNOWPACK documentation"
			"${DOC_PATH}/examples" "SNOWPACK examples"
			)
		ENDIF(INSTALL_DOC)
		
		#all other links: inishell, meteoio
		IF(CPACK_BUNDLE_ALL)
# 			SET(MENU_LINKS_INISHELL "bin/inishell.exe" "SNOWPACK's Inishell")
			SET(MENU_LINKS_METEOIO "https://meteoio.slf.ch" "MeteoIO Home page")
			IF(INSTALL_DOC)
				SET(MENU_LINKS_METEOIO
				"${METEOIO_DOC}/MeteoIO_doc.html" "MeteoIO documentation"
				"https://meteoio.slf.ch" "MeteoIO Home page"
				)
			ENDIF(INSTALL_DOC)
		ENDIF(CPACK_BUNDLE_ALL)
		
		SET(CPACK_NSIS_MENU_LINKS ${MENU_LINKS_SNOWPACK} ${MENU_LINKS_INISHELL} ${MENU_LINKS_NIVIZ} ${MENU_LINKS_METEOIO})
		
		#info for the add/uninstall software administration
		SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_NAME}-${VERSION}")
		SET(CPACK_NSIS_INSTALLED_ICON_NAME "doc\\\\icon.ico")
		SET(CPACK_NSIS_URL_INFO_ABOUT "http://www.slf.ch")
		SET(CPACK_NSIS_HELP_LINK "https://snowpack.slf.ch")
		SET(CPACK_NSIS_CONTACT "${CPACK_PACKAGE_CONTACT}")
		SET(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL "ON")
	ENDIF(WIN32)
ENDIF(UNIX)


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
