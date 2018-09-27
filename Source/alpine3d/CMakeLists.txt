CMAKE_MINIMUM_REQUIRED(VERSION 2.8.4)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/tools/cmake/")

PROJECT(alpine3d)
# Version of library
SET(VERSION_MAJOR 3 CACHE INTERNAL "")
SET(VERSION_MINOR 1 CACHE INTERNAL "")
SET(VERSION_PATCH 0 CACHE INTERNAL "")

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
INCLUDE("${CMAKE_SOURCE_DIR}/tools/cmake/SetCompilerOptions.cmake")
SET_COMPILER_OPTIONS()

#The following pre-processor variables are automatically defined:
# __linux, __GNUC__, __CYGWIN__, _WIN32, __MINGW32__, __APPLE__
IF(UNIX)
	SET(LIB_DEST "lib")
	#set a full rpath for both build & install tree
	SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
	SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
	# add the automatically determined parts of the RPATH
	SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
	# the RPATH to be used when installing, but only if it's not a system directory
	LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
	IF("${isSystemDir}" STREQUAL "-1")
		SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
	ENDIF("${isSystemDir}" STREQUAL "-1")

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
FIND_PACKAGE(Libsnowpack REQUIRED)

###########################################################
#plugins and applications
SET(ALPINE3D_STANDALONE ON CACHE BOOL "Compilation Alpine3D standalone ON or OFF")
SET(MPI OFF CACHE BOOL "Compile with MPI support ON or OFF")
SET(OPENMP OFF CACHE BOOL "Compile with OPENMP support ON or OFF")
SET(PETSC OFF CACHE BOOL "Compile with PETSc support ON or OFF")

IF(PETSC)
	SET(PETSC_FLAGS "-DENABLE_PETSC")
	UNSET(MPI CACHE)
	SET (MPI ON CACHE BOOL "Compile with MPI support ON or OFF")
ENDIF(PETSC)
IF(MPI)
	SET(MPI_FLAGS "-DENABLE_MPI")
ENDIF(MPI)
IF(OPENMP)
	SET(OPENMP_FLAGS "-fopenmp")
ENDIF(OPENMP)

###########################################################
#finally, SET compile flags
SET(CMAKE_CXX_FLAGS "${PETSC_FLAGS} ${MPI_FLAGS} ${OPENMP_FLAGS} ${_VERSION} ${ARCH} ${CFLAGS} ${EXTRA}" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_RELEASE "${OPTIM}" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_DEBUG "${DEBUG} ${WARNINGS} ${EXTRA_WARNINGS}" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_CXXFLAGS "$ENV{CXXFLAGS}" CACHE STRING "" FORCE)
MARK_AS_ADVANCED(FORCE CMAKE_CXX_FLAGS_CXXFLAGS)

#Add subdirectory, so they can inherit all variables defined above
ADD_SUBDIRECTORY(alpine3d)

IF(ALPINE3D_STANDALONE)
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
ENDIF(ALPINE3D_STANDALONE)

###########################################################
#for the install target
FILE(GLOB hfiles_gen "alpine3d/*.h")
FILE(GLOB hfiles_eb "alpine3d/ebalance/*.h")
FILE(GLOB hfiles_rot "alpine3d/runoff/*.h")
FILE(GLOB hfiles_drift "alpine3d/snowdrift/*.h")
INSTALL(FILES ${hfiles_gen} DESTINATION include/alpine3d COMPONENT headers)
INSTALL(FILES ${hfiles_eb} DESTINATION include/alpine3d/ebalance COMPONENT headers)
INSTALL(FILES ${hfiles_rot} DESTINATION include/alpine3d/runoff COMPONENT headers)
INSTALL(FILES ${hfiles_drift} DESTINATION include/alpine3d/snowdrift COMPONENT headers)

ADD_CUSTOM_TARGET(distclean make clean
	COMMAND cmake -E remove *~ alpine3d/*~ alpine3d/ebalance/*~ alpine3d/runoff/*~ alpine3d/snowdrift/*~
	COMMAND cmake -E remove CMakeCache.txt cmake_install.cmake
	COMMAND cmake -E remove alpine3d/*.cmake alpine3d/ebalance/*.cmake alpine3d/runoff/*.cmake alpine3d/snowdrift/*.cmake
	COMMAND cmake -E remove_directory CMakeFiles
	COMMAND cmake -E remove_directory alpine3d/CMakeFiles
	COMMAND cmake -E remove_directory alpine3d/ebalance/CMakeFiles
	COMMAND cmake -E remove_directory alpine3d/runoff/CMakeFiles
	COMMAND cmake -E remove_directory alpine3d/snowdrift/CMakeFiles
	COMMAND cmake -E remove bin/*
	COMMAND cmake -E remove lib/*
)

ADD_CUSTOM_TARGET(deploy
	DEPENDS alpine3d.app
	COMMAND cp -Rf "${METEOIO_ROOT}/lib/libmeteo*" bin/
	COMMAND cp -Rf "${LIBSNOWPACK_ROOT}/lib/libsnowpac*" bin/
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
IF(DOXYGEN_FOUND)
	ADD_CUSTOM_TARGET(doc "doxygen" ${PROJECT_SOURCE_DIR}/alpine3d/config.dox)
ENDIF(DOXYGEN_FOUND)
#Install the documentation if the user desires it
SET(INSTALL_DOC OFF CACHE BOOL "Install documentation ON or OFF")
IF(INSTALL_DOC)
	IF(WIN32)
		SET(DOC_PATH "doc")
	ELSE(WIN32)
		SET(DOC_PATH "share/doc/${CMAKE_PROJECT_NAME}")
	ENDIF(WIN32)
	SET(DOC_FILES License.txt Readme.txt alpine3d/icons/icon.ico)
	IF (DOXYGEN_FOUND)
		INSTALL(DIRECTORY doc/html DESTINATION ${DOC_PATH} COMPONENT documentation)
		LIST(APPEND DOC_FILES doc/Alpine3D_doc.html)
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
	#bundle MeteoIO
	FIND_PACKAGE(MeteoIO REQUIRED)
	IF(BUILD_SHARED_LIBS) #consider that if STATIC, we don't need the libs
		FILE(GLOB meteoio_libs "${METEOIO_ROOT}/lib/libmeteo*.${DLL_EXT}")
		IF(WIN32)
			INSTALL(FILES ${meteoio_libs} DESTINATION bin COMPONENT exe)
		ELSE(WIN32)
			INSTALL(FILES ${meteoio_libs} DESTINATION lib COMPONENT exe)
		ENDIF(WIN32)
	ENDIF(BUILD_SHARED_LIBS)
	#bundle Snwpack
	FIND_PACKAGE(Libsnowpack REQUIRED)
	IF(BUILD_SHARED_LIBS) #consider that if STATIC, we don't need the libs
		FILE(GLOB snowpack_libs "${LIBSNOWPACK_ROOT}/${LIB_DEST}/libsnowpack*.${DLL_EXT}")
		IF(WIN32)
			INSTALL(FILES ${snowpack_libs} DESTINATION bin COMPONENT exe)
		ELSE(WIN32)
			INSTALL(FILES ${snowpack_libs} DESTINATION lib COMPONENT exe)
		ENDIF(WIN32)
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
		#snowpack's documentation
		IF(WIN32)
			SET(SNOWPACK_INSTALLED_DOC "${LIBSNOWPACK_ROOT}/doc")
			SET(SNOWPACK_DOC ${DOC_PATH}/snowpack)
		ELSE(WIN32)
			SET(SNOWPACK_INSTALLED_DOC "${LIBSNOWPACK_ROOT}/share/doc/snowpack")
			SET(SNOWPACK_DOC ${DOC_PATH}/../snowpack)
		ENDIF(WIN32)
		IF(NOT IS_DIRECTORY "${SNOWPACK_INSTALLED_DOC}") #if we detected neighboring checked out version of Snowpack
			SET(SNOWPACK_INSTALLED_DOC "${LIBSNOWPACK_ROOT}/doc")
		ENDIF()
		INSTALL(DIRECTORY "${SNOWPACK_INSTALLED_DOC}/html" DESTINATION ${SNOWPACK_DOC} COMPONENT documentation)
		FILE(GLOB snowpack_docs "${SNOWPACK_INSTALLED_DOC}/*.*")
		INSTALL(FILES ${snowpack_docs} DESTINATION ${SNOWPACK_DOC} COMPONENT documentation)
	ENDIF(INSTALL_DOC)

	#bundle inishell
	FILE(GLOB inishell "${CMAKE_SOURCE_DIR}/tools/inish*.jar")
	IF("${inishell}" STREQUAL "")
		MESSAGE("Inishell not found! Please copy inishell.jar into the tools directory")
	ENDIF("${inishell}" STREQUAL "")
	INSTALL(FILES ${inishell} DESTINATION bin COMPONENT exe)
	
	#bundle the grid viewer
	FILE(GLOB grid_view "${CMAKE_SOURCE_DIR}/tools/Interface/view.*")
	IF("${grid_view}" STREQUAL "")
		MESSAGE("View not found! Please copy view.jar into the tools directory")
	ENDIF("${grid_view}" STREQUAL "")
	INSTALL(FILES ${grid_view} DESTINATION bin COMPONENT exe)

	#bundle niViz
	SET(CPACK_BUNDLE_NIVIZ ON CACHE BOOL "Should niViz also be bundled?")
	IF(CPACK_BUNDLE_NIVIZ)
		SET(NIVIZ_DIST_PATH "../niViz" CACHE PATH "Path to niViz 'dist' subdirectory")
		FILE(GLOB niviz_html "${NIVIZ_DIST_PATH}/*.html")
		IF("${niviz_html}" STREQUAL "")
			MESSAGE("niViz not found! Please provide a proper path to niViz 'dist' and don't forget to run 'make dist'")
		ENDIF("${niviz_html}" STREQUAL "")
		INSTALL(DIRECTORY "${NIVIZ_DIST_PATH}/assets" DESTINATION bin/niviz COMPONENT exe)
		INSTALL(FILES ${niviz_html} DESTINATION bin/niviz COMPONENT exe)
	ENDIF(CPACK_BUNDLE_NIVIZ)
ENDIF(CPACK_BUNDLE_ALL)

SET(CPACK_PACKAGE_DESCRIPTION "Alpine3D")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "This is the Alpine3D model distributed physical modeling of snow pack.")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/Welcome.txt")
SET(CPACK_PACKAGE_NAME "Alpine3D")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Alpine3D")
SET(CPACK_PACKAGE_CONTACT "Mathias Bavay <Mathias.Bavay@wsl.ch>")
SET(CPACK_PACKAGE_VENDOR "WSL Institute for Snow and Avalanche Research")
SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/alpine3d/icons\\\\icon.ico")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/License.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/Welcome.txt")
SET(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/Readme.txt")

SET(CPACK_PACKAGE_VERSION_MAJOR "${VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${VERSION_PATCH}")
SET(VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_ARCH "${CMAKE_SYSTEM_PROCESSOR}")
SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}${PACKAGE_EXTRA}-${VERSION}-${CPACK_ARCH}")

#installation by components for formats that support it
#SET(CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE 1)
SET(CPACK_ARCHIVE_COMPONENT_INSTALL ON)
SET(CPACK_COMPONENTS_ALL documentation examples exe)
SET(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "Alpine3D static library")
SET(CPACK_COMPONENT_LIBRARIES_DESCRIPTION "The Alpine3D static library")
SET(CPACK_COMPONENT_EXE_DISPLAY_NAME "Alpine3D application")
SET(CPACK_COMPONENT_EXE_DESCRIPTION "The Alpine3D standalone model")
SET(CPACK_COMPONENT_EXE_GROUP Runtime)
SET(CPACK_COMPONENT_RUNTIME_GROUP "Runtime")
SET(CPACK_COMPONENT_GROUP_RUNTIME_DESCRIPTION "What you need to run an application based on Alpine3D")
IF(INSTALL_DOC)
	SET(CPACK_COMPONENT_DOCUMENTATION_DISPLAY_NAME "Documentation")
	SET(CPACK_COMPONENT_DOCUMENTATION_DESCRIPTION "Documentation of the Alpine3D usage")
	SET(CPACK_COMPONENT_DOCUMENTATION_GROUP Docs)
	SET(CPACK_COMPONENT_EXAMPLES_DISPLAY_NAME "Examples")
	SET(CPACK_COMPONENT_EXAMPLES_DESCRIPTION "Various simulations with their example data")
	SET(CPACK_COMPONENT_EXAMPLES_GROUP Docs)
	#SET(CPACK_COMPONENT_EXAMPLES_DEPENDS libraries headers runtime)
	SET(CPACK_COMPONENT_DOCS_GROUP "Documentation")
	SET(CPACK_COMPONENT_GROUP_DOCS_DESCRIPTION "How to use Alpine3D")
ENDIF(INSTALL_DOC)

IF(UNIX)
	IF(APPLE)
		SET(CPACK_GENERATOR "DragNDrop")
		SET(CPACK_PACKAGING_INSTALL_PREFIX "/${CPACK_PACKAGE_NAME}")
		SET(CPACK_DMG_FORMAT "UDBZ")
		#SET(CPACK_GENERATOR "Bundle")
		SET(CPACK_BUNDLE_PLIST "${CMAKE_SOURCE_DIR}/tools/alpine3d.plist")
		SET(MACOSX_BUNDLE_INFO_STRING "${CPACK_PACKAGE_DESCRIPTION_SUMMARY}")
		SET(CPACK_BUNDLE_ICON "${CMAKE_SOURCE_DIR}/alpine3d/icons/icon.icns")
		SET(MACOSX_BUNDLE_ICON_FILE "${CPACK_BUNDLE_ICON}")
		SET(MACOSX_BUNDLE_GUI_IDENTIFIER "${CPACK_PACKAGE_DESCRIPTION}")
		SET(MACOSX_BUNDLE_LONG_VERSION_STRING "${VERSION}")
		SET(CPACK_BUNDLE_NAME "${CPACK_PACKAGE_NAME}")
		SET(MACOSX_BUNDLE_BUNDLE_NAME "${CPACK_PACKAGE_NAME}")
		SET(MACOSX_BUNDLE_SHORT_VERSION_STRING "${VERSION}")
		SET(MACOSX_BUNDLE_BUNDLE_VERSION "0")
		SET(MACOSX_BUNDLE_COPYRIGHT "${CPACK_PACKAGE_VENDOR}")
		#SET(CPACK_PACKAGING_INSTALL_PREFIX "/Applications/${CPACK_PACKAGE_NAME}")
		SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/alpine3d/icons/icon.icns")
		SET_SOURCE_FILES_PROPERTIES(${CPACK_BUNDLE_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
		INSTALL(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink bin/alpine3d /usr/bin/alpine3d)")
	ELSE(APPLE)
		SET(CPACK_ARCHIVE_COMPONENT_INSTALL OFF)
		SET(CPACK_GENERATOR "DEB;RPM;TGZ")
		SET(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
		SET(CPACK_DEBIAN_PACKAGE_NAME "alpine3d")
		SET(CPACK_DEBIAN_PACKAGE_SECTION "science")
		#objdump -p lib/libmeteoio.so | grep NEEDED | tr -s ' ' | cut -d' ' -f3 | tr '\n' ','
		IF(CPACK_BUNDLE_ALL OR BUILD_STATIC_LIBS)
			SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6, libc6")
		ELSE(CPACK_BUNDLE_ALL OR BUILD_STATIC_LIBS)
			SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++, meteoio,snowpack")
		ENDIF(CPACK_BUNDLE_ALL OR BUILD_STATIC_LIBS)
		#SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386") #dpkg --print-architecture
		IF("${CPACK_ARCH}" STREQUAL "x86_64")
			SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
		ELSE("${CPACK_ARCH}" STREQUAL "x86_64")
			SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "${CPACK_ARCH}")
		ENDIF("${CPACK_ARCH}" STREQUAL "x86_64")
		SET(CPACK_RPM_PACKAGE_NAME "alpine3d")
		SET(CPACK_RPM_PACKAGE_LICENSE "LPGLv3")
		SET(CPACK_RPM_PACKAGE_REQUIRES "libstdc++6, meteoio,snowpack")
		SET(CPACK_RPM_PACKAGE_ARCHITECTURE "${CPACK_ARCH}")
		SET(CPACK_RPM_PACKAGE_GROUP "Applications/Science")
	ENDIF(APPLE)
ELSE(UNIX)
	IF(WIN32)
		SET(CPACK_GENERATOR "NSIS")
		SET(CPACK_NSIS_MUI_ICON "${CPACK_PACKAGE_ICON}")
		SET(CPACK_PACKAGE_ICON "${CPACK_PACKAGE_ICON}")
		SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}")
		SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME}")
		#SET(CPACK_NSIS_EXTRA_INSTALL_COMMANDS "${EnvVarUpdate} $0 \"LIB\"  \"A\" \"HKLM\" \"C:\MyLib\"")
		SET(CPACK_NSIS_MODIFY_PATH ON)

		#Main links (they can get overwritten later)
		SET(MENU_LINKS_NIVIZ "http://niviz.org" "niViz online visualization" )
		SET(MENU_LINKS_ALPINE3D "https://models.slf.ch/p/alpine3d/" "Alpine3D Home page")
		IF(INSTALL_DOC)
			LIST(APPEND MENU_LINKS_ALPINE3D
			"${DOC_PATH}/Alpine3D_doc.html" "Alpine3D documentation"
			"${DOC_PATH}/examples" "Alpine3D examples"
			)
		ENDIF(INSTALL_DOC)

		#all other links: inishell, local niviz, meteoio
		IF(CPACK_BUNDLE_ALL)
			SET(MENU_LINKS_INISHELL "bin/inishell.jar" "Inishell")
			IF(CPACK_BUNDLE_NIVIZ)
				SET(MENU_LINKS_NIVIZ "bin/niviz/index.html" "NIVIZ visualization")
			ENDIF(CPACK_BUNDLE_NIVIZ)
			SET(MENU_LINKS_METEOIO "https://models.slf.ch/p/meteoio/" "MeteoIO Home page")
			IF(INSTALL_DOC)
				SET(MENU_LINKS_METEOIO
				"${METEOIO_DOC}/MeteoIO_doc.html" "MeteoIO documentation"
				"https://models.slf.ch/p/meteoio/" "MeteoIO Home page"
				)
			ENDIF(INSTALL_DOC)
		ENDIF(CPACK_BUNDLE_ALL)

		SET(CPACK_NSIS_MENU_LINKS ${MENU_LINKS_ALPINE3D} ${MENU_LINKS_INISHELL} ${MENU_LINKS_NIVIZ} ${MENU_LINKS_METEOIO})

		#info for the add/uninstall software administration
		SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_NAME}-${VERSION}")
		SET(CPACK_NSIS_INSTALLED_ICON_NAME "doc\\\\icon.ico")
		SET(CPACK_NSIS_URL_INFO_ABOUT "http://www.slf.ch")
		SET(CPACK_NSIS_HELP_LINK "https://models.slf.ch/p/alpine3d")
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
