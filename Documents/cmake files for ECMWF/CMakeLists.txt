cmake_minimum_required(VERSION 2.8)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/tools/cmake/")

# Version of library
SET(VERSION_MAJOR 3 CACHE INTERNAL "")
SET(VERSION_MINOR 3 CACHE INTERNAL "")
SET(VERSION_PATCH 0 CACHE INTERNAL "")

SET(LIBPREFIX "lib") #make sure ALL plugins are libXXX for ALL plateforms
SET(CMAKE_IMPORT_LIBRARY_PREFIX "lib") #for Visual C++ import libraries

PROJECT(snowpack)
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
#plugins and applications
SET(SNOWPACK_STANDALONE ON CACHE BOOL "Compilation SNOWPACK standalone ON or OFF")
SET(PROFEVAL OFF CACHE BOOL "Compilation PROFEVAL ON or OFF")
SET(INSTALL_DOC OFF CACHE BOOL "Install documentation ON or OFF")

SET(ENABLE_LAPACK OFF CACHE BOOL "Compile with the CLAPACK library?")
SET(PLUGIN_IMISIO OFF CACHE BOOL "Compilation IMISDBIO ON or OFF - only relevant for SLF")
SET(PLUGIN_CAAMLIO OFF CACHE BOOL "Compilation CAAMLIO ON or OFF to read CAAML profiles")

###########################################################
#platform specific. First, compiler specific flags, then system specific flags
SET(BUILD_SHARED_LIBS ON CACHE BOOL "Set to ON to build dynamic libraries")
IF(WIN32)
	SET(BUILD_STATIC_LIBS ON) #otherwise the lib files are not generated
ELSE(WIN32)
	SET(BUILD_STATIC_LIBS OFF CACHE BOOL "Set to ON to build static libraries")
ENDIF(WIN32)

INCLUDE("${CMAKE_MODULE_PATH}/BuildVersion.cmake")
BuildVersion()

IF(MSVC)
	IF(PLUGIN_IMISIO)
		SET(EXTRA "${EXTRA} /DIMISDBIO")
	ENDIF(PLUGIN_IMISIO)
	IF(PLUGIN_CAAMLIO)
		SET(EXTRA "${EXTRA} /DCAAMLIO")
	ENDIF(PLUGIN_CAAMLIO)
	#SET(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited configs"  FORCE)
	SET(WARNINGS "/W4")
	#SET(EXTRA_WARNINGS "/Wp64") #/Wall
	SET(OPTIM "/O2 /DNDEBUG /MD /DNOSAFECHECKS")
	SET(ARCH_OPTIM "/arch:SSE2")
	SET(ARCH_SAFE "")
	SET(DEBUG "/Z7 /Od /D__DEBUG /MDd")
	IF(DEBUG_ARITHM)
		SET(EXTRA "${EXTRA} /EHa")
	ENDIF(DEBUG_ARITHM)
	LIST(APPEND CFLAGS " /D_CRT_SECURE_NO_WARNINGS /D_USE_MATH_DEFINES /EHsc") #Za: strict ansi
	SET(_VERSION "/D_VERSION=\\\"${_versionString}\\\"")
ELSE(MSVC)
	#we consider that all other compilers support "-" options and silently ignore what they don't know
	IF(PLUGIN_IMISIO)
		SET(EXTRA "${EXTRA} -DIMISDBIO")
	ENDIF(PLUGIN_IMISIO)
	IF(PLUGIN_CAAMLIO)
		SET(EXTRA "${EXTRA} -DCAAMLIO")
	ENDIF(PLUGIN_CAAMLIO)
	IF(ENABLE_LAPACK)
		SET(EXTRA "${EXTRA} -DCLAPACK")
	ENDIF(ENABLE_LAPACK)
	SET(WARNINGS "-Wall")
	SET(DEEP_WARNINGS "-Wunused-value -Wshadow -Wpointer-arith -Wconversion -Winline -Wdisabled-optimization") #-Wfloat-equal -Wpadded
	SET(EXTRA_WARNINGS "-Wextra -pedantic -Weffc++ ${DEEP_WARNINGS}") #-Weffc++
	IF(WIN32)
		LIST(APPEND CFLAGS " -D_USE_MATH_DEFINES") #USE_MATH_DEFINES needed for Win32
	ENDIF(WIN32)
	SET(OPTIM "-g -O3 -DNDEBUG -DNOSAFECHECKS")
	IF("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "x86_64" OR "${CMAKE_SYSTEM_PROCESSOR}" MATCHES "AMD64")
		SET(ARCH_SAFE "")
	ENDIF()
	SET(DEBUG "-g3 -O0 -D__DEBUG")
	IF(DEBUG_ARITHM)
		SET(EXTRA "${EXTRA} -DDEBUG_ARITHM")
	ENDIF(DEBUG_ARITHM)
	SET(_VERSION "-D_VERSION=\\\"${_versionString}\\\"")
	
	IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		SET(PROFILING "-pg -fprofile-arcs") #add ${PROFILING} to the CFLAGS when necessary
		SET(EXTRA_WARNINGS "${EXTRA_WARNINGS} -Wunsafe-loop-optimizations -Wvector-operation-performance")
		IF(NOT ANDROID)
			SET(EXTRA_WARNINGS "${EXTRA_WARNINGS} -ansi")
			IF(WIN32) #for gcc on windows
				SET(CMAKE_SHARED_LINKER_FLAGS "--enable-auto-import")
				SET(CMAKE_EXE_LINKER_FLAGS "--enable-auto-import")
			ENDIF(WIN32)
		ENDIF(NOT ANDROID)
		EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
		IF(GCC_VERSION VERSION_GREATER 4.2 OR GCC_VERSION VERSION_EQUAL 4.2)
			SET(ARCH_OPTIM  "-march=native -mtune=native")
		ENDIF()
		IF(GCC_VERSION VERSION_GREATER 4.8 OR GCC_VERSION VERSION_EQUAL 4.8)
			SET(OPTIM "${OPTIM} -flto") #for gcc>4.5, but first implementations were slow, so it is safe to enforce 4.8
			# if set to ON, all binaries depending on the library have to be compiled the same way.
			#Then, do an "export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.4" and run with "ASAN_OPTIONS=symbolize=1"
			SET(LEAKS_CHECK OFF CACHE BOOL "Set to ON to dynamically check for memory corruption (and do the same for applications linked with MeteoIO)")
			IF (LEAKS_CHECK)
				SET(EXTRA "${EXTRA} -fsanitize=address -fno-omit-frame-pointer")
			ENDIF(LEAKS_CHECK)
		ENDIF()
	ENDIF()
	IF("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
		SET(PROFILING "-pg") #add ${PROFILING} to the CFLAGS when necessary
		SET(EXTRA "${EXTRA} -fcolor-diagnostics") #-fapple-pragma-pack does not seems necessary; -ftrapv should be replaced by sanitize=integer
		SET(OPTIM "${OPTIM} -flto")
		SET(WARNINGS "${WARNINGS} -Weverything -Wno-date-time -Wno-float-equal -Wno-documentation -Wno-documentation-unknown-command -Wno-old-style-cast -Wno-padded -Wno-missing-noreturn -Wno-weak-vtables -Wno-global-constructors -Wno-exit-time-destructors -Wno-sign-conversion -Wno-covered-switch-default") #obviously, we should try to fix the warnings!
		SET(DEEP_WARNINGS "${DEEP_WARNINGS} -Rpass=.*")
		SET(LEAKS_CHECK OFF CACHE BOOL "Set to ON to dynamically check for memory corruption (and do the same for applications linked with MeteoIO)")
			IF (LEAKS_CHECK)
				SET(EXTRA "${EXTRA} -ftrapv -fno-omit-frame-pointer") #-fsanitize=address,undefined,integer,undefined-trap but this is currently not supported by Apple
			ENDIF(LEAKS_CHECK)
		SET(ARCH_OPTIM  "-march=native")
	ENDIF()

ENDIF(MSVC)


#The following pre-processor variables are automatically defined:
# __linux, __GNUC__, __CYGWIN__, _WIN32, __MINGW32__, __APPLE__
IF(UNIX)
	IF(APPLE)
		SET(DLL_EXT "dylib")
		SET(STAT_EXT "a")
		#set(CMAKE_OSX_SYSROOT iphoneos2.2.1) #for IOS
		#set(CMAKE_OSX_ARCHITECTURES $(ARCHS_STANDARD_32_BIT)) #for IOS

		#set a full rpath for both build & install tree
		SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
		SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
		SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

		# add the automatically determined parts of the RPATH
		# which point to directories outside the build tree to the install RPATH
		SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

		# the RPATH to be used when installing, but only if it's not a system directory
		LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
		IF("${isSystemDir}" STREQUAL "-1")
			SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
		ENDIF("${isSystemDir}" STREQUAL "-1")
	ELSE(APPLE)
		SET(DLL_EXT "so*")
		SET(STAT_EXT "a")
	ENDIF(APPLE)
ELSE(UNIX)
	IF(WIN32)
		#SET(PLATFORM "-D_WIN32 -U__STRICT_ANSI__")
		ADD_DEFINITIONS(-U__STRICT_ANSI__)
		SET(DLL_EXT "dll")
		SET(STAT_EXT "a")
	ENDIF(WIN32)
ENDIF(UNIX)

IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE "release" CACHE STRING "Choose the type of build, options are: debug release " FORCE)
	SET(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE})
ENDIF(NOT CMAKE_BUILD_TYPE)

###########################################################
#targets providing sets of compiler options
if(NOT DEST)
  set(DEST "safe" CACHE STRING "Choose safe or optimized" FORCE)
endif(NOT DEST)

IF (DEST STREQUAL "safe")
	SET(ARCH  "${ARCH_SAFE}")
ENDIF(DEST STREQUAL "safe")

IF(DEST STREQUAL "optimized")
	SET(ARCH  "${ARCH_OPTIM}")
ENDIF(DEST STREQUAL "optimized")

###########################################################
#finally, SET compile flags
SET(CMAKE_CXX_FLAGS "${_VERSION} ${ARCH} -fPIC ${CFLAGS} ${EXTRA}" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_RELEASE "${OPTIM}" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_DEBUG "${DEBUG} ${WARNINGS} ${EXTRA_WARNINGS}" CACHE STRING "" FORCE)

MESSAGE(STATUS "DEBUG CXXFLAGS  : ${CMAKE_CXX_FLAGS_DEBUG}")
MESSAGE(STATUS "RELEASE CXXFLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

###########################################################
#for the install target
FILE(GLOB hfiles "snowpack/*.h")
FILE(GLOB hpluginfiles "snowpack/plugins/*.h")
FILE(GLOB hcorefiles "snowpack/snowpackCore/*.h")
INSTALL(FILES ${hfiles} DESTINATION include/snowpack COMPONENT headers)
INSTALL(FILES ${hpluginfiles} DESTINATION include/snowpack/plugins COMPONENT headers)
INSTALL(FILES ${hcorefiles} DESTINATION include/snowpack/snowpackCore COMPONENT headers)

ADD_CUSTOM_TARGET(distclean make clean
	COMMAND cmake -E remove *~ snowpack/*~ applications/snowpack/*~
	COMMAND cmake -E remove CMakeCache.txt cmake_install.cmake snowpack/cmake_install.cmake applications/snowpack/cmake_install.cmake tests/cmake_install.cmake tests/CTestTestfile.cmake
	COMMAND cmake -E remove_directory CMakeFiles
	COMMAND cmake -E remove_directory snowpack/CMakeFiles
	COMMAND cmake -E remove_directory applications/snowpack/CMakeFiles
	COMMAND cmake -E remove_directory tests/CMakeFiles
	COMMAND cmake -E remove tests/res1exp/*.cmake tests/res1exp/output/*
	COMMAND cmake -E remove_directory tests/res1exp/CMakeFiles
	COMMAND cmake -E remove_directory tests/res1exp/Testing
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
IF (DOXYGEN_FOUND)
	ADD_CUSTOM_TARGET(doc "doxygen" ${PROJECT_SOURCE_DIR}/snowpack/config.dox)
ENDIF (DOXYGEN_FOUND)
#Install the documentation if the user desires it
IF(INSTALL_DOC)
	IF(WIN32)
		SET(DOC_PATH "doc")
	ELSE(WIN32)
		SET(DOC_PATH "share/doc/${CMAKE_PROJECT_NAME}")
	ENDIF(WIN32)
	SET(DOC_FILES License.txt Readme.txt snowpack/icons/icon.ico)
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
	IF(INSTALL_DOC)
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

	#bundle inishell
	FILE(GLOB inishell "${CMAKE_SOURCE_DIR}/tools/inish*.jar")
	IF("${inishell}" STREQUAL "")
		MESSAGE("Inishell not found! Please copy inishell.jar into the tools directory")
	ENDIF("${inishell}" STREQUAL "")
	INSTALL(FILES ${inishell} DESTINATION bin COMPONENT exe)
ENDIF(CPACK_BUNDLE_ALL)

SET(CPACK_PACKAGE_DESCRIPTION "Snowpack Library")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "This is the SNOWPACK model ([Bartelt 2002], [Lehning 2002a], [Lehning 2002b]) for physical modeling of snow pack. This consists of a library and a program calling it so it can easily be integrated into other models.")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/Welcome.txt")
SET(CPACK_PACKAGE_NAME "Snowpack")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Snowpack")
SET(CPACK_PACKAGE_CONTACT "Mathias Bavay <Mathias.Bavay@wsl.ch>")
SET(CPACK_PACKAGE_VENDOR "WSL Institute for Snow and Avalanche Research")
SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/snowpack/icons\\\\icon.ico")
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
SET(CPACK_COMPONENTS_ALL documentation examples exe)
SET(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "Snowpack static library")
SET(CPACK_COMPONENT_LIBRARIES_DESCRIPTION "The SNOWPACK static library")
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
		SET(CPACK_GENERATOR "DEB;RPM;TGZ")
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
		SET(CPACK_RPM_PACKAGE_REQUIRES "libstdc++6, meteoio")
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
		IF(INSTALL_DOC)
			IF(CPACK_BUNDLE_ALL)
				SET(CPACK_NSIS_MENU_LINKS
				"https://models.slf.ch/snopviz/" "SNOPVIZ visualization"
				"bin/inishell.jar" "SNOWPACK inishell"
				"${DOC_PATH}/Snowpack_doc.html" "SNOWPACK documentation"
				"${DOC_PATH}/examples" "SNOWPACK examples"
				"https://models.slf.ch/p/snowpack/" "SNOWPACK Home page"
				"${METEOIO_DOC}/MeteoIO_doc.html" "MeteoIO documentation"
				"https://models.slf.ch/p/meteoio/" "MeteoIO Home page")
			ELSE(CPACK_BUNDLE_ALL)
				SET(CPACK_NSIS_MENU_LINKS
				"https://models.slf.ch/snopviz/" "SNOPVIZ visualization"
				"${DOC_PATH}/Snowpack_doc.html" "SNOWPACK documentation"
				"${DOC_PATH}/examples" "SNOWPACK examples"
				"https://models.slf.ch/p/snowpack/" "SNOWPACK Home page")
			ENDIF(CPACK_BUNDLE_ALL)
		ELSEIF(INSTALL_DOC)
			IF(CPACK_BUNDLE_ALL)
				SET(CPACK_NSIS_MENU_LINKS
				"https://models.slf.ch/snopviz/" "SNOPVIZ visualization"
				"bin/inishell.jar" "SNOWPACK inishell"
				"https://models.slf.ch/p/snowpack/" "SNOWPACK Home page"
				"https://models.slf.ch/p/meteoio/" "MeteoIO Home page")
			ELSE(CPACK_BUNDLE_ALL)
				SET(CPACK_NSIS_MENU_LINKS 
				"https://models.slf.ch/snopviz/" "SNOPVIZ visualization"
				"https://models.slf.ch/p/snowpack/" "SNOWPACK Home page")
			ENDIF(CPACK_BUNDLE_ALL)
		ENDIF(INSTALL_DOC)

		#info for the add/uninstall software administration
		SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_NAME}-${VERSION}")
		SET(CPACK_NSIS_INSTALLED_ICON_NAME "doc\\\\icon.ico")
		SET(CPACK_NSIS_URL_INFO_ABOUT "http://www.slf.ch")
		SET(CPACK_NSIS_HELP_LINK "https://models.slf.ch/p/snowpack")
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
## Special build type configuration for Coverage Testing.
# Note: Coverage testing only works with g++ (or some commercial tools)
OPTION(BUILD_TESTING "Build the QA tests" OFF)
IF(BUILD_TESTING)
	# Coverage statistic rign now only supportet with GCC (GNUGXX)
	IF(CMAKE_COMPILER_IS_GNUCXX)
		OPTION(BUILD_TESTING_WITH_COVERAGE "Use gcc flags to make covarage tests" OFF)
		IF(BUILD_TESTING_WITH_COVERAGE)
			# force the correct flags to make coverage trought a build type
			SET(CMAKE_BUILD_TYPE coverage)
			SET(CMAKE_CXX_FLAGS_COVERAGE "-Wall -O0 -ftest-coverage -fprofile-arcs" CACHE STRING "" FORCE)
			SET(CMAKE_C_FLAGS_COVERAGE "-Wall -O0 -ftest-coverage -fprofile-arcs" CACHE STRING "" FORCE)
			MARK_AS_ADVANCED(FORCE CMAKE_CXX_FLAGS_COVERAGE CMAKE_C_FLAGS_COVERAGE)
			# Force debug mode, so that no optimisation is made to get "Correct" results of coverage testing
			SET(CMAKE_BUILD_TYPE "coverage" CACHE STRING "Choose the type of build, options are: debug release coverage " FORCE)
			SET(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE})
		ENDIF(BUILD_TESTING_WITH_COVERAGE)
	ENDIF(CMAKE_COMPILER_IS_GNUCXX)

	INCLUDE(CTest) # This makes ENABLE_TESTING() and gives support for Dashboard
	ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTING)
