#SPDX-License-Identifier: LGPL-3.0-or-later
#Set different variables according to the detected compiler and processor
#based on $CMAKE_CXX_COMPILER_ID it sets the following variables:
# WARNINGS, EXTRA_WARNINGS, EXTRA, OPTIM, ARCH, DEBUG, _VERSION, PROFILING
# It can also edit CMAKE_SHARED_LINKER_FLAGS and CMAKE_EXE_LINKER_FLAGS

INCLUDE("${CMAKE_SOURCE_DIR}/tools/cmake/BuildVersion.cmake")
BuildVersionGIT()

#TODO: replace all defs such as /D__DEBUG or -DDEBUG_ARITHM by add_compile_definitions() once moving to cmake >=3.12 (released 11.2018)

MACRO (SET_COMPILER_OPTIONS)
	SET(USER_COMPILER_OPTIONS "" CACHE STRING "Provide some extra compiler options")
	MARK_AS_ADVANCED(FORCE USER_COMPILER_OPTIONS)
	SET(EXTRA "${EXTRA} ${USER_COMPILER_OPTIONS}")

	IF(NOT (CMAKE_CXX_COMPILER_ID MATCHES "MSVC"))
		#we consider that all other compilers support "-" options and silently ignore what they don't know
		IF(ENABLE_LAPACK)
			SET(EXTRA "${EXTRA} -DCLAPACK")
		ENDIF(ENABLE_LAPACK)

		IF(DEBUG_ARITHM)
			SET(EXTRA "${EXTRA} -DDEBUG_ARITHM")
		ENDIF(DEBUG_ARITHM)
	ENDIF()

	###########################################################
	IF(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
		SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)	#this is required for building libraries
		SET(EXTRA "${EXTRA} /D_USE_MATH_DEFINES")	#USE_MATH_DEFINES needed for VC++
		IF(DEBUG_ARITHM)
			SET(EXTRA "${EXTRA} /EHa")
		ENDIF(DEBUG_ARITHM)
		
		#SET(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited configs" FORCE)
		SET(WARNINGS "/W4 /D_CRT_SECURE_NO_WARNINGS /EHsc") #Za: strict ansi EHsc: handle c++ exceptions /w35045: inform about Spectre mitigation
		#SET(EXTRA_WARNINGS "/Wp64") #/Wall
		SET(WARNINGS "${WARNINGS} /experimental:external /external:I c:/Windows /external:W0")
		SET(OPTIM "/O2 /DNDEBUG /DEBUG:FASTLINK /MD /DNOSAFECHECKS")
		SET(ARCH_OPTIM "/arch:AVX2")
		SET(ARCH_SAFE "")
		IF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64")
			SET(ARCH_SAFE "/arch:SSE2")
		ENDIF()
		SET(DEBUG "/Z7 /Od /D__DEBUG /MDd")
		SET(_VERSION "/D_VERSION=${_versionString}")
		
	###########################################################
	ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL Intel)
		SET(WARNINGS_OFF "-Wno-long-long -Wno-unknown-pragmas -wd2015,11071")
		SET(WARNINGS "-Wall -Wswitch ${WARNINGS_OFF}")
		SET(DEEP_WARNINGS "-Wshadow -Wpointer-arith -Wconversion -Winline -Wdisabled-optimization") #-Wfloat-equal -Wpadded
		SET(EXTRA_WARNINGS "-Wextra -pedantic ${DEEP_WARNINGS}")
		SET(OPTIM "-g -O3 -DNDEBUG -DNOSAFECHECKS")
		IF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64")
			SET(ARCH_SAFE "-march=nehalem -mtune=skylake")
			SET(ARCH_OPTIM "-march=native -mtune=native")
		ENDIF()
		SET(DEBUG "-g3 -O0 -D__DEBUG")
		SET(_VERSION "-D_VERSION=${_versionString}")
		
	###########################################################
	ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL Cray)
		SET(WARNINGS "-hlist=m -h negmsgs -h msglevel_3 -h nomessage=870") #870: accept multibyte chars
		#SET(EXTRA_WARNINGS "-h msglevel_2")
		SET(OPTIM "-O3 -hfp3 -h msglevel_4 -DNDEBUG -DNOSAFECHECKS")
		IF($ENV{CRAY_CPU_TARGET} MATCHES "^$")
			IF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64")
				SET(ARCH_SAFE "-h cpu=x86-64")
				MESSAGE("No CRAY_CPU_TARGET set, setting it to x86-64; please consider loading the proper target module.")
			ELSE()
				MESSAGE("No CRAY_CPU_TARGET set; please consider loading the proper target module.")
			ENDIF()
		ENDIF()
		SET(DEBUG "-g -D__DEBUG")
		SET(_VERSION "-D_VERSION=${_versionString}")
	
	###########################################################
	ELSEIF(CMAKE_CXX_COMPILER_ID MATCHES "^GNU$")
		IF(WIN32)
			LIST(APPEND CFLAGS " -D_USE_MATH_DEFINES") #USE_MATH_DEFINES needed for Win32
		ENDIF(WIN32)
		
		SET(WARNINGS "-Wall -Wno-long-long -Wswitch -Wno-unknown-pragmas")
		SET(DEEP_WARNINGS "-Wunused-value -Wshadow -Wpointer-arith -Winline -Wdisabled-optimization -Wctor-dtor-privacy") #-Wfloat-equal -Wpadded -Wconversion
		SET(EXTRA_WARNINGS "-Wextra -pedantic ${DEEP_WARNINGS}") #-Weffc++
		SET(OPTIM "-g -O3 -DNDEBUG -DNOSAFECHECKS")
		IF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64")
			IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11.0)
				SET(ARCH_SAFE "-march=nehalem -mtune=skylake")
			ELSE()
				SET(ARCH_SAFE "-march=x86-64-v2 -mtune=core-avx2")
			ENDIF()
		ENDIF()
		SET(DEBUG "-g3 -O0 -D__DEBUG")
		SET(_VERSION "-D_VERSION=${_versionString}")
		
		SET(PROFILING "-pg -fprofile-arcs") #add ${PROFILING} to the CFLAGS when necessary
		SET(EXTRA_WARNINGS "${EXTRA_WARNINGS} -Wunsafe-loop-optimizations -Wwrite-strings")
		IF(NOT ANDROID)
			SET(EXTRA_WARNINGS "${EXTRA_WARNINGS} -ansi")
			IF(WIN32) #for gcc on windows
				SET(CMAKE_SHARED_LINKER_FLAGS "--enable-auto-import")
				SET(CMAKE_EXE_LINKER_FLAGS "--enable-auto-import")
			ENDIF(WIN32)
		ENDIF(NOT ANDROID)
		SET(ARCH_OPTIM "-march=native -mtune=native")
		SET(EXTRA_WARNINGS "${EXTRA_WARNINGS} -Wvector-operation-performance")
		IF(NOT WIN32)
			FIND_PROGRAM(CMAKE_GCC_AR NAMES ${_CMAKE_TOOLCHAIN_PREFIX}gcc-ar${_CMAKE_TOOLCHAIN_SUFFIX} HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
			FIND_PROGRAM(CMAKE_GCC_NM NAMES ${_CMAKE_TOOLCHAIN_PREFIX}gcc-nm HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
			FIND_PROGRAM(CMAKE_GCC_RANLIB NAMES ${_CMAKE_TOOLCHAIN_PREFIX}gcc-ranlib HINTS ${_CMAKE_TOOLCHAIN_LOCATION})

			IF(CMAKE_GCC_AR AND CMAKE_GCC_NM AND CMAKE_GCC_RANLIB)
				#for cmake>3.9: set CMAKE_INTERPROCEDURAL_OPTIMIZATION to TRUE
				SET(USE_LTO_OPTIMIZATIONS ON CACHE BOOL "Use Link Time Optmizations when compiling (memory heavy while compiling)")
				MARK_AS_ADVANCED(FORCE USE_LTO_OPTIMIZATIONS)
				IF(USE_LTO_OPTIMIZATIONS)
					SET(OPTIM "${OPTIM} -flto=auto -fno-fat-lto-objects")
				ENDIF()
				SET( CMAKE_AR "${CMAKE_GCC_AR}" )
				SET( CMAKE_NM "${CMAKE_GCC_NM}" )
				SET( CMAKE_RANLIB "${CMAKE_GCC_RANLIB}" )
			ELSE()
				MESSAGE( WARNING "GCC indicates LTO support, but binutils wrappers could not be found. Disabling LTO." )
			ENDIF()
		ENDIF(NOT WIN32)
		#if set to ON, all binaries depending on the library have to be compiled the same way.
		#Then, do an "export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.4" and run with "ASAN_OPTIONS=symbolize=1"
		SET(LEAKS_CHECK OFF CACHE BOOL "Set to ON to dynamically check for memory corruption (and do the same for applications linked with MeteoIO)")
		IF (LEAKS_CHECK)
			SET(EXTRA "${EXTRA} -fsanitize=address -fno-omit-frame-pointer")
		ENDIF(LEAKS_CHECK)

	###########################################################
	ELSEIF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
		IF(WIN32)
			LIST(APPEND CFLAGS " -D_USE_MATH_DEFINES") #USE_MATH_DEFINES needed for Win32
		ENDIF(WIN32)
		
		SET(WARNINGS_OFF "-Wno-long-long -Wno-float-equal -Wno-documentation -Wno-documentation-unknown-command -Wno-old-style-cast -Wno-padded -Wno-missing-noreturn -Wno-weak-vtables -Wno-switch-enum -Wno-covered-switch-default -Wno-global-constructors -Wno-exit-time-destructors -Wno-unknown-pragmas -Wno-format-nonliteral -Wno-date-time -Wno-unused-template -Wno-disabled-macro-expansion -Wno-c++98-compat -Wno-c++98-compat-pedantic")
		SET(WARNINGS "-Wall -Wswitch -Weverything ${WARNINGS_OFF}") #obviously, we should try to fix the warnings! Keeping in mind that some of these W are half buggy...
		SET(DEEP_WARNINGS "-Wunused-value -Wshadow -Wpointer-arith -Wconversion -Winline -Wdisabled-optimization -Wctor-dtor-privacy") #-Rpass=.* for static analysis
		SET(EXTRA_WARNINGS "-Wextra -pedantic -Weffc++ ${DEEP_WARNINGS}")
		SET(OPTIM "-g -O3 -DNDEBUG -DNOSAFECHECKS -flto")
		IF(APPLE)
			OPTION(BUILD_FAT_BINARIES "Compile fat binaries, for x86_64 and arm64" OFF)
			IF(BUILD_FAT_BINARIES)
				SET(CMAKE_OSX_ARCHITECTURES "x86_64;arm64")
			ENDIF()
		ELSE(APPLE)
			IF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64")
				IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.0)
					SET(ARCH_SAFE "-march=nehalem -mtune=skylake")
				ELSE()
					SET(ARCH_SAFE "-march=x86-64-v2 -mtune=core-avx2")
				ENDIF()
			ENDIF()
		ENDIF(APPLE)
		SET(DEBUG "-g3 -O0 -D__DEBUG")
		SET(_VERSION "-D_VERSION=${_versionString}")
		
		SET(PROFILING "-pg") #add ${PROFILING} to the CFLAGS when necessary
		SET(EXTRA "${EXTRA} -fcolor-diagnostics") #-fapple-pragma-pack does not seems necessary; -ftrapv should be replaced by sanitize=integer
		SET(LEAKS_CHECK OFF CACHE BOOL "Set to ON to dynamically check for memory corruption (and do the same for applications linked with our software)")
			IF (LEAKS_CHECK)
				SET(EXTRA "${EXTRA} -ftrapv -fno-omit-frame-pointer") #-fsanitize=address,undefined,integer,undefined-trap but this is currently not supported by Apple
			ENDIF(LEAKS_CHECK)

		#Considering that on all supported platforms, CLang has some versions that might no accept the march=native flag
		#An alternative solution would be to consider that only Apple's Clang 12 does not support march=native
		IF(CMAKE_SYSTEM_PROCESSOR MATCHES "ARM64")
			include(CheckCXXCompilerFlag)
			CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_NATIVE)
			IF(COMPILER_SUPPORTS_NATIVE)
				SET(ARCH_OPTIM "-march=native -mtune=native")
			ELSE(COMPILER_SUPPORTS_NATIVE)
				SET(ARCH_OPTIM ${ARCH_SAFE})
			ENDIF(COMPILER_SUPPORTS_NATIVE)
		ELSE()
			SET(ARCH_OPTIM "-march=native -mtune=native")
		ENDIF()
	ENDIF()
	
	###########################################################
	#targets providing SETs of compiler options
	IF(NOT DEST)
		SET(DEST "safe" CACHE STRING "Choose safe or optimized" FORCE)
	ENDIF(NOT DEST)

	IF (DEST STREQUAL "safe")
		SET(ARCH "${ARCH_SAFE}")
	ENDIF(DEST STREQUAL "safe")

	IF(DEST STREQUAL "optimized")
		SET(ARCH "${ARCH_OPTIM}")
	ENDIF(DEST STREQUAL "optimized")

ENDMACRO (SET_COMPILER_OPTIONS)
