#building version number in variable _versionString

MACRO (GETDATE TODAY)
	IF(CMAKE_VERSION VERSION_GREATER 2.8.11)
		STRING(TIMESTAMP TODAY "%Y%m%d")
	ELSE(CMAKE_VERSION VERSION_GREATER 2.8.11)
		IF (WIN32)
			EXECUTE_PROCESS(COMMAND "cmd" " /C date /T" OUTPUT_VARIABLE ${TODAY})
			STRING(REGEX REPLACE "(..)/(..)/(....) .*\n" "\\3\\2\\1" ${TODAY} ${${TODAY}}) #US format
			STRING(REGEX REPLACE "(..)-(..)-(....) .*\n" "\\3\\2\\1" ${TODAY} ${${TODAY}}) #UK format
			STRING(REGEX REPLACE "(..)\\.(..)\\.(....) .*\n" "\\3\\2\\1" ${TODAY} ${${TODAY}}) #CH format
		ELSEIF(UNIX)
			EXECUTE_PROCESS(COMMAND "date" "+%Y-%m-%d" OUTPUT_VARIABLE ${TODAY})
			string(REGEX REPLACE "(....)-(..)-(..).*" "\\1\\2\\3" ${TODAY} ${${TODAY}})
		ELSE (WIN32)
			MESSAGE(SEND_ERROR "date not implemented")
			SET(${TODAY} 000000)
		ENDIF (WIN32)
	ENDIF(CMAKE_VERSION VERSION_GREATER 2.8.11)
ENDMACRO (GETDATE)

MACRO(BuildVersionSVN)
	FIND_PACKAGE(Subversion)
	IF(Subversion_FOUND)
		SET(VERSION_FROM_SVN OFF CACHE BOOL "Retrieve software version from Subversion")
		IF(VERSION_FROM_SVN)
			Subversion_WC_INFO(${PROJECT_SOURCE_DIR} project) #HACK: if not an svn tree, it does not work
			GETDATE(TODAY)
			SET(_versionString "${TODAY}.${project_WC_REVISION}")
		ELSE(VERSION_FROM_SVN)
			SET(_versionString "${VERSION_MAJOR}.${VERSION_MINOR}${VERSION_PATCH}")
		ENDIF(VERSION_FROM_SVN)
	ELSE(Subversion_FOUND)
		SET(_versionString "${VERSION_MAJOR}.${VERSION_MINOR}${VERSION_PATCH}")
	ENDIF(Subversion_FOUND)
ENDMACRO(BuildVersionSVN)

MACRO(BuildVersionGIT)
	FIND_PACKAGE(Git)
	IF(GIT_FOUND)
		SET(VERSION_FROM_GIT OFF CACHE BOOL "Retrieve software version from Git")
		IF(VERSION_FROM_GIT)
			execute_process(
				COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
				WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
				OUTPUT_VARIABLE project_WC_REVISION
				ERROR_QUIET
				OUTPUT_STRIP_TRAILING_WHITESPACE
			)
			GETDATE(TODAY)
			SET(_versionString "${TODAY}.${project_WC_REVISION}")
		ELSE(VERSION_FROM_GIT)
			SET(_versionString "${VERSION_MAJOR}.${VERSION_MINOR}${VERSION_PATCH}")
		ENDIF(VERSION_FROM_GIT)
	ELSE(GIT_FOUND)
		SET(_versionString "${VERSION_MAJOR}.${VERSION_MINOR}${VERSION_PATCH}")
	ENDIF(GIT_FOUND)
ENDMACRO(BuildVersionGIT)
