# - Try to find PROJ 
# Once done this will define 
# 
#  PROJ_FOUND - system has PROJ 
#  PROJ_INCLUDE_DIRS - the PROJ include directory 
#  PROJ_LIBRARIES - Link these to use PROJ 
#  PROJ_DEFINITIONS - Compiler switches required for using PROJ 
# 

if (PROJ_LIBRARIES AND PROJ_INCLUDE_DIRS) 
  # in cache already 
  set(PROJ_FOUND TRUE) 
else (PROJ_LIBRARIES AND PROJ_INCLUDE_DIRS) 

  find_path(PROJ_INCLUDE_DIR 
    NAMES 
      proj.h 
    PATHS 
if(WIN32) 
      ${PROJ_DEV_PATH}/include 
endif(WIN32) 
        /usr/include 
        /usr/local/include 
        /opt/local/include 
        /sw/include 
        ${CMAKE_INSTALL_PREFIX}/include 
        ${CMAKE_SOURCE_DIR}/Win32/GDAL/include 
    PATH_SUFFIXES 
        proj 

  ) 
  #mark_as_advanced(PROJ_INCLUDE_DIR) 

  find_library(LIBPROJ_LIBRARY 
    NAMES 
        proj 
        libproj 
        proj_6_0 
        proj_6_1 
        proj_6_2 
    PATHS 
if(WIN32) 
      ${PROJ_DEV_PATH}/lib 
endif(WIN32) 
        /usr/lib 
        /usr/lib64 
        /usr/local/lib 
        /opt/local/lib 
        /sw/lib 
        ${CMAKE_INSTALL_PREFIX}/lib 
        ${CMAKE_SOURCE_DIR}/Win32/GDAL/lib 
  ) 
  #mark_as_advanced(LIBPROJ_LIBRARY) 

  if (LIBPROJ_LIBRARY) 
    set(LIBPROJ_FOUND TRUE) 
  else (LIBPROJ_LIBRARY) 
      message(FATAL_ERROR "Not found PROJ library") 
  endif (LIBPROJ_LIBRARY) 

  set(PROJ_INCLUDE_DIRS 
    ${PROJ_INCLUDE_DIR} 
  ) 

  if (LIBPROJ_FOUND) 
    set(PROJ_LIBRARIES 
      ${PROJ_LIBRARIES} 
      ${LIBPROJ_LIBRARY} 
    ) 
  endif (LIBPROJ_FOUND) 

  if (PROJ_INCLUDE_DIRS AND PROJ_LIBRARIES) 
     set(PROJ_FOUND TRUE) 
     set(PROJ_VERSION "6.0.0") 
  endif (PROJ_INCLUDE_DIRS AND PROJ_LIBRARIES) 

  if (PROJ_FOUND) 
    if (NOT PROJ_FIND_QUIETLY) 
      message(STATUS "Found PROJ: ${PROJ_LIBRARIES}") 
    endif (NOT PROJ_FIND_QUIETLY) 
  else (PROJ_FOUND) 
    if (PROJ_FIND_REQUIRED) 
      message(FATAL_ERROR "Could not find PROJ") 
    endif (PROJ_FIND_REQUIRED) 
  endif (PROJ_FOUND) 

  # show the PROJ_INCLUDE_DIRS and PROJ_LIBRARIES variables only in the advanced view 
  #mark_as_advanced(PROJ_INCLUDE_DIRS PROJ_LIBRARIES) 

endif (PROJ_LIBRARIES AND PROJ_INCLUDE_DIRS)
