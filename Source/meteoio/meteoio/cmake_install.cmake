# Install script for directory: /media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/meteoio/meteoio

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/nander/usr")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "libraries" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmeteoio.so.2.7.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmeteoio.so.2"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmeteoio.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "/home/nander/usr/lib")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/meteoio/lib/libmeteoio.so.2.7.0"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/meteoio/lib/libmeteoio.so.2"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/meteoio/lib/libmeteoio.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmeteoio.so.2.7.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmeteoio.so.2"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmeteoio.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "::::::::::::::::::::"
           NEW_RPATH "/home/nander/usr/lib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

