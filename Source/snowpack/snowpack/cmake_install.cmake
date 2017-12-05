# Install script for directory: /media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack

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

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "exe" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsnowpack.so.3.4.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsnowpack.so.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsnowpack.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "/home/nander/usr/lib:/usr/lib/x86_64-linux-gnu/liblapack.so:/usr/lib/x86_64-linux-gnu/libblas.so")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/lib/libsnowpack.so.3.4.1"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/lib/libsnowpack.so.3"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/lib/libsnowpack.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsnowpack.so.3.4.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsnowpack.so.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsnowpack.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/usr/lib/x86_64-linux-gnu/liblapack.so:/usr/lib/x86_64-linux-gnu/libblas.so:/home/nander/usr/lib:"
           NEW_RPATH "/home/nander/usr/lib:/usr/lib/x86_64-linux-gnu/liblapack.so:/usr/lib/x86_64-linux-gnu/libblas.so")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

