# Install script for directory: /media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack

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

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/snowpack" TYPE FILE FILES
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/Constants.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/DataClasses.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/Hazard.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/Laws_sn.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/MainPage.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/Meteo.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/Saltation.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/SnowDrift.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/SnowpackConfig.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/Stability.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/StabilityAlgorithms.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/Utils.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/libsnowpack.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/vanGenuchten.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/snowpack/plugins" TYPE FILE FILES
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/plugins/AsciiIO.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/plugins/CaaMLIO.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/plugins/ImisDBIO.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/plugins/SmetIO.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/plugins/SnowpackIO.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/plugins/SnowpackIOInterface.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/snowpack/snowpackCore" TYPE FILE FILES
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/snowpackCore/Aggregate.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/snowpackCore/Canopy.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/snowpackCore/Metamorphism.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/snowpackCore/PhaseChange.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/snowpackCore/ReSolver1d.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/snowpackCore/SeaIce.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/snowpackCore/Snowpack.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/snowpackCore/Solver.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/snowpackCore/VapourTransport.h"
    "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/snowpackCore/WaterTransport.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/snowpack/cmake_install.cmake")
  include("/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/applications/snowpack/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/snowpack/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
