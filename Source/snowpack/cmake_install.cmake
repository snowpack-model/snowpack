# Install script for directory: /home/nander/snowpack_greenland_new

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
    "/home/nander/snowpack_greenland_new/snowpack/Constants.h"
    "/home/nander/snowpack_greenland_new/snowpack/DataClasses.h"
    "/home/nander/snowpack_greenland_new/snowpack/Hazard.h"
    "/home/nander/snowpack_greenland_new/snowpack/Laws_sn.h"
    "/home/nander/snowpack_greenland_new/snowpack/MainPage.h"
    "/home/nander/snowpack_greenland_new/snowpack/Meteo.h"
    "/home/nander/snowpack_greenland_new/snowpack/Saltation.h"
    "/home/nander/snowpack_greenland_new/snowpack/SnowDrift.h"
    "/home/nander/snowpack_greenland_new/snowpack/SnowpackConfig.h"
    "/home/nander/snowpack_greenland_new/snowpack/Stability.h"
    "/home/nander/snowpack_greenland_new/snowpack/StabilityAlgorithms.h"
    "/home/nander/snowpack_greenland_new/snowpack/Utils.h"
    "/home/nander/snowpack_greenland_new/snowpack/libsnowpack.h"
    "/home/nander/snowpack_greenland_new/snowpack/vanGenuchten.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/snowpack/plugins" TYPE FILE FILES
    "/home/nander/snowpack_greenland_new/snowpack/plugins/AsciiIO.h"
    "/home/nander/snowpack_greenland_new/snowpack/plugins/CaaMLIO.h"
    "/home/nander/snowpack_greenland_new/snowpack/plugins/ImisDBIO.h"
    "/home/nander/snowpack_greenland_new/snowpack/plugins/SmetIO.h"
    "/home/nander/snowpack_greenland_new/snowpack/plugins/SnowpackIO.h"
    "/home/nander/snowpack_greenland_new/snowpack/plugins/SnowpackIOInterface.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/snowpack/snowpackCore" TYPE FILE FILES
    "/home/nander/snowpack_greenland_new/snowpack/snowpackCore/Aggregate.h"
    "/home/nander/snowpack_greenland_new/snowpack/snowpackCore/Canopy.h"
    "/home/nander/snowpack_greenland_new/snowpack/snowpackCore/Metamorphism.h"
    "/home/nander/snowpack_greenland_new/snowpack/snowpackCore/PhaseChange.h"
    "/home/nander/snowpack_greenland_new/snowpack/snowpackCore/ReSolver1d.h"
    "/home/nander/snowpack_greenland_new/snowpack/snowpackCore/SeaIce.h"
    "/home/nander/snowpack_greenland_new/snowpack/snowpackCore/Snowpack.h"
    "/home/nander/snowpack_greenland_new/snowpack/snowpackCore/Solver.h"
    "/home/nander/snowpack_greenland_new/snowpack/snowpackCore/VapourTransport.h"
    "/home/nander/snowpack_greenland_new/snowpack/snowpackCore/WaterTransport.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/nander/snowpack_greenland_new/snowpack/cmake_install.cmake")
  include("/home/nander/snowpack_greenland_new/applications/snowpack/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/nander/snowpack_greenland_new/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
