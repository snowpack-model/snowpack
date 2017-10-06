# Install script for directory: /home/nander/meteoio

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/meteoio" TYPE FILE FILES
    "/home/nander/meteoio/meteoio/Config.h"
    "/home/nander/meteoio/meteoio/DataCreator.h"
    "/home/nander/meteoio/meteoio/DataGenerator.h"
    "/home/nander/meteoio/meteoio/FileUtils.h"
    "/home/nander/meteoio/meteoio/Graphics.h"
    "/home/nander/meteoio/meteoio/GridsManager.h"
    "/home/nander/meteoio/meteoio/IOExceptions.h"
    "/home/nander/meteoio/meteoio/IOHandler.h"
    "/home/nander/meteoio/meteoio/IOInterface.h"
    "/home/nander/meteoio/meteoio/IOManager.h"
    "/home/nander/meteoio/meteoio/IOUtils.h"
    "/home/nander/meteoio/meteoio/MainPage.h"
    "/home/nander/meteoio/meteoio/MathOptim.h"
    "/home/nander/meteoio/meteoio/MessageBoxX11.h"
    "/home/nander/meteoio/meteoio/Meteo1DInterpolator.h"
    "/home/nander/meteoio/meteoio/Meteo2DInterpolator.h"
    "/home/nander/meteoio/meteoio/MeteoIO.h"
    "/home/nander/meteoio/meteoio/MeteoProcessor.h"
    "/home/nander/meteoio/meteoio/ResamplingAlgorithms2D.h"
    "/home/nander/meteoio/meteoio/TimeSeriesManager.h"
    "/home/nander/meteoio/meteoio/Timer.h"
    "/home/nander/meteoio/meteoio/exports.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/meteoio/dataClasses" TYPE FILE FILES
    "/home/nander/meteoio/meteoio/dataClasses/Array1D.h"
    "/home/nander/meteoio/meteoio/dataClasses/Array2D.h"
    "/home/nander/meteoio/meteoio/dataClasses/Array3D.h"
    "/home/nander/meteoio/meteoio/dataClasses/Array4D.h"
    "/home/nander/meteoio/meteoio/dataClasses/Buffer.h"
    "/home/nander/meteoio/meteoio/dataClasses/Coords.h"
    "/home/nander/meteoio/meteoio/dataClasses/CoordsAlgorithms.h"
    "/home/nander/meteoio/meteoio/dataClasses/DEMObject.h"
    "/home/nander/meteoio/meteoio/dataClasses/Date.h"
    "/home/nander/meteoio/meteoio/dataClasses/Grid2DObject.h"
    "/home/nander/meteoio/meteoio/dataClasses/Grid3DObject.h"
    "/home/nander/meteoio/meteoio/dataClasses/Matrix.h"
    "/home/nander/meteoio/meteoio/dataClasses/MeteoData.h"
    "/home/nander/meteoio/meteoio/dataClasses/StationData.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/meteoio/dataGenerators" TYPE FILE FILES
    "/home/nander/meteoio/meteoio/dataGenerators/AllSkyLWGenerator.h"
    "/home/nander/meteoio/meteoio/dataGenerators/AllSkySWGenerator.h"
    "/home/nander/meteoio/meteoio/dataGenerators/ClearSkyLWGenerator.h"
    "/home/nander/meteoio/meteoio/dataGenerators/ClearSkySWGenerator.h"
    "/home/nander/meteoio/meteoio/dataGenerators/ConstGenerator.h"
    "/home/nander/meteoio/meteoio/dataGenerators/ESOLIPGenerator.h"
    "/home/nander/meteoio/meteoio/dataGenerators/GeneratorAlgorithms.h"
    "/home/nander/meteoio/meteoio/dataGenerators/IswrAlbedoGenerator.h"
    "/home/nander/meteoio/meteoio/dataGenerators/PPHASEGenerator.h"
    "/home/nander/meteoio/meteoio/dataGenerators/PrecUnsplit.h"
    "/home/nander/meteoio/meteoio/dataGenerators/RelHumGenerator.h"
    "/home/nander/meteoio/meteoio/dataGenerators/SinGenerator.h"
    "/home/nander/meteoio/meteoio/dataGenerators/StdPressGenerator.h"
    "/home/nander/meteoio/meteoio/dataGenerators/TauCLDGenerator.h"
    "/home/nander/meteoio/meteoio/dataGenerators/TsGenerator.h"
    "/home/nander/meteoio/meteoio/dataGenerators/template.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/meteoio/meteoResampling" TYPE FILE FILES
    "/home/nander/meteoio/meteoio/meteoResampling/Accumulate.h"
    "/home/nander/meteoio/meteoio/meteoResampling/DailyAverageResampling.h"
    "/home/nander/meteoio/meteoio/meteoResampling/DailySolarResampling.h"
    "/home/nander/meteoio/meteoio/meteoResampling/LinearResampling.h"
    "/home/nander/meteoio/meteoio/meteoResampling/NearestNeighbour.h"
    "/home/nander/meteoio/meteoio/meteoResampling/NoResampling.h"
    "/home/nander/meteoio/meteoio/meteoResampling/ResamplingAlgorithms.h"
    "/home/nander/meteoio/meteoio/meteoResampling/SolarResampling.h"
    "/home/nander/meteoio/meteoio/meteoResampling/template.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/meteoio/spatialInterpolations" TYPE FILE FILES
    "/home/nander/meteoio/meteoio/spatialInterpolations/ALSScaleAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/AvgAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/AvgLapseAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/ConstAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/IDWAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/IDWLapseAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/IDWLapseLocalAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/ILWREpsAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/InterpolationAlgorithms.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/ListonWindAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/NearestNeighbourAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/NoneAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/ODKrigAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/ODKrigLapseAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/PPhaseAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/RHListonAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/RyanWindAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/SnowPsumAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/StdPressAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/SwRadAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/UserAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/WinstralAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/WinstralListonAlgorithm.h"
    "/home/nander/meteoio/meteoio/spatialInterpolations/template.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/meteoio/plugins" TYPE FILE FILES
    "/home/nander/meteoio/meteoio/plugins/A3DIO.h"
    "/home/nander/meteoio/meteoio/plugins/ALPUG.h"
    "/home/nander/meteoio/meteoio/plugins/ARCIO.h"
    "/home/nander/meteoio/meteoio/plugins/ARPSIO.h"
    "/home/nander/meteoio/meteoio/plugins/BormaIO.h"
    "/home/nander/meteoio/meteoio/plugins/CNRMIO.h"
    "/home/nander/meteoio/meteoio/plugins/CosmoXMLIO.h"
    "/home/nander/meteoio/meteoio/plugins/DBO.h"
    "/home/nander/meteoio/meteoio/plugins/GRIBIO.h"
    "/home/nander/meteoio/meteoio/plugins/GSNIO.h"
    "/home/nander/meteoio/meteoio/plugins/GeotopIO.h"
    "/home/nander/meteoio/meteoio/plugins/GrassIO.h"
    "/home/nander/meteoio/meteoio/plugins/ImisIO.h"
    "/home/nander/meteoio/meteoio/plugins/NetCDFIO.h"
    "/home/nander/meteoio/meteoio/plugins/OshdIO.h"
    "/home/nander/meteoio/meteoio/plugins/PGMIO.h"
    "/home/nander/meteoio/meteoio/plugins/PNGIO.h"
    "/home/nander/meteoio/meteoio/plugins/PSQLIO.h"
    "/home/nander/meteoio/meteoio/plugins/SASEIO.h"
    "/home/nander/meteoio/meteoio/plugins/SMETIO.h"
    "/home/nander/meteoio/meteoio/plugins/SNIO.h"
    "/home/nander/meteoio/meteoio/plugins/exports.h"
    "/home/nander/meteoio/meteoio/plugins/libMatioWrapper.h"
    "/home/nander/meteoio/meteoio/plugins/libncpp.h"
    "/home/nander/meteoio/meteoio/plugins/libsmet.h"
    "/home/nander/meteoio/meteoio/plugins/picojson.h"
    "/home/nander/meteoio/meteoio/plugins/template.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/meteoio/meteoLaws" TYPE FILE FILES
    "/home/nander/meteoio/meteoio/meteoLaws/Atmosphere.h"
    "/home/nander/meteoio/meteoio/meteoLaws/Meteoconst.h"
    "/home/nander/meteoio/meteoio/meteoLaws/Sun.h"
    "/home/nander/meteoio/meteoio/meteoLaws/Suntrajectory.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/meteoio/meteoFilters" TYPE FILE FILES
    "/home/nander/meteoio/meteoio/meteoFilters/FilterBlock.h"
    "/home/nander/meteoio/meteoio/meteoFilters/FilterDeGrass.h"
    "/home/nander/meteoio/meteoio/meteoFilters/FilterMAD.h"
    "/home/nander/meteoio/meteoio/meteoFilters/FilterMax.h"
    "/home/nander/meteoio/meteoio/meteoFilters/FilterMin.h"
    "/home/nander/meteoio/meteoio/meteoFilters/FilterMinMax.h"
    "/home/nander/meteoio/meteoio/meteoFilters/FilterNoChange.h"
    "/home/nander/meteoio/meteoio/meteoFilters/FilterPotentialSW.h"
    "/home/nander/meteoio/meteoio/meteoFilters/FilterRate.h"
    "/home/nander/meteoio/meteoio/meteoFilters/FilterStdDev.h"
    "/home/nander/meteoio/meteoio/meteoFilters/FilterSuppr.h"
    "/home/nander/meteoio/meteoio/meteoFilters/FilterTimeconsistency.h"
    "/home/nander/meteoio/meteoio/meteoFilters/FilterTukey.h"
    "/home/nander/meteoio/meteoio/meteoFilters/FilterUnheatedPSUM.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcAdd.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcAggregate.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcExpSmoothing.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcIIR.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcMult.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcPSUMDistribute.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcShade.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcUndercatch_Forland.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcUndercatch_Hamon.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcUndercatch_WMO.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcUnventilatedT.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcWMASmoothing.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcessingBlock.h"
    "/home/nander/meteoio/meteoio/meteoFilters/ProcessingStack.h"
    "/home/nander/meteoio/meteoio/meteoFilters/WindowedFilter.h"
    "/home/nander/meteoio/meteoio/meteoFilters/template.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/meteoio/meteoStats" TYPE FILE FILES
    "/home/nander/meteoio/meteoio/meteoStats/libfit1D.h"
    "/home/nander/meteoio/meteoio/meteoStats/libfit1DCore.h"
    "/home/nander/meteoio/meteoio/meteoStats/libinterpol1D.h"
    "/home/nander/meteoio/meteoio/meteoStats/libinterpol2D.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/nander/meteoio/meteoio/dataClasses/cmake_install.cmake")
  include("/home/nander/meteoio/meteoio/dataGenerators/cmake_install.cmake")
  include("/home/nander/meteoio/meteoio/plugins/cmake_install.cmake")
  include("/home/nander/meteoio/meteoio/meteoLaws/cmake_install.cmake")
  include("/home/nander/meteoio/meteoio/meteoFilters/cmake_install.cmake")
  include("/home/nander/meteoio/meteoio/meteoResampling/cmake_install.cmake")
  include("/home/nander/meteoio/meteoio/spatialInterpolations/cmake_install.cmake")
  include("/home/nander/meteoio/meteoio/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/nander/meteoio/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
