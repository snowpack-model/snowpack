#SPDX-License-Identifier: LGPL-3.0-or-later
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")

SET(dataGenerators_sources
	dataGenerators/GeneratorAlgorithms.cc
	dataGenerators/AllSkyLWGenerator.cc
	dataGenerators/AllSkySWGenerator.cc
	dataGenerators/ClearSkyLWGenerator.cc
	dataGenerators/ClearSkySWGenerator.cc
	dataGenerators/ConstGenerator.cc
	dataGenerators/ESOLIPGenerator.cc
	dataGenerators/HumidityGenerator.cc
	dataGenerators/IswrAlbedoGenerator.cc
	dataGenerators/MeteoIndexGenerator.cc
	dataGenerators/PrecSplitting.cc
	dataGenerators/RadiationComponents.cc
	dataGenerators/SinGenerator.cc
	dataGenerators/StdPressGenerator.cc
	dataGenerators/TauCLDGenerator.cc
	dataGenerators/TsGenerator.cc
	dataGenerators/WindComponents.cc
)
