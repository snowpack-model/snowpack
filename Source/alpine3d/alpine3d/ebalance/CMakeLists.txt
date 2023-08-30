INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}")

SET(ebalance_sources
	ebalance/EnergyBalance.cc
	ebalance/RadiationField.cc
	ebalance/TerrainRadiationAlgorithm.cc
	ebalance/TerrainRadiationSimple.cc
	ebalance/TerrainRadiationHelbig.cc
	ebalance/ViewFactors.cc
	ebalance/ViewFactorsHelbig.cc
	ebalance/ViewFactorsSectors.cc
	ebalance/ViewFactorsCluster.cc
	ebalance/SolarPanel.cc
	ebalance/TerrainRadiationComplex.cc
	ebalance/SnowBRDF.cc
)

IF(MPI)
	FIND_PACKAGE(MPI REQUIRED)
	INCLUDE_DIRECTORIES(SYSTEM ${MPI_INCLUDE_PATH})
ENDIF(MPI)
