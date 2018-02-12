INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}")

SET(ebalance_sources
	ebalance/EnergyBalance.cc
	ebalance/RadiationField.cc
	ebalance/TerrainRadiationAlgorithm.cc
	ebalance/TerrainRadiation.cc
	ebalance/TerrainRadiationSimple.cc
	ebalance/TerrainRadiationHelbig.cc
	ebalance/ViewFactors.cc
	ebalance/ViewFactorsHelbig.cc
	ebalance/ViewFactorsSectors.cc
	ebalance/ViewFactorsCluster.cc
)

IF(PETSC)
	FIND_PACKAGE(PETSc REQUIRED)
	INCLUDE_DIRECTORIES(${PETSC_INCLUDES}})
	LIST(APPEND ebalance_sources TerrainRadiationPETSc.cc)
ENDIF(PETSC)
IF(MPI)
	FIND_PACKAGE(MPI REQUIRED)
	INCLUDE_DIRECTORIES(SYSTEM ${MPI_INCLUDE_PATH})
ENDIF(MPI)
