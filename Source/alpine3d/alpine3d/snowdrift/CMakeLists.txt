INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}")

SET(snowdrift_sources
	snowdrift/SnowDrift.cc
	snowdrift/Saltation.cc
	snowdrift/Suspension.cc
	snowdrift/Sublimation.cc
	snowdrift/checksum.cc
	snowdrift/SnowDriftFEInit.cc
	snowdrift/SnowDriftFENumerics.cc
	snowdrift/SnowDriftFEControl.cc
	snowdrift/Cell.cc
)
