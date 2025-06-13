#SPDX-License-Identifier: LGPL-3.0-or-later
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")

SET(gridResampling_sources
	gridResampling/GridResamplingAlgorithms.cc
	gridResampling/GridNearestResampling.cc
	gridResampling/GridLinearResampling.cc
	gridResampling/GridTimeseriesResampling.cc
)
