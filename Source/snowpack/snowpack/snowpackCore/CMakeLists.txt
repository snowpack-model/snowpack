#SPDX-License-Identifier: LGPL-3.0-or-later
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")

SET(core_sources
	snowpackCore/Snowpack.cc
	snowpackCore/Aggregate.cc
	snowpackCore/WaterTransport.cc
	snowpackCore/VapourTransport.cc
	snowpackCore/ReSolver1d.cc
	snowpackCore/Canopy.cc
	snowpackCore/Metamorphism.cc
	snowpackCore/PhaseChange.cc
	snowpackCore/SeaIce.cc
	snowpackCore/Solver.cc
	snowpackCore/SalinityTransport.cc
)
