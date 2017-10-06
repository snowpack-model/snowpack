
## Test implicitsolver

FIND_PACKAGE(MeteoIO)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${METEOIO_INCLUDE_DIR})
SET(extra_libs ${extra_libs} ${METEOIO_LIBRARIES})


# include sources
SET(test_sources
	HeatEquationAnalytical.cc
	TestSnowpack.cc
	implicitSolverTest.cc
	
)


# generate executable
ADD_EXECUTABLE(implicitSolverTest ${test_sources})
TARGET_LINK_LIBRARIES(implicitSolverTest ${LIBRARIES} ${SHAREDNAME})

# add the tests
ADD_TEST(implicitsolver.smoke implicitsolver.sh)
SET_TESTS_PROPERTIES(implicitsolver.smoke PROPERTIES LABELS smoke)



