
## Test linearsolver

FIND_PACKAGE(MeteoIO)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${METEOIO_INCLUDE_DIR})
SET(extra_libs ${extra_libs} ${METEOIO_LIBRARIES})


# generate executable
ADD_EXECUTABLE(linearSolverTest linearSolverTest.cc)
TARGET_LINK_LIBRARIES(linearSolverTest ${LIBRARIES})

# add the tests
ADD_TEST(linearsolver.smoke linearsolver.sh)
SET_TESTS_PROPERTIES(linearsolver.smoke PROPERTIES LABELS smoke)



