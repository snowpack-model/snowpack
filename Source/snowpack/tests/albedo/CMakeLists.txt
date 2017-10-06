
## Test albedo

FIND_PACKAGE(MeteoIO)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${METEOIO_INCLUDE_DIR})
SET(extra_libs ${extra_libs} ${METEOIO_LIBRARIES})


# include sources
SET(test_sources
	albedoTest.cc
)


# generate executable
ADD_EXECUTABLE(albedoTest ${test_sources})
TARGET_LINK_LIBRARIES(albedoTest ${LIBRARIES} ${SHAREDNAME})

# add the tests
#ADD_TEST(albedo.smoke albedo.sh)
#SET_TESTS_PROPERTIES(albedo.smoke PROPERTIES LABELS smoke)



