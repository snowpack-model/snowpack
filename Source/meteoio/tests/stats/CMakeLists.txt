
## Test statistical functions
# generate executable
ADD_EXECUTABLE(stats stats.cc)
TARGET_LINK_LIBRARIES(stats ${METEOIO_LIBRARIES})

# add the tests
ADD_TEST(stats.smoke stats)
SET_TESTS_PROPERTIES(stats.smoke PROPERTIES LABELS smoke)
