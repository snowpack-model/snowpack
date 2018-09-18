
## Test Sun
# generate executable
ADD_EXECUTABLE(sun sun.cc)
TARGET_LINK_LIBRARIES(sun ${METEOIO_LIBRARIES})

# add the tests
ADD_TEST(sun.smoke sun.sh)
SET_TESTS_PROPERTIES(sun.smoke
                     PROPERTIES LABELS smoke
                     FAIL_REGULAR_EXPRESSION "error|differ")



