
## Test arrays
# generate executable
ADD_EXECUTABLE(arrays arrays.cc)
TARGET_LINK_LIBRARIES(arrays ${METEOIO_LIBRARIES})

# add the tests
ADD_TEST(arrays.smoke arrays)
SET_TESTS_PROPERTIES(arrays.smoke PROPERTIES LABELS smoke)




