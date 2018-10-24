
## Test Coords
# generate executable
ADD_EXECUTABLE(coordinates coordinates.cc)
TARGET_LINK_LIBRARIES(coordinates ${METEOIO_LIBRARIES})

# add the tests
ADD_TEST(coordinates.smoke coordinates)
SET_TESTS_PROPERTIES(coordinates.smoke PROPERTIES LABELS smoke)



