
## Test Meteo reading
# generate executable
ADD_EXECUTABLE(2D_interpolations 2D_interpolations.cc)
TARGET_LINK_LIBRARIES(2D_interpolations ${METEOIO_LIBRARIES})

# add the tests
ADD_TEST(2D_interpolations.smoke 2D_interpolations 2009-01-19T12:00)
SET_TESTS_PROPERTIES(2D_interpolations.smoke
                     PROPERTIES LABELS smoke)
