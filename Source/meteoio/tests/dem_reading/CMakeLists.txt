
## Test Meteo reading
# generate executable
ADD_EXECUTABLE(dem_reading dem_reading.cc)
TARGET_LINK_LIBRARIES(dem_reading ${METEOIO_LIBRARIES})

# add the tests
ADD_TEST(dem_reading.smoke dem_reading)
SET_TESTS_PROPERTIES(dem_reading.smoke 
					PROPERTIES LABELS smoke)




