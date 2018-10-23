
## Test Meteo reading
# generate executable
ADD_EXECUTABLE(meteo_reading_int meteo_reading.cc)
TARGET_LINK_LIBRARIES(meteo_reading_int ${METEOIO_LIBRARIES})

# add the tests
ADD_TEST(meteo_reading_interpol.smoke meteo_reading_int)
SET_TESTS_PROPERTIES(meteo_reading_interpol.smoke 
					PROPERTIES LABELS smoke)




