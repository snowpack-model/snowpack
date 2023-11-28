#SPDX-License-Identifier: LGPL-3.0-or-later
## Test Meteo fstream
# generate executable
ADD_EXECUTABLE(fstream ofstream_test.cc)
TARGET_LINK_LIBRARIES(fstream ${METEOIO_LIBRARIES})

# add the tests
ADD_TEST(fstream.smoke fstream)
SET_TESTS_PROPERTIES(fstream.smoke 
					PROPERTIES LABELS smoke)

