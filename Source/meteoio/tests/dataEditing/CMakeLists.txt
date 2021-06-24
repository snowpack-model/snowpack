#SPDX-License-Identifier: LGPL-3.0-or-later
## Test Meteo reading
# generate executable
ADD_EXECUTABLE(data_editing data_editing.cc)
TARGET_LINK_LIBRARIES(data_editing ${METEOIO_LIBRARIES})

# add the tests
ADD_TEST(data_editing.smoke data_editing)
SET_TESTS_PROPERTIES(data_editing.smoke 
                     PROPERTIES LABELS smoke)


