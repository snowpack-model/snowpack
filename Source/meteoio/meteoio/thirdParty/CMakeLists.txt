#SPDX-License-Identifier: LGPL-3.0-or-later
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")
INCLUDE_DIRECTORIES(SYSTEM "${PROJECT_SOURCE_DIR}/meteoio/thirdParty/Eigen/")

SET(thirdParty_sources
	thirdParty/tinyexpr.c
)

# Add the Eigen library as part of the install target
INSTALL(
	DIRECTORY "${PROJECT_SOURCE_DIR}/meteoio/thirdParty/Eigen/"
	DESTINATION include/meteoio/thirdParty/Eigen
        COMPONENT headers
)
