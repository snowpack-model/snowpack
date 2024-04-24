#SPDX-License-Identifier: LGPL-3.0-or-later
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")
INCLUDE_DIRECTORIES(SYSTEM "${PROJECT_SOURCE_DIR}/meteoio/thirdParty/Eigen/")
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/meteoio/thirdParty/ctsa/)

SET(thirdParty_sources
	thirdParty/tinyexpr.c
	thirdParty/ctsa/autoutils.c
	thirdParty/ctsa/boxcox.c
	thirdParty/ctsa/boxjenkins.c
	thirdParty/ctsa/brent.c
	thirdParty/ctsa/conjgrad.c
	thirdParty/ctsa/conv.c
	thirdParty/ctsa/ctsa.c
	thirdParty/ctsa/dist.c
	thirdParty/ctsa/emle.c
	thirdParty/ctsa/erfunc.c
	thirdParty/ctsa/filter.c
	thirdParty/ctsa/hsfft.c
	thirdParty/ctsa/initest.c
	thirdParty/ctsa/lls.c
	thirdParty/ctsa/lnsrchmp.c
	thirdParty/ctsa/matrix.c
	thirdParty/ctsa/neldermead.c
	thirdParty/ctsa/newtonmin.c
	thirdParty/ctsa/nls.c
	thirdParty/ctsa/optimc.c
	thirdParty/ctsa/pdist.c
	thirdParty/ctsa/polyroot.c
	thirdParty/ctsa/pred.c
	thirdParty/ctsa/real.c
	thirdParty/ctsa/regression.c
	thirdParty/ctsa/secant.c
	thirdParty/ctsa/spectrum.c
	thirdParty/ctsa/stats.c
	thirdParty/ctsa/stl.c
	thirdParty/ctsa/talg.c
	thirdParty/ctsa/errors.c
	thirdParty/ctsa/unitroot.c
	thirdParty/ctsa/seastest.c
	thirdParty/ctsa/wtmath.c
	thirdParty/ctsa/wavefilt.c
	thirdParty/ctsa/wavelib.c
)

# Add the Eigen library as part of the install target
INSTALL(
	DIRECTORY "${PROJECT_SOURCE_DIR}/meteoio/thirdParty/Eigen/"
	DESTINATION include/meteoio/thirdParty/Eigen
        COMPONENT headers
)
