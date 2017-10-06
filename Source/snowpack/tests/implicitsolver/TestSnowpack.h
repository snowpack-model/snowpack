/*
 * TestSnowpack.h
 *
 *  Created on: Jun 30, 2017
 *      Author: julien
 */

#ifndef TESTS_IMPLICITSOLVER_TESTSNOWPACK_H_
#define TESTS_IMPLICITSOLVER_TESTSNOWPACK_H_

#include <snowpack/snowpackCore/Snowpack.h>

class TestSnowpack : public Snowpack {

 public:

  TestSnowpack(const SnowpackConfig& i_cfg);

  bool testCompTempProfile(const CurrentMeteo& Mdata, SnowStation& Xdata,
                           BoundCond& Bdata, const bool& ThrowAtNoConvergence);

  void setBC(BoundaryCondition bc);
};

#endif /* TESTS_IMPLICITSOLVER_TESTSNOWPACK_H_ */
