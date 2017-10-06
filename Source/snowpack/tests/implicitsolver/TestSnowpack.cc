/*
 * TestSnowpack.cc
 *
 *  Created on: Jun 30, 2017
 *      Author: julien
 */

#include "TestSnowpack.h"

TestSnowpack::TestSnowpack(const SnowpackConfig& i_cfg)
    : Snowpack(i_cfg) {
  // TODO Auto-generated constructor stub
}



bool TestSnowpack::testCompTempProfile(const CurrentMeteo& Mdata,
                                       SnowStation& Xdata, BoundCond& Bdata,
                                       const bool& ThrowAtNoConvergence) {

  return this->compTemperatureProfile(Mdata, Xdata, Bdata, ThrowAtNoConvergence);
}

void TestSnowpack::setBC(BoundaryCondition bc) {

  surfaceCode = bc;

}
