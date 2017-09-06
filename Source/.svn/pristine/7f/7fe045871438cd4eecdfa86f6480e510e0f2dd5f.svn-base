#include <iostream>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

//This is the a basic example of spatial interpolations. It does not check any exceptions, it only tries to be as c-like as possible
//provide date as ISO formatted, for example 2008-12-01T15:35:00 and
//it will retrieve and interpolate the data for this date according to the io.ini configuration file
int main(int /*argc*/, char** argv) {
	Date d1;

	//initializing the io handlers according to the config file
	Config cfg("io.ini");
	IOManager io(cfg);

	//reading the dem (necessary for several spatial interpolations algoritms)
	DEMObject dem;
	io.readDEM(dem);

	//we assume that the time given on the command line is in TZ=+1
	IOUtils::convertString(d1,argv[1], 1.);

	//performing spatial interpolations
	Grid2DObject param;
	io.getMeteoData(d1, dem, MeteoData::TA, param);
	io.write2DGrid(param, MeteoGrids::TA, d1);
	io.getMeteoData(d1, dem, MeteoData::PSUM, param);
	io.write2DGrid(param, MeteoGrids::PSUM, d1);
	io.getMeteoData(d1, dem, MeteoData::RH, param);
	io.write2DGrid(param, MeteoGrids::RH, d1);

	return 0;
}
