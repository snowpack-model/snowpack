#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

//This is a basic example of using as dem: the dem is read, the grid coordinates of a point given by its (lat,long) are retrieved
//and a sub-dem is extracted starting at these coordinates and extending dist_x and dist_y and written out.
int main(void) {
	DEMObject dem;
	Config cfg("io.ini");
	IOManager io(cfg);

	//reading dem
	dem.setUpdatePpt(DEMObject::SLOPE);
	io.readDEM(dem);

	//writing some statistics about this dem
	//dem.grid2D.getMin() scans the DEM grid to get the min, while dem.min_altitude is cached and therefore very cheap
	//The raw content of the 2D grids can also be accessed, for example dem.grid2D.getMin(IOUtils::RAW_NODATA). In this case, there would be no interpretation of some values as nodata.
	std::cout << "DEM information: \n";
	std::cout << "\tmin=" << dem.grid2D.getMin() << " max=" << dem.grid2D.getMax() << " mean=" << dem.grid2D.getMean() << "\n";
	std::cout << "\tmin slope=" << dem.min_slope << " max slope=" << dem.max_slope << std::endl;

	io.write2DGrid(dem, MeteoGrids::DEM, Date(0.));

	Grid2DObject slope(dem.cellsize, dem.llcorner, dem.slope);
	io.write2DGrid(slope, MeteoGrids::SLOPE, Date(0.));
	Grid2DObject azi(dem.cellsize, dem.llcorner, dem.azi);
	io.write2DGrid(azi,"azi.png");

	return 0;
}
