#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

//This is an example of what is possible using grids. As usual, for more information, have a look at the html documentation
int main(void) {
	Grid2DObject grid1, grid2;
	Config cfg("io.ini");
	const double TZ = cfg.get("TIME_ZONE", "Input");
	IOManager io(cfg);

	//reading the 2D grid by filename
	io.read2DGrid(grid1, "sub_dem.asc");
	//reading the 2D grid by parameter and date
	io.read2DGrid(grid2, MeteoGrids::TA, Date(2008, 12, 10, 12, 30, 0, TZ));

	//debug style output
	std::cout << "Initial air temperatures grid: " << grid2.toString() << "\n";

	//simple arithmetic operations (the nodata values are preserved)
	grid2 -= 273.15;
	std::cout << "Air temperatures grid in celsius: " << grid2.toString() << "\n";

	//operations between grids
	Grid2DObject grid3(grid1.getNx(), grid1.getNy(), grid1.cellsize, grid1.llcorner);
	grid3 = grid1 * (grid2 / 100.); //the parenthesis are required because of constness
	//std::cout << "Grids multiplication: " << grid3.toString() << "\n";

	//locate the cell that contains a specific point in real world coordinates
	Coords point("CH1903", "");
	point.setXY(559500., 221500., 1050.);
	grid1.gridify(point); //computes the point position (i,j) in the grid
	std::cout << "Position " << point.toString(Coords::LATLON) << " is point (" << point.getGridI() << "," << point.getGridJ() << ") in the grid with elevation=" << grid1.grid2D(point.getGridI(), point.getGridJ()) << "\n";

	//now let's make a grid subset: from point (2,2) and size 5x5
	Grid2DObject subgrid(grid1, 2, 2, 5, 5);
	std::cout << "The subgrid of grid1 is: " << subgrid.toString() << "\n";

	//we can check if two grids are "compatible", having the same geolocalization
	if(subgrid.isSameGeolocalization(grid1))
		std::cout << "The grids have the same geolocalization\n\n";
	else
		std::cout << "The grids do NOT have the same geolocalization\n\n";

	//it is not possible to make grid operations if they don't have the same size
	std::cout << "Q: What happens when we try to divide two grids that don't have the same size?\n\n";
	try {
		subgrid /= grid1;
	} catch(const std::exception &e) {
		std::cout << e.what();
		std::cout << "\nA: this is not possible and throws an exception!\n";
	}

	return 0;
}
