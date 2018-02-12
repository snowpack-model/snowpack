/***********************************************************************************/
/*  Copyright 2009-2015 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of Alpine3D.
    Alpine3D is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Alpine3D is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Alpine3D.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
#include <meteoio/MeteoIO.h>
#include "Runoff.h"

using namespace mio; //The MeteoIO namespace is called mio
using namespace std;

int main(int /*argc*/, char** argv) {
	Config cfg("runoff.ini");
	IOManager io(cfg);
	Runoff runoff(io, 1.3); //set rain/snow threshold at 1.3 C

	//reading dem
	DEMObject dem;
	dem.setUpdatePpt(DEMObject::SLOPE);
	io.readDEM(dem);

	//getting start and end dates
	double TZ;
	cfg.getValue("TIME_ZONE", "Input", TZ);
	Date start_date, end_date;
	IOUtils::convertString(start_date,argv[1], TZ);
	const std::string end_date_str = string(argv[2]);
	if (end_date_str == "NOW") { //interpret user provided end date
		end_date.setFromSys();
		end_date.setTimeZone(TZ);
		end_date.rnd(1800, mio::Date::DOWN);
	} else {
		IOUtils::convertString(end_date, end_date_str, TZ);
	}
	cerr << "Runoff simulation from " << start_date.toString(Date::ISO) << " to " << end_date.toString(Date::ISO) << "\n";

	//feeding Runoff the data and running it
	runoff.initialize(dem, cfg);
	mio::Array2D<double> dummy;
	mio::Array2D<double> dummy_c;
	const double step = 1./24.;
	const int info_interval = floor( (end_date.getJulian(true) - start_date.getJulian(true)) / (step*10.) );
	int counter = 0;
	for (Date d1=start_date; d1<=end_date; d1 += step) {
		counter++;
		if (counter % info_interval == 0) cout << "\tAt " << d1.toString(Date::ISO) << "\n";
		Grid2DObject rot;
		try {
			io.read2DGrid(rot, MeteoGrids::ROT, d1);
		} catch (const IOException &e) {
			cerr << "Stopping simulation: " << e.what() << endl;
			break;
		}
		runoff.setRunoff(dummy, rot.grid2D, dummy_c, dummy, dummy);
		runoff.output(d1);
	}

	return 0;
}
