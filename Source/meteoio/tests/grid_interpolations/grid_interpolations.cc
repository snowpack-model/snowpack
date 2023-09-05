// SPDX-License-Identifier: LGPL-3.0-or-later
#include <iostream>
#include <stdexcept>
#include <vector>

#include <meteoio/MeteoIO.h>

int main(int /*argc*/, char** /*argv*/)
{
	mio::Config cfg("./io.ini");
	mio::IOManager io(cfg);

	mio::DEMObject dem;
	io.readDEM(dem);

	//dates for dummy grids (must be covered by ../input/meteo):
	const std::string begin_date( "2008-12-01T00:00" );
	const std::string end_date( "2008-12-05T00:00" );
	const double TZ = cfg.get("TIME_ZONE", "Input"); //get dates according to time zone
	mio::Date sdate, edate;
	mio::IOUtils::convertString(sdate, begin_date, TZ); //start date
	mio::IOUtils::convertString(edate, end_date, TZ); //end date

	//generate a couple of mockup grids to test grid resampling with:
	const bool gen_grids = false;
	if (gen_grids) {
		mio::Grid2DObject mock_grid(dem); //base grid is the provided DEM
		mio::Date dt(sdate);
		while (dt < edate) //grids for each day between cmdline args
		{
			int year, month, day, hour, min, sec;
			dt.getDate(year, month, day, hour, min, sec);
			mio::Date dt_hours(dt);
			static const size_t hours[3] = {8, 12, 16};
			for (int ii = 0; ii < 3; ++ii) { //a couple of grids per day
				dt_hours.setDate(year, month, day, hours[ii], 0, 0);
				mock_grid = ii; //set all coordinates to some dummy value
				io.write2DGrid(mock_grid, mio::MeteoGrids::TA, dt_hours);
			}
			dt = dt + 1.;
		}
		return 0;
	} //endif gen_grids


	/* VSTATIONS */
	std::vector< std::vector<mio::MeteoData> > mvec;
	io.getMeteoData(sdate, edate, mvec); //get meteo data from a VStation according to INI
	io.writeMeteoData(mvec);

	/* TEMPORAL GRID RESAMPLING */
	mio::Grid2DObject grid_ta;

	//currently we must manually switch from resampling to regridding:
	cfg.deleteKey("RESAMPLING_STRATEGY", "InputEditing");
	cfg.addKey("REGRIDDING_STRATEGY", "InputEditing", "GRID_1DINTERPOLATE");
	mio::IOManager io2(cfg); //reload strategies

	mio::Date inbetween(sdate);
	inbetween.setDate(2008, 12, 1, 10, 0, 0); //pick some date that is not there as raw grid
	io2.getMeteoData(inbetween, dem, mio::MeteoData::TA, grid_ta);
	io2.write2DGrid(grid_ta, mio::MeteoGrids::TA, inbetween);

	return EXIT_SUCCESS;
}
