#include <iostream>
#include <algorithm>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio
using namespace std;

const double grid_epsilon = 1e-3; //1e-4 is still too tight because of truncated results when writing data out for the ref.
const bool gen_ref = false; //wether to generate ref files or to test

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

	std::string date_str = d1.toString(Date::ISO);
	std::replace( date_str.begin(), date_str.end(), ':', '.');

	//performing spatial interpolations
	Grid2DObject param, ref;
	int status = EXIT_SUCCESS;

	io.getMeteoData(d1, dem, MeteoData::TA, param);
	if(gen_ref) io.write2DGrid(param, MeteoGrids::TA, d1);
	else {
		io.read2DGrid(ref, date_str+"_TA_ref.asc");
		if(ref.grid2D.checkEpsilonEquality(param.grid2D, grid_epsilon)==false) {
			cout << "TA grids don't match!\n"; status = EXIT_FAILURE;
		}
	}

	io.getMeteoData(d1, dem, MeteoData::TSS, param);
	if(gen_ref) io.write2DGrid(param, MeteoGrids::TSS, d1);
	else {
		io.read2DGrid(ref, date_str+"_TSS_ref.asc");
		if(ref.grid2D.checkEpsilonEquality(param.grid2D, grid_epsilon)==false) {
			cout << "TSS grids don't match!\n"; status = EXIT_FAILURE;
		}
	}

	io.getMeteoData(d1, dem, MeteoData::TSG, param);
	if(gen_ref) io.write2DGrid(param, MeteoGrids::TSG, d1);
	else {
		io.read2DGrid(ref, date_str+"_TSG_ref.asc");
		if(ref.grid2D.checkEpsilonEquality(param.grid2D, grid_epsilon)==false) {
			cout << "TSG grids don't match!\n"; status = EXIT_FAILURE;
		}
	}

	io.getMeteoData(d1, dem, MeteoData::PSUM, param);
	if(gen_ref) io.write2DGrid(param, MeteoGrids::PSUM, d1);
	else {
		io.read2DGrid(ref, date_str+"_PSUM_ref.asc");
		if(ref.grid2D.checkEpsilonEquality(param.grid2D, grid_epsilon)==false) {
			cout << "PSUM grids don't match!\n"; status = EXIT_FAILURE;
		}
	}

	io.getMeteoData(d1, dem, MeteoData::RH, param);
	if(gen_ref) io.write2DGrid(param, MeteoGrids::RH, d1);
	else {
		io.read2DGrid(ref, date_str+"_RH_ref.asc");
		if(ref.grid2D.checkEpsilonEquality(param.grid2D, grid_epsilon)==false) {
			cout << "RH grids don't match!\n"; status = EXIT_FAILURE;
		}
	}

	io.getMeteoData(d1, dem, MeteoData::VW, param);
	if(gen_ref) io.write2DGrid(param, MeteoGrids::VW, d1);
	else {
		io.read2DGrid(ref, date_str+"_VW_ref.asc");
		if(ref.grid2D.checkEpsilonEquality(param.grid2D, grid_epsilon)==false) {
			cout << "VW grids don't match!\n"; status = EXIT_FAILURE;
		}
	}

	io.getMeteoData(d1, dem, MeteoData::RSWR, param);
	if(gen_ref) io.write2DGrid(param, MeteoGrids::RSWR, d1);
	else {
		io.read2DGrid(ref, date_str+"_RSWR_ref.asc");
		if(ref.grid2D.checkEpsilonEquality(param.grid2D, grid_epsilon)==false) {
			cout << "RSWR grids don't match!\n"; status = EXIT_FAILURE;
		}
	}

	io.getMeteoData(d1, dem, MeteoData::P, param);
	if(gen_ref) io.write2DGrid(param, MeteoGrids::P, d1);
	else {
		io.read2DGrid(ref, date_str+"_P_ref.asc");
		if(ref.grid2D.checkEpsilonEquality(param.grid2D, grid_epsilon)==false) {
			cout << "P grids don't match!\n"; status = EXIT_FAILURE;
		}
	}

	if(!gen_ref) io.write2DGrid(param, MeteoGrids::P, d1); //trying to write one grid out

	return status;
}
