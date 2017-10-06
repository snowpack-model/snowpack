#include <stdio.h>
#include <fstream>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio
using namespace std;

// Constants
const double epsilon = 1.0e-7;

// --- File names ---
unsigned int n_files= 3;
string files[] ={"DEM.asc", "AZI.asc","SLOPE.asc"};
string prefix_ref("ref_");

// --- Results ---
// Header:				min		max			mean		s. min	s. max		[s. == slope]
double r_DEM[]		={	193.,	4204.81,	1302.38,	0.,		43.486};
double r_SUB_DEM[]	={	403.4,	4027.3,		1291.28,	0.,		40.0628};

// controll basic values on dem
bool simpleDEMcontroll(DEMObject& dem, double results[]) {
	bool status = true;

	if(!IOUtils::checkEpsilonEquality(dem.grid2D.getMin(), results[0], epsilon)){
		cerr << "error on Min : " << dem.grid2D.getMin() << " != " << results[0] << endl;
		status=false;
	}

	if(!IOUtils::checkEpsilonEquality(dem.grid2D.getMax(), results[1], epsilon)){
		cerr << "error on Max : " << dem.grid2D.getMax() << " != " << results[1]<< endl;
		status=false;
	}

	if(!IOUtils::checkEpsilonEquality(dem.grid2D.getMean(), results[2],  1.0e-2)){ // HACK adapt epsilon for test...
		cerr << "error on Mean : " << dem.grid2D.getMean() << " != " << results[2]<< endl;
		status=false;
	}


	if(!IOUtils::checkEpsilonEquality(dem.min_slope, results[3], epsilon)){
		cerr << "error on Slope Min : " << dem.min_slope << " != " << results[3]<< endl;
		status=false;
	}

	if(!IOUtils::checkEpsilonEquality(dem.max_slope, results[4],  1.0e-4)){ // HACK adapt epsilon for test ....
		cerr << "error on Slope Max : " << dem.max_slope << " != " << results[4]<< endl;
		status=false;
	}

	return status;
}

// Make output files
bool makeDEMfiles() {
	bool status=true;
	cout << " ----- Read DEM, make subfiles and some basic controll \n";

	cout << " ---- Init Values \n";
	DEMObject dem;
	Config cfg("io.ini");
	IOManager io(cfg);

	cout << " ---- If output files exist, empty them \n"; // HACK need to make this ???
	for(unsigned int i = 0; i < n_files; i++){
		ofstream ofs(files[i].c_str(), ofstream::out | ofstream::trunc);
		ofs << " " ;
		ofs.close();
	}

	cout << " ---- Read DEM \n";
	dem.setUpdatePpt(DEMObject::SLOPE);
	io.readDEM(dem);

	cout << " ---- Generate file with slope values \n";
	if (MeteoGrids::getParameterName(MeteoGrids::SLOPE)!="SLOPE") {
		std::cout << "getParameterName() returned " << MeteoGrids::getParameterName(MeteoGrids::SLOPE) << " when SLOPE was expected!\n";
		status=false;
	}
	Grid2DObject slope(dem.cellsize, dem.llcorner, dem.slope);
	io.write2DGrid(slope, MeteoGrids::SLOPE, Date(0.));

	cout << " ---- Generate file with azimut values \n";
	if (MeteoGrids::getParameterName(MeteoGrids::AZI)!="AZI") {
		std::cout << "getParameterName() returned " << MeteoGrids::getParameterName(MeteoGrids::AZI) << " when AZI was expected!\n";
		status=false;
	}
	Grid2DObject azi(dem.cellsize, dem.llcorner, dem.azi);
	io.write2DGrid(azi, MeteoGrids::AZI, Date(0.));

	cout << " ---- Gridify \n"; // HACK how to controll, wath are points possible to controll
	//retrieving grid coordinates of a real world point
	Coords point;
	point.copyProj(dem.llcorner); //we use the same projection parameters as the DEM
	point.setLatLon(46.232103, 7.362185, IOUtils::nodata);
	dem.gridify(point);

	cout << " ---- Get out a Sub dem and writh it in file \n";
	//computing grid distances from real world distances
	const double dist_x=70000., dist_y=120000.;
	const unsigned int ncols = (unsigned int)ceil(dist_x/dem.cellsize);
	const unsigned int nrows = (unsigned int)ceil(dist_y/dem.cellsize);
	DEMObject sub_dem(dem, point.getGridI(), point.getGridJ(), ncols, nrows);
	sub_dem.grid2D.setKeepNodata(true); //HACK this removes error, but should it not take it in constructor ???
	io.write2DGrid(sub_dem, MeteoGrids::DEM, Date(0.));

	cout << " --- Controll some basic data from DEM\n";
	if(!simpleDEMcontroll(dem, r_DEM)){
		cerr << " Error on controlling basic values of DEM object !!! \n";
		status=false;
	}

	cout << " --- Controll some basic data from SUB_DEM\n";
	if(!simpleDEMcontroll(sub_dem, r_SUB_DEM)){
		cerr << " Error on controlling basic values of DEM object !!! \n";
		status=false;
	}

	return status;
}

bool compareFiles() {
	bool status=true;
	for(unsigned int i = 0; i < n_files; i++) {
		// ------ Compare reference file with generated results ---------
		cout << " --- Comparing reference file with generated output file " << files[i] << endl;
		ifstream ifref((prefix_ref+files[i]).c_str());
		ifstream ifout(files[i].c_str());
		string l_ref, l_out;
		size_t lineNr = 0;

		while (!ifref.eof()) {
			if (ifout.eof()) {
				cerr << "Generated file has less lines than the reference for " << files[i] << endl;
				return false;
			}

			getline(ifref,l_ref);
			getline(ifout,l_out);
			lineNr++;
			if (l_ref!=l_out) {
				cerr << "Generated and reference files differ for " << files[i] << " at line " << lineNr << endl;
				cerr << "ref : \n " << l_ref << endl;
				cerr << "out : \n " << l_out << endl;
				status=false;
			}
		}

		if (!ifout.eof()) {
			cerr << "Generated file has more lines than the reference for " <<  files[i] << endl;
			status=false;
		}
		ifout.close();
		ifref.close();
	}

	return status;
}

int main(void) {

	if(!makeDEMfiles()) {
		exit(1);
	}

	if(!compareFiles()) {
		exit(1);
	}

	return 0;
}
