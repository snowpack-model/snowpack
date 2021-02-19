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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOB.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Alpine3D.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <alpine3d/ebalance/SnowBRDF.h>


#include <fstream> // file operations

using namespace mio;


SnowBRDF::SnowBRDF(const mio::Config& cfg){
	initialize(cfg);
}

SnowBRDF::~SnowBRDF() {}

// initialize BRDF-array from files
void SnowBRDF::initialize(const mio::Config& cfg){
	unsigned n_thi=33;
	unsigned n_thv=60;
	unsigned n_phi=120;
	for (unsigned i = 1; i < n_thi+1; ++i)
	{
		std::string path;
		cfg.getValue("BRDFPATH", "EBalance", path);
		std::string mid="/weighted_cos=";
		char theta[3];
		sprintf(theta, "%d", 3*i);
		std::string filename = path + mid + theta;
		std::ifstream file(filename);

		if(file.is_open()){
		for (unsigned j = 0; j < n_phi; ++j)
			{
				for (unsigned k = 0; k < n_thv; ++k)
				{
					file >> BRDF_data[i-1][j][k];
				}
			}
		}
		if(file.is_open()==0) std::cout<<"[e] in BRDF: file ’"<<filename<<"’not found.\n";
		if(file.is_open()==0) throw std::invalid_argument( "File Not Found\n" );
	}
}

// linearly interpolates BRDF_data and gives out corresponding value.
double SnowBRDF::get_RF(double cth_i, double cphi, double cth_v){

	int cthi_plus, cthi_minus, cthv_plus, cthv_minus, phi_plus, phi_minus;
	double d_i,d_v,d_phi; // difference to data point
	double phi;
	double R3D[2][2][2];
	double R2D[2][2];
	double R1D[2];
	double R;

	d_i=std::fmod(100*cth_i/3,1);
	cthi_minus=100*cth_i/3-d_i-1; //could just use (int)100*cth_i/3 b

	if(cthi_minus<0){
		cthi_minus=0;
		d_i=0;
	}

	cthi_plus=cthi_minus+1;

	if(cthi_minus==32) cthi_plus=32;

	d_v=std::fmod(59*cth_v,1);
	cthv_minus=59*cth_v-d_v;

	if(cthv_minus<0){
		cthv_minus=0;
		d_v=0;
	}

	cthv_plus=cthv_minus+1;
	if(cthv_minus==59) cthv_plus=59;

	if (cphi>1) cphi=1;
	phi=180-acos(cphi)*180/Cst::PI; // data starts at 180 deg
	d_phi=std::fmod(120*phi/180,1);
	phi_minus=120*phi/180-d_phi-1;

	if(phi_minus<0){
		phi_minus=0;
		d_phi=0;
	}

	phi_plus=phi_minus+1;
	if(phi_minus==119) phi_plus=119;

	R3D[0][0][0]=BRDF_data[cthi_minus][phi_minus][cthv_minus];
	R3D[1][0][0]=BRDF_data[cthi_plus][phi_minus][cthv_minus];
	R3D[0][1][0]=BRDF_data[cthi_minus][phi_plus][cthv_minus];
	R3D[0][0][1]=BRDF_data[cthi_minus][phi_minus][cthv_plus];
	R3D[1][1][0]=BRDF_data[cthi_plus][phi_plus][cthv_minus];
	R3D[1][0][1]=BRDF_data[cthi_plus][phi_minus][cthv_plus];
	R3D[0][1][1]=BRDF_data[cthi_minus][phi_plus][cthv_plus];
	R3D[1][1][1]=BRDF_data[cthi_plus][phi_plus][cthv_plus];

	for (int i = 0; i < 2; ++i){
		for (int j = 0; j < 2; ++j){
			R2D[i][j]=d_v*R3D[i][j][1]+(1-d_v)*R3D[i][j][0];
		}
	}

	for (int j = 0; j < 2; ++j){
		R1D[j]=d_phi*R2D[j][1]+(1-d_phi)*R2D[j][0];
	}

	R=d_i*R1D[1]+(1-d_i)*R1D[0];

	return R;
}
