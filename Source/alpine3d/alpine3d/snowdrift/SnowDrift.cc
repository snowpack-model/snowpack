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
/*------------------------------------------------------------------------------------------+
 |  This module is contains the driving routines for a three dimensional numerical model of |
 |  snow drift. It was started in January 00 when Michael had a broken arm and was thinking |
 |  about Betty's surprise birthday party with Fondue which was coming up                   |
 +------------------------------------------------------------------------------------------*/
/********************************************************************************************/

#include <assert.h>
#include <vector>

#include <alpine3d/snowdrift/SnowDrift.h>
#include <alpine3d/AlpineMain.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"

using namespace mio;
using namespace std;

/************************************************************************************/

const double SnowDriftA3D::kinematicViscosityAir = 1.3E-5;	//Returns the kinematic viscosity of air (m2 s-1), actually a function of air temperature, The value here is taken from Liston and Sturm (1998)
const double SnowDriftA3D::USTAR=0.87; //ustar for boundary layer
const double SnowDriftA3D::molecularWeightofWater = 18.015E-3; //[kg mol-1]
const double SnowDriftA3D::thermalConductivityofAtm = 0.024;//[J m-1 s-1 K-1] Value taken from Liston and Sturm (1998)
const double SnowDriftA3D::c_red = 1.0; //Defines red. of saltation conc. in suspension sol.
const double SnowDriftA3D::grain_size = 0.000680;
const double SnowDriftA3D::tau_thresh = 0.094; //Original Value by Judith: 0.094
const double SnowDriftA3D::z0 = 0.01; //Wind Field Z0 - includes larger surface features
const bool SnowDriftA3D::thresh_snow = true;//Flag to determine whether ustar_thresh is calculated from the Snowpack properties


SnowDriftA3D::SnowDriftA3D(const DEMObject& dem, const mio::Config& cfg) 
                        : saltation_obj(cfg), auxLayerHeight(0.02), io(cfg), snowpack(NULL), eb(NULL), 
                        cH(dem, IOUtils::nodata), sp(dem, IOUtils::nodata), rg(dem, IOUtils::nodata), N3(dem, IOUtils::nodata), rb(dem, IOUtils::nodata),
                        nx(0), ny(0), nz(0), vw(dem, IOUtils::nodata), rh(dem, IOUtils::nodata), ta(dem, IOUtils::nodata), p(dem, IOUtils::nodata), 
                        psum(dem, IOUtils::nodata), psum_ph(dem, IOUtils::nodata), STATIONARY(true)
{
	const string wind_field_string = cfg.get("WINDFIELDS", "Input");

	vector<string> TA_interpol;
	cfg.getValues("TA::algorithms", "Interpolations2D", TA_interpol);
	if (TA_interpol.empty())
		throw InvalidArgumentException("No 2D interpolation algorithm defined for TA", AT);
	for (size_t ii=0; ii<TA_interpol.size(); ii++) {
		const string current = TA_interpol[ii];
		if (current!="CST" && current!="AVG" && current!="AVG_LAPSE") {
			if (MPIControl::instance().master()) 
				cout << "[W] for Snowdrift with sublimation, it is recommended to use CST or AVG or AVG_LAPSE as 2D interpolations for TA\n";
			break;
		}
	}

	buildWindFieldsTable(wind_field_string);
	Initialize();
	InitializeFEData();
}

std::string SnowDriftA3D::getGridsRequirements() const
{
	return "HS SP RG N3 RB STORE SWE"; //see setSnowSurfaceData() and saltation
}

void SnowDriftA3D::buildWindFieldsTable(const std::string& wind_field_string)
{//parses the given string containing the wind fields and number of steps, for example: "NW 3 SE 5 N 6"
//the result is written into a vector of wind_field structures, wind_fields
	int previous_time_step = 0, steps, fields_count=0;
	std::string word;
	std::istringstream iss(wind_field_string, std::istringstream::in);
	struct WIND_FIELD curr_field;

	while ( iss >> word ) {
		curr_field.wind = word;
		if (!(iss>>word)) {
			throw InvalidArgumentException("Missing number of time steps for wind field "+curr_field.wind, AT);
		}
		IOUtils::convertString(steps, word);
		curr_field.start_step = previous_time_step;
		previous_time_step += steps;
		fields_count ++;
		wind_fields.push_back(curr_field);
	}
	wind_field_index = -1;
	cout << "[I] Snowdrift will use " << fields_count << " wind situations" << endl;
}

bool SnowDriftA3D::isNewWindField(const unsigned int current_step)
{//this function returns true if a new wind field should be read
	const unsigned int next_index = wind_field_index+1;

	if (wind_field_index==-1) {
		return true;
	}
	if ((next_index<wind_fields.size()) && (current_step >= wind_fields[next_index].start_step)) {
		return true;
	} else {
		return false;
	}
}

/**
 * @brief Initialize for snowdrift
 * Allocate saltation variable, snow height and new snow mass per bottom element
 */
void SnowDriftA3D::Initialize()
{
	Grid3DObject z_readMatr;
	io.read3DGrid(z_readMatr, wind_fields[0].wind+":DEM");
	z_readMatr.llcorner = ta.llcorner; //most probably the user did not properly specify the ARPS_XCOORDS and ARPS_YCOORDS
	
	nx = z_readMatr.getNx();
	ny = z_readMatr.getNy();
	nz = z_readMatr.getNz();

	/* Allocate some memory and Initialize */
	/* Before entering the time integration loop allocate the saltation variable,
	the snow height and new snow mass per bottom element */
	saltation.resize(nx,ny, 0.);
	c_salt.resize(nx,ny);
	
	mns.set(nx, ny, ta.cellsize, ta.llcorner, 0.);
	mns_subl.resize(nx,ny);
	mns_nosubl.resize(nx,ny);
	dif_mns_subl.resize(nx,ny);
	flux_x.resize( nx*ny );
	flux_y.resize( nx*ny );
	flux_z.resize( nx*ny );
	flux_x_subl.resize( nx*ny );
	flux_y_subl.resize( nx*ny );
	flux_z_subl.resize( nx*ny );

	InitializeNodes(z_readMatr);
	ConstructElements();
}

void SnowDriftA3D::Destroy() {}

double SnowDriftA3D::getTiming() const
{
	return timer.getElapsed();
}

SnowDriftA3D::~SnowDriftA3D()
{
	Destroy();
}

/**CompleteNodes
 * @brief Initialize concentration according to the rain rate (uniform)
 * Now calculate the missing parameters
 * Slope in direction of the wind
 */
void SnowDriftA3D::CompleteNodes()
{
	if ( new_wind_status == true ) {

		/* The turbulent settling velocity for the whole domain is now set to WS0 (formerly it
			contained the */
		/* Look at subscale turbulence effect */
		nodes_wstar.grid3D.resize(nx,ny,nz, WS0);

		// LH_CHANGE BEGIN: re-definition of the arps-mesh; Adapted by ML on August 12 2006
		cout << "[i] Snowdrift: updating the mesh..."<<endl;

		// the second arps layer (topography) becomes the bottom (ie boundary)
		// layer of the suspension grid, calculate the slope first
		for (size_t ii=0; ii<nx; ii++){
			for (size_t jj=0; jj<ny; jj++){
				nodes_slope(ii,jj,0) = (nodes_v(ii,jj,1)*atan(-nodes_sy(ii,jj,1)) + nodes_u(ii,jj,1)*atan(-nodes_sx(ii,jj,1)))
				/sqrt(Optim::pow2(nodes_v(ii,jj,1)) + Optim::pow2(nodes_u(ii,jj,1)));
				nodes_slope(ii,jj,1) = nodes_slope(ii,jj,0);
				nodes_u(ii,jj,0) = 0.;
				nodes_v(ii,jj,0) = 0.;
				nodes_w(ii,jj,0) = 0.;
				nodes_wstar(ii,jj,0) = nodes_wstar(ii,jj,1);
				nodes_e(ii,jj,0) = nodes_e(ii,jj,1);
			}
		}

		//an additional suspension grid layer of nodes is added between the
		//first and second arps layer

		//the new layer is located between the first and second
		//arps layer, adjusted by the factor auxLayerHeight. The x,y
		//coordinates remain the same as long as the mesh remains
		//regular in these directions

		for (size_t ii=0; ii<nx; ii++){
			for (size_t jj=0; jj<ny; jj++){
				const double salt_height = nodes_z(ii,jj,1) - nodes_z(ii,jj,0);

				//the direction of the wind field of a node in the new layer is
				//assumed to be equal to the node lying in the next higher layer above.
				//and the magnitude is scaled by a factor according to a logarithmic wind profile
				const double fac = log( salt_height/z0 ) / log( salt_height/(auxLayerHeight * z0) );

				nodes_u(ii,jj,1) = fac*nodes_u(ii,jj,2);
				nodes_v(ii,jj,1) = fac*nodes_v(ii,jj,2);
				nodes_w(ii,jj,1) = fac*nodes_w(ii,jj,2);

				//the turbulent kinetic energy of the new layer is given by a
				//weighted average
				nodes_e(ii,jj,1) = (1.-auxLayerHeight) *nodes_e(ii,jj,0) + auxLayerHeight * nodes_e(ii,jj,2);
				//nodes_slope(ii,jj,1) = nodes_slope(ii,jj,2);
			}
		}
	} // End if new wind
}

/**
 * @brief Initialize nodes
 * Each (3D) node receives its (x,y,z), slope, aspect, normals
 */
void SnowDriftA3D::InitializeNodes(const mio::Grid3DObject& z_readMatr)
{
	//we create a DEMObject using the Corripio slope algorithm
	const double cellsize = z_readMatr.cellsize;

	nodes_z = z_readMatr; //full copy
	nodes_x.set(z_readMatr, 0.);
	nodes_y.set(z_readMatr, 0.);
	nodes_sx.set(z_readMatr, 0.);
	nodes_sy.set(z_readMatr, 0.);
	nodes_e.set(z_readMatr, 0.);
	nodes_c.set(z_readMatr, 0.);
	nodes_Tair.set(z_readMatr, 0.);
	nodes_slope.set(z_readMatr, 0.);

	if (SUBLIMATION) {
		nodes_RH.set(z_readMatr, 0.); //relative humidity
		nodes_q.set(z_readMatr, 0.);
		nodes_Subl.set(z_readMatr, 0.);

		nodes_tmp_c.set(z_readMatr, 0.);
		nodes_Tair_ini.set(z_readMatr, 0.);
		nodes_q_ini.set(z_readMatr, 0.);
		nodes_Subl_ini.set(z_readMatr, 0.);
	}
	
	for (unsigned int kk=0; kk<nz; kk++){
		for (unsigned int jj=0;jj<ny;jj++){
			for (unsigned int ii=0;ii<nx;ii++){
			nodes_x(ii,jj,kk) = ii*cellsize;
			nodes_y(ii,jj,kk) = jj*cellsize;
			}
		}
	}

	DEMObject curr_layer(ta, false, DEMObject::CORR); //take the geolocalization from ta
	for (unsigned int kk=0; kk<4; kk++){
		nodes_z.extractLayer(kk, curr_layer);
		curr_layer.update(); //compute the slopes, normals, min/max, etc only for the bottom 3 layers
		for (unsigned int jj=0;jj<ny;jj++){
			for (unsigned int ii=0;ii<nx;ii++){
				nodes_sx(ii,jj,kk) = curr_layer.Nx(ii,jj);
				nodes_sy(ii,jj,kk) = curr_layer.Ny(ii,jj);
			}
		}
	}

	// LH_CHANGE BEGIN: re-definition of the arps-mesh
	cout <<"[i] Snowdrift adding artificial layer..."<<endl;

	// the second arps layer (topography) becomes the bottom (ie boundary)
	// layer of the suspension grid
	for (unsigned int ii=0;ii<nx;ii++){
		for (unsigned int jj=0;jj<ny;jj++){
			nodes_z(ii,jj,0) = z_readMatr(ii,jj,1);
			nodes_x(ii,jj,0) = nodes_x(ii,jj,1);
			nodes_y(ii,jj,0) = nodes_y(ii,jj,1);
			nodes_sx(ii,jj,0) = nodes_sx(ii,jj,1);
			nodes_sy(ii,jj,0) = nodes_sy(ii,jj,1);
		}
	}
	//an additional suspension grid layer of nodes is added between the
	//first and second arps layer
	//the new layer is located between the first and second
	//arps layer, adjusted by the factor auxLayerHeight. The x,y
	//coordinates remain the same as long as the mesh remains
	//regular in these directions

	for (unsigned int ii=0;ii<nx;ii++){
	    for (unsigned int jj=0;jj<ny;jj++){
		const double salt_height = auxLayerHeight*(nodes_z(ii,jj,2) - nodes_z(ii,jj,0));
		if (salt_height<0.) std::cout << "[E] Invalid height ARPS data at (" << ii << "," << jj << ") = " << salt_height << "\n";
		nodes_z(ii,jj,1) = nodes_z(ii,jj,0)+salt_height;
		nodes_x(ii,jj,1) = nodes_x(ii,jj,2);
		nodes_y(ii,jj,1) = nodes_y(ii,jj,2);
		nodes_sx(ii,jj,1) = nodes_sx(ii,jj,2);
		nodes_sy(ii,jj,1) = nodes_sy(ii,jj,2);
	    }
	}
}

void SnowDriftA3D::ConstructElements()
{
	//arrays for mapping local node indices onto global ones
	dofMap.resize( (nx-1)*(ny-1)*(nz-1) + 1, 8);
	/*
	* the dofMap matrix is used together with c and Psi, i.e. with vectors
	* that contain only the degrees of freedom; for building up dofMap we take
	* all elements of the domain but the application does only make
	* sense for nodes that are on the interior of the domain; we go from
	* 1 to (nx-1)*(ny-1)*(nz-1) elements and assign them their nodes
	* (from 1 to (nx-2)*(ny-2)*(nz-2) nodes)
	*/
	for (unsigned int kk = 0; kk < (nz-1); kk++ ) {
		for (unsigned int jj = 0; jj < (ny-1); jj++ ) {
			for (unsigned int ii = 1; ii <= (nx-1); ii++ ) {
			const int elem = kk*(nx-1)*(ny-1) + jj*(nx-1) + ii;

			dofMap(elem,0) = (jj-1)*(nx-2) + ii-1 + (kk-1)*(nx-2)*(ny-2);
			dofMap(elem,1) = (jj-1)*(nx-2) + (ii+0) + (kk-1)*(nx-2)*(ny-2);
			dofMap(elem,2) = (jj)*(nx-2) + (ii+0) + (kk-1)*(nx-2)*(ny-2);
			dofMap(elem,3) = (jj)*(nx-2) + ii-1 + (kk-1)*(nx-2)*(ny-2);
			dofMap(elem,4) = (jj-1)*(nx-2) + ii-1 + kk*(nx-2)*(ny-2);
			dofMap(elem,5) = (jj-1)*(nx-2) + (ii+0) + kk*(nx-2)*(ny-2);
			dofMap(elem,6) = (jj)*(nx-2) + (ii+0) + kk*(nx-2)*(ny-2);
			dofMap(elem,7) = (jj)*(nx-2) + ii-1 + kk*(nx-2)*(ny-2);
			}
		}
	}


	/*
	* nodeMap
	* this is the nodeMap matrix with all the elements from 1 to
	* (nx-1)*(ny-1)*(nz-1) and a mapping on all nodes of the domain note
	* that this approach is good for the enumeration of the real ARPS
	* mesh with the underground layer and is also good shit for the mesh
	* used for the calculation, i.e. without underground but with the
	* artificial nodes layer
	*/
	nodeMap.resize( (nx-1)*(ny-1)*(nz-1), 8);
	int element = 0;
	for (unsigned int iz=0; iz<(nz-1); iz++) {
		for (unsigned int iy=0; iy<(ny-1); iy++) {
			for (unsigned int ix=0; ix<(nx-1); ix++) {
				nodeMap(element,0) = ix + iy*nx + iz*nx*ny;
				nodeMap(element,1) = nodeMap(element,0) + 1;
				nodeMap(element,2) = nodeMap(element,0) + nx+1;
				nodeMap(element,3) = nodeMap(element,0) + nx;
				nodeMap(element,4) = nodeMap(element,0) + nx*ny;
				nodeMap(element,5) = nodeMap(element,4) + 1;
				nodeMap(element,6) = nodeMap(element,4) + nx+1;
				nodeMap(element,7) = nodeMap(element,4) + nx;
				element++;
			}
		}
	}


	elems.resize((nx-1)*(ny-1)*(nz-1), 8);
	int nelems=0;
	for (unsigned int iz=0; iz<(nz-1); iz++) {
		for (unsigned int iy=0; iy<(ny-1); iy++)  {
			for (unsigned int ix=0; ix<(nx-1); ix++) {
				elems(nelems,0) = ix + iy*nx + iz*nx*ny;
				elems(nelems,4) = elems[nelems][0] + nx*ny;
				elems(nelems,1) = elems[nelems][0] + 1;
				elems(nelems,3) = elems[nelems][0] + nx;
				elems(nelems,2) = elems[nelems][0] + nx+1;
				elems(nelems,5) = elems[nelems][4] + 1;
				elems(nelems,7) = elems[nelems][4] + nx;
				elems(nelems,6) = elems[nelems][4] + nx+1;
				++nelems;
			}
		}
	}
}


void SnowDriftA3D::SnowMassChange(bool setbound, const mio::Date& calcDate)
{
	//LH_DEBUG BEGIN
	static int timestepp=0;
	timestepp++;
	
	//loop over all interior nodes, boundary nodes are set to zero
	//subsequently
	for (unsigned int iy = 1; iy < ny-1; iy++) {
		for (unsigned int ix = 1; ix < nx-1; ix++) {
			dif_mns_subl(ix,iy)=0.; //reset
			if (!SALTATION) {
				mns(ix,iy) = 0.;
				continue; //saltation is not computed, only initiliazing 
			}

			// second ARPS layer is topography
			// Now first layer is on the ground and has a ground parallel velocity as Henning has tested
			// n = iy*nx + ix;

			// Enumerate the elements and nodes used to calculate saltation;
			// units of saltation[i][j] are kg/m/s
			double salt[9];
			salt[0] = salt[4] = saltation(ix-1,iy);
			salt[1] = saltation(ix,iy-1);
			salt[2] = saltation(ix+1,iy);
			salt[3] = saltation(ix,iy+1);
			salt[5] = saltation(ix-1,iy-1);
			salt[6] = saltation(ix+1,iy-1);
			salt[7] = saltation(ix+1,iy+1);
			salt[8] = saltation(ix-1,iy+1);

			Array2D<int> grid_ixiy(10,2);//tested Array2D<int> grid_ixiy(2,10);

			grid_ixiy(0,0)=ix-1;
			grid_ixiy(0,1)=iy;

			grid_ixiy(1,0)=ix;
			grid_ixiy(1,1)=iy-1;

			grid_ixiy(2,0)=ix+1;
			grid_ixiy(2,1)=iy;

			grid_ixiy(3,0)=ix;
			grid_ixiy(3,1)=iy+1;

			grid_ixiy(4,0)=ix-1;
			grid_ixiy(4,1)=iy;

			grid_ixiy(5,0)=ix;
			grid_ixiy(5,1)=iy-1;

			grid_ixiy(6,0)=ix-1;
			grid_ixiy(6,1)=iy-1;

			grid_ixiy(7,0)=ix+1;
			grid_ixiy(7,1)=iy-1;

			grid_ixiy(8,0)=ix+1;
			grid_ixiy(8,1)=iy+1;

			grid_ixiy(9,0)=ix-1;
			grid_ixiy(9,1)=iy+1;

			const double DX=nodes_x.grid3D(grid_ixiy(2,0),grid_ixiy(2,1),1)-nodes_x.grid3D(grid_ixiy(1,0),grid_ixiy(1,1),1);
			const double DY=nodes_y.grid3D(grid_ixiy(2,0),grid_ixiy(2,1),1)-nodes_y.grid3D(grid_ixiy(1,0),grid_ixiy(1,1),1);
			const double DZ=nodes_z.grid3D(grid_ixiy(2,0),grid_ixiy(2,1),1)-nodes_z.grid3D(grid_ixiy(1,0),grid_ixiy(1,1),1);

			const double DX1=nodes_x.grid3D(grid_ixiy(0,0),grid_ixiy(0,1),1)-nodes_x.grid3D(grid_ixiy(1,0),grid_ixiy(1,1),1);
			const double DY1=nodes_y.grid3D(grid_ixiy(0,0),grid_ixiy(0,1),1)-nodes_y.grid3D(grid_ixiy(1,0),grid_ixiy(1,1),1);
			const double DZ1=nodes_z.grid3D(grid_ixiy(0,0),grid_ixiy(0,1),1)-nodes_z.grid3D(grid_ixiy(1,0),grid_ixiy(1,1),1);

			double f_cross[3];
			f_cross[0] = DY*DZ1-DZ*DY1;  /* cross product of DX(vector) and DX1(vector) */
			f_cross[1] = DZ*DX1-DX*DZ1;
			f_cross[2] = DX*DY1-DY*DX1;

			const double area = sqrt( Optim::pow2(f_cross[0]) + Optim::pow2(f_cross[1]) + Optim::pow2(f_cross[2]) );

			/* area of the element is the norm of the cross product */
			double salt_flux = 0.0; /* The four contributions to erosion / deposition at a node */

			//now the edge loop in order to calculate the contribution
			//of each edge to the final salt_flux of the element
			double salt_flux_frac[4];
			
			double fac1[4];
			for (size_t i = 0; i < 4; i++) {
				//4ptdiv
				const double u = 0.25 * ( nodes_u.grid3D(ix,iy,1) + nodes_u.grid3D(grid_ixiy(i,0),grid_ixiy(i,1),1) + nodes_u.grid3D(grid_ixiy(i+1,0),grid_ixiy(i+1,1),1) +  nodes_u.grid3D(grid_ixiy(i+6,0),grid_ixiy(i+6,1),1));
				const double v = 0.25 * ( nodes_v.grid3D(ix,iy,1) + nodes_v.grid3D(grid_ixiy(i,0),grid_ixiy(i,1),1) + nodes_v.grid3D(grid_ixiy(i+1,0),grid_ixiy(i+1,1),1) +  nodes_v.grid3D(grid_ixiy(i+6,0),grid_ixiy(i+6,1),1));
				const double w = 0.25 * ( nodes_w.grid3D(ix,iy,1) + nodes_w.grid3D(grid_ixiy(i,0),grid_ixiy(i,1),1) + nodes_w.grid3D(grid_ixiy(i+1,0),grid_ixiy(i+1,1),1) +  nodes_w.grid3D(grid_ixiy(i+6,0),grid_ixiy(i+6,1),1));

				// vector from node i to i+1
				const double dx = ( nodes_x.grid3D(grid_ixiy(i+1,0),grid_ixiy(i+1,1),1) - nodes_x.grid3D(grid_ixiy(i,0),grid_ixiy(i,1),1));
				const double dy = ( nodes_y.grid3D(grid_ixiy(i+1,0),grid_ixiy(i+1,1),1) - nodes_y.grid3D(grid_ixiy(i,0),grid_ixiy(i,1),1));
				const double dz = ( nodes_z.grid3D(grid_ixiy(i+1,0),grid_ixiy(i+1,1),1) - nodes_z.grid3D(grid_ixiy(i,0),grid_ixiy(i,1),1));

				// only valid if grid is rectangular, take vector from
				// node i+1 to node i+2 as the normal vector on (dx,dy,dz)
				const double dx1 = ( nodes_x.grid3D(grid_ixiy(i+2,0),grid_ixiy(i+2,1),1) - nodes_x.grid3D(grid_ixiy(i+1,0),grid_ixiy(i+1,1),1) );
				const double dy1 = ( nodes_y.grid3D(grid_ixiy(i+2,0),grid_ixiy(i+2,1),1) - nodes_y.grid3D(grid_ixiy(i+1,0),grid_ixiy(i+1,1),1) );
				const double dz1 = ( nodes_z.grid3D(grid_ixiy(i+2,0),grid_ixiy(i+2,1),1) - nodes_z.grid3D(grid_ixiy(i+1,0),grid_ixiy(i+1,1),1) );

				// if there is no wind, set fac1 to zero
				// LH_REMARK: testing floats on equality? dangerous!
				if ( (u == 0.0) && (v == 0.0) && (w==0.0) ) {
					printf("Scheiss\n");
					    fac1[i] = 0;
				} else {
					//fac1=velocity(unitvec.) * lineunitvector that is perp.to
					//dx(vec)*abs(dx(vector))/area of the element;
					//the units of fac1 are 1/m, or equivalently: *m/m2 makes more physical sense
					fac1[i] = (u * dx1 + v * dy1 + w * dz1)
						/ sqrt( Optim::pow2(u) + Optim::pow2(v) + Optim::pow2(w) )
						/ sqrt( Optim::pow2(dx1) + Optim::pow2(dy1) + Optim::pow2(dz1) )
						/ area * sqrt( Optim::pow2(dx) + Optim::pow2(dy) + Optim::pow2(dz) );
				}
				salt_flux_frac[i] = fac1[i] * 1/4* (salt[i]+salt[i+1]+salt[i+5]+saltation(ix,iy));  /*flux over 1 of the four edges*/

				if (salt_flux_frac[i] > 0 && (salt[i+5] == 0)){
					// point at other side doesn't have saltation, block incoming flux
					salt_flux_frac[i] = 0;
				}
				if (salt_flux_frac[i] < 0 && (saltation(ix,iy) == 0)){
					// erosion necessary but there was no snow, block outgoing flux
					salt_flux_frac[i] = 0;
				}

				salt_flux+=salt_flux_frac[i];/*sum over the four edges*/

				// here we calculate the saltation flux on edge i by multiplication
				// of the saltation by the lenght of the effective edge and by division
				// of  the area of the element, therefore we obtain a flux in units of kg/m2/s

			} // end edge loop


			// scalar product of flux and normal unit vector to the surface, -ve
			// sign because the unit vector is in +ve z direction; units of flux_i[j]
			// is kg/m2/s , therefore also for diff_flux because we multiplicate
			// by unit vector f that is dimensionless ;
			/* The four contributions to erosion / deposition at a node */
			double diff_flux = -( flux_x[nx*iy + ix] * f_cross[0]
				+ flux_y[nx*iy + ix] * f_cross[1]
				+ flux_z[nx*iy + ix] * f_cross[2] ) / area;

			const double diff_flux_subl = (SUBLIMATION)? -( flux_x_subl[nx*iy + ix] * f_cross[0]
					+ flux_y_subl[nx*iy + ix] * f_cross[1] + flux_z_subl[nx*iy + ix] * f_cross[2] ) / area : 0.;

			mns(ix,iy) = dt_main * ( salt_flux + diff_flux );
			if (SUBLIMATION){
				mns_subl(ix,iy) = dt_main * ( salt_flux + diff_flux_subl);
				//cout<<"mns (nosubl)= " <<mns(ix,iy)<<" mns (subl)=" << mns_subl(ix,iy)<<endl;
			}
			/*the total flux is the sum of saltation and diff_flux*/

			//this test to limit crazy drift is normally done in SnowInterface, do this here to make sure to be able to compare sims with and without subl
			if (mns(ix,iy) > 37.) {
				//printf("Crazy drift: ix %d,  iy %d  mm/h %f\n", ix, iy, mns(ix,iy));
				mns(ix,iy) = 37.;
			}
			if (mns(ix,iy) < -37.) {
				//printf("Crazy drift: ix %d,  iy %d  mm/h %f\n", ix, iy, mns(ix,iy));
				mns(ix,iy) = -37.;
			}

			if (SUBLIMATION){
				if (mns_subl(ix,iy) > 37.) {
					printf("Crazy drift (subl): ix %d,  iy %d  mm/h %f\n", ix, iy, mns_subl(ix,iy));
					mns_subl(ix,iy) = 37.;
				}
				if (mns_subl(ix,iy) < -37.) {
					printf("Crazy drift (subl): ix %d,  iy %d  mm/h %f\n", ix, iy, mns_subl(ix,iy));
					mns_subl(ix,iy) = -37.;
				}

				dif_mns_subl(ix,iy)=mns(ix,iy)-mns_subl(ix,iy);
			}else{
				dif_mns_subl(ix,iy)=0.;
			}

			if (SUBLIMATION_OUTPUT){
				mns_nosubl(ix,iy)=mns(ix,iy);//just to save this in file
			}

			if (SUBLIMATION){
				//take value with sublimation as standard
				mns(ix,iy)=mns_subl(ix,iy);
				diff_flux=diff_flux_subl;
			}

 //			if ( (diff_flux < 0.0) && (saltation(ix,iy) == 0.0) && (salt_flux == 0.0)) {
// 				mns(ix,iy) = std::max(0.0,mns(ix,iy));
// 			}

#if WRITE_DRIFT_FLUXES
			const string fname = "../output/" + calcDate.toString(Date::NUM) + ".dep";
			FILE *difFile = fopen(fname.c_str(), "a");
			if (difFile == NULL)
				throw AccessException("Can not open file '"+fname+"'", AT);
			fprintf(difFile, "%s %d %d %f %f %f\n", calcDate.toString(Date::ISO).c_str(), ix, iy ,dt_main*salt_flux, dt_main*diff_flux, mns(ix,iy));
			fclose( difFile );
#endif
			if (fabs(mns(ix,iy)) > 1000.) {
				cout<<"Halt Salt:"<<ix<<" "<<iy<<" "<<salt_flux<<" "<<diff_flux<<endl;
			}
		}
	}/* End node loop */


	if (SUBLIMATION_OUTPUT){
		//Allows to obtain sublimation amount in 1 simulation > used for season simulation. NOTE: needs to be corrected for slope when comparing to SWE-grid
		const string fname= string("../results/")+ calcDate.toString(Date::NUM)+ string(".dif_subl");
		FILE *difFile = fopen(fname.c_str(), "w");
		if (difFile == NULL) {
			cout <<"Cannot open file "<<fname.c_str()<<endl;
			return;
		}
		for (unsigned int iy = 0; iy < ny; iy++){
			for (unsigned int ix = 0;ix < nx; ix++){
				fprintf(difFile, " %f", dif_mns_subl(ix,iy));
			}
			fprintf(difFile,"\n");
		}
		fclose(difFile);

		string fname2= string("../results/")+ calcDate.toString(Date::NUM)+ string(".mns_nosubl");
		FILE *mnsFile;
		mnsFile=fopen(fname2.c_str(), "w");
		if (mnsFile == NULL) {
			cout <<"Cannot open file "<<fname2.c_str()<<endl;
			return;
		}

		for (unsigned int iy = 0; iy < ny; iy++){
			for (unsigned int ix = 0;ix < nx; ix++){
				fprintf(mnsFile, " %f", mns_nosubl(ix,iy));
			}
			fprintf(mnsFile,"\n");
		}
		fclose(mnsFile);
	}

	// Now set border values to zero
	if ( setbound ) {
		for (unsigned int iy = 0 ; iy < ny; iy++ ) {
			mns(0,iy) = 0.0;
			mns(nx-1,iy) = 0.0;
			mns(1,iy) = 0.0;
			mns(nx-2,iy) = 0.0;
		}
		for (unsigned int ix = 0; ix < nx; ix++ ) {
			mns(ix,0) = 0.0;
			mns(ix,ny-1) = 0.0;
			mns(ix,1) = 0.0;
			mns(ix,ny-2) = 0.0;
		}
	}
}

/**
 * @brief Main: Calls the essential routines
 * @param calcDate date of current time step
 */
void SnowDriftA3D::Compute(const Date& calcDate)
{
	timer.restart();
	
	const double max_hs = cH.grid2D.getMax();
	const double max_psum = psum.grid2D.getMax();
	const double min_psum_ph = psum_ph.grid2D.getMin();
	if (max_hs==0. && (max_psum==0. || min_psum_ph>=1.)) {
		cout<<"[i] SnowDrift: no snow on the ground & no solid precipitation\n";
		mns=0.;
		if (snowpack!=NULL) {
			snowpack->setSnowMassChange(mns, calcDate);
		}
		return;
	}

	/* Calculate the Saltation Fluxes */
	cout<<"[i] SnowDrift starting Saltation\n";
	if (SALTATION) compSaltation(true);

	/* Calculate the suspension drift for the wind field until stationary */
	cout<<"[i] SnowDrift starting Suspension\n";
	//  DEBUG("***Checksum BEFORE Suspension: nodes=%lf, nodes[c]=%lf", checksum(nodes), checksum_c(nodes_c));
	Suspension();
	//  DEBUG("***Checksum after Suspension: nodes=%lf, nodes[c]=%lf", checksum(nodes), checksum_c(nodes_c));

	if (FIELD3D_OUTPUT){
 		debugOutputs(calcDate, string("../results/"), OUT_CONC);
 	}

	SnowMassChange(true, calcDate);

	std::cout << "[i] SnowDrift simulation done for " << calcDate.toString(Date::ISO) << "\n";

	if (snowpack!=NULL) {
		snowpack->setSnowMassChange(mns, calcDate);
	}

	//DEBUG("Checksum SnowDrift: nodes=%lf, saltation=%lf, c_salt=%lf", checksum_c(nodes_c), checksum(saltation), checksum(c_salt));
	timer.stop();

}  /* End SnowDrift */

void SnowDriftA3D::GetTResults(double  outtime_v[15], double outtime_tau[15], double outtime_salt[15], double outtime_diff[15])
{
	const size_t sz=15*sizeof(double);
	memcpy(outtime_v, time_v, sz);
	memcpy(outtime_tau, time_tau, sz);
	memcpy(outtime_salt, time_salt, sz);
	memcpy(outtime_diff, time_diff, sz);
}

void SnowDriftA3D::setSnowPack(SnowpackInterface &mysnowpack)
{
	snowpack=&mysnowpack;
}

void SnowDriftA3D::setEnergyBalance(EnergyBalance &myeb)
{
	eb=&myeb;
}

void SnowDriftA3D::setSnowSurfaceData(const mio::Grid2DObject& cH_in, const mio::Grid2DObject& sp_in, const mio::Grid2DObject& rg_in,
                                      const mio::Grid2DObject& N3_in, const mio::Grid2DObject& rb_in)
{
	cH = cH_in;
	sp = sp_in;
	rg = rg_in;
	N3 = N3_in;
	rb = rb_in;
}

void SnowDriftA3D::debugOutputs(const Date& calcDate, const std::string& outpath, const DRIFT_OUTPUT& output_type)
{
	const std::string ext = (output_type==OUT_CONC)? ".ctn" : ".sub"; 
	const std::string fname = outpath + calcDate.toString(Date::NUM) + ext;
	
	writeOutput(fname); //write output file of snowdrift
}


//HACK: this should be done by MeteoIO
/**
 * @brief Write output
 * Writes the values of several nodes-fields
 */
void SnowDriftA3D::writeOutput(const std::string& fname)
{
	FILE *logfile = fopen(fname.c_str(), "w");
	if (logfile == NULL) {
		cout<<"Cannot open file"<<fname.c_str()<<endl;
		return;
	}
	fprintf(logfile, "#1:j 2:ix 3:iy 4:iz \t 5:z 6:c");
	if (SUBLIMATION_OUTPUT) fprintf(logfile, " \t 7:sublimation 8:RH 9:q 10:Ta");
	fprintf(logfile, " \t 11:u 12:v 13:w \n");
	
	for (unsigned int iy=0; iy<ny; iy++) {
		for (unsigned int ix=0; ix<nx; ix++) {
		    //compute int dz c
		    double intc=0;
		    for (unsigned int iz=0; iz<nz-1; iz++) {
			const int j = iz*nx*ny + iy*nx + ix;
			    intc += c[j]*(nodes_z.grid3D(ix,iy,iz+1)-nodes_z.grid3D(ix,iy,iz));
		    }
		    intc /= (nodes_z.grid3D(ix,iy,nz-1)-nodes_z.grid3D(ix,iy,0));
		    
		    for (unsigned int iz=0; iz<nz; iz++) {
			const int j = iz*nx*ny + iy*nx + ix;
			fprintf(logfile, "%d %d %d %d \t %f %f", j, ix, iy, iz, nodes_z.grid3D(ix,iy,iz), c[j]*1000);
			
			if (SUBLIMATION_OUTPUT)
			fprintf(logfile, " \t %f %f %f %f",
				nodes_Subl.grid3D(ix,iy,iz)*1000, nodes_RH.grid3D(ix,iy,iz), nodes_q.grid3D(ix,iy,iz)*1000, nodes_Tair.grid3D(ix,iy,iz));
			
			fprintf(logfile, " \t %f %f %f\n", nodes_u.grid3D(ix,iy,iz), nodes_v.grid3D(ix,iy,iz), nodes_w.grid3D(ix,iy,iz));
		    }
		}
	}
	fclose (logfile);
}

/**
* @brief Sets the required meteo fields
*/
void SnowDriftA3D::setMeteo (const unsigned int& steps, const Grid2DObject& new_psum, const mio::Grid2DObject& new_psum_ph, const Grid2DObject& new_p, const Grid2DObject& /*new_vw*/,
                          const Grid2DObject& new_rh, const Grid2DObject& new_ta, const Grid2DObject& new_tsg, const Grid2DObject& new_ilwr, const mio::Date& calcDate,
                          const std::vector<mio::MeteoData>& vecMeteo)
{
	if (vecMeteo.empty())
		throw mio::NoDataException("No meteo data!", AT);
	
	// find the first (often the only one) MeteoData having iswr, ta and ilwr
	// it should exist, otherwise the checkInputsRequirements method would have failed.
	bool found = false;
	for (size_t i = 0; i < vecMeteo.size() ; ++i) {
		ta_1D = vecMeteo[i](MeteoData::TA);
		station_altitude = vecMeteo[i].meta.position.getAltitude();
		if (ta_1D != IOUtils::nodata && station_altitude!=IOUtils::nodata) {
			found = true;
			break;
		}
	}
	if (!found)
		throw mio::NoDataException("No TA station could be found at "+vecMeteo[0].date.toString(mio::Date::ISO), AT);
	
	//Common meteo fields
	psum = new_psum;
	psum_ph = new_psum_ph;
	rh = new_rh;
	ta = new_ta;
	tsg = new_tsg;
	p = new_p;
	new_wind_status = isNewWindField(steps);

	if (new_wind_status) {
		wind_field_index++;
		const std::string filename( wind_fields[wind_field_index].wind );
		io.read3DGrid(nodes_u, filename+":u");
		io.read3DGrid(nodes_v, filename+":v");
		io.read3DGrid(nodes_w, filename+":w");
		if (READK) io.read3DGrid(nodes_K, filename+":kmh");
		cout <<"[i] Snowdrift: ARPS wind field successfully read"<<endl;
	}

	mio::Grid2DObject dw(vw, IOUtils::nodata);	// dw field with vw as template for dimensions
	//adjust the meteo fields that depend on the 3D wind field
	for (size_t jj=0; jj<ny; jj++) {
		for (size_t ii=0; ii<nx; ii++) {
			vw.grid2D(ii,jj) = Optim::fastSqrt_Q3( Optim::pow2(nodes_u.grid3D(ii,jj,2)) + Optim::pow2(nodes_v.grid3D(ii,jj,2)) ); //Third layer is first layer in the air
			dw.grid2D(ii,jj) = (atan2(nodes_u.grid3D(ii,jj,2), nodes_v.grid3D(ii,jj,2))) * mio::Cst::to_deg;
		}
	}
	
	CompleteNodes();
	if (SUBLIMATION) 	initializeTRH();
	
	//TODO: feedback mecanism: make it more general!
	if (snowpack!=NULL) snowpack->setMeteo(psum, psum_ph, vw, dw, rh, ta, tsg, calcDate);
	if (eb!=NULL) eb->setMeteo(new_ilwr, ta, rh, p, calcDate);
}

/**
* @brief Initialize RH, T and q
* Initialize fields of humidity and temperature, necessary for sublimation calculations
*/
void SnowDriftA3D::initializeTRH()
{
	//start with first level
	for (unsigned int ix=0; ix<nx; ix++){
		for (unsigned int iy=0; iy<ny; iy++){
			nodes_RH.grid3D(ix,iy,0) = rh.grid2D(ix,iy); //constant field, single measurement at Wan3
			nodes_Tair.grid3D(ix,iy,0) = ta_1D;
			nodes_q.grid3D(ix,iy,0) = Atmosphere::relToSpecHumidity(nodes_z.grid3D(ix,iy,0), nodes_Tair.grid3D(ix,iy,0),nodes_RH.grid3D(ix,iy,0));
			nodes_Tair.grid3D(ix,iy,0) = ta_1D-(nodes_z.grid3D(ix,iy,0) - station_altitude)*Cst::dry_adiabatique_lapse_rate;
			nodes_RH.grid3D(ix,iy,0) = RH_from_q(nodes_Tair.grid3D(ix,iy,0),nodes_q.grid3D(ix,iy,0), nodes_z.grid3D(ix,iy,0));
		}
	}

	//adjust rest of domain
	for (unsigned int iz=1; iz<nz; iz++){
		for (unsigned int iy=0; iy<ny; iy++){
			for (unsigned int ix=0;ix<nx; ix++){
			nodes_q.grid3D(ix,iy,iz)=nodes_q.grid3D(ix,iy,iz-1);
			nodes_Tair.grid3D(ix,iy,iz)=ta_1D-(nodes_z.grid3D(ix,iy,iz) - station_altitude)*Cst::dry_adiabatique_lapse_rate;
			nodes_RH.grid3D(ix,iy,iz)=RH_from_q(nodes_Tair.grid3D(ix,iy,iz),nodes_q.grid3D(ix,iy,iz), nodes_z.grid3D(ix,iy,iz));
			}
		}
	}
	nodes_Tair_ini=nodes_Tair;
	nodes_q_ini=nodes_q;
}

#pragma GCC diagnostic pop
