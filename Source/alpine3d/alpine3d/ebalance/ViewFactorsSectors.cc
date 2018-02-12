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
#include <alpine3d/ebalance/ViewFactorsSectors.h>

ViewFactorsSectors::ViewFactorsSectors(const mio::Config& cfg, const mio::DEMObject &dem_in) : dem(dem_in)
{
        cellsize = dem.cellsize;
        dimx = dem.getNx();
        dimy = dem.getNy();

        sky_vf.resize(dimx, dimy, 1);

        cfg.getValue("sw_radius", "EBalance", sw_radius);

        hSections = 60;
        vSections = 30;

        schrittweitenLimit = 1;
        //max_shade_distance = (dem.grid2D.getMax() - dem.grid2D.getMin()) / tan(5.*M_PI/180.);

        max_shade_distance = std::numeric_limits<double>::max();

        calcHorizonField();
        fill_vf_map();
}

double ViewFactorsSectors::getSkyViewFactor(const int &i, const int &j)
{
        return sky_vf(i,j);
}

void ViewFactorsSectors::getSkyViewFactor(mio::Array2D<double> &o_sky_vf) const {
	o_sky_vf = sky_vf;
}

void ViewFactorsSectors::calcHorizonField()
{
	std::cout << "[i] computing sector view factors" << std::endl;

	//horizonMap.resize(dem.getNx(), dem.getNy(), hSections);

	const double alpha_factor = (2.*M_PI/hSections);

	viewCells.resize(dem.getNx(),dem.getNy(),hSections,vSections);

	for (unsigned int iy1=0; iy1<dem.getNy(); iy1++) {
		for (unsigned int ix1=0; ix1<dem.getNx(); ix1++) {

			double cell_alt = dem.grid2D(ix1, iy1);
			if (cell_alt==mio::IOUtils::nodata) continue;

			for (int z=0; z<hSections;z++) {
				std::vector<unsigned int> viewCellsTemp(vSections);
				getHorizonForRay(ix1,iy1, (z+0.5)*alpha_factor, viewCellsTemp);
				for (int i=0;i<vSections; i++) {
					viewCells(ix1,iy1,z,i) = viewCellsTemp[i];
				}
			}
		}
	}
}

void ViewFactorsSectors::calcSchrittweite(const double& altitude, const double& dH, const double& distance,
                                          const double& horizon_tan_angle, unsigned int& schrittweite) const
{
	const double dalpha = dH/(distance*schrittweite);

	if (!(dalpha == 0))
		schrittweite = (int)((horizon_tan_angle-altitude)/dalpha);
	else
		schrittweite = (int)((horizon_tan_angle-altitude)/2.);
	if (schrittweite == 0)
		schrittweite = 1;
	else if (schrittweite > schrittweitenLimit)
		schrittweite = schrittweitenLimit;
}


double ViewFactorsSectors::getHorizonForRay(const unsigned int& ix1, const unsigned int& iy1, const double& alpha,
                                            std::vector<unsigned int>& i_viewCells)
{
	bool horizon_found = false;
	int nb_cells = 1;
	const double cell_alt = dem.grid2D(ix1, iy1);
	double horizon_tan_angle = 0;
	double min_horizon_angle;
	unsigned int schrittweite = 1;
	double old_altitude = cell_alt;
	int old_sectorLatitudeAlpha=0;

	//update horizon_tan_angle
	if (dem.Nz(ix1,iy1) == 0) {
		horizon_tan_angle = -999.;
		std::cout << "fehler" << std::endl;
	} else {
		horizon_tan_angle = -1 * (dem.Nx(ix1,iy1)*sin(alpha)+dem.Ny(ix1,iy1)*cos(alpha))/dem.Nz(ix1,iy1);
	}
	min_horizon_angle = atan(horizon_tan_angle);

	if (ix1==0 || (int)ix1==dimx-1 || iy1==0 || (int)iy1==dimy-1) return 0.; //a border cell is not shadded

	int horzionHeightsCounter = 0;

	//if (ix1==35 && iy1 ==90)
	//       horzionMinHorizonAngle(horzionHeightsAngle,0) = min_horizon_angle;

	while (!horizon_found) {
		nb_cells+=schrittweite;
		const int ix2 = ix1 + (int)mio::Optim::round( ((double)(nb_cells))*sin(alpha) ); //alpha is a bearing
		const int iy2 = iy1 + (int)mio::Optim::round( ((double)(nb_cells))*cos(alpha) ); //alpha is a bearing

		if (ix2<=0 || ix2>=(int)dimx-1 || iy2<=0 || iy2>=(int)dimy-1) break; //we are out of the dem

		const double new_altitude = dem.grid2D(ix2, iy2);
		if (new_altitude==mio::IOUtils::nodata) {
			schrittweite = 1;
			continue; //we continue at nodata cells
		}

		const double distance = sqrt( (double)( (ix2-ix1)*(ix2-ix1) + (iy2-iy1)*(iy2-iy1)) ) * cellsize;
		calcSchrittweite(new_altitude, new_altitude - old_altitude, distance, horizon_tan_angle, schrittweite);
		old_altitude = new_altitude;


		const double tan_angle = (new_altitude-cell_alt)/distance;
		if (tan_angle>horizon_tan_angle) {
			horizon_tan_angle = tan_angle;

			//int new_sectorLatitudeAlpha = (int)(((double)((atan(tan_angle)-min_horizon_angle)*(2*vSections)))/M_PI);
			const int new_sectorLatitudeAlpha = (int)(((double)((atan(tan_angle)-min_horizon_angle)*(2*vSections)))/M_PI);

			for (int sectorLatitudeAlpha=old_sectorLatitudeAlpha;sectorLatitudeAlpha<new_sectorLatitudeAlpha;sectorLatitudeAlpha++) {
				if (sectorLatitudeAlpha<vSections) {
					if (sectorLatitudeAlpha >= 0) {
						i_viewCells[sectorLatitudeAlpha] = ix2+iy2*dimx;
					} else {
						std::cout << "Fehler: viewCell Index kleiner 0" << std::endl;
					}
				}
			}

			old_sectorLatitudeAlpha = new_sectorLatitudeAlpha;

		}

		horzionHeightsCounter++;
		if (distance>max_shade_distance) horizon_found=true; //maximum lookup distance reached
	}

	return horizon_tan_angle;
}

void ViewFactorsSectors::fill_vf_map()
{
	const double sw_radius2 = sw_radius * sw_radius;
	for ( int i_shoot = 0; i_shoot < dimx; i_shoot++ ) {
		std::cout << i_shoot * 100 / dimx << "%" << std::endl;
		for ( int j_shoot = 0; j_shoot < dimy; j_shoot++ ) {

			const double z_shoot=dem.grid2D(i_shoot,j_shoot);
			for ( int v = 0; v < vSections; v++ ) {
				const double theta_min = M_PI * (double)v / (double)(2 * vSections);
				const double theta_max = M_PI * (double)(v+1)/ (double)(2 * vSections);
				//const double psi = (2. * M_PI)/(double)hSections;

				for ( int h = 0; h < hSections; h++ ) {
					const unsigned int cell = viewCells(i_shoot, j_shoot, h, v);
					if (cell == 0) continue;

					const unsigned int i = cell % dimx;
					const unsigned int j = (int) (cell/dimx);

					const double bx = (i_shoot - i) * cellsize;	// bx = dx (m) going from (i_shoot,j_shoot) to (i,j)
					const double by = (j_shoot - j) * cellsize;	// by = dy (m) going from (i_shoot,j_shoot) to (i,j)
					const double bz = z_shoot - dem.grid2D(i,j);	// bz = dz (m) going from (i_shoot,j_shoot) to (i,j)

					// distance between the surfaces
					const double dist = sqrt( (bx * bx) + (by * by) + (bz * bz) );
					const double dist2 = (bx * bx) + (by * by) + (bz * bz);

					if (!( dist2 <= sw_radius2 && dist2>0.)) continue;

					const double viewFactor = getPrecViewFactor(cos(theta_min), cos(theta_max), dist);

					const unsigned int shootID = j_shoot+i_shoot*dimy;
					const unsigned int ID = j+i*dimy + shootID*dimx*dimy;

					sky_vf(i_shoot, j_shoot) -= viewFactor;

					if (vf.find(ID) == vf.end())
						vf.insert(std::pair<int,double>(ID, viewFactor));
					else
						vf.at(ID) += viewFactor;


				}
			}
		}
	}
}

double ViewFactorsSectors::GetViewfactor(const int i, const int j, const int a, const int b)
{
	const unsigned int shootID = j+i*dimy;
	const unsigned int ID = b+a*dimy + shootID*dimx*dimy;

	if (vf.find(ID) == vf.end())
		return 0.;
	return vf.at(ID);
}

int ViewFactorsSectors::getViewCells(const int i, const int j, const int h, const int v)
{
	return viewCells(i,j,h,v);
}


/**
        Opimizted viewfactor calculation
*/
double ViewFactorsSectors::getPrecViewFactor(const double &cos_theta_min, const double &/*cos_theta_max*/, const double &radius)
{

	double viewFactor=0;

	int sizeI = 1;
	int sizeJ = 1;

	if (radius > cellsize * 10) {
		sizeI = 1;
		sizeJ = 1;
	}


	double x = cellsize * -0.5 + cellsize / (double) sizeI * 0.5;

	for (double i=0; i<sizeI; i++) {
		double y = cellsize * -0.5 + cellsize / (double) sizeJ * 0.5;
		for (double j=0; j<sizeJ; j++) {

			const double a1 = (cos_theta_min*radius-x)*(cos_theta_min*radius-x);
			const double a2 = y*y;
			const double a3 = (a1+a2)*(a1+a2);
			const double a4 = radius*radius*(1-cos_theta_min*cos_theta_min);
			const double cos_theta_min_h =  sqrt((a1+a2)/(a3+a4));

			//const double b1 = (cos_theta_max*radius-x)*(cos_theta_max*radius-x);
			//const double b2 = y*y;
			//const double b3 = (a1+a2)*(a1+a2);
			//const double b4 = radius*radius*(1-cos_theta_max*cos_theta_max);
			const double cos_theta_max_h =  sqrt((a1+a2)/(a3+a4));

			const double psi = (2. * M_PI)/(double)hSections;
			const double psi_h = atan((radius*sin(psi/2)-y)/(radius*cos(psi/2)-x)) + atan((radius*sin(psi/2)+y)/(radius*cos(psi/2)-x));

			viewFactor += (cos(2.*acos(cos_theta_min_h))-cos(2.*acos(cos_theta_max_h)))*psi_h/(4.*M_PI)*cellsize*cellsize/(double)(sizeI*sizeJ);

			y +=  cellsize / (double)sizeJ;
		}
		x += cellsize / (double)sizeI;
	}
	return viewFactor/(cellsize*cellsize);
}
