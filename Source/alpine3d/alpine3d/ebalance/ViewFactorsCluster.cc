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
#include <alpine3d/ebalance/ViewFactorsCluster.h>

ViewFactorsCluster::ViewFactorsCluster(const mio::Config& cfg, const mio::DEMObject &dem_in) : dem(dem_in)
{
	cellsize = dem.cellsize;
	dimx = dem.getNx();
	dimy = dem.getNy();

	sky_vf.resize(dimx, dimy, 1);

	cfg.getValue("sw_radius", "EBalance", sw_radius);
	cfg.getValue("sub_crit", "EBalance", sub_crit);

	hSections = 60;
	vSections = 30;

	schrittweitenLimit = 1;
	//max_shade_distance = (dem.grid2D.getMax() - dem.grid2D.getMin()) / tan(5.*M_PI/180.);

	max_shade_distance = std::numeric_limits<double>::max();

	calcVF_cluster();
	fill_vf_map();
}

double ViewFactorsCluster::getSkyViewFactor(const int &i, const int &j)
{
	return sky_vf(i,j);
}

void ViewFactorsCluster::getSkyViewFactor(mio::Array2D<double> &o_sky_vf) const {
	o_sky_vf = sky_vf;
}

void ViewFactorsCluster::calcSchrittweite( const double &altitude, const double &dH, const double &distance, const double &horizon_tan_angle, unsigned int &schrittweite) const
{
	const double dalpha = dH/(distance*schrittweite);

	if (!(dalpha == 0)) {
		schrittweite = (int)((horizon_tan_angle-altitude)/dalpha);
	} else {
		schrittweite = (int)((horizon_tan_angle-altitude)/2.);
	}

	if (schrittweite == 0) {
		schrittweite = 1;
	} else if (schrittweite > schrittweitenLimit) {
		schrittweite = schrittweitenLimit;
	}
}


void ViewFactorsCluster::calcVF_cluster()
{
	std::cout << "[i] computing cluster view factors" << std::endl;

	vf_cluster.resize(dem.getNx(),dem.getNy(),hSections,vSections);
	vc_cluster.resize(dem.getNx(),dem.getNy(),hSections,vSections);

	for (int x=0;x<dimx; x++) {
		std::cout << x*100/dimx << "%" << std::endl;
		for (int y=0; y<dimy; y++) {
			//std::cout << "  " << y*100/dimy << "%" << std::endl;

			mio::Array1D<double> tan_h(hSections);

			for (unsigned int h=0; h<hSections; h++)
				tan_h(h) = -1 * (dem.Nx(x,y)*sin(h*2.*M_PI/hSections)+dem.Ny(x,y)*cos(h*2.*M_PI/hSections))/dem.Nz(x,y);

			bool cont = true;
			int r=1;
			while (cont) {
				cont = false;
				int i = x-r;
				int j = y-r;
				while (i<=x+r) {
					cont = VF_calc(x,y,i,j,tan_h) || cont;
					i++;
				}
				i--;
				while (j<y+r) {
					j++;
					cont = VF_calc(x,y,i,j,tan_h) || cont;
				}
				while (i>x-r) {
					i--;
					cont = VF_calc(x,y,i,j,tan_h) || cont;
				}
				while (j>y-r+1) {
					j--;
					cont = VF_calc(x,y,i,j,tan_h) || cont;
				}
				r++;
			}
		}
	}
}

bool ViewFactorsCluster::VF_calc(const unsigned int ix1, const unsigned int iy1, const int ix2, const int iy2, mio::Array1D<double> &tan_h)
{
	//std::cout << ix1 << " " << iy1 << " " <<  ix2 << " " << iy2 << std::endl;
	if (ix2<0 || ix2>dimx-1 || iy2<0 || iy2>dimy-1) return false;

	const double cell_alt = dem.grid2D(ix1,iy1);
	const double new_altitude = dem.grid2D(ix2,iy2);
	if (new_altitude==mio::IOUtils::nodata) return true;
	const double distance = sqrt( (double)( (ix2-ix1)*(ix2-ix1) + (iy2-iy1)*(iy2-iy1)) ) * cellsize;

	const double tan_angle = (new_altitude-cell_alt)/distance;

	double alpha = atan2((iy2-(int)iy1),(ix2-(int)ix1));
	if (alpha < 0)
		alpha += 2.*M_PI;
	const double dalpha = atan((cellsize/2./distance));

	//std::cout << "alpha: " << alpha << std::endl;

	unsigned int bottom = (unsigned int)round(alpha*hSections/(2.*M_PI));
	if (bottom>=hSections)
		bottom = hSections -1 ;
	const unsigned int top = (unsigned int)round((alpha+dalpha)*hSections/(2.*M_PI));

	//std::cout << "bottom: " << bottom << std::endl;

	const double horizon_tan_angle = (dem.Nx(ix1,iy1)*sin(alpha)+dem.Ny(ix1,iy1)*cos(alpha))/dem.Nz(ix1,iy1);
	const double horizon_angle = atan(horizon_tan_angle);


	if (tan_angle>tan_h(bottom)) {
		unsigned int j = (unsigned int)round((atan(tan_angle)+horizon_angle)*vSections/(M_PI/2.));

		if (j >= vSections)
			j = vSections-1;

		vf_cluster(ix1, iy1, bottom, j) += GetSymetricPartOfViewFactor(ix1,iy1,ix2,iy2) / (EuclidianDistance(dem.Nx(ix1,iy1), dem.Ny(ix1,iy1), 1.) * cellsize * cellsize);
		vc_cluster(ix1, iy1, bottom, j) = ix2*dimx+iy2;
	}
	for (unsigned int i=bottom; i<top; i++) {
		if (tan_angle>tan_h(i)) tan_h(i) = tan_angle;
	}

	return true;
}

void ViewFactorsCluster::fill_vf_map()
{
	std::cout << "[i] computing cluster view factors.. fill map" << std::endl;
	for (int i=0;i<dimx;i++) {
		for (int j=0;j<dimy;j++) {
			for (unsigned int h=0; h<hSections; h++) {
				for (unsigned int v=0; v<vSections; v++) {
					const unsigned int shootID = j+i*dimy;
					const unsigned int ID = vc_cluster(i,j,h,v) + shootID*dimx*dimy;

					const double viewFactor = vf_cluster(i,j,h,v);

					sky_vf(i, j) -= viewFactor;

					vf.insert(std::pair<int,double>(ID, viewFactor));
				}
			}
		}
	}
}


double ViewFactorsCluster::GetViewfactor(const int i, const int j, const int a, const int b)
{
	const unsigned int shootID = j+i*dimy;
	const unsigned int ID = b+a*dimy + shootID*dimx*dimy;

	if (vf.find(ID) == vf.end())
		return 0.;
	return vf.at(ID);
}

double ViewFactorsCluster::GetSymetricPartOfViewFactor(const int i, const int j, const int a, const int b)
{
	//Return the symetric part of the view factor (the view factor whitout the division of the area)
	//int splitij;                   // splitting factors for the respective grid cell dimension
	//int splitab;                   // splitting factors for the respective grid cell dimension
	//double area_uv, area_kl;       // small surface areas because of splitting uv in ij-plane and kl in ab-plane
	//double x_uv, y_uv, z_uv;       // vector components of the small cells in the ij-plane
	//double x_kl, y_kl, z_kl;       // vector components of the small cells in the ab-plane
	double cosangle_ij;            // angle between normal to [i,j] and beam from [i,j]
	double cosangle_ab;            // angle between normal to [a,b] and beam from [i,j]
	double add_vf;                 // for summing up the small view factors in the substructured grid cell
	double dist;                   // actual centercenter distance between every small pair of areas

	// x and y values of the normal vector and z value of a cell
	double sx_source = dem.Nx(i,j);
	double sy_source = dem.Ny(i,j);
	double z_source = dem.grid2D(i,j);

	// z-values of the edge vectors at the forced to be rectangular area ij.
	// The "official" z value is considered at the center of the cell, just indicate
	// the part of the cellsize must be move in x and y.
	// They are forced on the ij-plane around a centered coordinate.

	// z-values of the vectors in the center of the 3 areas around grid point (i,j)
	const double iza = GetZEdgeVector(z_source, sx_source, sy_source, -0.5, -0.5);
	const double izb = GetZEdgeVector(z_source, sx_source, sy_source, 0.5, -0.5);
	const double izc = GetZEdgeVector(z_source, sx_source, sy_source, -0.5, 0.5);
	//const double izd = GetZEdgeVector(z_source, sx_source, sy_source, 0.5, 0.5);

	// z-values of the vectors in the center of the 3 areas around grid point (a,b)
	const double aza = GetZEdgeVector(dem.grid2D(a,b), dem.Nx(a,b), dem.Ny(a,b), -0.5, -0.5);
	const double azb = GetZEdgeVector(dem.grid2D(a,b), dem.Nx(a,b), dem.Ny(a,b), 0.5, -0.5);
	const double azc = GetZEdgeVector(dem.grid2D(a,b), dem.Nx(a,b), dem.Ny(a,b), -0.5, 0.5);
	//const double azd = GetZEdgeVector(dem.grid2D(a,b), dem.Nx(a,b), dem.Ny(a,b), 0.5, 0.5);

	double lengthj;
	double lengthi;
	double lengtha;
	double lengthb;

	// actual grid cell dimensions of (i,j) considering slope angle
	if (iza != mio::IOUtils::nodata  && izb != mio::IOUtils::nodata)
		lengthi = EuclidianDistance(cellsize, 0, (izb - iza));
	else
		lengthi = cellsize;
	if (iza != mio::IOUtils::nodata  && izc != mio::IOUtils::nodata)
		lengthj = EuclidianDistance(0, cellsize, (izc - iza));
	else
		lengthj = cellsize;


	// actual grid cell dimensions of (a,b) considering slope angle
	if (azb != mio::IOUtils::nodata  && aza != mio::IOUtils::nodata)
		lengtha = EuclidianDistance(cellsize,0,(azb - aza));
	else
		lengtha = cellsize;
	if (azc != mio::IOUtils::nodata  && aza != mio::IOUtils::nodata)
		lengthb = EuclidianDistance(0,cellsize,(azc - aza));
	else
		lengthb = cellsize;

	// surface area = |n|, assumed that the 3D-area's are parallelograms:
	// |n| represents the inclined area for "1 unit" of x and y.
	// Multiply with the flat area give the inclined area for a cell.
	// large unsplitted surface areas
	const double area_ij = EuclidianDistance(sx_source, sy_source, 1.) * cellsize * cellsize;
	const double area_ab = EuclidianDistance(dem.Nx(a,b), dem.Ny(a,b), 1.) * cellsize * cellsize;

	// the following is valid also for inclined surfaces as the planes are enlarged
	// by forcing the corner z-values on the plane with sx, sy, sz at the lower left
	// corner of the grid cell
	// distances between the viewers large and the viewed large surface
	const double bx = (a - i) * cellsize;		// bx = dx (m) going from [a,b] to [i,j]
	const double by = (b - j) * cellsize;		// by = dy (m) going from [a,b] to [i,j]
	const double bz = dem.grid2D(a,b) - z_source;	// bz = dz (m) going from [a,b] to [i,j]

	// center distance between the large surfaces
	dist = EuclidianDistance(bx,by,bz);

	// substructuring if the side length of area ij divided by the distance exceeds
	// the preassumed value 'sub_crit'
	// McCluney (1994) assumes that the maximum side length should be lower than 10%
	// of the distance as a rough criteria until when a finite source can still be
	// treated as a point source (no substructuring is required)
	if ( lengthi / dist >= sub_crit || lengthj / dist >= sub_crit ||
		lengtha / dist >= sub_crit || lengthb / dist >= sub_crit ) {
		// computing the necessary splitting factors for all border lengths
		// and find and store the largest splitting factor
		const int splitij = std::max(GetSplitFactor(lengthi, dist), GetSplitFactor(lengthj, dist));
		const int splitab = std::max(GetSplitFactor(lengtha, dist),GetSplitFactor(lengthb, dist));

		// small surface areas in the subdivided areas ij and ab
		const double area_uv = area_ij / (splitij * splitij);
		const double area_kl = area_ab / (splitab * splitab);

		add_vf = 0.;                  // every small view factor needs to be calculated and summed

		//std::cout << splitij << " " << splitab << std::endl;
		//splitij = 5;
		//splitab = 5;

		for ( int u = 0; u < splitij; u++ ) {

			for ( int v = 0; v < splitij; v++ ) {
				// find the coordinate values of the small cell uv in the plane ij
				// vec(a) + vec(b)-vec(a) * (u + 0.5) / spliti
				//        + vec(c)-vec(a) * (v + 0.5) / splitj

				const double x_uv = GetCoordinateForSmallCell(i, i+1, i, u, v, splitij);
				const double y_uv = GetCoordinateForSmallCell(j, j, j+1, u, v, splitij);
				const double z_uv = GetCoordinateForSmallCell(iza, izb, izc, u, v, splitij);

				for ( int k = 0; k < splitab; k++ ) {

					//std::cout << k << std::endl;
					for ( int l = 0; l < splitab; l++ ) {

						// find the coordinate values of the small cell kl in the plane ab
						// vec(a) + vec(b)-vec(a) * (k + 0.5) / splita
						//        + vec(c)-vec(a) * (l + 0.5) / splitb
						const double x_kl = GetCoordinateForSmallCell(a, a+1, a, k, l, splitab);
						const double y_kl = GetCoordinateForSmallCell(b, b, b+1, k, l, splitab);
						const double z_kl = GetCoordinateForSmallCell(aza, azb, azc, k, l, splitab);

						dist = EuclidianDistance((x_kl - x_uv) * cellsize,(y_kl - y_uv)*cellsize,(z_kl - z_uv));

						// cosangle_ij: n * (b) / (|n| * |b|)
						// the angle between normal to area_ij and beam from area_uv
						cosangle_ij = CosAngleBetween2Vectors(sx_source, sy_source, 1,(x_kl - x_uv) * cellsize, (y_kl - y_uv) * cellsize, (z_kl - z_uv));

						// cosangle_ab: n * (-b) / (|n| * |b|)
						// the angle between normal to area_ab and beam from area_uv
						cosangle_ab = CosAngleBetween2Vectors(dem.Nx(a,b), dem.Ny(a,b), 1, (x_uv - x_kl) * cellsize, (y_uv - y_kl) * cellsize, (z_uv - z_kl));


						if ( cosangle_ij > 0. && cosangle_ab > 0. ) {
							add_vf += (cosangle_ij * cosangle_ab * area_uv * area_kl) / (M_PI * dist * dist);
						}
					}
				}
			}
		}
		//std::cout << "add_vf: " << add_vf << std::endl;
		return add_vf;
		// end of if condition : if a length of an area > 10 % of the distance
	} else {
		add_vf = 0;

		// cosangle_ij: n * (b) / (|n| * |b|)
		// the angle between normal to area_ij and beam from area_ij
		cosangle_ij = CosAngleBetween2Vectors(sx_source, sy_source, 1, bx, by, bz);

		// cosangle_ab: n * (-b) / (|n| * |b|)
		// the angle between normal to area_ab and beam from area_ij
		cosangle_ab = CosAngleBetween2Vectors(dem.Nx(a,b), dem.Ny(a,b), 1, -bx, -by, -bz);


		if ( cosangle_ij > 0. && cosangle_ab > 0. ) {
			add_vf = (cosangle_ij * cosangle_ab * area_ab * area_ij) /
				(M_PI * dist * dist);
			return add_vf;
		}
	}

	return 0;
}

double ViewFactorsCluster::GetZEdgeVector(const double z, const double nx, const double ny, const double deltaX, const double deltaY)
{
	//Return the z value of a cell's edge. The z value is considered at the center of the cell.
	return (z + nx * (-deltaX) * cellsize + ny * (-deltaY) * cellsize);
}

double ViewFactorsCluster::EuclidianDistance(const double x, const double y, const double z)
{
	//Return the euclidian distance between 2 points in 3D
	return sqrt(x*x + y*y + z*z);
}

int ViewFactorsCluster::GetSplitFactor(const double sideLength, const double dist)
{
	/*
	  Return the split factor (of the area) for the sum (numerical approximation of the integral)
	*/
	int split = 1;
	if (sideLength / dist >= sub_crit)
		{
			split = (int)((sideLength / (dist * sub_crit)) + 0.5);
			if (split < 1)
				{
					split = 1;
				}
		}
	return split;
}

double ViewFactorsCluster::GetCoordinateForSmallCell(const double veca, const double vecb, const double vecc, const int i, const int j, const int splitFactor)
{
	/*
	  find the coordinate values of the small cell uv in the plane ij
	  vec(a) + vec(b)-vec(a) * (u + 0.5) / spliti
       + vec(c)-vec(a) * (v + 0.5) / splitj
	  Note that spliti = splitj
	*/
	return veca + (vecb - veca) * (i + 0.5) * (1. / splitFactor) + (vecc - veca) * (j + 0.5) * (1. / splitFactor);
}
double ViewFactorsCluster::CosAngleBetween2Vectors(const double nx, const double ny, const double nz, const double bx, const double by, const double bz)
{
	// Return the cos value of an angle between 2 vectors like this :
	// n * (b) / (|n| * |b|)

	return (nx * bx + ny * by + nz * bz) / (EuclidianDistance(nx, ny, nz) * EuclidianDistance(bx, by, bz));
}
