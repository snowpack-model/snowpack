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
#include <alpine3d/ebalance/ViewFactorsHelbig.h>

#include <ctime>
#include <cstdio>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"

//void writeLogFile(std::string filename, mio::Array2D<double> value);

const double ViewFactorsHelbig::to_rad = M_PI/180.;

const double ViewFactorsHelbig::vf_thresh = 1e-4; //below this threshold, view factor is forced to 0

ViewFactorsHelbig::ViewFactorsHelbig(const mio::Config& cfg, const mio::DEMObject &dem_in) : io(cfg), dem(dem_in)
{
	double lw_radius;
	double sw_radius;
	cfg.getValue("vf_in_ram", "EBalance", vf_in_ram);
	cfg.getValue("sub_crit", "EBalance", sub_crit);
	cfg.getValue("sw_radius", "EBalance", sw_radius);
	cfg.getValue("lw_radius", "EBalance", lw_radius);
	cfg.getValue("vf_file", "Input", vf_file_in, mio::IOUtils::nothrow);
	cfg.getValue("tvfarea", "Input", tvfarea_file_in, mio::IOUtils::nothrow);
	cfg.getValue("vf_file", "Output", vf_file_out, mio::IOUtils::nothrow);
	cfg.getValue("tvfarea", "Output", tvfarea_file_out, mio::IOUtils::nothrow);

	cellsize = dem.cellsize;
	dimx = dem.getNx();
	dimy = dem.getNy();

	LW_distance_index = (int)ceil(lw_radius / cellsize);
	SW_distance_index = (int)ceil(sw_radius / cellsize);

	sky_vf.resize(dimx, dimy);
	vf_t.resize(dimx, dimy);

	// switch for the storage of view factors in Progressive Refinement iteration
	if (vf_in_ram) {
		const int N = dimx * dimy;
		vf.resize ( N, N ); //required because storage is done in 1D hash table...
	}

	InitializeViewFactor();

	//writeLogFile("skyNora", sky_vf);
}

double ViewFactorsHelbig::getSkyViewFactor(const int &i, const int &j) {
        return sky_vf(i,j);
}

void ViewFactorsHelbig::getSkyViewFactor(mio::Array2D<double> &o_sky_vf) const {
	o_sky_vf = sky_vf;
}

double ViewFactorsHelbig::EuclidianDistance(const double x, const double y, const double z)
{
//Return the euclidian distance between 2 points in 3D
	return sqrt(x*x + y*y + z*z);
}

double ViewFactorsHelbig::DoubleRounding(const double aNumber, const double precision)
{
//Round a double, the goal is : |DoubleRounding(-X)| = |DoubleRounding(X)| (X > 0)
//precision is <1 (ex: 1e-4)
	const double temp = aNumber / precision;
	if ( temp >= 0. ) {
		return ( floor( temp + 0.5 ) * precision );
	} else {
		return ( ceil( temp - 0.5 ) * precision );
	}
}

int ViewFactorsHelbig::GetStep(const int a, const int b)
{
/*
Return the number of cell to move for one step to go from A to B.
It's the same method for row and column.
Possible values are 0 (same row or column), 1 or -1.
*/
	const double difference = a-b;

	if (difference < 0.) {
		return 1;
	} else if (difference > 0.) {
		return -1;
	}

	return 0;
}

double ViewFactorsHelbig::GetMax(const int step, const double azi, const bool x_Mode)
{
/*
Return the length components in x or y that are necessary
to leave the current cell in direction of the destination.
This is in 2D representation of the "top view" of the simulated area.
*/
	if (step!=0) {
		if (x_Mode) {
			return (0.5 * cellsize) / sin( azi );
		} else {
			return (-0.5 * cellsize) / cos( azi );
		}
	}
	return 0.;
}

double ViewFactorsHelbig::GetDelta(const int step, const double azi, const bool x_Mode)
{
/*
Return the length component in y or x to cross the width or the height
of one grid cell.
This is in 2D representation of the "top view" of the simulated area.
*/
	if (step!=0) {
		if (x_Mode) {
			return (cellsize / sin( azi ));
		} else {
			return ((-1.)*cellsize / cos(azi));
		}
	}
	return 0.;
}

double ViewFactorsHelbig::GetPseudoAzimuth(const int i, const int j, const int m, const int t)
{
/*
Return the angle between the "ray" (from the origin to the dest) and the "south"
(vector -1;-1) in counterclockwise.
This is in 2D representation of the "top view" of the grid.
*/

	const double flatdist = EuclidianDistance((i-m),(j-t),0);
	const double azi = acos(-(double)(t-j) / flatdist);

	if ( m-i > 0) {
		return azi;
	} else {
		return (2. * M_PI - azi);
	}
}

void ViewFactorsHelbig::GetNextCell(int &x, int &y, double &tMaxX, double &tMaxY, const int stepsX, const int stepsY, const double tDeltaX, const double tDeltaY)
{//Compute and set the next intermediate cell that is on the "ray"
//(form the origin to the dest)

	if ((fabs(tMaxX) < fabs(tMaxY) && stepsX != 0) || tMaxY==0) {
		x = x + stepsX;
		tMaxX = tMaxX + tDeltaX;
	} else {
		if ((fabs(tMaxY) < fabs(tMaxX) && stepsY != 0) || tMaxX == 0) {
			y = y + stepsY;
			tMaxY = tMaxY + tDeltaY;
		} else {
			x = x + stepsX;
			tMaxX = tMaxX + tDeltaX;
			y = y + stepsY;
			tMaxY = tMaxY + tDeltaY;
		}
	}
}

bool ViewFactorsHelbig::CellOutsideGrid(const int a, const int b)
{//Return is the cell is outside the grid
	return  (a < 0 || a > dimx - 1 || b < 0 || b > dimy - 1);

}

double ViewFactorsHelbig::DistanceBetween2Cells(const int i, const int j, const int a, const int b, const double bz)
{
/*
Deprecated
It's the basic implementation of the distance between 2 cells
*/
	double bx = (a - i) * cellsize;
	double by = (b - j) * cellsize;

	return EuclidianDistance(bx,by,bz);
}

double ViewFactorsHelbig::DistanceOnTheBeamMiddleCell(const int i, const int j, const int a, const int b, const int m, const int t, const double bz)
{
/*
Deprecated
Return the distance between cell ij and the intermerdiate cell ab.
The distance is on the beam between cell ij and cell mt, in the middle of the
beam's segment of cell ab.
*/
	double bx = (a-i) * cellsize;
	double by = (b-j) * cellsize;

	//Compute the horizontal distance between source and intermediate cell
	double hyp = EuclidianDistance(bx,by,0);

	//Compute the vector between source and destination cell
	double cx = (m-i)*cellsize;
	double cy = (t-j)*cellsize;

	//Compute the cosangle between the source - intermediate and the source -destination vector
	double cosangle = CosAngleBetween2Vectors(bx, by, 0, cx, cy, 0);

	return EuclidianDistance(cosangle*hyp,bz,0);
}

double ViewFactorsHelbig::DistanceOnTheBeamBorderCell(const int i, const int j, const int a, const int b, const int m, const int t, const double bz)
{
/*
Return the distance between cell ij and the intermerdiate cell ab.
This is the distance between the center of the cell ij and the border
of the cell ab that the beam, between center cell ij and center cell mt, meets.
*/
	double p1x, p1y, v1x/*, v1y*/; //Starting point and vector for horizontal border.
	double p2x, p2y/*, v2x*//*, v2y*/; //Starting point and vector for vertical border.
	double p3x, p3y, v3x, v3y; //Starting point and vector for the beam
	double r1x, r1y, bvar1, avar1; //Intersection point and lambdas for the horizontal border
	double r2x, r2y, bvar2/*, avar2*/; //Intersection point and lambdas for the vertical border

	//This is the horizontal vector
	v1x = cellsize;
	//v1y = 0;

	//This is the vertical vector
	//v2x = 0;
	//v2y = cellsize;

	//This is the beam vector
	v3x= (m-i)*cellsize;
	v3y= (t-j)*cellsize;

	//This is the starting point of the vertical border
	p2y = b*cellsize;
	p2x = a*cellsize;

	//This is the starting point of the horizontal border
	p1x = p2x;
	p1y = p2y;

	//This is the sarting point of the beam (it's the center of the cell)
	p3x = (i+0.5)*cellsize;
	p3y = (j+0.5)*cellsize;

	if ((i > m && bz > 0) || (i < m && bz < 0)) {
		//The vertical border from the source is the right border
		p2x += cellsize;

	}

	if ((j > t && bz > 0) || (j < t && bz<0)) {
		//The horizontal border from the source is the top border
		p1y += cellsize;
	}

	if (i == m) {
		//The beam has an infinite slope.
		//The only response possible is for the horizontal border
		r1y = p1y;
		bvar1 = (p1y-p3y)/v3y;
		r1x = bvar1*v3x + p3x;
		return EuclidianDistance(r1x-p3x, r1y-p3y, bz);
	}

	if (j == t) {
		//The beam has no slope
		//The only response possible is for the vertical border
		r2x = p2x;
		bvar1 = (p2x-p3x)/v3x;
		r2y = bvar1*v3y+p3y;
		return EuclidianDistance(r2x-p3x, r2y-p3y, bz);
	}

	//The 2 answers are possible, compute the 2 points.
	r1y = p1y;
	bvar1 = (p1y-p3y)/v3y;
	r1x = bvar1*v3x + p3x;
	avar1= (r1x - p1x) / v1x;

	r2x = p2x;
	bvar2 = (p2x-p3x)/v3x;
	r2y = bvar2*v3y+p3y;
	//avar2 = (r2y-p2y)/v2y;

	//If the lambda for the horizontal is between 0 and 1 the intersection
	//is on the horizontal border it's not necessary to control the other param
	//because it's not possible to have 2 points on the cell (or it's the same)
	if (avar1 <= 1 && avar1 >= 0) {
		return EuclidianDistance(r1x-p3x, r1y-p3y, bz);
	} else {
		return EuclidianDistance(r2x-p3x, r2y-p3y, bz);
	}
}

bool ViewFactorsHelbig::IsAngleHigher(const double view_angleP, const int i, const int j, int &m, int &t)
{
/*
Return if the angle in parameter is higher than the all intermediate angle cell
(angle between the intermerdiate and the source cells).
If this angle not the highest, the destination cell values is changed
to the cell that has the first higher intermerdiate angle.
The visibility test according to Amanatides(1987) for ray tracing
*/
	const double zSource = dem.grid2D(i,j);
	bool visible;
	double tMaxX, tMaxY;
	double tDeltaX, tDeltaY;
	int stepsx, stepsy;
	double view_angle, testview_angle;
	double azi;
	int a, b;
	double bz;
	double dist;

	//Round value to avoid that internal rounding or inaccuracies
	// prevent that the comparison between them go wrong
	view_angle = DoubleRounding(view_angleP,1.e-4);

	//Get the angle from the ray and the south
	azi = GetPseudoAzimuth(i, j, m, t);

	//Get the step in x and y
	stepsx = GetStep(i, m);
	stepsy = GetStep(j, t);

	//Get Max in x and y
	tMaxX= GetMax(stepsx, azi, true);
	tMaxY= GetMax(stepsy, azi, false);

	//Get Delta in x and y
	tDeltaX = GetDelta(stepsx, azi, true);
	tDeltaY = GetDelta(stepsy, azi, false);

	// Round values to avoid that internal rounding or inaccuracies
	// prevent that the comparison between them go wrong

	tMaxX = DoubleRounding(tMaxX,1.e-4);
	tDeltaX = DoubleRounding(tDeltaX,1.e-4);
	tMaxY = DoubleRounding(tMaxY,1.e-4);
	tDeltaY = DoubleRounding(tDeltaY,1.e-4);

	// Set the initial intermerdiate cell (start on the source cell)
	a = i;
	b = j;

	visible = true;

	// do while the source and the destination seems to bo visible
	do
	{
		// Get the next intermerdiate cells
		GetNextCell(a, b, tMaxX, tMaxY, stepsx, stepsy, tDeltaX, tDeltaY);

		if (CellOutsideGrid(a,b)) {
			// The new cell isn't inside the grid, the 2 cells is
			// visible
			break;
		}

		// If the intermerdiate cell is the destination cell, the 2
		// cells are visible
		if ( a == m && b == t ) {
			break;
		}

		// Get the distance on the beam
		bz = dem.grid2D(a,b) - zSource;

		//There is 3 differents technics
		//double DistanceOnTheBeamBorderCell, DistanceOnTheBeamMiddleCell and DistanceBetween2Cells
		dist = DistanceOnTheBeamBorderCell(i, j, a, b, m, t, bz);

		testview_angle = bz/dist;

		// Round value to avoid that internal rounding or inaccuracies
		// prevent that the comparison between them go wrong
		testview_angle = DoubleRounding(testview_angle,1.e-4);

		// If the intermerdiante angle is higher, the 2 cells aren't visible
		if ( testview_angle > view_angle ) {
			visible = false;
			// Set the destination cell to the first
			// higher intermerdiate cell (and not the highest!)
			m = a;
			t = b;
		}

	} while ( visible );

	return visible;
}

bool ViewFactorsHelbig::Is2CellsVisible(int i, int j, int m, int t)
{//Return if 2 cells are visible

	// If the source and the origin are the same, return false
	if (i==m && j == t) {
		return false;
	}
	if (i>m || ((i==m) && (j>t))) {
		//we swap (i,j) and (m,t) so that the source cell is always the most on the lower right
		//This solves the problem that numerically it is hard to guarantee A->B exactly equivalent to B->A
		//(rounding errors leading to slightly different paths)
		const int tempm = m; m = i; i = tempm; //swap i
		const int tempt = t; t = j; j = tempt; //swap j
	}

	// Compte the angle between the source and the destination
	const double bx = (double)(m - i) * cellsize;
	const double by = (double)(t - j) * cellsize;
	const double bz = dem.grid2D(m,t) - dem.grid2D(i,j);
	const double dist = EuclidianDistance(bx,by,bz);
	const double view_angle = bz / dist;

	//Return if the angle is higher than intermerdiate angle
	return IsAngleHigher(view_angle, i, j, m, t);

}
void ViewFactorsHelbig::ApplySunBorderTreatment(const int i, const int j, int &m, int &t)
{
// special border treatment of ij to assure that for border cells the horizon is, by looking along the
// border, the farthest that grid cell can view
// this is done by giving that cell the appropriate edge cell to determine the horizon for that
	if ( (j == 0 && t == 0) || (j == dimy - 1 && t == dimy - 1) ) {
			if ( m < i ) {
				switch ( j ) {
					case 0:
						m = 0;
						t = 0;
						break;
					default:
						m = 0;
						t = dimy - 1;
						break;
				}
			} else {
				switch ( j ) {
					case 0:
						m = dimx - 1;
						t = 0;
						break;
					default:
						m = dimx - 1;
						t = dimy - 1;
						break;
				}
			}
	} // end of if condition : upper and lower border

	if ( (i == 0 && m == 0) || (i == dimx - 1 && m == dimx -1) ) {
		if ( t < j ) {
			switch ( i ) {
				case 0:
					m = 0;
					t = 0;
					break;
				default:
					m = dimx - 1;
					t = 0;
					break;
			}
		} else {
			switch ( i ) {
				case 0:
					m = 0;
					t = dimy - 1;
					break;
				default:
					m = dimx - 1;
					t = dimy - 1;
					break;
			}
		}
	} // end of if condition : left and right border
}

//solar_elev is in degrees
int ViewFactorsHelbig::ComputeFirstObstacle(const double solar_elev, const int i, const int j, const int mt, const int tt,
                                        int& horizon_x, int& horizon_y)
{/* Compute and set the first obstacle on the ray between the source and
    the destination that has a higher angle value.*/
	horizon_x = mt;
	horizon_y = tt;

	// If the source and the destination are the same, the horizon is this cell
	if (horizon_x!=i || horizon_y!=j) {
		// Apply eventual border treatement
		ApplySunBorderTreatment(i, j, horizon_x, horizon_y);

		// Variable that will contains the first higher angle
		// (or the destination if no higher angle)
		double view_angle = sin( solar_elev*to_rad );

		// Detect is the angle is the highest. Set the "horizon" to the first higher obstacle
		IsAngleHigher(view_angle, i, j, horizon_x, horizon_y);
	}

	return EXIT_SUCCESS;
}

double ViewFactorsHelbig::GetZEdgeVector(const double z, const double nx, const double ny, const double deltaX, const double deltaY)
{
//Return the z value of a cell's edge. The z value is considered at the center of the cell.
	return (z + nx * (-deltaX) * cellsize + ny * (-deltaY) * cellsize);
}

int ViewFactorsHelbig::GetSplitFactor(const double sideLength, const double dist)
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
double ViewFactorsHelbig::GetCoordinateForSmallCell(const double veca, const double vecb, const double vecc, const int i, const int j, const int splitFactor)
{
/*
find the coordinate values of the small cell uv in the plane ij
vec(a) + vec(b)-vec(a) * (u + 0.5) / spliti
       + vec(c)-vec(a) * (v + 0.5) / splitj
Note that spliti = splitj
*/
	return veca + (vecb - veca) * (i + 0.5) * (1. / splitFactor) + (vecc - veca) * (j + 0.5) * (1. / splitFactor);
}
double ViewFactorsHelbig::CosAngleBetween2Vectors(const double nx, const double ny, const double nz, const double bx, const double by, const double bz) {
// Return the cos value of an angle between 2 vectors like this :
// n * (b) / (|n| * |b|)

	return (nx * bx + ny * by + nz * bz) / (EuclidianDistance(nx, ny, nz) * EuclidianDistance(bx, by, bz));
}

double ViewFactorsHelbig::GetSymetricPartOfViewFactor(const int i, const int j, const int a, const int b)
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

double ViewFactorsHelbig::GetViewfactor(const int i, const int j, const int a, const int b )
{
/*
Return the view factor between 2 cells.
If the storage of the view factor is enabled the symetric part of the view factor
is read on memory
else it's computed
*/


        //The method GetViewFactor control if the view factor is stored
        if (vf_in_ram) {
                const unsigned int ij = j * dimx + i;
                const unsigned int ab = a + b*dimx;	//1D addressing of the shooting cell for the view factors matrix
                const double sx_0 = dem.Nx(i,j);
                const double sy_0 = dem.Ny(i,j);
                const double vf_0 = vf(ij, ab);
                //view factor is symetric part (from vf() ) divided by the "area" of the current cell
                return vf_0 / ( sqrt( sx_0*sx_0 + sy_0*sy_0 + 1. ) * cellsize * cellsize );
        } else {

                double symetric_part;

                // If the view factor is stored return the value, or check if the 2 cells are
                // visible and compute and return the value
                if (vf_in_ram) {
                        const int ij = j * dimx + i;
                        const int ab = b * dimx + a;
                        symetric_part = vf(ij, ab);
                } else {
                        if (Is2CellsVisible(i, j, a, b)) {
                                symetric_part = GetSymetricPartOfViewFactor(i, j, a, b);
                        } else {
                                return 0.;
                        }
                }

                //Divide the symetric part of the view factor to obtain the view factor
                // x and y values of the normal vector and z value of a cell
                const double sx_source = dem.Nx(i,j);
                const double sy_source = dem.Ny(i,j);
                //const double z_source = dem.grid2D(i,j);
                const double area_ij = EuclidianDistance(sx_source, sy_source, 1.) * cellsize * cellsize;
                return symetric_part / area_ij;
        }

}

/**
* @brief Computes the view factors for each point of the grid
* @return (int) EXIT_SUCCESS
*/
int ViewFactorsHelbig::InitGridViewFactors()
{
	double t1, t0;		// clock variables for measuring the time needed in the main ComputeHorizon
	double delay;		// for the delay of calling the clock function
	long int count_vf = 0;	// counts number of View factors(i,j,a,b) > 0
	const int maxDistIdx = std::max(LW_distance_index, SW_distance_index); // Maximal radius distance

	t0 = clock();
	t1 = clock();
	delay = t1 - t0;

	//For each cell of the grid
	for (int i = 0; i < dimx ; i++) {
		for (int j = 0; j < dimy; j++) {
			//Compute the area and the indice of the external cell

			const int ij = j * dimx + i;

			//For each cell of the grid
			const int min_x = std::max(0, i-maxDistIdx);
			const int max_x = std::min(i+maxDistIdx, (int)dimx);
			const int min_y = std::max(0, j-maxDistIdx);
			const int max_y = std::min(j+maxDistIdx, (int)dimy);

			for (int m = min_x; m < max_x; m++) {
				for (int t= min_y; t < max_y ; t++) {

					//Compute the indice of the internal cell
					const int ab = t * dimx + m;

					// Check on the indices of 2 cells to compute only one symetric part
					if (ab <= ij) {
						if (Is2CellsVisible(i, j, m, t)) {
							// Get the symetric part
							const double temp = GetSymetricPartOfViewFactor ( i, j, m, t);
							if (temp > vf_thresh) {
								//we only account for vf large enough
								count_vf++;

								// Add the view factory value to the sum of view factor cells
								// The division by the area is done after, only for the sky view_factor
								vf_t(i,j) += temp;
								vf_t(m,t) += temp;
								if (vf_in_ram) {
									// Only if storage is enabled
									vf.setElement(ij,ab,temp);
								}
							}


						}
					}
				}
			}
		}
	}

	t1 = clock();
	std::cout << "[i] " << count_vf << " view factors computed in " << (t1 - t0 - delay) / CLOCKS_PER_SEC << " seconds" << std::endl;

	return(EXIT_SUCCESS);
}


/**
* @brief Computes the sky view factors of the whole grid
* @return (int) EXIT_SUCCESS
*/
int ViewFactorsHelbig::InitSkyViewFactors()
{
	for (int m = 0; m < dimx; m++) {
		for (int t= 0; t <dimy ; t++) {
			const double tempvf_t = vf_t(m,t) / (EuclidianDistance(dem.Nx(m,t), dem.Ny(m,t), 1.) * cellsize * cellsize);
			sky_vf(m,t) = 1.0 - tempvf_t;
		}
	}

	return(EXIT_SUCCESS);
}

bool ViewFactorsHelbig::InitializeViewFactor ( )
{
	if (vf_file_in.empty() || tvfarea_file_in.empty()) {
		mio::Timer timer_nora;
		timer_nora.start();
		//computing the view factors
		try {
			std::cout << "[i] computing grid view factors" << std::endl;
			InitGridViewFactors();
		} catch (std::bad_alloc&) {
			std::cout << "[E] ebalance : Exceeded memory space " << std::endl;
			fflush( stdout );
			exit(1);
		}

		//computing sky and terrain view factors
		InitSkyViewFactors();
		
		timer_nora.stop();
		std::cout << "noras viewfactor calculation: " << timer_nora.getElapsed() << std::endl;
		
		//writing view factors out
		if (!vf_file_out.empty()) {
			const mio::Grid2DObject tmp(dem.cellsize, dem.llcorner, sky_vf);
			io.write2DGrid(tmp, vf_file_out);
		}
		if (!tvfarea_file_out.empty()) {
			const mio::Grid2DObject tmp(dem.cellsize, dem.llcorner, vf_t);
			io.write2DGrid(tmp, tvfarea_file_out);
		}
	} else {
		mio::Grid2DObject tmp(dem, mio::IOUtils::nodata);
		io.read2DGrid(tmp, vf_file_in);
		sky_vf = tmp.grid2D;
		io.read2DGrid(tmp, tvfarea_file_in);
		vf_t = tmp.grid2D;
	}

	// initialization min surface patches
	min_area = cellsize * cellsize / cos( 70. * M_PI / 180. );
	// initialization min terrain view factor sum
	min_vterr = 0.999;

	// at model start the min, max values of various geometric parameters are computed since required for the radiation exchange
	for ( int m = 0; m < dimx; m++ ) {
		for ( int t = 0; t < dimy; t++ ) {
			// derivation of minimum (inclined) grid patches for the stopping criteria of the radiosity eqn.
			if ( (sqrt( dem.Nx(m,t) * dem.Nx(m,t) + dem.Ny(m,t) * dem.Ny(m,t) + 1. ) * cellsize * cellsize) < min_area ) {
				min_area = sqrt( dem.Nx(m,t) * dem.Nx(m,t) + dem.Ny(m,t) * dem.Ny(m,t) + 1. ) * cellsize * cellsize;
			}

			// derivation of minimum terrain view factors
			if ( (1.0 - sky_vf(m,t)) < min_vterr && (1.0 - sky_vf(m,t)) > 0.001 ) {
			// min vterr has to be larger zero, this is
			// possible as these patches will not take part in the iterating process anyway
				min_vterr = 1.0 - sky_vf(m,t);
			}
		}
	}

	return EXIT_SUCCESS;
} // end of InitializeViewFactor

void ViewFactorsHelbig::setVF_IN_RAM(bool v)
{
        vf_in_ram = v;
}

double ViewFactorsHelbig::getSymetricTerrainViewFactor(const int &i, const int &j)
{
        return vf_t(i,j);
}

#pragma GCC diagnostic pop
