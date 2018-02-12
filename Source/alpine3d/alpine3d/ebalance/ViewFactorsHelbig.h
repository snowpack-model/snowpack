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
#ifndef VIEWFACTORSHELBIG_H
#define VIEWFACTORSHELBIG_H

#include <meteoio/MeteoIO.h>
#include <alpine3d/ebalance/VFSymetricMatrix.h>
#include <alpine3d/ebalance/ViewFactorsAlgorithm.h>

#include <string>

class ViewFactorsHelbig : public ViewFactorsAlgorithm {
	public:
		ViewFactorsHelbig(const mio::Config& cfg, const mio::DEMObject &dem_in);
		double getSkyViewFactor(const int &i, const int &j);
		void getSkyViewFactor(mio::Array2D<double> &o_sky_vf) const;
		double GetViewfactor(const int i, const int j, const int a, const int b);
		double getSymetricTerrainViewFactor(const int &i, const int &j);
		
		double min_area;             // min (inclined) surface area
		double min_vterr;            // min terrain view factor
		bool vf_in_ram;

	private:
		mio::IOManager io;
		mio::Array2D<double> sky_vf, vf_t;
		double cellsize;
		int dimx, dimy;
		mio::DEMObject dem;
		std::string vf_file_in, tvfarea_file_in; //where to read the sky view factors and the terrain view factor x surface
		std::string vf_file_out, tvfarea_file_out; //where to write the sky view factors and the terrain view factor x surface
		double sub_crit;
		VFSymetricMatrix<float, double> vf; // view factor matrix with dynamic dimension
		int LW_distance_index, SW_distance_index;

		const static double to_rad;

		const static double vf_thresh; //below this threshold, view factor is forced to 0

		double EuclidianDistance(const double x, const double y, const double z);
		double DoubleRounding(const double aNumber, const double precision);
		int GetStep(const int a, const int b);
		double GetMax(const int step, const double azi, const bool x_Mode);
		double GetDelta(const int step, const double azi, const bool x_Mode);
		double GetPseudoAzimuth(const int i, const int j, const int m, const int t);
		void GetNextCell(int& x, int& y, double& tMaxX, double& tMaxY, 
		                 const int stepsX, const int stepsY, const double tDeltaX, const double tDeltaY);
		bool CellOutsideGrid(const int a, const int b);
		double DistanceBetween2Cells(const int i, const int j, const int a, const int b, const double bz);
		double DistanceOnTheBeamMiddleCell(const int i, const int j, const int a, const int b,
		                                   const int m, const int t, const double bz);
		double DistanceOnTheBeamBorderCell(const int i, const int j, const int a, const int b,
		                                   const int m, const int t, const double bz);
		bool IsAngleHigher(const double view_angleP, const int i, const int j, int &m, int &t);
		bool Is2CellsVisible(int i, int j, int m, int t);
		void ApplySunBorderTreatment(const int i, const int j, int &m, int &t);
		int ComputeFirstObstacle(const double solar_elev, const int i, const int j, const int mt, const int tt,
		                         int& horizon_x, int& horizon_y);
		double GetZEdgeVector(const double z, const double nx, const double ny, const double deltaX, const double deltaY);
		int GetSplitFactor(const double sideLength, const double dist);
		double GetCoordinateForSmallCell(const double veca, const double vecb, const double vecc, 
		                                 const int i, const int j, const int splitFactor);
		double CosAngleBetween2Vectors(const double nx, const double ny, const double nz, const double bx, 
		                               const double by, const double bz);
		double GetSymetricPartOfViewFactor(const int i, const int j, const int a, const int b);
		int InitGridViewFactors();
		int InitSkyViewFactors();
		bool InitializeViewFactor();
		
		void setVF_IN_RAM(bool);
};

#endif
