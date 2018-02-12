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
#ifndef VIEWFACTORSCLUSTER_H
#define VIEWFACTORSCLUSTER_H

#include <meteoio/MeteoIO.h>
#include <alpine3d/ebalance/ViewFactorsAlgorithm.h>
#include <map>

class ViewFactorsCluster : public ViewFactorsAlgorithm {

	public:
		ViewFactorsCluster(const mio::Config& cfg, const mio::DEMObject &dem_in);
		double getSkyViewFactor(const int &i, const int &j);
		void getSkyViewFactor(mio::Array2D<double> &o_sky_vf) const;
		double GetViewfactor(const int i, const int j, const int a, const int b);

	private:
		mio::Array4D<double> vf_cluster;
		mio::Array4D<unsigned int> vc_cluster;
		mio::Array2D<double> sky_vf;
		std::map<int, double> vf;

		mio::DEMObject dem;
		double cellsize;
		double sw_radius;
		double max_shade_distance;
		double sub_crit;

		int dimx, dimy;
		unsigned int hSections, vSections;
		unsigned int schrittweitenLimit;

		void fill_vf_map();
		void calcHorizonField();
		void calcSchrittweite(const double& altitude, const double& dH, const double& distance,
		                      const double& horizon_tan_angle, unsigned int& schrittweite) const;
		double getHorizonForRay(const unsigned int& ix1, const unsigned int& iy1, const double& alpha,
		                        std::vector<unsigned int>& viewCells, std::vector<unsigned int>& vf_clusterT);
		double GetSymetricPartOfViewFactor(const int i, const int j, const int a, const int b);
		double GetZEdgeVector(const double z, const double nx, const double ny, const double deltaX, const double deltaY);
		double EuclidianDistance(const double x, const double y, const double z);
		int GetSplitFactor(const double sideLength, const double dist);
		double GetCoordinateForSmallCell(const double veca, const double vecb, const double vecc,
		                                 const int i, const int j, const int splitFactor);
		double CosAngleBetween2Vectors(const double nx, const double ny, const double nz, const double bx,
		                               const double by, const double bz);

		void calcVF_cluster();
		bool VF_calc(const unsigned int ix1, const unsigned int iy1, const int ix2, const int iy2, mio::Array1D<double> &tan_h);
};

#endif
