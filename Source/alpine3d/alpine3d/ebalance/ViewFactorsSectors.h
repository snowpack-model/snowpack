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
#ifndef VIEWFACTORSSECTORS_H
#define VIEWFACTORSSECTORS_H

#include <meteoio/MeteoIO.h>
#include <alpine3d/ebalance/ViewFactorsAlgorithm.h>
#include <map>

class ViewFactorsSectors : public ViewFactorsAlgorithm {
	public:
		ViewFactorsSectors(const mio::Config& cfg, const mio::DEMObject &dem_in);
		double getSkyViewFactor(const int &i, const int &j);
		void getSkyViewFactor(mio::Array2D<double> &o_sky_vf) const;
		double GetViewfactor(const int i, const int j, const int a, const int b);
		int getViewCells(const int i, const int j, const int h, const int v);

	private:
		std::map<int, double> vf;
		mio::Array4D<unsigned int> viewCells;
		mio::Array2D<double> sky_vf;
		mio::DEMObject dem;

		double cellsize;
		double max_shade_distance;
		double sw_radius;

		int dimx, dimy;
		int hSections, vSections;
		unsigned int schrittweitenLimit;

		void fill_vf_map();
		void calcHorizonField();
		void calcSchrittweite(const double& altitude, const double& dH, const double& distance,
		                      const double& horizon_tan_angle, unsigned int& schrittweite) const;
		double getHorizonForRay(const unsigned int& ix1, const unsigned int& iy1, const double& alpha,
		                        std::vector<unsigned int>& viewCells);
		double getPrecViewFactor(const double& cos_theta_min, const double& cos_theta_max, const double& radius);
};

#endif
