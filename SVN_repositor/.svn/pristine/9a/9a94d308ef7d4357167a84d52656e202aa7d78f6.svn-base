/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef COORDSALGORITHMS_H
#define COORDSALGORITHMS_H

#include <string>

namespace mio {
/**
 * @class CoordsAlgorithms
 * @brief A static class to handle geographic algorithms.
 * This class offers methods to handle lat/lon coordinates, to transform lat/lon coordinates, to compute distances on a sphere or
 * to convert between coordinate systems.
 *
 * @ingroup data_str
 * @author Mathias Bavay
 * @date   2015-09-22
*/

class CoordsAlgorithms {
public:
	//Handling of lat/lon
	static double dms_to_decimal(const std::string& dms);
	static std::string decimal_to_dms(const double& decimal);
	static void parseLatLon(const std::string& coordinates, double& lat, double& lon);
	static std::string printLatLon(const double& latitude, const double& longitude);
	static double lat_degree_lenght(const double& latitude);
	static double lon_degree_lenght(const double& latitude);
	
	//EPSG helper methods
	static short int str_to_EPSG(const std::string& coordsystem, const std::string& coordparam);
	static void EPSG_to_str(const int& epsg, std::string& coordsystem, std::string& coordparam);

	//handling of rotated coordinates
	static void rotatedToTrueLatLon(const double& lat_N, const double& lon_N, const double& lat_rot, const double& lon_rot, double &lat_true, double &lon_true);
	static void trueLatLonToRotated(const double& lat_N, const double& lon_N, const double& lat_true, const double& lon_true, double &lat_rot, double &lon_rot);

	//handling of distances on a sphere
	static double cosineDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha);
	static void cosineInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon);
	static double VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha);
	static void VincentyInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon);

	//Coordinates conversions
	static void WGS84_to_CH1903(const double& lat_in, const double& long_in, double& east_out, double& north_out);
	static void CH1903_to_WGS84(const double& east_in, const double& north_in, double& lat_out, double& long_out);
	static void WGS84_to_UTM(const double& lat_in, double long_in, const std::string& coordparam, double& east_out, double& north_out);
	static void UTM_to_WGS84(double east_in, double north_in, const std::string& coordparam, double& lat_out, double& long_out);
	static void WGS84_to_UPS(const double& lat_in, const double& long_in, const std::string& coordparam, double& east_out, double& north_out);
	static void UPS_to_WGS84(const double& east_in, const double& north_in, const std::string& coordparam, double& lat_out, double& long_out);
	static void WGS84_to_PROJ4(const double& lat_in, const double& long_in, const std::string& coordparam, double& east_out, double& north_out);
	static void PROJ4_to_WGS84(const double& east_in, const double& north_in, const std::string& coordparam, double& lat_out, double& long_out);

 private:
	static int getUTMZone(const double& latitude, const double& longitude, std::string& zone_out);
	static void parseUTMZone(const std::string& zone_info, char& zoneLetter, short int& zoneNumber);

 private:
	///Keywords for selecting an ellipsoid to use
	enum ELLIPSOIDS_NAMES {
		E_WGS84, ///<Globally useable WGS84 ellipsoid
		E_GRS80, ///<GRS80 ellispoid, equivalent to WGS84 but used by America and Australia
		E_AIRY, ///<Airy ellispoid, good fit for the UK
		E_INTL1924, ///<International 1924 ellispoid, good for most of Europe
		E_CLARKE1880, ///<Clarke 1880, good for Africa
		E_GRS67 ///<GRS67 ellispoid, good for South America
	};
	struct ELLIPSOID {
		double a;
		double b;
	};
	static const struct ELLIPSOID ellipsoids[6];
};
} //end namespace

#endif
