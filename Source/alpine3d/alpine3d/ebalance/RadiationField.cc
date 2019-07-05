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
#include <cmath>
#include <limits>

#include <alpine3d/ebalance/RadiationField.h>
#include <alpine3d/MPIControl.h>

RadiationField::RadiationField()
              : date(), dem(), dem_band(), direct(), diffuse(), Sun(),
                vecMeta(), vecMd(), vecCorr(),
                dem_mean_altitude(0.), cellsize(0.), dem_dimx(0), band_dimx(0), dimy(0), startx(0),
                day(true), night(false) {}

RadiationField::RadiationField(const mio::DEMObject& in_dem, const size_t& in_startx, const size_t& in_nx)
              : date(), dem(), dem_band(), direct(), diffuse(), Sun(),
                vecMeta(), vecMd(), vecCorr(),
                dem_mean_altitude(0.), cellsize(0.), dem_dimx(0), band_dimx(0), dimy(0), startx(0),
                day(true), night(false)
{
	setDEM(in_dem, in_startx, in_nx);
}

void RadiationField::setDEM(const mio::DEMObject& in_dem)
{
	setDEM(in_dem, (unsigned)0, mio::IOUtils::unodata);
}

void RadiationField::setDEM(const mio::DEMObject& in_dem, const size_t& in_startx, const size_t& in_nx)
{
	dem = in_dem;
	dem_dimx=dem.getNx();
	dimy=dem.getNy();
	cellsize=dem.cellsize;

	startx = in_startx;
	if (in_nx==mio::IOUtils::unodata)
		band_dimx = dem_dimx;
	else
		band_dimx = in_nx;
	
	//set the Sun to the midle of the dem
	mio::Coords dem_cntr( dem.llcorner );
	const double sun_easting = static_cast<double>( mio::Optim::round( dem.cellsize*static_cast<double>(dem.getNx())/2. ) );
	const double sun_northing = static_cast<double>( mio::Optim::round( dem.cellsize*static_cast<double>(dem.getNy())/2. ) );
	dem_cntr.moveByXY(  sun_easting, sun_northing );
	dem_mean_altitude = dem.grid2D.getMean();
	Sun.setLatLon( dem_cntr.getLat(), dem_cntr.getLon(), dem_mean_altitude );
}

void RadiationField::setStations(const std::vector<mio::MeteoData>& vecMeteo, const mio::Grid2DObject& albedo)
{
	if (vecMeteo.empty())
		throw mio::NoDataException("No meteo data provided!", AT);
	if (dem.empty()) //with the dem, the Sun has been set
		throw mio::NoDataException("The DEM has not been set", AT);
	
	//reset the state variables
	vecMeta.clear();
	vecMd.clear();
	vecCorr.clear();
	day = true;
	night = true;
	Sun.resetAltitude( dem_mean_altitude ); //it has been reset when computing the cells
	Sun.setDate(vecMeteo.front().date.getJulian(false), vecMeteo.front().date.getTimeZone()); //we have at least one station
	
	const double domain_alb = albedo.grid2D.getMean();
	if (domain_alb==mio::IOUtils::nodata)
			throw mio::IOException("[E] All cells have nodata albedo!! This can happen if all the cells are nodata/undefined in the land use file", AT);
	
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		mio::Coords location( vecMeteo[ii].meta.position );
		const bool has_meta = (location.getLat()!=mio::IOUtils::nodata) && (location.getLon()!=mio::IOUtils::nodata) && (location.getAltitude()!=mio::IOUtils::nodata);
		const bool has_meteo = (vecMeteo[ii](mio::MeteoData::ISWR)!=mio::IOUtils::nodata) 
		                                  && (vecMeteo[ii](mio::MeteoData::TA)!=mio::IOUtils::nodata) 
		                                  && (vecMeteo[ii](mio::MeteoData::RH)!=mio::IOUtils::nodata);
		if (has_meta && has_meteo) {
			const bool in_grid = albedo.gridify(location);
			
			double local_albedo( domain_alb );
			if (!in_grid) {
				const double HS = vecMeteo[ii](mio::MeteoData::HS);
				if (HS!=mio::IOUtils::nodata) //no big deal if we can not adapt the albedo
					local_albedo = (HS>=0.1)? 0.85 : 0.23; //snow or grass albedo
			} else {
				const int i_station = location.getGridI();
				const int j_station = location.getGridJ();
				const double tmp_albedo = albedo(i_station, j_station);
				if (tmp_albedo!=mio::IOUtils::nodata)
					local_albedo = tmp_albedo;
			}
			
			Sun.calculateRadiation(vecMeteo[ii](mio::MeteoData::TA), vecMeteo[ii](mio::MeteoData::RH), vecMeteo[ii](mio::MeteoData::P), local_albedo);
			bool local_day, local_night;
			double Md;
			const double corr = Sun.getCorrectionFactor(vecMeteo[ii](mio::MeteoData::ISWR), Md, local_day, local_night);
			vecCorr.push_back( corr );
			vecMd.push_back( Md );
			vecMeta.push_back( vecMeteo[ii].meta );
			day = local_day && day;
			night = local_night && night;
		}
	}
	
	if (vecMd.empty())
		throw mio::NoDataException("No suitable radiation station at "+vecMeteo.front().date.toString(mio::Date::ISO), AT);
}

void RadiationField::setMeteo(const mio::Grid2DObject& in_ta, const mio::Grid2DObject& in_rh,
                              const mio::Grid2DObject& in_p, const mio::Grid2DObject& in_albedo)
{
	//check that we can proceed
	if (dem.empty())
		throw mio::InvalidArgumentException("[E] Please set DEM before setting meteo grids!", AT);
	if (vecMd.empty())
		throw mio::InvalidArgumentException("[E] Please set radiation station before setting meteo grids!", AT);
	const size_t in_nx = in_ta.getNx(), in_ny = in_ta.getNy();
	if (band_dimx!=in_nx || dimy!=in_ny) //we only check TA in order to keep checks cheaper
		throw mio::InvalidArgumentException("[E] Given DEM and TA grid don't match!", AT);
	
	if (dem_band.empty()) dem_band.set(in_ta, 0.); //we just need to have the proper geolocalization for the IDW
	direct.set(in_ta, 0.); //reset to a band size (the "night" case might have set to the full dem)
	diffuse.set(in_ta, 0.); //reset to a band size (the "night" case might have set to the full dem)
	if (night) return; //no iswr at night
	
	mio::Grid2DObject Md;
	mio::Interpol2D::IDW(vecMd, vecMeta, dem_band, Md, 1000., 1.); //fixed scaling parameter of 1km
	mio::Grid2DObject corr_glob;
	mio::Interpol2D::IDW(vecCorr, vecMeta, dem_band, corr_glob, 1000., 1.); //fixed scaling parameter of 1km
	
	//get solar position for shading
	double solarAzimuth, solarElevation;
	Sun.position.getHorizontalCoordinates(solarAzimuth, solarElevation);
	const double tan_sun_elev = tan(solarElevation*mio::Cst::to_rad);
	
	for (size_t jj = 0; jj < dimy; jj++ ) {
		for (size_t i_dem = startx; i_dem < (startx+band_dimx); i_dem++ ) {
			const size_t i_band = i_dem - startx;
			if (in_albedo(i_band,jj)==mio::IOUtils::nodata || dem(i_dem,jj)==mio::IOUtils::nodata) {
				diffuse(i_band,jj) = mio::IOUtils::nodata;
				direct(i_band,jj) = mio::IOUtils::nodata;
				continue;
			}
			
			Sun.resetAltitude( dem(i_dem, jj) );
			Sun.calculateRadiation(in_ta(i_band,jj), in_rh(i_band,jj), in_p(i_band,jj), in_albedo(i_band,jj));
			double cell_toa, cell_direct, cell_diffuse;
			Sun.getHorizontalRadiation(cell_toa, cell_direct, cell_diffuse);
			
			if (day) {
				const double tan_horizon = mio::DEMAlgorithms::getHorizon(dem, i_dem, jj, solarAzimuth);
				const double global = cell_direct + cell_diffuse; //redo the splitting according to the interpolated Md
				
				if ( tan_sun_elev<tan_horizon ) { //cell is shaded
					cell_direct = 0.;
					cell_diffuse = global*Md(i_band, jj);
				} else {
					const double slope_azi=dem.azi(i_dem,jj);
					const double slope_angle=dem.slope(i_dem,jj);
					cell_diffuse = global*Md(i_band, jj);
					cell_direct = global*(1.-Md(i_band, jj));
					cell_direct = mio::SunTrajectory::projectHorizontalToSlope( solarAzimuth, solarElevation, slope_azi, slope_angle, cell_direct );
				}
			} else { //this is either dawn or dusk
				cell_diffuse += cell_direct;
				cell_direct=0.;
			}

			diffuse(i_band,jj) = corr_glob(i_band, jj)*cell_diffuse;
			direct(i_band,jj) = corr_glob(i_band, jj)*cell_direct;
		}
	}
}

void RadiationField::getPositionSun(double& o_solarAzimuth, double& o_solarElevation) const
{
	if (vecMd.empty())
		throw mio::InvalidArgumentException("[E] Please set radiation station before getting Sun's position!", AT);

	Sun.position.getHorizontalCoordinates(o_solarAzimuth, o_solarElevation);
}

void RadiationField::getBandOffsets(size_t& o_startx, size_t& o_stopx) const
{
	o_startx = startx;
	o_stopx = band_dimx;
}

void RadiationField::getRadiation(mio::Array2D<double>& o_direct, mio::Array2D<double>& o_diffuse) const
{
	o_direct = direct.grid2D;
	o_diffuse = diffuse.grid2D;
}
