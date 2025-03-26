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
#include "MeteoObj.h"
#include <alpine3d/SnowpackInterfaceWorker.h>
#include <alpine3d/AlpineMain.h>

using namespace std;
using namespace mio;

//This is a function to retrieve soil temperature at a given depth (measured from
//the soil surface, not from the snow surface) for output purposes.
inline double getSoilTemperature(const SnowStation& pixel, const double& depth)
{
	if(pixel.getNumberOfNodes() == 0) {
		return IOUtils::nodata;
	} else if(depth < 1e-5) {//we want the temperature at the soil surface
		return pixel.Ndata[pixel.SoilNode].T;
	} else if(depth >= pixel.Ground) {//we want the temperature below the lowest node
		return pixel.Ndata.front().T;
	} else {
		const double height( pixel.Ground - depth );

		//Looking for the node located directly above the specified depth
		size_t iNode(0);
		while (pixel.Ndata[iNode].z < height)
			++iNode;

		const NodeData& aboveNode( pixel.Ndata[iNode] );
		if(aboveNode.z == height) {//exact match
			return aboveNode.T;
		} else {//compute linear interpolation between the two nodes
			const NodeData& belowNode( pixel.Ndata[iNode-1] );
			const double a = (aboveNode.T - belowNode.T)/(aboveNode.z - belowNode.z);
			return a*(height - belowNode.z) + belowNode.T;
		}
	}
}

//This is a function to retrieve snow density from the surface to a given depth for output purposes.
inline double getSnowDensityDepth(const SnowStation& pixel, const double& depth)
{
	if(pixel.getNumberOfNodes() == 0) {
		return IOUtils::nodata;
	} else if(depth < 1e-5) {//we want the density of the surface element
		return pixel.Edata[pixel.getNumberOfElements()-1].Rho;
	} else {
		size_t i = pixel.getNumberOfElements();
		double H = 0.;
		double M = 0.;
		while (i-- > 0) {
			const double dz = pixel.Edata[i].L;
			H += dz;
			M += dz * pixel.Edata[i].Rho;
			if (H > depth) {
				return (H==0)?(IOUtils::nodata):(M/H);
			}
		}
		return (H==0)?(IOUtils::nodata):(M/H);
	}
}

//This is a function to retrieve soil runoff at a given depth (measured from
//the soil surface, not from the snow surface) for output purposes.
inline double getSoilRunoff(const SnowStation& pixel, const SurfaceFluxes& surfaceFlux, const double& depth)
{
	if(pixel.getNumberOfNodes() == 0) {
		return IOUtils::nodata;
	} else if(depth < 1e-5) {//we want the water flux at the soil surface
		return surfaceFlux.mass[SurfaceFluxes::MS_SURFACE_MASS_FLUX];
	} else if(depth >= pixel.Ground) {//we want the water flux below the lowest node
		return pixel.Ndata.front().water_flux;
	} else {
		const double height( pixel.Ground - depth );

		//Looking for the node located directly above the specified depth
		size_t iNode(0);
		while (pixel.Ndata[iNode].z < height)
			++iNode;

		const NodeData& aboveNode( pixel.Ndata[iNode] );
		if(aboveNode.z == height) {//exact match
			return aboveNode.water_flux;
		} else {//compute linear interpolation between the two nodes
			const NodeData& belowNode( pixel.Ndata[iNode-1] );
			const double a = (aboveNode.water_flux - belowNode.water_flux)/(aboveNode.z - belowNode.z);
			return a*(height - belowNode.z) + belowNode.water_flux;
		}
	}
}

//This is a function to retrieve a given parameter at a given depth in the snow
//for output purposes.
inline double getValueAtDepth(const SnowStation& pixel, const double& depth)
{
	const size_t nNodes = pixel.getNumberOfNodes();
	if (nNodes==0) return mio::IOUtils::nodata; //no nodes -> nodata

	const size_t max_ii = nNodes - 1;
	const double offset = pixel.Ground;
	const double snow_height = pixel.Ndata[max_ii].z - offset;
	const double height = snow_height - depth;

	if (height<0.) return mio::IOUtils::nodata; //not deep enough -> nodata
	if (depth==0.) return pixel.Ndata[max_ii].T; //if depth==0 -> surface

	//looking for the first node less than depth.
	//At this point, we know that we have >= than depth snow
	size_t ii = 0;
	while ( (pixel.Ndata[ii].z-offset) <= height ) {
		ii++;
	}

	if ( (pixel.Ndata[ii].z-offset) == height) return pixel.Ndata[ii].T; //we found it exactly

	//compute linear interpolation between the two points y = az+b
	const double a = (pixel.Ndata[ii].T-pixel.Ndata[ii-1].T) / (pixel.Ndata[ii].z-pixel.Ndata[ii-1].z);
	const double b = pixel.Ndata[ii-1].T;
	const double z = height - (pixel.Ndata[ii-1].z-offset);
	return a*z + b;
}

// Returns the average snow temperature in the top x cm, or, when snow depth is less than x cm, the average snow temperature.
inline double getAvgAtDepth(const SnowStation& pixel, const double& depth, const bool& is_temp)
{

	double thickness = 0.;
	double sum = 0.;
	if (pixel.getNumberOfElements() > pixel.SoilNode) {
		size_t i = pixel.getNumberOfElements();
		while (i-- > pixel.SoilNode) {
			const double val = (is_temp)? pixel.Edata[i].Te : pixel.Edata[i].Rho;
			if (thickness < depth) {
				sum += pixel.Edata[i].L * val;
				thickness += pixel.Edata[i].L;
			} else { //we found the first element that is too deep
				if (thickness!=depth) { //just in case, we might have found an element exactly at "depth"
					const double thick_contrib = depth - thickness;
					sum += thick_contrib * val;
					thickness += thick_contrib;
				}
				break;
			}
		}
	}

	const double value = (thickness == 0.) ? (IOUtils::nodata) : (sum / thickness);
	return value;
}

/******************************************************************************
 * Constructor / Destructor
 ******************************************************************************/

/**
 * @brief Constructs and initialise Snowpack Interface Worker. Create one worker
 * and init values for slice which the worker need to simulate it slice.
 * @param[in] io_cfg configuration to pass to Snowpack
 * @param[in] dem_in gives the demographic Data. Also tetermines size and position of the slice
 * @param[in] landuse_in gives the landuse Data. Also tetermines size and position of the landuse for slice
 * @param[in] pts_in gives the spezial points. For this points more output is done then for the others. Calcualtion is the same.
 * @param[in] snow_stations gives a vector of pointers to the SnowStation objects relevant for this thread
 * @param[in] snow_stations_coord provide the (ii,jj) coordinates of each snow station. This is necessary because the slices might be irregular
 * @param[in] grids_not_computed_in_worker vector of grids that are required but will be computed outside the workers
 */
SnowpackInterfaceWorker::SnowpackInterfaceWorker(const mio::Config& io_cfg,
                                                 const mio::DEMObject& dem_in,
                                                 const mio::Grid2DObject& landuse_in,
                                                 const std::vector< std::pair<size_t,size_t> >& pts_in,
                                                 const std::vector<SnowStation*>& snow_stations,
                                                 const std::vector<std::pair<size_t,size_t> >& snow_stations_coord,
                                                 const std::vector<std::string>& grids_not_computed_in_worker)
 : SnowStationsCoord(snow_stations_coord), sn_cfg(io_cfg), sn(sn_cfg), meteo(sn_cfg), stability(sn_cfg, false), sn_techsnow(sn_cfg), dem(dem_in),
   dimx(dem.getNx()), dimy(dem.getNy()), SnowStations(snow_stations),
   isSpecialPoint(snow_stations.size(), false), landuse(landuse_in), store(dem_in, 0.), erodedmass(dem_in, 0.), grids(), snow_pixel(), meteo_pixel(),
   surface_flux(), soil_temp_depths(), soil_runoff_depths(), snow_density_depths(), calculation_step_length(0.), height_of_wind_value(0.),
   snow_temp_depth(IOUtils::nodata), snow_avg_temp_depth(IOUtils::nodata), snow_avg_rho_depth(IOUtils::nodata),
   enable_simple_snow_drift(false), useDrift(false), useEBalance(false), useCanopy(false)
{

	sn_cfg.getValue("CALCULATION_STEP_LENGTH", "Snowpack", calculation_step_length);
	sn_cfg.getValue("HEIGHT_OF_WIND_VALUE", "Snowpack", height_of_wind_value); //currently unused
	sn_cfg.getValue("CANOPY", "Snowpack", useCanopy);
	io_cfg.getValue("SNOW_TEMPERATURE_DEPTH", "Output", snow_temp_depth, IOUtils::nothrow);
	io_cfg.getValue("SNOW_AVG_TEMPERATURE_DEPTH", "Output", snow_avg_temp_depth, IOUtils::nothrow);
	io_cfg.getValue("SNOW_AVG_DENSITY_DEPTH", "Output", snow_avg_rho_depth, IOUtils::nothrow);

	//check if simple snow drift is enabled
	enable_simple_snow_drift = false;
	sn_cfg.getValue("SIMPLE_SNOW_DRIFT", "Alpine3D", enable_simple_snow_drift, IOUtils::nothrow);

	//create the vector of output grids
	std::vector<std::string> params = sn_cfg.get("GRIDS_PARAMETERS", "output");
	if (snow_temp_depth!=IOUtils::nodata) params.push_back("TSNOW");
	if (snow_avg_temp_depth!=IOUtils::nodata) params.push_back("TSNOW_AVG");
	if (snow_avg_rho_depth!=IOUtils::nodata) params.push_back("RHOSNOW_AVG");

	//handle the soil temperatures, runoff and snow densities
	io_cfg.getValue("SOIL_TEMPERATURE_DEPTHS", "Output", soil_temp_depths, IOUtils::nothrow);
	static const size_t n_tsoil_depths = SnGrids::TSOIL_MAX - SnGrids::TSOIL1 + 1;
	for (size_t ii=0; ii<std::min(soil_temp_depths.size(), n_tsoil_depths); ii++) {
		if(ii==n_tsoil_depths-1) {
			params.push_back("TSOIL_MAX");
		} else {
			params.push_back( "TSOIL"+mio::IOUtils::toString(ii+1) );
		}
	}
	sn_cfg.getValue("SOIL_RUNOFF_DEPTHS", "Output", soil_runoff_depths, IOUtils::nothrow);
	static const size_t n_runoffdepths = SnGrids::SOIL_RUNOFF_MAX - SnGrids::SOIL_RUNOFF1 + 1;
	for (size_t ii=0; ii<std::min(soil_runoff_depths.size(), n_runoffdepths); ii++) {
		if(ii==n_runoffdepths-1) {
			params.push_back("SOIL_RUNOFF_MAX");
		} else {
			params.push_back( "SOIL_RUNOFF"+mio::IOUtils::toString(ii+1) );
		}
	}
	io_cfg.getValue("SNOW_DENSITY_DEPTHS", "Output", snow_density_depths, IOUtils::nothrow);
	for (size_t ii=0; ii<snow_density_depths.size(); ii++) {
		params.push_back( "RHO"+mio::IOUtils::toString(ii+1) );
	}
	uniqueOutputGrids(params);
	initGrids(params, grids_not_computed_in_worker);

	const CurrentMeteo meteoPixel; //this is only necessary in order to have something for fillGrids()
	const SurfaceFluxes surfaceFlux; //this is only necessary in order to have something for fillGrids()

	for(size_t ii = 0; ii < SnowStationsCoord.size(); ++ii) {
		size_t ix = SnowStationsCoord.at(ii).first;
		size_t iy = SnowStationsCoord.at(ii).second;
		if (SnowpackInterfaceWorker::skipThisCell(landuse(ix,iy), dem(ix,iy))) { //skip nodata cells as well as water bodies, etc
			if (!pts_in.empty() && is_special(pts_in, ix, iy)){
				std::cout << "[W] POI (" << ix << "," << iy << ") will be skipped (nodata, water body, etc)" << std::endl;
			}
			continue;
		}

		if (SnowStations[ii]==NULL) continue; //for safety: skip cells initialized with NULL
		const SnowStation &snowPixel = *SnowStations[ii];
		fillGrids(ix, iy, meteoPixel, snowPixel, surfaceFlux);
		if (!pts_in.empty()) { //is it a special point?
			isSpecialPoint[ii] = is_special(pts_in, ix, iy);
		}
	}
}

SnowpackInterfaceWorker::~SnowpackInterfaceWorker()
{
	//sorry for this cryptic syntax, this is to guarantee execution order with the "," operator
	while (!SnowStations.empty())
		static_cast<void>((SnowStations.back()!=NULL)? delete SnowStations.back() : (void)0) , SnowStations.pop_back();
}

/******************************************************************************
 * Methods that have to do with output
 ******************************************************************************/

/** @brief Make sure all requested grids only appear once
 * @param output_grids vector of requeste grids to sort and filter
 */
void SnowpackInterfaceWorker::uniqueOutputGrids(std::vector<std::string>& output_grids)
{
	for (size_t ii = 0; ii<output_grids.size(); ++ii)
		IOUtils::toUpper(output_grids[ii]);

	std::sort(output_grids.begin(), output_grids.end());
	const std::vector<std::string>::iterator it = std::unique(output_grids.begin(), output_grids.end());
	output_grids.resize(std::distance(output_grids.begin(),it));
}

/** @brief Initialize and add to the grid map the requested grids
 * @param[in] params string representation of the grids to add
 * @param[in] grids_not_computed_in_worker vector of grids that are required but will be computed outside the workers
 */
void SnowpackInterfaceWorker::initGrids(std::vector<std::string> params,
                                        const std::vector<std::string>& grids_not_computed_in_worker)
{
	for (size_t ii = 0; ii<params.size(); ++ii) {
		IOUtils::toUpper(params[ii]); //make sure all parameters are upper case
		const size_t param_idx = SnGrids::getParameterIndex( params[ii] );
		const auto position = std::find(grids_not_computed_in_worker.begin(), grids_not_computed_in_worker.end(), params[ii]);
		if(position<grids_not_computed_in_worker.end()) continue;
		
		if (param_idx==IOUtils::npos)
			throw UnknownValueException("Unknown meteo grid '"+params[ii]+"' selected for gridded output", AT);

		const std::map< SnGrids::Parameters, mio::Grid2DObject >::const_iterator it = grids.find( static_cast<SnGrids::Parameters>(param_idx) );
		if (it==grids.end()) { //the parameter did not already exist, adding it
			Grid2DObject tmp_grid(dem, IOUtils::nodata);
			grids[ static_cast<SnGrids::Parameters>(param_idx) ] = tmp_grid;
		}
	}
}

void SnowpackInterfaceWorker::getOutputSpecialPoints(std::vector<SnowStation*>& ptr_snow_pixel,
                              std::vector<CurrentMeteo*>& ptr_meteo_pixel,
                              std::vector<SurfaceFluxes*>& ptr_surface_flux)
{
	for (size_t ii=0; ii<snow_pixel.size(); ii++) {
		ptr_snow_pixel.push_back( &snow_pixel[ii] );
		ptr_meteo_pixel.push_back( &meteo_pixel[ii] );
		ptr_surface_flux.push_back( &surface_flux[ii] );
	}
}

void SnowpackInterfaceWorker::clearSpecialPointsData()
{
	snow_pixel.clear();
	meteo_pixel.clear();
	surface_flux.clear();
}

/**
 * @brief method that returns SnowCover files for specific date.
 * @param snow_station vector that will be filled with the SnowStation data for each pixels
 */
void SnowpackInterfaceWorker::getOutputSNO(std::vector<SnowStation*>& snow_station) const
{
	for(size_t i = 0; i < SnowStationsCoord.size(); ++i) {
		size_t ix = SnowStationsCoord.at(i).first;
		size_t iy = SnowStationsCoord.at(i).second;
		if (SnowpackInterfaceWorker::skipThisCell(landuse(ix,iy), dem(ix,iy))) continue; //skip nodata cells as well as water bodies, etc
		const size_t index_SnowStation = i;
		if (SnowStations[index_SnowStation]==NULL) continue; //for safety: skipped cells were initialized with NULL
		#pragma omp critical (snow_station_lock)
		snow_station.push_back( SnowStations[index_SnowStation] );
	}
}

/****************************************************************************
 * Getters and Setters
 ****************************************************************************/

/**
 * @brief Method that the Master can search the neded data (in grids) from Worker (Pull from client)
 * @param param says which grid param the Master wants to have
 * @return the 2D output grid, which gives back the data to the master
 */
mio::Grid2DObject SnowpackInterfaceWorker::getGrid(const SnGrids::Parameters& param) const
{
	const std::map< SnGrids::Parameters, mio::Grid2DObject >::const_iterator it = grids.find(param);
	if (it==grids.end()) {
		throw UnknownValueException("Undeclared grid '"+SnGrids::getParameterName(param)+"' requested from SnowpackInterfaceWorker", AT);
	}

	return it->second;
}

/**
 * @brief Retrieve one point (ii,jj) from the specified grid.
 * @param param says which grid param the Master wants to have
 * @param ii ii index
 * @param jj jj index
 * @return grid value at point (ii,jj) for parameter param
 */
double& SnowpackInterfaceWorker::getGridPoint(const SnGrids::Parameters& param, const size_t& ii, const size_t& jj)
{
	const std::map< SnGrids::Parameters, mio::Grid2DObject >::const_iterator it( grids.find(param) );
	if (it==grids.end())
		throw UnknownValueException("Undeclared grid '"+SnGrids::getParameterName(param)+"' requested from SnowpackInterfaceWorker", AT);

	return grids[ param ](ii,jj);
}

/**
 * @brief Fill all the grids stored in the **grids** map of 2D grids with all the necessary values for the provided pixel.
 * Before calling this method, make sure that snowPixel is not NULL (as it could be for pixels that should be skipped)
 * @param ii ii index
 * @param jj jj index
 * @param meteoPixel meteorological forcing for the pixel
 * @param snowPixel canopy/snow/soil stratigraphy information
 * @param surfaceFlux surface fluxes for the pixel
 */
void SnowpackInterfaceWorker::fillGrids(const size_t& ii, const size_t& jj, const CurrentMeteo& meteoPixel, const SnowStation& snowPixel, const SurfaceFluxes& surfaceFlux)
{
	std::map< SnGrids::Parameters, mio::Grid2DObject >::iterator it;
	for (it=grids.begin(); it!=grids.end(); ++it) {
		double value = IOUtils::nodata;
		switch (it->first) {
			case SnGrids::ISWR_BELOW_CAN:
				value = (useEBalance && useCanopy)? meteoPixel.iswr : IOUtils::nodata; break;
			case SnGrids::HS:
				value = (snowPixel.cH - snowPixel.Ground) /  snowPixel.cos_sl; break; //slope2horiz
			case SnGrids::TSS:
				if (!useCanopy || snowPixel.Cdata->zdispl < 0.) {
					value = snowPixel.Ndata.back().T;
				} else { // radiative surface temperature above canopy
					value =  pow(snowPixel.Cdata->rlwrac/Constants::stefan_boltzmann, 0.25);
				}
				break;
			case SnGrids::TSG:
				value = snowPixel.Ndata[snowPixel.SoilNode].T; break;
			case SnGrids::TS0:
				value = (snowPixel.SoilNode>0)? snowPixel.Ndata[snowPixel.SoilNode].T : IOUtils::nodata; break;
			case SnGrids::TSNOW:
				value = getValueAtDepth(snowPixel, snow_temp_depth); break;
			case SnGrids::TSNOW_AVG:
				value = getAvgAtDepth(snowPixel, snow_temp_depth, true); break;
			case SnGrids::RHOSNOW_AVG:
				value = getAvgAtDepth(snowPixel, snow_temp_depth, false); break;
			case SnGrids::SWE:
				value = surfaceFlux.mass[SurfaceFluxes::MS_TOTALMASS] / snowPixel.cos_sl; break; //slope2horiz
			case SnGrids::RSNO: {
				const double hs = (snowPixel.cH - snowPixel.Ground);
				value = (hs>0.)? snowPixel.mass_sum / hs : IOUtils::nodata;
				break; }
			case SnGrids::TOP_ALB:
				if (!useCanopy || snowPixel.Cdata->zdispl < 0.) {
					value = snowPixel.Albedo;
				} else { // above canopy albedo
					value =  snowPixel.Cdata->totalalb;
				}
				break;
			case SnGrids::SURF_ALB:
				value = snowPixel.Albedo; break;
			case SnGrids::SP:
				value = (!snowPixel.Edata.empty())? snowPixel.Edata.back().sp : IOUtils::nodata; break;
			case SnGrids::RB:
				value = (!snowPixel.Edata.empty())? snowPixel.Edata.back().rb : IOUtils::nodata; break;
			case SnGrids::RG:
				value = (!snowPixel.Edata.empty())? snowPixel.Edata.back().rg : IOUtils::nodata; break;
			case SnGrids::N3:
				value = (!snowPixel.Edata.empty())? snowPixel.Edata.back().N3 : IOUtils::nodata; break;
			case SnGrids::MS_SNOWPACK_RUNOFF:
				value = surfaceFlux.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] / snowPixel.cos_sl; break;
			case SnGrids::MS_SURFACE_MASS_FLUX:
				value = surfaceFlux.mass[SurfaceFluxes::MS_SURFACE_MASS_FLUX] / snowPixel.cos_sl; break;
			case SnGrids::MS_SOIL_RUNOFF:
				value = surfaceFlux.mass[SurfaceFluxes::MS_SOIL_RUNOFF] / snowPixel.cos_sl; break;
			case SnGrids::MS_RAIN:
				value = surfaceFlux.mass[SurfaceFluxes::MS_RAIN]; break;
			case SnGrids::MS_HNW:
				value = surfaceFlux.mass[SurfaceFluxes::MS_HNW]; break;
			case SnGrids::MS_WIND:
				value = surfaceFlux.mass[SurfaceFluxes::MS_WIND] / snowPixel.cos_sl; break;
			case SnGrids::MS_WATER:
				value = surfaceFlux.mass[SurfaceFluxes::MS_WATER] / snowPixel.cos_sl; break;
			case SnGrids::MS_WATER_SOIL:
				value = surfaceFlux.mass[SurfaceFluxes::MS_WATER_SOIL] / snowPixel.cos_sl; break;
			case SnGrids::MS_ICE_SOIL:
				value = surfaceFlux.mass[SurfaceFluxes::MS_ICE_SOIL] / snowPixel.cos_sl; break;
			case SnGrids::SFC_SUBL:
				value = -surfaceFlux.mass[SurfaceFluxes::MS_SUBLIMATION] / snowPixel.cos_sl; break; //slope2horiz
			case SnGrids::STORE:
				value = store(ii,jj); break;
			case SnGrids::ERODEDMASS:
				value = erodedmass(ii,jj) / snowPixel.cos_sl; break;
			case SnGrids::MS_SNOW_DHS:
				value = (surfaceFlux.mass[SurfaceFluxes::MS_SNOW_DHS] != IOUtils::nodata) ? (M_TO_MM(surfaceFlux.mass[SurfaceFluxes::MS_SNOW_DHS])/snowPixel.cos_sl) : (IOUtils::nodata); break;
			case SnGrids::MS_SUBL_DHS:
				value = (surfaceFlux.mass[SurfaceFluxes::MS_SUBL_DHS] != IOUtils::nodata) ? (M_TO_MM(surfaceFlux.mass[SurfaceFluxes::MS_SUBL_DHS])/snowPixel.cos_sl) : (IOUtils::nodata); break;
			case SnGrids::MS_SETTLING_DHS:
				value = (surfaceFlux.mass[SurfaceFluxes::MS_SETTLING_DHS] != IOUtils::nodata) ? (M_TO_MM(surfaceFlux.mass[SurfaceFluxes::MS_SETTLING_DHS])/snowPixel.cos_sl) : (IOUtils::nodata); break;
			case SnGrids::MS_EROSION_DHS:
				value = (surfaceFlux.mass[SurfaceFluxes::MS_EROSION_DHS] != IOUtils::nodata) ? (M_TO_MM(surfaceFlux.mass[SurfaceFluxes::MS_EROSION_DHS])/snowPixel.cos_sl) : (IOUtils::nodata); break;
			case SnGrids::GLACIER:
				value = (!snowPixel.isGlacier(true))? 1. : IOUtils::nodata; break; //glaciated pixels receive IOUtils::nodata
			case SnGrids::GLACIER_EXPOSED:
				value = (!snowPixel.isGlacier(false))? 1. : IOUtils::nodata; break; //glaciated pixels receive IOUtils::nodata
			case SnGrids::ET:
				value = -(surfaceFlux.mass[SurfaceFluxes::MS_SUBLIMATION]+surfaceFlux.mass[SurfaceFluxes::MS_EVAPORATION])/snowPixel.cos_sl; //slope2horiz
				// Add part from Canopy
				value += useCanopy?(snowPixel.Cdata->transp+snowPixel.Cdata->intevap)/snowPixel.cos_sl:0; //slope2horiz
				break;
			default:
				if (it->first>=SnGrids::TSOIL1 && it->first<=SnGrids::TSOIL_MAX) //dealing with soil temperatures
				{
					value = (soil_temp_depths.empty())? IOUtils::nodata : getSoilTemperature(snowPixel,
                                                                       soil_temp_depths[ it->first - SnGrids::TSOIL1 ]);
				}
				else if (it->first>=SnGrids::SOIL_RUNOFF1 && it->first<=SnGrids::SOIL_RUNOFF_MAX) //dealing with soil runoff
				{
					value = (soil_runoff_depths.empty())? IOUtils::nodata : getSoilRunoff(snowPixel, surfaceFlux,
                                                                  soil_runoff_depths[ it->first - SnGrids::SOIL_RUNOFF1 ])/ snowPixel.cos_sl;
				}
				else if (it->first>=SnGrids::RHO1 && it->first<=SnGrids::RHO5) //dealing with snow densities
				{
					value = (snow_density_depths.empty())? IOUtils::nodata : getSnowDensityDepth(snowPixel, snow_density_depths[ it->first - SnGrids::RHO1 ]);
				}
				else
				{
					throw InvalidArgumentException("Invalid parameter requested " + SnGrids::getParameterName( it->first ), AT);
				}
		}
		it->second(ii,jj) = value;
	}
}

/**
 * @brief method which prepares all data for simulation and then access correctly the
 * Snowpack model interfaces.
 * @param date current simulation time step
 * @param psum precipitation grid (kg m-2)
 * @param psum_ph precipitation phase grid (between 0 and 1)
 * @param psum_tech technical precipitation grid (kg m-2)
 * @param rh relative humidity grid (% or 1)
 * @param ta air temperature grid (K)
 * @param vw wind velocity grid (m s-1)
 * @param dw wind direction grid (degrees north)
 * @param mns map of the Precipitation (mm/h) HACK get this map only if per pull from Master if Drift is used !!
 * @param shortwave incoming shortwave radiation grid (W m-2)
 * @param diffuse diffuse radiation from the sky grid (W m-2)
 * @param longwave incoming longwave grid (W m-2)
 * @param solarElevation solar elevation (in dec)
 */
void SnowpackInterfaceWorker::runModel(const mio::Date &date,
                                       const mio::Grid2DObject &psum,
                                       const mio::Grid2DObject &psum_ph,
                                       const mio::Grid2DObject &psum_tech,
                                       const mio::Grid2DObject &rh,
                                       const mio::Grid2DObject &ta,
                                       const mio::Grid2DObject &tsg,
                                       const mio::Grid2DObject &vw,
                                       const mio::Grid2DObject &vw_drift,
                                       const mio::Grid2DObject &dw,
                                       const mio::Grid2DObject &mns,
                                       const mio::Grid2DObject &shortwave,
                                       const mio::Grid2DObject &diffuse,
                                       const mio::Grid2DObject &longwave,
                                       const double solarElevation)
{
	const Meteo::ATM_STABILITY USER_STABILITY = meteo.getStability();
	const std::string bcu_watertransportmodel_snow = sn_cfg.get("WATERTRANSPORTMODEL_SNOW", "SnowpackAdvanced");
	const std::string bcu_reduce_n_elements = sn_cfg.get("REDUCE_N_ELEMENTS", "SnowpackAdvanced");
	const std::string bcu_adjust_height_of_meteo= sn_cfg.get("ADJUST_HEIGHT_OF_METEO_VALUES", "SnowpackAdvanced");
	const std::string bcu_adjust_height_of_wind = sn_cfg.get("ADJUST_HEIGHT_OF_WIND_VALUE", "SnowpackAdvanced");

	CurrentMeteo meteoPixel(sn_cfg);
	meteoPixel.date = date;
	meteoPixel.elev = solarElevation*Cst::to_rad; //HACK: Snowpack uses RAD !!!!!

	if (enable_simple_snow_drift) {
		// reset eroded mass grid
		erodedmass.set(erodedmass, 0.);
	}

	// make SN calculations....
	for(size_t i = 0; i < SnowStationsCoord.size(); ++i) {
		size_t ix = SnowStationsCoord.at(i).first;
		size_t iy = SnowStationsCoord.at(i).second;
		if (SnowpackInterfaceWorker::skipThisCell(landuse(ix,iy), dem(ix,iy))) continue; //skip nodata cells as well as water bodies, etc
		const size_t index_SnowStation = i;
		if (SnowStations[index_SnowStation]==NULL) continue; //for safety: skipped cells were initialized with NULL
		SnowStation &snowPixel = *SnowStations[index_SnowStation];

		const bool isGlacier = snowPixel.isGlacier(false);

		//In case of ice and firn pixels, use BUCKET model for water transport:
		const int land = (round_landuse(landuse(ix,iy)) - 10000) / 100;
		if (land==13 || land==14) {
			sn_cfg.addKey("WATERTRANSPORTMODEL_SNOW", "SnowpackAdvanced", "BUCKET");
		}
		// In case of glacier pixel, remove meteo height correction and try to merge elemnts
		if (land==14) {
			sn_cfg.addKey("REDUCE_N_ELEMENTS", "SnowpackAdvanced", "TRUE");
			sn_cfg.addKey("ADJUST_HEIGHT_OF_METEO_VALUES","SnowpackAdvanced","FALSE");
			sn_cfg.addKey("ADJUST_HEIGHT_OF_WIND_VALUE","SnowpackAdvanced","FALSE");
		}
		// Set curent meteo variables from 2D fields to single pixel
		const double previous_albedo = getGridPoint(SnGrids::TOP_ALB, ix, iy);
		meteoPixel.rh = rh(ix,iy);
		meteoPixel.ta = ta(ix,iy);
		meteoPixel.vw = vw(ix,iy);
		meteoPixel.vw_drift = (vw_drift(ix,iy) > 0.) ? (vw_drift(ix,iy)) : (0.);	//negative vw_drift is used by the Simple Snow Drift to store the positive sx values (sheltered pixels)
		meteoPixel.dw = dw(ix,iy);
		meteoPixel.iswr = shortwave(ix,iy);
		meteoPixel.rswr = previous_albedo*meteoPixel.iswr;
		meteoPixel.tss = snowPixel.Ndata[snowPixel.getNumberOfElements()].T; //we use previous timestep value
		meteoPixel.ts0 = (snowPixel.SoilNode>0)? snowPixel.Ndata[snowPixel.SoilNode].T : tsg(ix,iy); //we use previous timestep value
		meteoPixel.ea = Atmosphere::blkBody_Emissivity(longwave(ix,iy), meteoPixel.ta); //to be consistent with Snowpack
		meteoPixel.psum_ph = psum_ph(ix,iy);
		meteoPixel.psum_tech = psum_tech(ix, iy);
		meteoPixel.hs = IOUtils::nodata;
		if ( (meteoPixel.tss != IOUtils::nodata && meteoPixel.tss<=100) || (meteoPixel.ts0 != IOUtils::nodata && meteoPixel.ts0<=100) ) {
			cout << "[E] Pixel (" << ix << "," << iy << ") too cold! tss=" << meteoPixel.tss << " ts0=" << meteoPixel.ts0 << std::endl;
		}

		// Now determine perpendicular to slope mass balance from drift or precip
		// Note that psum must be in mm/h and SNOWPACK will converted it to mm/CALCULTION_STEP_LENGTH
		if ( useDrift && mns(ix,iy)!=IOUtils::nodata ) {
			double drift_mass = mns(ix,iy);
			// For extreme terrain, the drift might go crazy: limit deposition / erosion
			if (fabs(drift_mass)>37.) {
				cout << "crazy drift at (" << ix << "," << iy << ") = " << drift_mass << " mm/h\n";
				drift_mass = (drift_mass>0.)? 37. : -37.; //set to deposition limit
			}

			meteoPixel.psum = drift_mass;
		} else {
			meteoPixel.psum = psum(ix,iy) * snowPixel.cos_sl; //horiz2slope
		}

		if (useEBalance) meteoPixel.diff = diffuse(ix,iy);
		// Reset Canopy Surface Data to 0 before calling Meteo and Canopy module
		if (useCanopy) snowPixel.Cdata->initializeSurfaceExchangeData();
		//if the current pixel contains soil layers, then activate soil modeling in Snowpack
		sn.setUseSoilLayers( snowPixel.hasSoilLayers() ); //TODO: when Snowpack runs, it should check it

		// exposed glacier special case
		if (isGlacier) {
			const std::string tmp_sw_mode = sn_cfg.get("SW_MODE", "Snowpack");
			if (tmp_sw_mode == "BOTH") {
				//switch to glacier albedo (when sw_mode != BOTH, the calculation internally relies on SnLaws and should not be overwritten here!)
				const std::string tmp_variant = sn_cfg.get("VARIANT", "SnowpackAdvanced");
				snowPixel.Albedo = ((tmp_variant == "POLAR" || tmp_variant == "ANTARCTICA") ? (Constants::blueice_albedo) : (Constants::glacier_albedo));
			}
			if (meteoPixel.ta>IOUtils::C_TO_K(5.)) {
				//switch to STABLE atmosphere on glacier if TA>5°C
				meteo.setStability(Meteo::MO_HOLTSLAG);
			}
		}
		bool adjust_height_of_wind_value;
		sn_cfg.getValue("ADJUST_HEIGHT_OF_WIND_VALUE", "SnowpackAdvanced", adjust_height_of_wind_value);
		//compute ustar, psi_s, z0
		meteo.compMeteo(meteoPixel, snowPixel, true, adjust_height_of_wind_value);
		SurfaceFluxes surfaceFlux;
		snowPixel.reset_water_fluxes();
		// run snowpack model itself
		double dIntEnergy = 0.; //accumulate the dIntEnergy over the snowsteps
		const unsigned int nr_snowsteps = (unsigned int)(dt_main/M_TO_S(calculation_step_length));

		for (unsigned int snowsteps = 0; snowsteps < nr_snowsteps; snowsteps++) {
			/* Update the store variable */
			/* david: why += ? isnt store reset every timestep ? */
			/* Michi: This is exactly the point, store is only set to zero when enough precipitation
			has accumulated to create at least one new (finite element) layer */

			if (snowsteps == 0) meteoPixel.psum /= nr_snowsteps;
			store(ix,iy) += meteoPixel.psum;

			// Reset fluxes, etc.
			snowPixel.ErosionMass = 0.;
			snowPixel.hn = 0.;
			snowPixel.rho_hn = 0.;

			try {
				BoundCond Bdata;
				sn.runSnowpackModel(meteoPixel, snowPixel, store(ix,iy), Bdata, surfaceFlux);
				surfaceFlux.collectSurfaceFluxes(Bdata, snowPixel, meteoPixel);
				surfaceFlux.mass[SurfaceFluxes::MS_HNW] += (snowPixel.hn * snowPixel.rho_hn) / snowPixel.cos_sl;
				surfaceFlux.mass[SurfaceFluxes::MS_WIND] += snowPixel.ErosionMass;
				dIntEnergy += snowPixel.dIntEnergy; //it is reset at every new call to runSnowpackModel
				meteoPixel.hs = snowPixel.cH - snowPixel.Ground; //do not reproject here, otherwise Snowpack outputs would get messed up
				if (enable_simple_snow_drift) erodedmass(ix,iy) += snowPixel.ErosionMass; //store the eroded mass
			} catch (const std::bad_alloc&) { //don't try anything fancy when running low on memory
				const int lus =(int)floor( landuse(ix,iy));
				const double slope2horiz = (1. / snowPixel.cos_sl);
				const double snow = (snowPixel.cH - snowPixel.Ground) * slope2horiz;
				cout << "[E] Could not allocate memory in Snowpack for cell (" << ix << "," << iy << ") LUS=" << lus << " ";
				cout << "with " << std::fixed << std::setprecision(4) <<  snow << " m snow in " << snowPixel.getNumberOfElements()-snowPixel.SoilNode << "/" << snowPixel.getNumberOfElements() << " elements.\n";
				fflush( stdout );
				throw;
			} catch (const std::exception& e) {
				const int lus =(int)floor( landuse(ix,iy));
				const double slope2horiz = (1. / snowPixel.cos_sl);
				const double snow = (snowPixel.cH - snowPixel.Ground) * slope2horiz;
				if (nr_snowsteps > 1) {
					surfaceFlux.multiplyFluxes(1./nr_snowsteps);
					if (useCanopy)
						snowPixel.Cdata->multiplyFluxes(1./nr_snowsteps);
					meteoPixel.psum *= nr_snowsteps;
					snowPixel.dIntEnergy = dIntEnergy;
				}

				ostringstream ss;
				ss << "[E] Snowpack exception: " << e.what() << "\n";
				ss << "[E] at cell (" << ix << "," << iy << ") LUS=" << lus << " ";
				ss << "with " << std::fixed << std::setprecision(4) << snow << " m snow in " << snowPixel.getNumberOfElements()-snowPixel.SoilNode << "/" << snowPixel.getNumberOfElements() << " elements.\n";
				if (isSpecialPoint[index_SnowStation])
					gatherSpecialPoints(meteoPixel, snowPixel, surfaceFlux); //gather special point data, in order to get as much information as possible
				throw IOException(ss.str(), AT);
			}
		}

		//some variables are now wrong if we ran multiple Snowpack steps -> recompute them!
		if (nr_snowsteps > 1) {
			surfaceFlux.multiplyFluxes(1./nr_snowsteps);
			if (useCanopy)
				snowPixel.Cdata->multiplyFluxes(1./nr_snowsteps);
			meteoPixel.psum *= nr_snowsteps;
			snowPixel.dIntEnergy = dIntEnergy;
		}

		//switch stability back to normal if it was changed
		if (meteo.getStability()!=USER_STABILITY) meteo.setStability(USER_STABILITY);
		//if the glacier is still exposed, force the albedo back to glacier albedo
		if (isGlacier) {
			const std::string tmp_sw_mode = sn_cfg.get("SW_MODE", "Snowpack");
			if (tmp_sw_mode == "BOTH") {
				//switch to glacier albedo (when sw_mode != BOTH, the calculation internally relies on SnLaws and should not be overwritten here!)
				const std::string tmp_variant = sn_cfg.get("VARIANT", "SnowpackAdvanced");
				surfaceFlux.pAlbedo = snowPixel.Albedo = ((tmp_variant == "POLAR" || tmp_variant == "ANTARCTICA") ? (Constants::blueice_albedo) : (Constants::glacier_albedo));
			}
		}
		if (!std::isfinite( getGridPoint(SnGrids::TOP_ALB, ix, iy) )) {
			//if the albedo is nan, infinity, etc reset it to its previous
			//value to try to rescue the pixel...
			cerr << "[E] pixel (" << ix << "," << iy << ") found with a nan/infinit albedo ["<<  getGridPoint(SnGrids::TOP_ALB, ix, iy) <<"]; reseting to " << previous_albedo << std::endl;
			getGridPoint(SnGrids::TOP_ALB, ix, iy) = previous_albedo;
		}

		try{
			stability.checkStability(meteoPixel, snowPixel);
		} catch(...) {
#ifndef SNOWPACK_CORE
			snowPixel.S_4 = Stability::max_stability; //nothing else to do...
#endif
		}

		surfaceFlux.mass[SurfaceFluxes::MS_TOTALMASS] = 0.0;
		const std::vector<ElementData> EMS( snowPixel.Edata );
		for (size_t e=0; e<EMS.size(); e++) {
			if (EMS[e].theta[SOIL] <= 0.) {
				surfaceFlux.mass[SurfaceFluxes::MS_TOTALMASS] += EMS[e].M;
			}
		}

		// Output special points and grids
		if (isSpecialPoint[index_SnowStation]) gatherSpecialPoints(meteoPixel, snowPixel, surfaceFlux);
		fillGrids(ix, iy, meteoPixel, snowPixel, surfaceFlux);

		//Restore original water transport scheme that has been changed for ice & firn
		if (land==13 || land==14) {
			sn_cfg.addKey("WATERTRANSPORTMODEL_SNOW", "SnowpackAdvanced", bcu_watertransportmodel_snow);
		}
		//Restore original keys that were modified for glacier pixels
		if (land==14) {
			sn_cfg.addKey("REDUCE_N_ELEMENTS", "SnowpackAdvanced", bcu_reduce_n_elements);
			sn_cfg.addKey("ADJUST_HEIGHT_OF_METEO_VALUES", "SnowpackAdvanced", bcu_adjust_height_of_meteo);
			sn_cfg.addKey("ADJUST_HEIGHT_OF_WIND_VALUE", "SnowpackAdvanced", bcu_adjust_height_of_wind);
		}
	}
}

void SnowpackInterfaceWorker::grooming(const mio::Date &current_date, const mio::Grid2DObject &grooming_map)
{
	for(size_t ii = 0; ii < SnowStationsCoord.size(); ++ii) {
		const size_t ix = SnowStationsCoord.at(ii).first;
		const size_t iy = SnowStationsCoord.at(ii).second;
		if (SnowpackInterfaceWorker::skipThisCell(landuse(ix,iy), dem(ix,iy))) continue; //skip nodata cells as well as water bodies, etc
		if (SnowStations[ii]==NULL) continue; //for safety: skipped cells were initialized with NULL

		if (grooming_map(ix, iy)==IOUtils::nodata || grooming_map(ix, iy)==0) continue;

		if (sn_techsnow.prepare(current_date)) sn_techsnow.preparation(*SnowStations[ii]);
	}
}

void SnowpackInterfaceWorker::gatherSpecialPoints(const CurrentMeteo& meteoPixel, const SnowStation& snowPixel, const SurfaceFluxes& surfaceFlux)
{
	meteo_pixel.push_back( meteoPixel );
	snow_pixel.push_back( snowPixel );
	surface_flux.push_back( surfaceFlux );
}

/**
 * @brief optimised way to round landuse
 * @param landuse_dbl is the landuse to round
 */
int SnowpackInterfaceWorker::round_landuse(const double& landuse_dbl)
{
	return (int)(landuse_dbl + 0.0001);
}

/**
 * @brief check if a cell should be simulated or skipped
 * @param landuse_val land use parameter for this pixel
 * @param dem_val dem altitude for this pixel
 */
bool SnowpackInterfaceWorker::skipThisCell(const double& landuse_val, const double& dem_val)
{
	//determines from the landuse and dem if this cell has to be included in the computation
	//landuse codes are PREVAH codes
	if (landuse_val==IOUtils::nodata) return true;

	const int land = (SnowpackInterfaceWorker::round_landuse(landuse_val) - 10000) / 100;
	if (land==9 || land==10 || land==12 || land==16 || land==17 || land>=30) return true; //undefined
	if (land==1) return true;//water

	if (dem_val==IOUtils::nodata) return true; //no DEM data

	return false;
}

bool SnowpackInterfaceWorker::is_special(const std::vector< std::pair<size_t,size_t> >& pts_in, const size_t& ix, const size_t& iy)
{
	for (size_t ii=0; ii<pts_in.size(); ii++) {
		if ((pts_in[ii].first == ix) && (pts_in[ii].second == iy)) return true;
	}

	return false;
}

/**
 * @brief method called by SnowpackInterface to retrieve all snow_pixels, to read the lateral flow variable
 * @param ptr_snow_pixel to vector of SnowStations, where the variable SlopeParFlux contains the lateral flow
 */
void SnowpackInterfaceWorker::getLateralFlow(std::vector<SnowStation*>& ptr_snow_pixel)
{
	for(size_t ii = 0; ii < SnowStationsCoord.size(); ++ii) {
		#pragma omp critical (snow_station_lock)
		ptr_snow_pixel.push_back( SnowStations[ii] );
	}
}

/**
 * @brief method called by SnowpackInterface to send back the source/sink term for treating lateral flow
 * @param ptr_snow_pixel to vector of SnowStations, in which the lwc_source variable is set to contain the lateral flow.
 */
void SnowpackInterfaceWorker::setLateralFlow(const std::vector<SnowStation*>& ptr_snow_pixel)
{
	for(size_t ii = 0; ii < SnowStationsCoord.size(); ++ii) {
		if (SnowStations[ii]!=NULL) {
			for(size_t n=0; n<SnowStations[ii]->getNumberOfElements(); n++) {
				SnowStations[ii]->Edata[n].lwc_source = ptr_snow_pixel[ii]->Edata[n].lwc_source;
			}
		}
	}
}
