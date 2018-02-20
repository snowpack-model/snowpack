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
#include <alpine3d/SnowpackInterfaceWorker.h>
#include <alpine3d/AlpineMain.h>

using namespace std;
using namespace mio;

/************************************************************
 * static section                                           *
 ************************************************************/
const size_t SnGrids::nrOfParameters =  SnGrids::lastparam - SnGrids::firstparam + 1;
std::vector<std::string> SnGrids::paramname;
const bool SnGrids::__init = SnGrids::initStaticData();

bool SnGrids::initStaticData()
{
	//the order must be the same as in the enum
	paramname.push_back("TA");
	paramname.push_back("RH");
	paramname.push_back("VW");
	paramname.push_back("ISWR");
	paramname.push_back("ISWR_DIFF");
	paramname.push_back("ISWR_DIR");
	paramname.push_back("ILWR");
	paramname.push_back("HS");
	paramname.push_back("PSUM");
	paramname.push_back("PSUM_PH");
	paramname.push_back("TSG");
	paramname.push_back("TSS");
	paramname.push_back("TS0");
	paramname.push_back("TSNOW");
	paramname.push_back("TSNOW_AVG");
	paramname.push_back("RHOSNOW_AVG");
	paramname.push_back("TSOIL");
	paramname.push_back("SWE");
	paramname.push_back("RSNO");
	paramname.push_back("TOP_ALB");
	paramname.push_back("SURF_ALB");
	paramname.push_back("SP");
	paramname.push_back("RB");
	paramname.push_back("RG");
	paramname.push_back("N3");
	paramname.push_back("MS_SNOWPACK_RUNOFF");
	paramname.push_back("MS_SOIL_RUNOFF");
	paramname.push_back("SFC_SUBL");
	paramname.push_back("STORE");
	paramname.push_back("GLACIER");
	paramname.push_back("GLACIER_EXPOSED");
	
	if (paramname.size()!=(SnGrids::lastparam+1))
		throw IOException("Wrong number of string representations for the SnGrids parameters! You forgot to update the \"paramname\" vector.", AT);

	return true;
}

const std::string& SnGrids::getParameterName(const size_t& parindex)
{
	if (parindex >= SnGrids::nrOfParameters)
		throw IndexOutOfBoundsException("Trying to get name for parameter that does not exist", AT);

	return paramname[parindex];
}

size_t SnGrids::getParameterIndex(const std::string& parname)
{
	for (size_t ii=0; ii<SnGrids::nrOfParameters; ii++) {
		if (paramname[ii] == parname) return ii;
	}

	return IOUtils::npos; //parameter not a part of SnGrids
}


/************************************************************
 * SnowpackInterfaceWorker                                           *
 ************************************************************/
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
 * @param io_cfg configuration to pass to Snowpack
 * @param dem_in gives the demographic Data. Also tetermines size and position of the slice
 * @param landuse_in gives the landuse Data. Also tetermines size and position of the landuse for slice
 * @param pts_in gives the spezial points. For this points more output is done then for the others. Calcualtion is the same.
 * @param snow_stations gives a vector of pointers to the SnowStation objects relevant for this thread
 * @param offset_in gives the offsett on X for this slice (needed to read data and error messages)
 */
SnowpackInterfaceWorker::SnowpackInterfaceWorker(const mio::Config& io_cfg,
                                                 const mio::DEMObject& dem_in,
                                                 const mio::Grid2DObject& landuse_in,
                                                 const std::vector< std::pair<size_t,size_t> >& pts_in,
                                                 const std::vector<SnowStation*>& snow_stations,
                                                 const size_t offset_in)
 : sn_cfg(io_cfg), sn(sn_cfg), meteo(sn_cfg), stability(sn_cfg, false), dem(dem_in),
   dimx(dem.getNx()), dimy(dem.getNy()),  offset(offset_in), SnowStations(snow_stations), isSpecialPoint(snow_stations.size(), false), 
   landuse(landuse_in), store(dem_in, 0.), grids(), snow_pixel(), meteo_pixel(), surface_flux(), 
   calculation_step_length(0.), height_of_wind_value(0.),
   soil_temp_depth(IOUtils::nodata), snow_temp_depth(IOUtils::nodata), snow_avg_temp_depth(IOUtils::nodata), snow_avg_rho_depth(IOUtils::nodata),
   useDrift(false), useEBalance(false), useCanopy(false)
{
	if (SnowStations.size() != (dimx*dimy)) throw IOException("Initial snow data does not match the provided grid size", AT);
	
	sn_cfg.getValue("CALCULATION_STEP_LENGTH", "Snowpack", calculation_step_length);
	sn_cfg.getValue("HEIGHT_OF_WIND_VALUE", "Snowpack", height_of_wind_value); //currently unused
	sn_cfg.getValue("CANOPY", "Snowpack", useCanopy);
	io_cfg.getValue("SOIL_TEMPERATURE_DEPTH", "Output", soil_temp_depth, IOUtils::nothrow);
	io_cfg.getValue("SNOW_TEMPERATURE_DEPTH", "Output", snow_temp_depth, IOUtils::nothrow);
	io_cfg.getValue("SNOW_AVG_TEMPERATURE_DEPTH", "Output", snow_avg_temp_depth, IOUtils::nothrow);
	io_cfg.getValue("SNOW_AVG_DENSITY_DEPTH", "Output", snow_avg_rho_depth, IOUtils::nothrow);
	
	//create the vector of output grids
	std::vector<std::string> params = sn_cfg.get("GRIDS_PARAMETERS", "output");
	if (soil_temp_depth!=IOUtils::nodata) params.push_back("TSOIL");
	if (snow_temp_depth!=IOUtils::nodata) params.push_back("TSNOW");
	if (snow_avg_temp_depth!=IOUtils::nodata) params.push_back("TSNOW_AVG");
	if (snow_avg_rho_depth!=IOUtils::nodata) params.push_back("RHOSNOW_AVG");
	initGrids(params);
	
	const CurrentMeteo meteoPixel; //this is only necessary in order to have something for fillGrids()
	const SurfaceFluxes surfaceFlux; //this is only necessary in order to have something for fillGrids()
	
	for (size_t iy = 0; iy < dimy; iy++) {
		for (size_t ix = 0; ix < dimx; ix++) {
			if (SnowpackInterfaceWorker::skipThisCell(landuse(ix,iy), dem(ix,iy))) { //skip nodata cells as well as water bodies, etc
				if (!pts_in.empty() && is_special(pts_in, ix, iy))
					std::cout << "[W] POI (" << ix+offset << "," << iy << ") will be skipped (nodatat, water body, etc)" << std::endl;
				continue;
			}
			
			const size_t index_SnowStation = ix + dem.getNx()*iy;
			if (SnowStations[index_SnowStation]==NULL) continue; //for safety: skip cells initialized with NULL
			const SnowStation &snowPixel = *SnowStations[index_SnowStation];
			fillGrids(ix, iy, meteoPixel, snowPixel, surfaceFlux);
			
			if (!pts_in.empty()) { //is it a special point?
				isSpecialPoint[index_SnowStation] = is_special(pts_in, ix, iy);
			}
		}
	}
}

SnowpackInterfaceWorker::~SnowpackInterfaceWorker()
{
	//sorry for this cryptic syntax, this is to guarantee execution order with the "," operator
	while (!SnowStations.empty()) 
		((SnowStations.back()!=NULL)? delete SnowStations.back() : (void)0) , SnowStations.pop_back();
}

/******************************************************************************
 * Methods that have to do with output
 ******************************************************************************/

/** @brief Initialize and add to the grid map the requested grids
 * @param params string representation of the grids to add
 */
void SnowpackInterfaceWorker::initGrids(std::vector<std::string>& params)
{
	for (size_t ii = 0; ii<params.size(); ++ii) {
		IOUtils::toUpper(params[ii]); //make sure all parameters are upper case
		
		const size_t param_idx = SnGrids::getParameterIndex( params[ii] );
		if (param_idx==IOUtils::npos)
			throw UnknownValueException("Unknow meteo grid '"+params[ii]+"' selected for gridded output", AT);
		
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
	for (size_t iy=0;iy<dimy;iy++) {
		for (size_t ix=0;ix<dimx;ix++) {
			if (SnowpackInterfaceWorker::skipThisCell(landuse(ix,iy), dem(ix,iy))) continue; //skip nodata cells as well as water bodies, etc
			const size_t index_SnowStation = dem.getNy()*ix + iy;

			if (SnowStations[index_SnowStation]==NULL) continue; //for safety: skipped cells were initialized with NULL
			#pragma omp critical (snow_station_lock)
			snow_station.push_back( SnowStations[index_SnowStation] );
		}
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
			case SnGrids::TA:
				value = meteoPixel.ta; break;
			case SnGrids::RH:
				value = meteoPixel.rh; break;
			case SnGrids::VW:
				value = meteoPixel.vw; break;
			case SnGrids::ISWR:
				value = meteoPixel.iswr; break;
			case SnGrids::ISWR_DIFF:
				value = (useEBalance)? meteoPixel.diff : IOUtils::nodata; break;
			case SnGrids::ISWR_DIR:
				value = (useEBalance)? meteoPixel.iswr - meteoPixel.diff : IOUtils::nodata; break;
			case SnGrids::ILWR:
				value = Atmosphere::blkBody_Radiation(meteoPixel.ea, meteoPixel.ta); break;
			case SnGrids::HS:
				value = (snowPixel.cH - snowPixel.Ground) /  snowPixel.cos_sl; break; //slope2horiz
			case SnGrids::PSUM:
				value = meteoPixel.psum; break;
			case SnGrids::PSUM_PH:
				value = meteoPixel.psum_ph; break;
			case SnGrids::TSS:
				if (!useCanopy || snowPixel.Cdata.zdispl < 0.) {
					value = snowPixel.Ndata.back().T;
				} else { // radiative surface temperature above canopy
					value =  pow(snowPixel.Cdata.rlwrac/Constants::stefan_boltzmann, 0.25);
				}
				break;
			case SnGrids::TS0:
				value = (snowPixel.SoilNode>0)? snowPixel.Ndata[snowPixel.SoilNode].T : IOUtils::nodata; break;
			case SnGrids::TSNOW:
				value = getValueAtDepth(snowPixel, snow_temp_depth); break;
			case SnGrids::TSNOW_AVG:
				value = getAvgAtDepth(snowPixel, snow_temp_depth, true); break;
			case SnGrids::RHOSNOW_AVG:
				value = getAvgAtDepth(snowPixel, snow_temp_depth, false); break;
			case SnGrids::TSOIL:
				value = (soil_temp_depth != IOUtils::nodata)? getSoilTemperature(snowPixel, soil_temp_depth) : IOUtils::nodata; break;
			case SnGrids::SWE:
				value = surfaceFlux.mass[SurfaceFluxes::MS_TOTALMASS] /  snowPixel.cos_sl; break; //slope2horiz
			case SnGrids::RSNO: {
				const double hs = (snowPixel.cH - snowPixel.Ground);
				value = (hs>0.)? snowPixel.mass_sum / hs : IOUtils::nodata; 
				break; }
			case SnGrids::TOP_ALB:
				if (!useCanopy || snowPixel.Cdata.zdispl < 0.) {
					value = snowPixel.Albedo;
				} else { // above canopy albedo
					value =  snowPixel.Cdata.totalalb;
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
				value = surfaceFlux.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF]; break;
			case SnGrids::MS_SOIL_RUNOFF:
				value =surfaceFlux.mass[SurfaceFluxes::MS_SOIL_RUNOFF]; break;
			case SnGrids::SFC_SUBL:
				value = -surfaceFlux.mass[SurfaceFluxes::MS_SUBLIMATION] /  snowPixel.cos_sl; break; //slope2horiz
			case SnGrids::STORE:
				value = store(ii,jj); break;
			case SnGrids::GLACIER:
				value = (!snowPixel.isGlacier(true))? 1. : IOUtils::nodata; break; //glaciated pixels receive IOUtils::nodata
			case SnGrids::GLACIER_EXPOSED:
				value = (!snowPixel.isGlacier(false))? 1. : IOUtils::nodata; break; //glaciated pixels receive IOUtils::nodata
			default:
				throw  UnknownValueException("Unknown/unsupported grid requested", AT);
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
 * @param rh relative humidity grid (% or 1)
 * @param ta air temperature grid (K)
 * @param vw wind velocity grid (m s-1)
 * @param mns map of the Precipitation (mm/h) HACK get this map only if per pull from Master if Drift is used !!
 * @param shortwave incoming shortwave radiation grid (W m-2)
 * @param diffuse diffuse radiation from the sky grid (W m-2)
 * @param longwave incoming longwave grid (W m-2)
 * @param solarElevation solar elevation (in dec)
 */
void SnowpackInterfaceWorker::runModel(const mio::Date &date,
                                       const mio::Grid2DObject &psum,
                                       const mio::Grid2DObject &psum_ph,
                                       const mio::Grid2DObject &rh,
                                       const mio::Grid2DObject &ta,
                                       const mio::Grid2DObject &vw,
                                       const mio::Grid2DObject &mns,
                                       const mio::Grid2DObject &shortwave,
                                       const mio::Grid2DObject &diffuse,
                                       const mio::Grid2DObject &longwave,
                                       const double solarElevation)
{
	const Meteo::ATM_STABILITY USER_STABILITY = meteo.getStability();
	const std::string bcu_watertransportmodel_snow = sn_cfg.get("WATERTRANSPORTMODEL_SNOW", "SnowpackAdvanced");
	const std::string bcu_watertransportmodel_soil = sn_cfg.get("WATERTRANSPORTMODEL_SOIL", "SnowpackAdvanced");
	
	CurrentMeteo meteoPixel(sn_cfg);
	meteoPixel.date = date;
	meteoPixel.elev = solarElevation*Cst::to_rad; //HACK: Snowpack uses RAD !!!!!
	if (soil_temp_depth!=IOUtils::nodata) {
		meteoPixel.zv_ts.push_back( soil_temp_depth );
		meteoPixel.ts.push_back( IOUtils::nodata );
	}

	// make SN calculations....
	for (size_t iy=0; iy<dimy; iy++) {
		for (size_t ix=0; ix<dimx; ix++) {
			if (SnowpackInterfaceWorker::skipThisCell(landuse(ix,iy), dem(ix,iy))) continue; //skip nodata cells as well as water bodies, etc
			const size_t index_SnowStation = ix + dem.getNx()*iy;
			if (SnowStations[index_SnowStation]==NULL) continue; //for safety: skipped cells were initialized with NULL
			SnowStation &snowPixel = *SnowStations[index_SnowStation];
			const bool isGlacier = snowPixel.isGlacier(false);

			//In case of ice and firn pixels, use BUCKET model for water transport:
			const int land = (round_landuse(landuse(ix,iy)) - 10000) / 100;
			if (land==13 || land==14) {
				sn_cfg.addKey("WATERTRANSPORTMODEL_SNOW", "SnowpackAdvanced", "BUCKET");
				sn_cfg.addKey("WATERTRANSPORTMODEL_SOIL", "SnowpackAdvanced", "BUCKET");
			}

			// Set curent meteo variables from 2D fields to single pixel
			const double previous_albedo = getGridPoint(SnGrids::TOP_ALB, ix, iy);
			meteoPixel.rh = rh(ix,iy);
			meteoPixel.ta = ta(ix,iy);
			meteoPixel.vw = vw(ix,iy);
			meteoPixel.iswr = shortwave(ix,iy);
			meteoPixel.rswr = previous_albedo*meteoPixel.iswr;
			meteoPixel.tss = snowPixel.Ndata[snowPixel.getNumberOfElements()].T; //we use previous timestep value
			meteoPixel.ts0 = (snowPixel.SoilNode>0)? snowPixel.Ndata[snowPixel.SoilNode].T : meteoPixel.tss; //we use previous timestep value
			if (soil_temp_depth!=IOUtils::nodata) meteoPixel.ts[0] = getGridPoint(SnGrids::TSOIL, ix,iy); //for consistency: previous time step value
			meteoPixel.ea = Atmosphere::blkBody_Emissivity(longwave(ix,iy), meteoPixel.ta); //to be consistent with Snowpack
			meteoPixel.psum_ph = psum_ph(ix,iy);
			meteoPixel.hs = IOUtils::nodata;
			if (meteoPixel.tss<=100 || meteoPixel.ts0<=100) {
				cout << "[E] Pixel (" << ix+offset << "," << iy << ") too cold! tss=" << meteoPixel.tss << " ts0=" << meteoPixel.ts0 << std::endl;
			}

			// Now determine perpendicular to slope mass balance from drift or precip
			// Note that psum must be in mm/h and SNOWPACK will converted it to mm/CALCULTION_STEP_LENGTH
			if ( useDrift && mns(ix,iy)!=IOUtils::nodata ) {
				double drift_mass = mns(ix,iy);
				// For extreme terrain, the drift might go crazy: limit deposition / erosion
				if (fabs(drift_mass)>37.) {
					cout << "crazy drift at (" << ix+offset << "," << iy << ") = " << drift_mass << " mm/h\n";
					drift_mass = (drift_mass>0.)? 37. : -37.; //set to deposition limit
				}

				meteoPixel.psum = drift_mass;
			} else {
				meteoPixel.psum= psum(ix,iy) * snowPixel.cos_sl; //horiz2slope
			}

			if (useEBalance) meteoPixel.diff = diffuse(ix,iy);
			// Reset Canopy Surface Data to 0 before calling Meteo and Canopy module
			if (useCanopy) snowPixel.Cdata.initializeSurfaceExchangeData();
			//if the current pixel contains soil layers, then activate soil modeling in Snowpack
			sn.setUseSoilLayers( snowPixel.hasSoilLayers() ); //TODO: when Snowpack runs, it should check it

			// exposed glacier special case
			if (isGlacier) {
				//switch to glacier albedo
				snowPixel.Albedo = Constants::glacier_albedo;
				if (meteoPixel.ta>IOUtils::C_TO_K(5.)) {
					//switch to STABLE atmosphere on glacier if TA>5°C
					meteo.setStability(Meteo::MO_HOLTSLAG);
				}
			}

			//compute ustar, psi_s, z0
			meteo.compMeteo(meteoPixel, snowPixel, true);
			SurfaceFluxes surfaceFlux;
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

				try {
					BoundCond Bdata;
					sn.runSnowpackModel(meteoPixel, snowPixel, store(ix,iy), Bdata, surfaceFlux);
					surfaceFlux.collectSurfaceFluxes(Bdata, snowPixel, meteoPixel);
					dIntEnergy += snowPixel.dIntEnergy; //it is reset at every new call to runSnowpackModel
					meteoPixel.hs = snowPixel.cH - snowPixel.Ground; //do not reproject here, otherwise Snowpack outputs would get messed up
				} catch (const std::bad_alloc&) { //don't try anything fancy when running low on memory
					const int lus =(int)floor( landuse(ix,iy));
					const double slope2horiz = (1. / snowPixel.cos_sl);
					const double snow = (snowPixel.cH - snowPixel.Ground) * slope2horiz;
					cout << "[E] Could not allocate memory in Snowpack for cell (" << ix+offset << "," << iy << ") LUS=" << lus << " ";
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
							snowPixel.Cdata.multiplyFluxes(1./nr_snowsteps);
						meteoPixel.psum *= nr_snowsteps;
						snowPixel.dIntEnergy = dIntEnergy;
					}

					ostringstream ss;
					ss << "[E] Snowpack exception: " << e.what() << "\n";
					ss << "[E] at cell (" << ix+offset << "," << iy << ") LUS=" << lus << " ";
					ss << "with " << std::fixed << std::setprecision(4) << snow << " m snow in " << snowPixel.getNumberOfElements()-snowPixel.SoilNode << "/" << snowPixel.getNumberOfElements() << " elements.\n";
					if (isSpecialPoint[index_SnowStation])
						gatherSpecialPoints(meteoPixel, snowPixel, surfaceFlux); //gather special point data, in order to get as much information as possible
					throw IOException(ss.str(), AT);
				}

				// Fill the surfaceFlux.mass variable for output
				surfaceFlux.mass[SurfaceFluxes::MS_HNW] += meteoPixel.psum;
			}
			
			 //some variables are now wrong if we ran multiple Snowpack steps -> recompute them!
			if (nr_snowsteps > 1) {
				surfaceFlux.multiplyFluxes(1./nr_snowsteps);
				if (useCanopy)
					snowPixel.Cdata.multiplyFluxes(1./nr_snowsteps);
				meteoPixel.psum *= nr_snowsteps;
				snowPixel.dIntEnergy = dIntEnergy;
			}

			//switch stability back to normal if it was changed
			if (meteo.getStability()!=USER_STABILITY) meteo.setStability(USER_STABILITY);
			//if the glacier is still exposed, force the albedo back to glacier albedo
			if (isGlacier) {
				snowPixel.Albedo = Constants::glacier_albedo;
				surfaceFlux.pAlbedo = Constants::glacier_albedo;
			}
			if (!std::isfinite( getGridPoint(SnGrids::TOP_ALB, ix, iy) )) { 
				//if the albedo is nan, infinity, etc reset it to its previous
				//value to try to rescue the pixel...
				cerr << "[E] pixel (" << ix+offset << "," << iy << ") found with a nan/infinit albedo ["<<  getGridPoint(SnGrids::TOP_ALB, ix, iy) <<"]; reseting to " << previous_albedo << std::endl;
				getGridPoint(SnGrids::TOP_ALB, ix, iy) = previous_albedo;
			}

			try{
				stability.checkStability(meteoPixel, snowPixel);
			} catch(...) {
				snowPixel.S_4 = Stability::max_stability; //nothing else to do...
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
				sn_cfg.addKey("WATERTRANSPORTMODEL_SOIL", "SnowpackAdvanced", bcu_watertransportmodel_soil);
			}
		}
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
 * @param pointer to vector of SnowStations, where the variable SlopeParFlux contains the lateral flow
 */
void SnowpackInterfaceWorker::getLateralFlow(std::vector<SnowStation*>& ptr_snow_pixel)
{
	for (size_t iy=0;iy<dimy;iy++) {
		for (size_t ix=0;ix<dimx;ix++) {
			const size_t index_SnowStation = dem.getNy()*ix + iy;
			#pragma omp critical (snow_station_lock)
			ptr_snow_pixel.push_back( SnowStations[index_SnowStation] );
		}
	}
}

/**
 * @brief method called by SnowpackInterface to send back the source/sink term for treating lateral flow
 * @param pointer to vector of SnowStations, in which the lwc_source variable is set to contain the lateral flow.
 */
void SnowpackInterfaceWorker::setLateralFlow(const std::vector<SnowStation*>& ptr_snow_pixel)
{
	for (size_t iy=0;iy<dimy;iy++) {
		for (size_t ix=0;ix<dimx;ix++) {
			const size_t index_SnowStation = dem.getNy()*ix + iy;
			if (SnowStations[index_SnowStation]!=NULL) {
				for(size_t n=0; n<SnowStations[index_SnowStation]->getNumberOfElements(); n++) {
					SnowStations[index_SnowStation]->Edata[n].lwc_source = ptr_snow_pixel[index_SnowStation]->Edata[n].lwc_source;
				}
			}
		}
	}
}
