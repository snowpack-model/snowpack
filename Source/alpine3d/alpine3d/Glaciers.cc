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
#include <alpine3d/Glaciers.h>
#include <alpine3d/MPIControl.h>

using namespace std;
using namespace mio;

Glaciers::Glaciers(const mio::Config& cfg) 
              : dem(), isGlacier(), flowpath(), src_altitude(), src_distance(), KBL(IOUtils::nodata), K(0.), scale(1.)
{
	init(cfg);
}

Glaciers::Glaciers(const mio::Config& cfg, const mio::DEMObject& in_dem) 
              : dem(in_dem), isGlacier(), flowpath(), src_altitude(), src_distance(), KBL(IOUtils::nodata), K(0.), scale(1.)
{
	init(cfg);
}

void Glaciers::init(const mio::Config& cfg)
{
	cfg.getValue("KATABATIC_LAYER_HEIGHT", "Snowpack", KBL, IOUtils::nothrow); //Katabatic Boundary Layer (KBL) thickness, 17m in GB
	cfg.getValue("KATABATIC_SCALING", "Snowpack", scale, IOUtils::nothrow); //a value of .25 seems to work quite well
	cfg.getValue("KATABATIC_K_COEFFICIENT", "Snowpack", K, IOUtils::nothrow); //K=0 for original GB model
	
	if (MPIControl::instance().master()) 
		std::cout << "[i] Enabled katabatic flows TA corrections\n";
}

void Glaciers::setDEM(const mio::DEMObject& in_dem)
{
	dem = in_dem;
	if (!isGlacier.empty()) {
		//apply dem mask
		for (size_t ii=0; ii<isGlacier.size(); ii++) {
			if (dem(ii)==IOUtils::nodata) 
				isGlacier(ii) = IOUtils::nodata;
		}
		hillslope_flow(isGlacier); //compute the flow paths
	}
}

/**
 * @brief Set the mask of where the glaciers are.
 * @param[in] glacierMask grid containing IOUtils::nodata at glaciated pixels
 */
void Glaciers::setGlacierMap(const Grid2DObject& glacierMask)
{
	//convert the mask that remove the glaciers to a mask that keeps the glaciers
	isGlacier.set(glacierMask, 0.); //non glacier pixels tagged with 0
	for (size_t idx=0; idx<glacierMask.size(); idx++) { //tag glacier pixels with 1
		if (glacierMask(idx)==IOUtils::nodata)
			isGlacier(idx) = 1.;
	}

	if (!dem.empty()) {
		//apply dem mask
		for (size_t ii=0; ii<isGlacier.size(); ii++) {
			if (dem(ii)==IOUtils::nodata) 
				isGlacier(ii) = IOUtils::nodata;
		}
		hillslope_flow(isGlacier); //compute the flow paths
	}
}

void Glaciers::getGrids(Grid2DObject &alt, Grid2DObject &dist) const
{
	alt = src_altitude;
	dist = src_distance;
}

/**
 * @brief Based on the air and surface temperatures as well as the fraction of the domain that is snow covered, 
 * enable or disable the katabatic flows (the fractions are currently set at 20%).
 * @param[in] hs snow height grid
 * @param[in] tss surface temperature grid
 * @param[in] ta air temperature grid
 * @param[in] isGlacier map of glaciated pixels
 * @return should katabatic flows be computed?
 */
bool Glaciers::enableKatabatikFlows(const mio::Grid2DObject& hs, const mio::Grid2DObject& tss, const mio::Grid2DObject& ta, const mio::Grid2DObject& isGlacier)
{
	if (!ta.isSameGeolocalization(hs))
		throw IOException("The temperature grid and the provided snow height grids don't have the same geolocalization!", AT);
	if (!ta.isSameGeolocalization(tss))
		throw IOException("The temperature grid and the provided dem don't have the same geolocalization!", AT);
	
	size_t glacier_pixels = 0;
	size_t non_glacier_pixels = 0;
	size_t n_snowfree = 0; //count cells that should experience katabatic winds
	size_t n_glacier_katabatic = 0; //count glacier cells that should experience katabatic winds
	
	for (size_t idx=0; idx<isGlacier.size(); idx++) {
		if (isGlacier(idx)==IOUtils::nodata) continue; //not in the simulation
		
		if (isGlacier(idx)==0) { //out of the glacier
			non_glacier_pixels++;
			if (hs(idx)<=0.) n_snowfree++;
		} else { //on the glacier
			glacier_pixels++;
			if (ta(idx)>tss(idx)) n_glacier_katabatic++;
		}
	}
	
	if (n_snowfree>(non_glacier_pixels/5) && n_glacier_katabatic>(glacier_pixels/5))
		return true;
	
	return false;
}

const mio::Grid2DObject Glaciers::correctTemperatures(const mio::Grid2DObject& hs, const mio::Grid2DObject& tss, const mio::Grid2DObject& ta) const
{
	Grid2DObject new_ta(ta);
	correctTemperatures(hs, tss, new_ta);
	return new_ta;
}

//TODO see KLAM_21 for cold air flows

/**
 * @brief  Greuell and Bohm katabatic flow temperature correction.
 * This corrects the air temperature over the glaciated areas according to Greuell, W. and Bohm, R.,
 * <i>"2 m temperatures along melting mid-latitude glaciers, and implications for the sensitivity of the mass
 * balance to variations in temperature"</i>, J. Glaciol., 44(146), 1998.
 * The implementation here follows the overview given by Petersen, L., Pellicciotti, F., Juszak, I., Carenzo, M.
 * and Brock, B. <i>"Suitability of a constant air temperature lapse rate over an Alpine glacier:
 * testing the Greuell and Bohm model as an alternative"</i>, Annals of Glaciology, 54(63), 2013.
 * 
 * The improvements suggested by Ayala, A., Pellicciotti, F., Shea, J. M., <i>"Modeling 2m air temperatures over
 * mountain glaciers: Exploring the influence of katabatic cooling and external warming"</i>, J. Geophys. Res. Atmos., <b>120</b>,
 * 2015, doi:10.1002/2015JDO23137.
 * 
 * The katabatic flows will only be considered active (ie computed and corrected for) if more than 20% of the non-glaciated
 * pixels are snow-free and if the air temperature is higher than the surface temperature at glaciated pixels.
 *
 * @param hs grid containing the snow heights (used to decide if the katabatic flows are active)
 * @param tss grid containing the surface temperatures (used to decide if the katabatic flows are active)
 * @param ta grid containing the air temperatures that will be corrected and returned.
 */
void Glaciers::correctTemperatures(const mio::Grid2DObject& hs, const mio::Grid2DObject& tss, mio::Grid2DObject& ta) const
{
	if (!dem.isSameGeolocalization(ta))
		throw IOException("The temperature grid and the provided dem don't have the same geolocalization!", AT);

	if (src_distance.empty())
		throw IOException("Please provide both dem and glacier map before computing temperature corrections!", AT);
	
	if (enableKatabatikFlows(hs, tss, ta, isGlacier))
		correctTemperatures(ta);
	else
		return;
	
	if (MPIControl::instance().master()) 
		std::cout << "[i] TA field corrected for glacier katabatic flow\n";
}

/**
 * @brief  Greuell and Bohm katabatic flow temperature correction.
 * This method forces the air temperature correction, without checking if it is appropriate or not. This is only
 * useful for calibration.
 * @param ta grid containing the air temperatures that will be corrected and returned.
 */
void Glaciers::correctTemperatures(mio::Grid2DObject& ta) const
{
	static const double CH = 0.002; //bulk heat transfer coefficient
	static const double Gamma = -mio::Cst::dry_adiabatique_lapse_rate;

	std::vector<double> X, Y;
	for (size_t idx=0; idx<isGlacier.size(); idx++) {
		if (isGlacier(idx)==1) {
			X.push_back( dem(idx) );
			Y.push_back( ta(idx) );
		}
	}
	
	const mio::Fit1D temp_fit(mio::Fit1D::NOISY_LINEAR, X, Y);

	for (size_t idx=0; idx<dem.size(); idx++) {
		if (dem(idx)==IOUtils::nodata) continue;
		
		const double slope = dem.slope(idx);
		const double downhill_dist = src_distance(idx)*scale;
		const double z0 = (src_altitude(idx)-dem(idx))*scale + dem(idx); //altitude where the air parcel enters the KBL
		//const double z0 = src_altitude(idx);
		
		const double T0 = temp_fit(z0) - mio::Cst::t_water_freezing_pt; //convert to celsius

		if (downhill_dist==IOUtils::nodata || T0==IOUtils::nodata || slope==IOUtils::nodata)
			continue;

		if (downhill_dist<=0.) continue;
		const double tan_slope = (z0-dem(idx)) / downhill_dist;
		const double cos_slope = Optim::invSqrt(1.+tan_slope*tan_slope); //1+tan^2 = 1/cos^2 and this is 30% faster than computing cos()
		const double Lr =KBL / CH * cos_slope; //KBL is H in the original paper
		const double b = Gamma * tan_slope;
		const double Teq = b * Lr;
		
		if (Lr>0.) {
			const double corr = ((T0 - Teq) * exp(-downhill_dist/Lr) + Teq + K*(downhill_dist/Lr) + mio::Cst::t_water_freezing_pt) - ta(idx); //in the paper: -(x-x0)/Lr
			ta(idx) = ta(idx) + 1.*corr;
		}
	}
}


/**
 * @brief This distributes the flow from one cell to its neighbours
 * @param[in] dem DEM of the domain
 * @param[in] A upsload area of the current cell
 * @param[in] ii first coordinate of the current cell
 * @param[in] jj second coordinate of the current cell
 * @param[out] flow grid that will be updated with the flow from the current cell
 * @param[out] src_altitude altitude of the source of the flow reaching each cell
 * @param[out] src_distance distance to the source of the flow reaching each cell
 * @return true if it could be distributed, false if some surrounding cells were higher
 */
bool Glaciers::hillslope_distribute_cell(const Grid2DObject& dem, const Grid2DObject& glacier_mask, const double& A, const size_t ii, const size_t jj, Grid2DObject &flow, Grid2DObject &src_altitude, Grid2DObject &src_distance)
{
	//set box around the current cell (ii,jj)
	const size_t jjmin = (jj>0)? jj-1 : 0;
	const size_t jjmax = (jj<dem.getNy()-1)? jj+1 : dem.getNy()-1;
	const size_t iimin = (ii>0)? ii-1 : 0;
	const size_t iimax = (ii<dem.getNx()-1)? ii+1 : dem.getNx()-1;

	//check that this is a top cell and get normalization factor
	double sum = 0.;
	for (size_t ll=jjmin; ll<=jjmax; ll++) {
		for (size_t mm=iimin; mm<=iimax; mm++) {
			if (glacier_mask(mm,ll)==1 && dem(mm,ll) > dem(ii,jj)) //some cells are higher and must be processed before
				return false;

			if (dem(ii,jj) > dem(mm,ll) && dem(mm,ll)!=IOUtils::nodata) { //(ii,jj) cell is naturally excluded
				sum += (dem(ii,jj) - dem(mm,ll)); //tan(beta) = delta_alt / L so tan(beta)*L = delta_alt
			}
		}
	}

	//initializa cell altitude and distance if necessary
	if (src_altitude(ii,jj) == IOUtils::nodata) {
		//ie cell is a top cell
		src_altitude(ii,jj) = dem(ii,jj);
		src_distance(ii,jj) = 0.;
	}

	//flow into the lower cells
	if (sum==0.) return true; //current cell is flat terrain
	const double C = A / sum;
	for (size_t ll=jjmin; ll<=jjmax; ll++) {
		for (size_t mm=iimin; mm<=iimax; mm++) {
			if (glacier_mask(mm,ll)!=1) continue; //not a glacier cell

			if (dem(ii,jj) > dem(mm,ll)) {
				const double flow_contrib = C * (dem(ii,jj) - dem(mm,ll));
				const double weight = flow_contrib/(flow(mm,ll)+flow_contrib);

				flow(mm,ll) += flow_contrib;

				//set altitude of the flow source
				if (src_altitude(mm,ll)==IOUtils::nodata)
					src_altitude(mm,ll) = src_altitude(ii,jj);
				else
					src_altitude(mm,ll) = weight * src_altitude(ii,jj) + (1.-weight)*src_altitude(mm,ll);

				//set distance of the flow source
				if (src_distance(mm,ll)==IOUtils::nodata)
					src_distance(mm,ll) = src_distance(ii,jj)+1.;
				else
					src_distance(mm,ll) = weight * (src_distance(ii,jj)+1.) + (1.-weight)*src_distance(mm,ll);
			}
		}
	}

	return true;
}

/**
 * @brief This is an implementation of the multiple flow direction algorithm.
 * This associates to each cell an upslope area that is based on the number of cells
 * that flow into it. This follows <i>"The prediction of hillslope flow paths for distributed
 * hydrological modelling using digital terrain models"</i>, Quinn P., Chevallier P., Planchon O.,
 * hydrological processes, <b>5</b>, 1991, pp 59-79.
 *
 * In this implementation, the cells that are higher than their neighbours are computed first,
 * and their elevation set to nodata. This is performed on the whole grid until all cells have
 * nodata as their elevation.
 * @param[in] glacier_mask pixels that are not glaciated are marked as IOUtils::nodata
 */
void Glaciers::hillslope_flow(Grid2DObject glacier_mask)
{
	if (!dem.isSameGeolocalization(glacier_mask))
		throw IOException("The glacier mask and the provided dem don't have the same geolocalization!", AT);
	
	flowpath.set(dem, 1.);
	src_altitude.set(dem, IOUtils::nodata);
	src_distance.set(dem, IOUtils::nodata);

	unsigned int nr_cells_left;
	const size_t ncols = dem.getNx(), nrows = dem.getNy();
	do {
		nr_cells_left = 0;
		for (size_t jj=0; jj<nrows; jj++) {
			for (size_t ii=0; ii<ncols; ii++) {
				if (glacier_mask(ii,jj)!=1) continue; //this also skips out of domain cells

				const double A = flowpath(ii,jj);
				const bool distributed = hillslope_distribute_cell(dem, glacier_mask, A, ii, jj, flowpath, src_altitude, src_distance);
				if (distributed)
					glacier_mask(ii,jj) = IOUtils::nodata; //mark the cell as done
				else
					nr_cells_left++;
			}
		}
	} while (nr_cells_left>0);

	src_distance *= dem.cellsize;
}
