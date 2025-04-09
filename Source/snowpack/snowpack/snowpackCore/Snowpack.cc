/*
 *  SNOWPACK stand-alone
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
*/
/*  This file is part of Snowpack.
    Snowpack is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snowpack is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snowpack.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * @file Snowpack.cc
 * @version 11.06
 * @bug     -
 * @brief This module contains the driving routines for the 1d snowpack model
 */

#include <snowpack/snowpackCore/Snowpack.h>
#include <snowpack/snowpackCore/Solver.h>
#include <snowpack/Meteo.h>
#include <snowpack/Constants.h>
#include <snowpack/Utils.h>
#include <snowpack/Laws_sn.h>
#include <snowpack/snowpackCore/WaterTransport.h>
#include <snowpack/snowpackCore/VapourTransport.h>
#include <snowpack/TechnicalSnow.h>
#include <snowpack/snowpackCore/Metamorphism.h>
#include <snowpack/snowpackCore/PhaseChange.h>

#include <assert.h>
#include <sstream>
#include <errno.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/
//Minimum allowed snowpack time step for solving the heat equation (automatic time stepping is applied when equation does not converge)
const double Snowpack::min_allowed_sn_dt = 0.01;

//Uses an empirically determined size of deposited hydrometeors as new snow grain size (mm)
const bool Snowpack::hydrometeor = false;

//Warning is issued if depth of snowfall is larger than this amount (m)
const double Snowpack::snowfall_warning = 0.5;

const unsigned int Snowpack::new_snow_marker = 0;

/// Min volumetric ice content allowed
const double Snowpack::min_ice_content = SnLaws::min_hn_density / Constants::density_ice;

/// @brief Define the assembly macro
void Snowpack::EL_INCID(const int &e, int Ie[]) {
	Ie[0] = e;
	Ie[1] = e+1;
}

/// @brief Define the node to element temperature macro
void Snowpack::EL_TEMP( const int Ie[], double Te0[], double Tei[], const std::vector<NodeData> &T0, const double Ti[] ) {
	Te0[ 0 ] = T0[ Ie[ 0 ] ].T;
	Te0[ 1 ] = T0[ Ie[ 1 ] ].T;
	Tei[ 0 ] = Ti[ Ie[ 0 ] ];
	Tei[ 1 ] = Ti[ Ie[ 1 ] ];
}

/// @brief Element right-hand side macro
void Snowpack::EL_RGT_ASSEM(double F[], const int Ie[], const double Fe[]) {
	F[Ie[0]] += Fe[0];
	F[Ie[1]] += Fe[1];
}

/************************************************************
 * non-static section                                       *
 ************************************************************/

Snowpack::Snowpack(const SnowpackConfig& i_cfg)
          : surfaceCode(), cfg(i_cfg), techsnow(i_cfg),
            variant(), forcing(), viscosity_model(), watertransportmodel_snow("BUCKET"), watertransportmodel_soil("BUCKET"),
            hn_density(), hn_density_parameterization(), sw_mode(), snow_albedo(), albedo_parameterization(), albedo_average_schmucki(), sw_absorption_scheme(),
            atm_stability_model(), albedo_NIED_av(0.75), albedo_fixedValue(Constants::glacier_albedo), hn_density_fixedValue(SnLaws::min_hn_density),
            meteo_step_length(0.), thresh_change_bc(-1.0), geo_heat(Constants::undefined), height_of_meteo_values(0.),
            height_new_elem(0.), sn_dt(0.), t_crazy_min(0.), t_crazy_max(0.), thresh_rh(0.), thresh_dtempAirSnow(0.),
            new_snow_dd(0.), new_snow_sp(0.), new_snow_dd_wind(0.), new_snow_sp_wind(0.), rh_lowlim(0.), bond_factor_rh(0.),
            new_snow_grain_size(0.), new_snow_bond_size(0.), hoar_density_buried(0.), hoar_density_surf(0.), hoar_min_size_buried(0.),
            minimum_l_element(0.), comb_thresh_l(IOUtils::nodata), t_surf(0.),
            allow_adaptive_timestepping(false), research_mode(false), useCanopyModel(false), enforce_measured_snow_heights(false), detect_grass(false),
            soil_flux(false), useSoilLayers(false), coupled_phase_changes(false), combine_elements(false), reduce_n_elements(0), force_add_snowfall(false), max_simulated_hs(-1.),
            change_bc(false), meas_tss(false), vw_dendricity(false),
            enhanced_wind_slab(false), snow_erosion("NONE"), redeposit_keep_age(false), alpine3d(false), ageAlbedo(true), soot_ppmv(0.), adjust_height_of_meteo_values(true),
            adjust_height_of_wind_value(false), advective_heat(false), heat_begin(0.), heat_end(0.),
            temp_index_degree_day(0.), temp_index_swr_factor(0.), forestfloor_alb(false), rime_index(false), newsnow_lwc(false), read_dsm(false), soil_evaporation(), soil_thermal_conductivity()
{
	cfg.getValue("FORCING", "Snowpack", forcing);
	cfg.getValue("ALPINE3D", "SnowpackAdvanced", alpine3d);
	cfg.getValue("VARIANT", "SnowpackAdvanced", variant);
	cfg.getValue("COUPLEDPHASECHANGES", "SnowpackAdvanced", coupled_phase_changes);
	if (variant=="SEAICE") coupled_phase_changes = true;	// to better deal with variable freezing point due to salinity

	//Define keys for new snow density computation
	cfg.getValue("HN_DENSITY", "SnowpackAdvanced", hn_density);
	cfg.getValue("TEMP_INDEX_DEGREE_DAY", "SnowpackAdvanced", temp_index_degree_day, IOUtils::nothrow);
	cfg.getValue("TEMP_INDEX_SWR_FACTOR", "SnowpackAdvanced", temp_index_swr_factor, IOUtils::nothrow);
	cfg.getValue("HN_DENSITY_PARAMETERIZATION", "SnowpackAdvanced", hn_density_parameterization);
	cfg.getValue("HN_DENSITY_FIXEDVALUE", "SnowpackAdvanced", hn_density_fixedValue);

	//Define keys for new snow information
	cfg.getValue("RIME_INDEX", "SnowpackAdvanced", rime_index);
	cfg.getValue("NEWSNOW_LWC", "SnowpackAdvanced", newsnow_lwc);
	cfg.getValue("READ_DSM", "SnowpackAdvanced", read_dsm);

	//Define keys for snow albedo computation
	cfg.getValue("SNOW_ALBEDO", "SnowpackAdvanced", snow_albedo);
	cfg.getValue("ALBEDO_PARAMETERIZATION", "SnowpackAdvanced", albedo_parameterization);
	cfg.getValue("ALBEDO_AVERAGE_SCHMUCKI", "SnowpackAdvanced", albedo_average_schmucki);
	if (albedo_parameterization=="NIED")
		cfg.getValue("ALBEDO_NIED_AV", "SnowpackAdvanced", albedo_NIED_av);
	else
		albedo_NIED_av=Constants::undefined;
	cfg.getValue("ALBEDO_FIXEDVALUE", "SnowpackAdvanced", albedo_fixedValue);
	cfg.getValue("ALBEDO_AGING", "SnowpackAdvanced", ageAlbedo);

	//Defines whether a multiband model is used for short wave radiation absorption
	cfg.getValue("SW_ABSORPTION_SCHEME", "SnowpackAdvanced", sw_absorption_scheme);

	// Defines whether soil layers are used
	cfg.getValue("SNP_SOIL", "Snowpack", useSoilLayers);
	/* Defines the management of the bottom boundary conditions with soil layers
	 * - 0 ==> Dirichlet, i.e fixed Temperature
	 * - 1 ==> Neumann, fixed geothermal heat flux GEO_HEAT */
	cfg.getValue("SOIL_FLUX", "Snowpack", soil_flux, IOUtils::nothrow);
	if (soil_flux || variant == "SEAICE") {
		// For sea ice, geo_heat is ocean heat flux
		cfg.getValue("GEO_HEAT", "Snowpack", geo_heat); //Constant geothermal heat flux at (great) depth (W m-2)
	} else {
		geo_heat = Constants::undefined;
	}

	/* Defines the management of the surface boundary conditions
	 * - 0: Neumann boundary conditions throughout
	 * - 1: Dirichlet if Tss < THRESH_CHANGE_BC, Neumann else */
	cfg.getValue("CHANGE_BC", "Snowpack", change_bc);
	if (change_bc)
		cfg.getValue("THRESH_CHANGE_BC", "Snowpack", thresh_change_bc);

	//Should be NODATA for data-sets which do not provide measured surface temperatures
	cfg.getValue("MEAS_TSS", "Snowpack", meas_tss);

	/*
	 * Defines how the height of snow is going to be handled
	 * - 0: Depth of snowfall is determined from the water equivalent of snowfall (PSUM)
	 * - 1: The measured height of snow is used to determine whether new snow has been deposited.
	 *      This setting MUST be chosen in operational mode. \n
	 *      This procedure has the disadvantage that if the snowpack settles too strongly
	 *      extra mass is added to the snowpack. \n
	 * New snow density is needed in both cases, either parameterized, measured, or fixed.
	 * Also check whether growing grass should be detected
	 */
	cfg.getValue("ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack", enforce_measured_snow_heights);
	cfg.getValue("DETECT_GRASS", "SnowpackAdvanced", detect_grass);

	/* Defines whether the canopy model is used
	 * NOTE: OUT_CANOPY must also be set to dump canopy parameters to file; see Constants_local.h
	 */
	cfg.getValue("CANOPY", "Snowpack", useCanopyModel);

	/* Define the heights of the meteo measurements above ground (m)
	 * Required for surface energy exchange computation and for drifting and blowing snow.
	 */
	cfg.getValue("HEIGHT_OF_METEO_VALUES", "Snowpack", height_of_meteo_values);

	/* Defines what shortwave radiation flux(es) to use
	 * - "INCOMING" (downwelling) SW radiation is used
	 * - "REFLECTED" SW radiation is used
	 * - "BOTH" downward and reflected SW radiation is used
	 * @note { If SW_MODE == "BOTH", the input must hold both fluxes! } */
	cfg.getValue("SW_MODE", "Snowpack", sw_mode);

	// Defines used atmospheric stability, used for determining if dynamic time steps may be required
	const std::string atm_stability_string = cfg.get("ATMOSPHERIC_STABILITY", "Snowpack");
	atm_stability_model = Meteo::getStability(atm_stability_string);

	// Allow dynamic time stepping in case of unstable atmospheric stratification
	cfg.getValue("ALLOW_ADAPTIVE_TIMESTEPPING", "SnowpackAdvanced", allow_adaptive_timestepping);

	/* Height of new snow element (m) [NOT read from CONSTANTS_User.INI] \n
	 * Controls the addition of new snow layers. Set in qr_ReadParameters() \n
	 * The value depends on ENFORCE_MEASURED_SNOW_HEIGHTS:
	 * - 0: 2.0*MINIMUM_L_ELEMENT (value depends on VARIANT)
	 * - 1: 0.02 */
	cfg.getValue("HEIGHT_NEW_ELEM", "SnowpackAdvanced", height_new_elem);
	cfg.getValue("MINIMUM_L_ELEMENT", "SnowpackAdvanced", minimum_l_element);
	if(minimum_l_element<=0.) throw IOException("MINIMUM_L_ELEMENT must be >0! Please fix your ini file.", AT);
	cfg.getValue("FORCE_ADD_SNOWFALL", "SnowpackAdvanced", force_add_snowfall);

	cfg.getValue("RESEARCH", "SnowpackAdvanced", research_mode);

	cfg.getValue("VISCOSITY_MODEL", "SnowpackAdvanced", viscosity_model);

	/* Precipitation only for humidity above and temperature difference within threshold (1)
	 * - thresh rh (default): 0.50
	 * 	- 2007-12-01: set THRESH_RH to 0.70 to be consistent with data range of ZWART new snow density model
	 * 	- 2008-01-21: set back THRESH_RH to 0.50 (IMIS sensor problem in operational mode)
	 * 	- Antarctica: 0.70
	 * - thresh dtempAirSnow: 3.0 */
	cfg.getValue("THRESH_RH", "SnowpackAdvanced", thresh_rh);
	cfg.getValue("THRESH_DTEMP_AIR_SNOW", "SnowpackAdvanced", thresh_dtempAirSnow);

	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	const double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	sn_dt = M_TO_S(calculation_step_length);
	meteo_step_length = cfg.get("METEO_STEP_LENGTH", "Snowpack");

	//Defines whether joining elements will be considered at all
	cfg.getValue("COMBINE_ELEMENTS", "SnowpackAdvanced", combine_elements);
	//Activates algorithm to reduce the number of elements deeper in the snowpack AND to split elements again when they come back to the surface
	//Only works when COMBINE_ELEMENTS == TRUE.
	cfg.getValue("REDUCE_N_ELEMENTS", "SnowpackAdvanced", reduce_n_elements);
	cfg.getValue("COMB_THRESH_L", "SnowpackAdvanced", comb_thresh_l, IOUtils::nothrow);
	if(comb_thresh_l == IOUtils::nodata) comb_thresh_l = SnowStation::comb_thresh_l_ratio * height_new_elem;	// If no comb_thresh_l specified, use the default one (i.e., a fixed ratio from height_new_elem)
	cfg.getValue("MAX_SIMULATED_HS", "SnowpackAdvanced", max_simulated_hs);
	if(max_simulated_hs > 0. && useSoilLayers) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Inconsistent setting: you cannot use MAX_SIMULATED_HS > 0 for simulations with soil.");
		throw IOException("Runtime Error in Snowpack::Snowpack()", AT);
	}

	//Warning is issued if snow tempeartures are out of bonds, that is, crazy
	cfg.getValue("T_CRAZY_MIN", "SnowpackAdvanced", t_crazy_min);
	cfg.getValue("T_CRAZY_MAX", "SnowpackAdvanced", t_crazy_max);
	cfg.getValue("FORESTFLOOR_ALB", "SnowpackAdvanced", forestfloor_alb);
	cfg.getValue("ENHANCED_WIND_SLAB", "SnowpackAdvanced", enhanced_wind_slab);

	/* Initial new snow parameters, see computeSnowFall()
	* - that rg and rb are equal to 0.5*gsz and 0.5*bsz, respectively. Both given in millimetres
	* - If VW_DENDRICITY is set, new snow dendricity is f(vw)
	* - BOND_FACTOR_RH new snow bonds get stronger for average winds >= SnLaws::event_wind_lowlim and
	*   mean relative humidity >= rh_lowlim */
	if (variant == "ANTARCTICA" || variant == "POLAR") {
		if (variant == "ANTARCTICA") {
			new_snow_dd = 0.5;
			new_snow_sp = 0.75;
			new_snow_dd_wind = 0.15;
			new_snow_sp_wind = 1.0;
		}
		if (variant == "POLAR") {
			// average between ANTARCTICA and DEFAULT
			new_snow_dd = 0.75;
			new_snow_sp = 0.625;
			new_snow_dd_wind = 0.325;
			new_snow_sp_wind = 0.875;
		}
		vw_dendricity = false;
		rh_lowlim = 0.7;
		bond_factor_rh = 3.0;
		enhanced_wind_slab = true;
		ageAlbedo = false;
	} else {
		new_snow_dd = 1.0;
		new_snow_sp = 0.5;
		new_snow_dd_wind = 0.5;
		new_snow_sp_wind = 0.75;
		vw_dendricity = true;
		rh_lowlim = 1.0;
		bond_factor_rh = 1.0; 
		// enhanced_wind_slab = true; //true; // set in namelist now
	}

	cfg.getValue("SNOW_EROSION", "SnowpackAdvanced", snow_erosion);
	std::transform(snow_erosion.begin(), snow_erosion.end(), snow_erosion.begin(), ::toupper);	// Force upper case
	cfg.getValue("REDEPOSIT_KEEP_AGE", "SnowpackAdvanced", redeposit_keep_age);

	cfg.getValue("NEW_SNOW_GRAIN_SIZE", "SnowpackAdvanced", new_snow_grain_size);
	new_snow_bond_size = 0.25 * new_snow_grain_size;

	/* Thresholds for surface hoar formation and burial
	 * NOTE that the value of the parameter ROUGHNESS_LENGTH in CONSTANTS_User.INI is critical for surface hoar formation,
	 * particularly for Dirichlet boundary conditions. Value should be < 1 mm. Other considerations favor larger values.
	 * - 0.0007 m : original calibration with the 98/99 data set
	 * - 0.002  m : favored operational value with Dirichlet bc */
	//Density of BURIED surface hoar (kg m-3), default: 125./ Antarctica: 200.
	cfg.getValue("HOAR_DENSITY_BURIED", "SnowpackAdvanced", hoar_density_buried);
	if ( ! (hoar_density_buried > Constants::eps) )
		throw IOException("HOAR_DENSITY_BURIED must be >0! Please fix your ini file.", AT);

	//Density of surface hoar (-> hoar index of surface node) (kg m-3)
	cfg.getValue("HOAR_DENSITY_SURF", "SnowpackAdvanced", hoar_density_surf);
	if ( ! (hoar_density_surf > Constants::eps) )
		throw IOException("HOAR_DENSITY_SURF must be >0! Please fix your ini file.", AT);

	//Minimum surface hoar size to be buried (mm). Increased by 50% for Dirichlet bc.
	cfg.getValue("HOAR_MIN_SIZE_BURIED", "SnowpackAdvanced", hoar_min_size_buried);

	//Watertransport models
	cfg.getValue("WATERTRANSPORTMODEL_SNOW", "SnowpackAdvanced", watertransportmodel_snow);
	cfg.getValue("WATERTRANSPORTMODEL_SOIL", "SnowpackAdvanced", watertransportmodel_soil);
	if (variant!="SEAICE" && coupled_phase_changes && useSoilLayers && watertransportmodel_soil=="RICHARDSEQUATION")
		throw IOException("COUPLEDPHASECHANGES cannot be set to true when WATERTRANSPORTMODEL_SOIL == RICHARDSEQUATION! Please correct your ini file.", AT);

	//Indicate if the meteo values can be considered at constant height above the snow surface (e.g., Col de Porte measurement method)
	cfg.getValue("ADJUST_HEIGHT_OF_METEO_VALUES", "SnowpackAdvanced", adjust_height_of_meteo_values);
	cfg.getValue("ADJUST_HEIGHT_OF_WIND_VALUE", "SnowpackAdvanced", adjust_height_of_wind_value);

	// Allow for the effect of a known advective heat flux
	cfg.getValue("ADVECTIVE_HEAT", "SnowpackAdvanced", advective_heat, IOUtils::nothrow);
	cfg.getValue("HEAT_BEGIN", "SnowpackAdvanced", heat_begin, IOUtils::nothrow);
	cfg.getValue("HEAT_END", "SnowpackAdvanced", heat_end, IOUtils::nothrow);

	/* Get the soil evaporation model to be used
	*  - EVAP_RESISTANCE: Resistance Approach, see Laws_sn.c␊
	*  - RELATIVE_HUMIDITY: Relative Humidity Approach, see Snowpack.cc
	*  - NONE: none, assume saturation pressure and no extra resistance */
	cfg.getValue("SOIL_EVAP_MODEL", "SnowpackAdvanced", soil_evaporation);
	/* Get the soil thermal conductivity model to be used
	*  - FITTED: Use fit values for soil thermal conductivity, see snLaws::compSoilThermalConductivity()
	*  - COSENZA2003: Use the soil thermal conductivity model by Cosenza et al. (2003)
	*  - RAW: Use simply Edata.soil[SOIL_K] + Edata.theta[WATER] * SnLaws::conductivity_water(Edata.Te)
	                    + Edata.theta[ICE] * SnLaws::conductivity_ice(Edata.Te) */
	cfg.getValue("SOIL_THERMAL_CONDUCTIVITY", "SnowpackAdvanced", soil_thermal_conductivity);

	// Soot/impurity in ppmv for albedo caclulations
	cfg.getValue("SOOT_PPMV", "SnowpackAdvanced", soot_ppmv);
}

void Snowpack::setUseSoilLayers(const bool& value) { //NOTE is this really needed?
	useSoilLayers = value;
}

/**
 * @brief Snow creep
 * -# The Thing ain't settling any more in case of ice, soil or water only
 * -# Enhanced densification for wind slabs of Metamorphism::wind_slab_depth (m); see also mm_Metamorphism()
 *    dry snow AND strong wind AND near the surface => enhance densification \n
 * -# Normal densification
 * -# Empirism for surface hoar. Two different rates are used for either large or small SH.
 *    Implemented by Sascha Bellaire on 28.11.2006.
 * @todo name parameters for the computation of CDot
 * @param Xdata
 * @param Mdata
 */
void Snowpack::compSnowCreep(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata)
{
	const bool prn_WRN = false;
	const size_t nN = Xdata.getNumberOfNodes();
	if (nN == (Xdata.SoilNode + 1))
		return;

	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;
	const size_t nE = Xdata.getNumberOfElements();
	double SigC = 0.; // Cauchy stress
	const double SigC_fac = Constants::g * Xdata.cos_sl;
	for(size_t e = nE; e --> 0; ) {
		const double oldStress = EMS[e].C;
		const double age = std::max(0., Mdata.date.getJulian() - EMS[e].depositionDate.getJulian());
		if (e < nE-1)
			SigC -= (EMS[e+1].M / 2.) * SigC_fac;
		SigC -= (EMS[e].M / 2.) * SigC_fac;
		EMS[e].C = SigC;
		assert(EMS[e].C<0.);
		if (EMS[e].CDot / SigC > 0.05) {
			EMS[e].CDot *= exp(-0.037 * S_TO_D(sn_dt));
		} else {
			EMS[e].CDot = 0.;
		}
		if ((e < nE-1) && (age > Constants::eps)) {
			if ((SigC - oldStress) < 0.)
				EMS[e].CDot += (SigC - oldStress);
		}
	}

	for (size_t e = Xdata.SoilNode; e < nE; e++) {
		double eta = SnLaws::smallest_viscosity; // snow viscosity
		if (EMS[e].Rho > 910. ||  EMS[e].theta[SOIL] > 0. || EMS[e].theta[ICE] < Constants::eps) {
			EMS[e].k[SETTLEMENT] = eta = 1.0e99;
		} else {
			EMS[e].k[SETTLEMENT] = eta = SnLaws::compSnowViscosity(variant, viscosity_model, watertransportmodel_snow, EMS[e], Mdata.date);
			if (!(eta > 0.01 * SnLaws::smallest_viscosity && eta <= 1.e11 * SnLaws::smallest_viscosity)
			        && (EMS[e].theta[ICE] > 2. * Snowpack::min_ice_content) && (EMS[e].theta[ICE] < 0.6)) {
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date,
				          "Viscosity=%e out of range! e=%d nE=%d rg=%lf rb=%lf dd=%lf sp=%lf theta_i=%lf theta_w=%lf",
				            eta, e, nE, EMS[e].rg, EMS[e].rb, EMS[e].dd, EMS[e].sp,
				              EMS[e].theta[ICE], EMS[e].theta[WATER]);
			}
			if (eta < SnLaws::smallest_viscosity) {
				if (prn_WRN) //HACK
					prn_msg(__FILE__, __LINE__, "wrn", Mdata.date,
					        "Viscosity=%e reset to SMALLEST_VISCOSITY! e=%d nE=%d", eta, e, nE);
				EMS[e].k[SETTLEMENT] = eta = SnLaws::smallest_viscosity;
			}
		}
		const double Sig0 = SnLaws::compLoadingRateStress(viscosity_model, EMS[e], Mdata.date); // "Sintering" stress
		const double L0 = EMS[e].L;
		double dL;

		if (EMS[e].mk%100 != 3) { //ALL except SH
			double wind_slab=1.;
			const double dz = NDS[nE].z - NDS[e].z;
			const double z_ref_vw = 3.;	// See p. 336 in Groot Zwaaftink et al. (doi: https://doi.org/10.5194/tc-7-333-2013)
			const double vw_ref = Meteo::windspeedProfile(Mdata, z_ref_vw);
			const double dv = vw_ref - Metamorphism::wind_slab_vw;
			// if (snow_erosion == "REDEPOSIT") {  //Why should this not be used with REDEPOSIT? 
			// 	wind_slab = 1.;
			// } else {
			{	if ((EMS[e].theta[WATER] < SnowStation::thresh_moist_snow)
				      && (vw_ref > Metamorphism::wind_slab_vw)
					&& ((dz < Metamorphism::wind_slab_depth) || (e == nE-1))) {
					if (Snowpack::enhanced_wind_slab) { //NOTE tested with Antarctic variant: effects heavily low density snow
						// fits original parameterization at Metamorphism::wind_slab_vw + 0.6 m/s
						wind_slab += 2.7 * Metamorphism::wind_slab_enhance
							         * Optim::pow3(dv) * (1. - dz / (1.25 * Metamorphism::wind_slab_depth));
					} else {
						// original parameterization by Lehning
						wind_slab += Metamorphism::wind_slab_enhance * dv;
					}
				}
			}
			EMS[e].Eps_vDot = wind_slab * (EMS[e].C + Sig0) / eta;
			dL = L0 * sn_dt * EMS[e].Eps_vDot;

			// Make sure settling is not larger than the space that is available (basically settling can at most reduce theta[AIR] to 0).
			// We also leave some room in case all liquid water freezes and thereby expands.
			const double MaxSettlingFactor = (watertransportmodel_snow=="RICHARDSEQUATION") ? (0.9) : (1. - Constants::eps); // An additional maximum settling factor, between 0 and 1. 1: allow maximize possible settling, 0: no settling allowed.
			dL = std::max(dL, std::min(0., -1.*MaxSettlingFactor*L0*(EMS[e].theta[AIR]-((Constants::density_water/Constants::density_ice)-1.)*(EMS[e].theta[WATER]+EMS[e].theta[WATER_PREF]))));

			// Limit dL when the element length drops below minimum_l_element. This element will be merged in WaterTransport::mergingElements later on.
			if ((L0 + dL) < (1.-Constants::eps)*minimum_l_element)
				dL = std::min(0., (1.-Constants::eps)*minimum_l_element - L0);	// Make sure the element length gets smaller than minimum_l_element.
		} else { //SH
			if (NDS[e+1].hoar > 0.006) { // TODO Large initial size, i.e., deposited hoar mass/HOAR_DENSITY_BURIED ??
				if ((Mdata.date.getJulian() - EMS[e].depositionDate.getJulian()) < 21.)
					dL = MM_TO_M(-0.391 * S_TO_D(sn_dt));
				else
					dL = MM_TO_M(-0.0807 * S_TO_D(sn_dt));
			} else {
				dL = MM_TO_M(-0.1107 * S_TO_D(sn_dt));
			}
			EMS[e].Eps_vDot = dL / (L0 * sn_dt);
			if ((L0 + dL) < 0. )
				dL = 0.;
		}

		if (variant == "CALIBRATION") {
			NDS[e].f = (Sig0 - EMS[e].Eps_Dot) / eta; // sigMetamo
			EMS[e].Eps_Dot /= eta;                    // sigReac
			NDS[e].udot = EMS[e].C / eta;          // Deformation due to load alone
			EMS[e].S = EMS[e].CDot / EMS[e].C;     // ratio loadrate to load addLoad to load
		}

		EMS[e].theta[WATER] *= L0 / (L0 + dL);
		EMS[e].theta[WATER_PREF] *= L0 / (L0 + dL);
		EMS[e].theta[ICE]   *= L0 / (L0 + dL);
		EMS[e].theta_i_reservoir   *= L0 / (L0 + dL);
		EMS[e].theta_i_reservoir_cumul   *= L0 / (L0 + dL);
		EMS[e].L0 = EMS[e].L = (L0 + dL);
		NDS[e+1].z = NDS[e].z + EMS[e].L;
		EMS[e].theta[AIR] = 1.0 - EMS[e].theta[WATER] - EMS[e].theta[WATER_PREF] - EMS[e].theta[ICE] - EMS[e].theta[SOIL];
		EMS[e].updDensity();
		Sdata.mass[SurfaceFluxes::MS_SETTLING_DHS] += dL;	// Update snow height change due to settling
		if (EMS[e].Rho <= Constants::eps || EMS[e].theta[AIR] < 0.  ) {
			prn_msg(__FILE__, __LINE__, "err", Date(),
			          "Volume contents: e=%d nE=%d rho=%lf ice=%lf wat=%lf wat_pref=%lf air=%le",
			            e, nE, EMS[e].Rho, EMS[e].theta[ICE], EMS[e].theta[WATER], EMS[e].theta[WATER_PREF], EMS[e].theta[AIR]);
			throw IOException("Runtime Error in compSnowCreep()", AT);
		}
		EMS[e].gradT = (NDS[e+1].T - NDS[e].T) / EMS[e].L;
	}
	// Update computed snow depth
	Xdata.cH = NDS[nN-1].z + NDS[nN-1].u;
}

/**
 * @brief Computes the element stiffness matrix and right hand side vector for a fully implicit time integration scheme \n
 * The matrices that must be evaluated are : \n
 * - [Se] = [Ce]/dt + [Ke] :  [Ce] is the element heat capacity matrix and [Ke] is the
 *                            element conductivity matrix.
 * - {Fe} = {Q} - [Ke]*{T0} : {Fe} is the element right-hand side vector containing the element heat flux. \n
 *                            {Q} arises from shortwave radiation -- the other heat fluxes are treated separately
 *                            when determining the meteo parameters.
 * @param Edata
 * @param dt Calculation time step length (s)
 * @param dvdz Wind pumping velocity gradient (s-1)
 * @param T0 Initial nodal temperatures (K)
 * @param Se Element heat capacitity (stiffness) matrix
 * @param Fe Element right hand side vector
 * @param VaporEnhance Vapor transport enhancement factor
 * @return false on error, true if no error occurred
 */
bool Snowpack::sn_ElementKtMatrix(ElementData &Edata, double dt, const double dvdz, double T0[ N_OF_INCIDENCES ], double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ], double Fe[ N_OF_INCIDENCES ], const double VaporEnhance)
{
	if (Edata.L < 0.0) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Negative length L=%e", Edata.L);
		return false;
	}

	//reset Fe_0 and Fe_1
	Fe[0] = Fe[1] = 0.0;

	// Find the conductivity of the element TODO: check thresholds
	double Keff;    // the effective thermal conductivity
	if (Edata.theta[SOIL] > 0.0) {
		Keff = SnLaws::compSoilThermalConductivity(Edata, dvdz, soil_thermal_conductivity);
	} else if (Edata.theta[ICE] > 0.55 || Edata.theta[ICE] < min_ice_content) {
		// Note: no soil when inside this if-block.
		Keff = Edata.theta[AIR] * Constants::conductivity_air + Edata.theta[ICE] * Constants::conductivity_ice +
		           (Edata.theta[WATER]+Edata.theta[WATER_PREF]) * Constants::conductivity_water;
	} else {
		Keff = SnLaws::compSnowThermalConductivity(Edata, dvdz, !alpine3d); //do not show the warning for Alpine3D
	}

	// mimics effect of vapour transport if liquid water present in snowpack
	Keff *= VaporEnhance;
	Edata.k[TEMPERATURE] = Keff;
	const double k = Keff/Edata.L;   //Conductivity. Divide by the length to save from doing it during the matrix operations

	// Evaluate the stiffness matrix
	Se[0][0] = Se[1][1] = k;
	Se[0][1] = Se[1][0] = -k;
	// Heat the element via short-wave radiation
	if (Edata.sw_abs < 0.0 && !advective_heat) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "NEGATIVE Shortwave Radiation %e", Edata.sw_abs);
		return false;
	}
	Fe[1] += Edata.sw_abs;

	// Add the implicit time integration term to the right hand side
	Fe[0] -= (Se[0][0] * T0[0] + Se[0][1] * T0[1]);
	Fe[1] -= (Se[1][0] * T0[0] + Se[1][1] * T0[1]);

	// Now add the heat capacitity matrix
	Edata.heatCapacity();
	const double c = Edata.c[TEMPERATURE] * Edata.L * Edata.Rho / (6. * dt); //heat capacity
	Se[0][0] += 2. * c;
	Se[1][1] += 2. * c;
	Se[0][1] += c;
	Se[1][0] += c;

	// Add the source/sink term resulting from phase changes
	Fe[1] += Edata.Qph_up * 0.5 * Edata.L;
	Fe[0] += Edata.Qph_down * 0.5 * Edata.L;

	// Add the source/sink term resulting from phase changes (due to water vapor transport)
	Fe[1] += Edata.Qmm * 1.0 * Edata.L;
	Fe[0] += Edata.Qmm * 1.0 * Edata.L;

	//Se[1][1] += Edata.Qmm*Edata.L/Edata.Te;
	//Fe[1] += Edata.Qmm*Edata.L;

	return true;
}

/**
* @brief Update all boundary energy fluxes but solar irradiance \n
* The fluxes are first updated before entering compTemperatureProfile
* and used to impose Neumann boundary conditions when required.
* That way we ensure that these initial fluxes are used for each iteration in the
* semi-explicit solution while they are not altered by the implicit solution.
*
* - qs = alpha*(Tair - Tss) : Sensible heat transfer \n
* 	- alpha is a coefficient based on the famous Monin - Obukhov theory.
*     It is strongly dependent on the wind speed. Ta and Tss are
*     air and surface temperatures (K), respectively
* - ql : Latent heat transfer \n
* 	- eA and eS are the vapor pressures of air and snow, respectively.
*     Includes consideration of soil (one active element).
* - qr = gamma*(Tair - Tss) :  Energy convected by rain \n
* - lw_net = ES*SB*(e*Tair^4 - Tss^4) : Long wave radiation heat transfer \n
* 	- For the implicit solution,this will be treated as a convective
*     boundary condition, similar to the sensible heat exchange.
*     ES = emissivity_snow, SB = stefan-boltzmann constant.
* - qg : geothermal heat flux or heat flux at lower boundary \n
*
* The fluxes are then updated after compTemperatureProfile to be able to compute a correct energy balance
* @param Mdata
* @param Xdata
*/
void Snowpack::updateBoundHeatFluxes(BoundCond& Bdata, SnowStation& Xdata, const CurrentMeteo& Mdata)
{
	// Determine actual height of meteo values above Xdata.SoilNode:
	double actual_height_of_meteo_values;	// Height with reference Xdata.SoilNode
	if(!adjust_height_of_meteo_values) {
		// Case of fixed height above snow surface (e.g., weather model)
		actual_height_of_meteo_values = height_of_meteo_values + Xdata.cH - Xdata.Ground + ( (Xdata.findMarkedReferenceLayer() == Constants::undefined) ? (0.) : (Xdata.findMarkedReferenceLayer())  - Xdata.Ground);
	} else {
		// Case of fixed height above ground surface (e.g., weather station)
		actual_height_of_meteo_values = height_of_meteo_values;
	}

	const double alpha = SnLaws::compSensibleHeatCoefficient(Mdata, Xdata, actual_height_of_meteo_values) * Constants::density_air * Constants::specific_heat_air;
	const double Tair = Mdata.ta;
	const double Tss = Xdata.Ndata[Xdata.getNumberOfNodes()-1].T;

	assert(Tair>=t_crazy_min && Tair<=t_crazy_max);
	assert(Tss>=t_crazy_min && Tss<=t_crazy_max);

	if (forcing == "ATMOS") {
		// For atmospheric forcing
		Bdata.qs = alpha * (Tair - Tss);

		Bdata.ql = SnLaws::compLatentHeat_Rh(soil_evaporation, Mdata, Xdata, actual_height_of_meteo_values);

		if (Xdata.getNumberOfElements() > 0) {
		  	// Limit fluxes in case of explicit treatment of boundary conditions
			const double theta_r = ((watertransportmodel_snow=="RICHARDSEQUATION" && Xdata.getNumberOfElements()>Xdata.SoilNode) || (watertransportmodel_soil=="RICHARDSEQUATION" && Xdata.getNumberOfElements()==Xdata.SoilNode)) ? (PhaseChange::RE_theta_threshold) : (PhaseChange::theta_r);
			const double max_ice = ((watertransportmodel_snow=="RICHARDSEQUATION" && Xdata.getNumberOfElements()>Xdata.SoilNode) || (watertransportmodel_soil=="RICHARDSEQUATION" && Xdata.getNumberOfElements()==Xdata.SoilNode)) ? (ReSolver1d::max_theta_ice * (1. - Constants::eps)) : (1.);
			if (Xdata.Edata[Xdata.getNumberOfElements()-1].theta[WATER] > theta_r + Constants::eps		// Water and ice ...
			    && Xdata.Edata[Xdata.getNumberOfElements()-1].theta[ICE] > Constants::eps			// ... coexisting
			    && Xdata.Edata[Xdata.getNumberOfElements()-1].theta[ICE] < max_ice) {
				Bdata.qs = std::min(350., std::max(-350., Bdata.qs));
				Bdata.ql = std::min(250., std::max(-250., Bdata.ql));
			}
		}
	} else {
		// For mass balance forcing
		Bdata.qs = alpha * (Tair - Tss);
		const double theta_r = ((watertransportmodel_snow=="RICHARDSEQUATION" && Xdata.getNumberOfElements()>Xdata.SoilNode) || (watertransportmodel_soil=="RICHARDSEQUATION" && Xdata.getNumberOfElements()==Xdata.SoilNode)) ? (PhaseChange::RE_theta_r) : (PhaseChange::theta_r);
		if (Mdata.sublim != IOUtils::nodata) {
			if (Xdata.getNumberOfElements() > 0 && Xdata.Edata[Xdata.getNumberOfElements()-1].theta[WATER] > theta_r) {
				// Case of evaporation
				Bdata.ql = (Mdata.sublim * Constants::lh_vaporization) / sn_dt;
			} else {
				// Case of sublimation
				Bdata.ql = (Mdata.sublim * Constants::lh_sublimation) / sn_dt;
			}
		} else {
			Bdata.ql = 0.;
		}
	}

	if (Mdata.psum>0. && Mdata.psum_ph>0.) { //there is some rain
		const double gamma = ((Mdata.psum * Mdata.psum_ph) / sn_dt) * Constants::specific_heat_water;
		Bdata.qr = gamma * (Tair - Tss);
	} else {
		Bdata.qr = 0.;
	}

	const double emmisivity = (Xdata.getNumberOfElements() > Xdata.SoilNode) ? Constants::emissivity_snow : Xdata.SoilEmissivity;

	const double lw_in  = emmisivity * Atmosphere::blkBody_Radiation(Mdata.ea, Tair);
	Bdata.lw_out = emmisivity * Constants::stefan_boltzmann * Optim::pow4(Tss);
	if (Mdata.lw_net == IOUtils::nodata) {
		// Default
		Bdata.lw_net = lw_in - Bdata.lw_out;
	} else {
		// NET_LW provided
		Bdata.lw_net = Mdata.lw_net;
	}

	if (Mdata.geo_heat != IOUtils::nodata) {
		// If geo_heat is provided in CurrentMeteo,  use it.
		Bdata.qg = Mdata.geo_heat;
	} else if (geo_heat != Constants::undefined) {
		// Otherwise check if geo_heat is defined
		Bdata.qg = geo_heat;
	} else {
		Bdata.qg = 0.;
	}
}

/**
 * @brief Imposes Neumann boundary conditions at the surface. \n
 * The fluxes are first updated in updateBoundHeatFluxes. \n
 * This splitting ensures that these initial fluxes are used for each iteration in the
 * semi-explicit solution while they are not altered by the implicit solution.
 * Note that long wave radiation heat transfer is treated as a convection boundary condition
 * for the implicit solution, similar to the sensible heat exchange: \n
 *            lw_net = delta*(e*Ta - T_iter) \n
 * where delta is a function of both Ta and T_iter and is computed by SnLaws::compLWRadCoefficient
 * @param Mdata
 * @param Bdata
 * @param Xdata
 * @param T_snow Initial (snow) surface temperature (K)
 * @param T_iter Iterated (snow) surface temperature (K)
 * @param Se Element stiffness matrix
 * @param Fe Element right hand side vector
 */
void Snowpack::neumannBoundaryConditions(const CurrentMeteo& Mdata, BoundCond& Bdata, const SnowStation& Xdata,
                                         const double& T_snow, const double& T_iter,
                                         double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ],
                                         double Fe[ N_OF_INCIDENCES ])
{
	// First zero out the interiour node contribution
	Se[0][0] = Se[0][1] = Se[1][0] = Se[1][1] = Fe[0] = Fe[1] = 0.0;

	// Special case for MASSBAL forcing, surf_melt is the energy flux.
	if(forcing=="MASSBAL") {
		if (Mdata.surf_melt != IOUtils::nodata) Fe[1] += Mdata.surf_melt * Constants::lh_fusion / sn_dt;
		return;
	}

	// Determine actual height of meteo values above Xdata.SoilNode:
	double actual_height_of_meteo_values;	// Height with reference Xdata.SoilNode
	if(!adjust_height_of_meteo_values) {
		// Case of fixed height above snow surface (e.g., weather model)
		actual_height_of_meteo_values = height_of_meteo_values + Xdata.cH - Xdata.Ground;
	} else {
		// Case of fixed height above ground surface (e.g., weather station)
		actual_height_of_meteo_values = height_of_meteo_values;
	}

	const double T_air = Mdata.ta;
	const size_t nE = Xdata.getNumberOfElements();
	const double T_s = Xdata.Edata[nE-1].Te;

	// Now branch between phase change cases (semi-explicit treatment) and
	// dry snowpack dynamics/ice-free soil dynamics (implicit treatment)
	const double theta_r = ((watertransportmodel_snow=="RICHARDSEQUATION" && Xdata.getNumberOfElements()>Xdata.SoilNode) || (watertransportmodel_soil=="RICHARDSEQUATION" && Xdata.getNumberOfElements()==Xdata.SoilNode)) ? (PhaseChange::RE_theta_threshold) : (PhaseChange::theta_r);
	const double max_ice = ((watertransportmodel_snow=="RICHARDSEQUATION" && Xdata.getNumberOfElements()>Xdata.SoilNode) || (watertransportmodel_soil=="RICHARDSEQUATION" && Xdata.getNumberOfElements()==Xdata.SoilNode)) ? (ReSolver1d::max_theta_ice * (1. - Constants::eps)) : (1.);

	if ((Xdata.Edata[nE-1].theta[WATER] > theta_r + Constants::eps		// Water and ice ...
	     && Xdata.Edata[nE-1].theta[ICE] > Constants::eps			// ... coexisting
	     && Xdata.Edata[nE-1].theta[ICE] < max_ice)
	     && variant != "SEAICE"
	     && (T_iter != T_snow)) {
		// Explicit
		// Now allow a temperature index method if desired by the user
		if ( (temp_index_degree_day > 0.) && (T_air > T_s)) {
			Fe[1] += temp_index_degree_day*(T_air - T_s) + temp_index_swr_factor*(1. - Xdata.Albedo)*Mdata.iswr;
		} else {
			Fe[1] += Bdata.ql + Bdata.lw_net + Bdata.qs + Bdata.qr;
		}
	} else { // Implicit
		// Sensible heat transfer: linear dependence on snow surface temperature
		const double alpha = SnLaws::compSensibleHeatCoefficient(Mdata, Xdata, actual_height_of_meteo_values) * Constants::density_air * Constants::specific_heat_air;
		Se[1][1] += alpha;
		Fe[1] += alpha * T_air;
		// Latent heat transfer: NON-linear dependence on snow surface temperature
		//NOTE: should it not be linearized then?
		Fe[1] += Bdata.ql;
		// Advected rain energy: linear dependence on snow surface temperature
		if (Mdata.psum > 0. && Mdata.psum_ph>0.) { //there is some rain
			const double gamma = ((Mdata.psum * Mdata.psum_ph) / sn_dt) * Constants::specific_heat_water;
			Se[1][1] += gamma;
			Fe[1] += gamma * T_air;
		}

		// Net longwave radiation: NON-linear dependence on snow surface temperature
		if (Mdata.lw_net == IOUtils::nodata) {
			const double delta = SnLaws::compLWRadCoefficient( T_iter, T_air, Mdata.ea);
			Se[1][1] += delta;
			Fe[1] += delta * pow( Mdata.ea, 0.25 ) * T_air;
		} else {
			// When NET_LW is provided, we force explicitly
			Fe[1] += Bdata.lw_net;
		}

		// Because of the implicit time integration, must subtract this term from the flux ....
		Fe[1] -= Se[1][1] * T_snow;
	} // end else
}

/**
 * @brief Imposes the Neumann boundary conditions at the bottom node. \n
 * Note that this can only be used if you have a deep enough soil,
 * because the heat flux is currently constant. For Ethan, we have had already a
 * version with upper and lower heat flux on a small snow sample and that worked perfectly.
 * @param flux Ground heat flux (W m-2)
 * @param T_snow Initial upper node temperature of bottom element (K)
 * @param Se Element stiffness matrix
 * @param Fe Element right hand side vector
 */
void Snowpack::neumannBoundaryConditionsSoil(const double& flux, const double& T_snow,
                                             double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ],
                                             double Fe[ N_OF_INCIDENCES ])
{
	// First zero out the interiour node contribution
	Se[0][0] = Se[0][1] = Se[1][0] = Se[1][1] = Fe[0] = Fe[1] = 0.0;

	// Use the numerical trick of an assumed temperature difference of 1 K over the boundary
	const double T_pseudo = T_snow - 1.;

	// For the implicit solution, we need to define a pseudo-exchange coefficient
	const double alpha = flux / (T_pseudo - T_snow); // The heat exchange coefficients
	Se[1][1] += alpha;
	Fe[1] += alpha * T_pseudo;

	// Because of the implicit time integration, must subtract this term from the flux ....
	Fe[1] -= Se[1][1] * T_snow;
}

double Snowpack::getParameterizedAlbedo(const SnowStation& Xdata, const CurrentMeteo& Mdata) const
{
	//please keep in mind that the radiation might have been tweaked in Meteo::compRadiation()
	const vector<NodeData>& NDS = Xdata.Ndata;
	const vector<ElementData>& EMS = Xdata.Edata;
	const size_t nN = Xdata.getNumberOfNodes();
	const size_t nE = Xdata.getNumberOfElements();

	double Albedo = Xdata.SoilAlb; //pure soil profile will remain with soil albedo

	// Parameterized albedo (statistical model) including correct treatment of PLASTIC and WATER_LAYER
	if (nE > Xdata.SoilNode) { //there are some non-soil layers
		size_t eAlbedo = nE-1;
		size_t marker = EMS[eAlbedo].mk % 10;

		while ((marker==8 || marker==9) && eAlbedo > Xdata.SoilNode && !Xdata.isGlacier(false)){ //If Water or ice layer (but not glacier), go one layer down
			eAlbedo--;
			marker = EMS[eAlbedo].mk % 10;
		}

		if (eAlbedo > Xdata.SoilNode && (EMS[eAlbedo].theta[SOIL] < Constants::eps2)) { // Snow, or glacier ice
			Albedo = SnLaws::parameterizedSnowAlbedo(snow_albedo, albedo_parameterization, albedo_average_schmucki, albedo_NIED_av, albedo_fixedValue, EMS[eAlbedo], NDS[eAlbedo+1].T, Mdata, Xdata, ageAlbedo);
			if (useCanopyModel && (Xdata.Cdata.height > 3.5)) { //forest floor albedo
				const double age = (forestfloor_alb) ? std::max(0., Mdata.date.getJulian() - Xdata.Edata[eAlbedo].depositionDate.getJulian()) : 0.; // day
				Albedo = (Albedo -.3)* exp(-age/7.) + 0.3;
			}
		} else { // PLASTIC, or soil
			Albedo = Xdata.SoilAlb;
		}
	}

	//enforce albedo range
	if (useCanopyModel && (Xdata.Cdata.height > 3.5)) { //forest floor albedo
		Albedo = std::max(0.05, std::min(0.95, Albedo));
	} else {
		if (research_mode) { // Treatment of "No Snow" on the ground in research mode
			const bool use_hs_meas = enforce_measured_snow_heights && (Xdata.meta.getSlopeAngle() <= Constants::min_slope_angle);
			const double hs = (use_hs_meas)? Xdata.mH - Xdata.Ground : Xdata.cH - Xdata.Ground;
			const bool snow_free_ground = (hs < 0.02) || (NDS[nN-1].T > IOUtils::C_TO_K(3.5)) || ((hs < 0.05) && (NDS[nN-1].T > IOUtils::C_TO_K(1.7)));
			if (snow_free_ground)
				Albedo = Xdata.SoilAlb;
		}

		if (!alpine3d) //for Alpine3D, the radiation has been differently computed
			Albedo = std::max(Albedo, Mdata.rswr / Constants::solcon);

		if (nE > Xdata.SoilNode) {
			// For snow
			Albedo = std::max(Constants::min_albedo, std::min(Constants::max_albedo, Albedo));
		} else {
			// For soil
			Albedo = std::max(0.05, std::min(0.95, Albedo));
		}
	}

	return Albedo;
}

double Snowpack::getModelAlbedo(const SnowStation& Xdata, CurrentMeteo& Mdata) const
{
	//please keep in mind that the radiation might have been tweaked in Meteo::compRadiation()
	const double pAlbedo = Xdata.pAlbedo;

	// Assign iswr and rswr correct values according to switch value
	if (sw_mode == "INCOMING") { // use incoming SW flux only
		Mdata.rswr = Mdata.iswr * pAlbedo;
	} else if (sw_mode == "REFLECTED") {// use reflected SW flux only
		Mdata.iswr = Mdata.rswr / pAlbedo;
	} else if (sw_mode == "BOTH") { // use measured albedo ...
		// ... while the ground is still snow covered according to HS measurements
		if (Mdata.mAlbedo != Constants::undefined) {
			if ( (!( (Mdata.mAlbedo < 2.*Xdata.SoilAlb) && ((Xdata.cH - Xdata.Ground) > 0.05)) ) && Mdata.mAlbedo <= 0.95)
				return Mdata.mAlbedo; //we have a measured albedo
			else
				Mdata.rswr = Mdata.iswr * pAlbedo;
		} else {
			// When mAlbedo is undefined, either rswr or iswr is undefined. Then, use parameterization of albedo. Note: in Main.cc, the rswr and iswr are brought in agreement when either one is missing. This is crucial!
			Mdata.rswr = Mdata.iswr * pAlbedo;
		}
	} else {
		prn_msg(__FILE__, __LINE__, "err", Mdata.date, "sw_mode = %s not implemented yet!", sw_mode.c_str());
		exit(EXIT_FAILURE);
	}

	return pAlbedo; //we do not have a measured albedo -> use parametrized
}

/**
 * @brief Computes the snow temperatures which are given by the following formula: \n
 * @par
 *             dT            d^2T
 *   rho*c(T)*----  -  k(T)*------   =  Q
 *             dt            dz^2
 * \n
 * with initial and Dirichlet and/or Neumann boundary conditions.
 *
 * T(z,t) = temperature (K);
 * rho = snow density (kg m-3); c = specific heat capacity (J K-1 kg-1); k = heat conductivity (W K-1 m-1); t = time (s);
 * \n
 * Note:  The equations are solved with a fully implicit time-integration scheme and the
 * system of finite element matrices are solved using a sparse matrix solver.
 * @param Xdata Snow profile data
 * @param[in] Mdata Meteorological forcing
 * @param Bdata Boundary conditions
 * @param[in] ThrowAtNoConvergence	If true, throw exception when temperature equation does not converge; if false, function will return false after non convergence and true otherwise.
 * @return true when temperature equation converged, false if it did not.
 */
bool Snowpack::compTemperatureProfile(const CurrentMeteo& Mdata, SnowStation& Xdata, BoundCond& Bdata, const bool& ThrowAtNoConvergence)
{
	int Ie[N_OF_INCIDENCES];                     // Element incidences
	double T0[N_OF_INCIDENCES];                  // Element nodal temperatures at t0
	double TN[N_OF_INCIDENCES];                  // Iterated element nodal temperatures
	double Se[N_OF_INCIDENCES][N_OF_INCIDENCES]; // Element stiffnes matrix
	double Fe[N_OF_INCIDENCES];                  // Element right hand side vector

	double *U=NULL, *dU=NULL, *ddU=NULL;         // Solution vectors

	// Dereference the pointers
	void *Kt = Xdata.Kt;
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;

	const size_t nN = Xdata.getNumberOfNodes();
	const size_t nE = Xdata.getNumberOfElements();

	double I0 = Mdata.iswr - Mdata.rswr; // Net irradiance perpendicular to slope

	const double theta_r = ((watertransportmodel_snow=="RICHARDSEQUATION" && Xdata.getNumberOfElements()>Xdata.SoilNode) || (watertransportmodel_soil=="RICHARDSEQUATION" && Xdata.getNumberOfElements()==Xdata.SoilNode)) ? (PhaseChange::RE_theta_threshold) : (PhaseChange::theta_r);
	if ( (temp_index_degree_day > 0.) && (nE > 0) && (EMS[nE-1].theta[WATER] > theta_r + Constants::eps)		// Water and ice ...
	   && (EMS[nE-1].theta[ICE] > Constants::eps) && (Mdata.ta > EMS[nE-1].Te))
		I0 = 0.;
	if (I0 < 0.) {
		prn_msg(__FILE__, __LINE__, "err", Mdata.date, " iswr:%lf  rswr:%lf  Albedo:%lf", Mdata.iswr, Mdata.rswr, Xdata.Albedo);
		exit(EXIT_FAILURE);
	}

	if (forcing=="MASSBAL") I0 = 0.; // no shortwave radiation absorption within snowpack with MASSBAL forcing

	// ABSORPTION OF SOLAR RADIATION WITHIN THE SNOWPACK
	// Simple treatment of radiation absorption in snow: Beer-Lambert extinction (single or multiband).
	try {
		SnLaws::compShortWaveAbsorption(sw_absorption_scheme, Xdata, I0);
	} catch(const exception&){
		prn_msg(__FILE__, __LINE__, "err", Mdata.date, "Runtime error in compTemperatureProfile");
		throw;
	}

	// TREAT AN ASSUMED ADVECTIVE HEAT SOURCE
	// Simple treatment of constant heating rate between two depths.
	if(advective_heat) {
		SnLaws::compAdvectiveHeat(Xdata, Mdata.adv_heat, heat_begin, heat_end);
	}

	// Take care of uppermost soil element
	if ((nE > Xdata.SoilNode+1) && (EMS[Xdata.SoilNode].sw_abs > EMS[Xdata.SoilNode+1].sw_abs)) {
		EMS[nE-1].sw_abs += (EMS[Xdata.SoilNode].sw_abs - EMS[Xdata.SoilNode+1].sw_abs);
		EMS[Xdata.SoilNode].sw_abs = EMS[Xdata.SoilNode+1].sw_abs;
	}

	// Set bare ground surface temperature with no soil and return
	if (nN == 1) {
		if ((Mdata.ts0 > Constants::meltfreeze_tk) && ((Mdata.ts0 - Mdata.ta) > 10.))
			NDS[0].T = (Mdata.ts0 + Mdata.ta) / 2.;
		else
			NDS[0].T = (Mdata.ts0 == IOUtils::nodata) ? (Mdata.ta) : (Mdata.ts0);
		return true;
	}

	if (Kt != NULL)
		ds_Solve(ReleaseMatrixData, (SD_MATRIX_DATA*)Kt, 0);
	ds_Initialize(static_cast<int>(nN), (SD_MATRIX_DATA**)&Kt);
	/*
	 * Define the structure of the matrix, i.e. its connectivity. For each element
	 * we compute the element incidences and pass the incidences to the solver.
	 * The solver assumes that the element incidences build a crique, i.e. the
	 * equations specified by the incidence set are all connected to each other.
	 * Initialize element data.
	*/
	for (int e = 0; e < static_cast<int>(nE); e++) {
		int Nodes[2] = {e, e+1};
		ds_DefineConnectivity( (SD_MATRIX_DATA*)Kt, 2, Nodes , 1, 0 );
	}

	/*
	 * Perform the symbolic factorization. By specifying the element incidences, we
	 * have simply declared which coefficients of the global matrix are not zero.
	 * However, when we factorize the matrix in a LU form there is some fill-in.
	 * Coefficients that were zero prior to start the factorization process will
	 * have a value different from zero thereafter. At this step the solver compute
	 * exactly how many memory is required to solve the problem and allocate this
	 * memory in order to store the numerical matrix. Then reallocate all the
	 * solution vectors.
	*/
	ds_Solve(SymbolicFactorize, (SD_MATRIX_DATA*)Kt, 0);

	// Make sure that these vectors are always available for use ....
	errno=0;
	U=(double *) realloc(U, nN*sizeof(double));
	if (errno != 0 || U==NULL) {
		free(U);
		prn_msg(__FILE__, __LINE__, "err", Date(), "%s (allocating  solution vector U)", strerror(errno));
		throw IOException("Runtime error in compTemperatureProfile", AT);
	}
	dU=(double *) realloc(dU, nN*sizeof(double));
	if (errno != 0 || dU==NULL) {
		free(U); free(dU);
		prn_msg(__FILE__, __LINE__, "err", Date(), "%s (allocating  solution vector dU)", strerror(errno));
		throw IOException("Runtime error in compTemperatureProfile", AT);
	}
	ddU=(double *) realloc(ddU, nN*sizeof(double));
	if (errno != 0 || ddU==NULL) {
		free(U); free(dU); free(ddU);
		prn_msg(__FILE__, __LINE__, "err", Date(), "%s (allocating  solution vector ddU)", strerror(errno));
		throw IOException("Runtime error in compTemperatureProfile", AT);
	}

	// Make sure that the global data structures know where the pointers are for the next integration step after the reallocation ....
	Xdata.Kt = Kt;

	// Set the temperature at the snowpack base to the prescribed value.
	// This only in case the soil_flux is not used.
	if (!(soil_flux)) {
		if ((EMS[0].theta[ICE] >= min_ice_content)) {
			// NOTE if there is water and ice in the base element, then the base temperature MUST be meltfreeze_tk
			if ((EMS[0].theta[WATER] > SnowStation::thresh_moist_snow)) {
				NDS[0].T = EMS[0].meltfreeze_tk;
			} else if (!useSoilLayers) {
				// To avoid temperatures above freezing while snow covered
				NDS[0].T = std::min(Mdata.ts0, EMS[0].meltfreeze_tk);
			} else {
				NDS[0].T = Mdata.ts0;
			}
		} else {
			NDS[0].T = Mdata.ts0;
		}
	}
	// Now treat sea ice variant, in which ocean heat flux is used already at this point to build or destroy sea ice based on the net energy balance, so just set the temperature of the lowest node to melting.
	if (variant == "SEAICE") {
		NDS[0].T = Xdata.Seaice->calculateMeltingTemperature(SeaIce::OceanSalinity);
	}

	// Copy Temperature at time0 into First Iteration
	for (size_t n = 0; n < nN; n++) {
		if(n==nN-1 && coupled_phase_changes && surfaceCode != DIRICHLET_BC) {
			//Correct the upper node, as it may have been forced to melting temperature for assessing the energy balance
			U[n] = NDS[n].T = 2. * Xdata.Edata[n-1].Te - NDS[n-1].T;
		} else {
			U[n] = NDS[n].T;
		}

		dU[n] = 0.0;
		ddU[n] = 0.0;
		if (!(U[n] > t_crazy_min && U[n] < t_crazy_max)) {
			if (alpine3d) {
				const double Tnode_orig = U[n];
				const double T_mean_down = (n>=1)? 0.5*(NDS[n].T+NDS[n-1].T) : IOUtils::nodata;
				const double T_mean_up = (n<(nN-1))? 0.5*(NDS[n].T+NDS[n+1].T) : IOUtils::nodata;
				if (T_mean_down>t_crazy_min && T_mean_down<t_crazy_max) U[n] = T_mean_down;
				else if (T_mean_up>t_crazy_min && T_mean_up<t_crazy_max) U[n] = T_mean_up;
				if (U[n] <= t_crazy_min) U[n] = .5*( IOUtils::C_TO_K(0.) + t_crazy_min); //too cold -> reset to avg(Tmin, 0C)
				if (U[n] >= t_crazy_max) U[n] = .5*( IOUtils::C_TO_K(0.) + t_crazy_max); //too hot -> reset to avg(Tmax, 0C)

				prn_msg(__FILE__, __LINE__, "err", Mdata.date, "Temperature out of bound at beginning of iteration for node %d / %d (soil node=%d)! Reset from %.2lf to %.2lf", n, nN, Xdata.SoilNode, Tnode_orig, U[n]);
				for(size_t ii=0; ii<nE; ++ii) {
					std::cout << "   " << ii << "\t" << U[ii] << " " << EMS[ii].Te << " " << U[ii+1] << " " << " " << EMS[ii].L << " " << EMS[ii].Rho << "\n";
					if (ii==Xdata.SoilNode)
						std::cout << "---- end soil ----\n";
				}
			} else {
				prn_msg(__FILE__, __LINE__, "err", Mdata.date, "Temperature out of bound at beginning of iteration!");
				prn_msg(__FILE__, __LINE__, "msg", Date(), "At node n=%d (nN=%d, SoilNode=%d): T=%.2lf", n, nN, Xdata.SoilNode, U[n]);

				free(U); free(dU); free(ddU);
				throw IOException("Runtime error in compTemperatureProfile", AT);
			}
		}
	}

	// Set the iteration counters, as well as the phase change boolean values
	unsigned int iteration = 0;   // iteration counter (not really required)
	bool NotConverged = true;     // true if iteration did not converge
	// Set the default solution routine convergence parameters
	unsigned int MaxItnTemp = 40; // maximum 40 iterations for temperature field
	double ControlTemp = 0.01;    // solution convergence to within 0.01 degC
	// Determine the displacement depth d_pump and the wind pumping speed at the surface
	const double d_pump = SnLaws::compWindPumpingDisplacement(Xdata);
	double v_pump = (nE > Xdata.SoilNode || SnLaws::wind_pump_soil)? SnLaws::compWindPumpingVelocity(Mdata, d_pump) : 0.0;

	// The temperature equation was found to show slow convergence with only 1 or 2 elements left.
	// Likely, the reason is that the LW-radiation is only approximated as linear, but in reality it is not. When only 1 or 2 elements
	// are left, their temperature gets very sensitive to energy input and during the iterations, the temperature gets out of the
	// validity range for the linearization. Therefore, we increase the MaxItnTemp for these cases:
	if (nN==3) MaxItnTemp = 200;
	if (nN==2) MaxItnTemp = 400;
	if (nN==1 || coupled_phase_changes) MaxItnTemp = 2000;

	for(size_t e = nE; e -->0; ) Xdata.Edata[e].Qph_up = Xdata.Edata[e].Qph_down = 0.;	// Reset the energy flux to the adjecent nodes due to phase changes in the element
	std::vector<double> dth_i_up(nE, 0.);							// Initialize theta[ICE] change due to phase changes at the upper adjacent node
	std::vector<double> dth_i_down(nE, 0.);							// Initialize theta[ICE] change due to phase changes at the lower adjacent node

	// IMPLICIT INTEGRATION LOOP
	bool TempEqConverged = true;	// Return value of this function compTemperatureProfile(...)
	do {
		iteration++;
		// Reset the matrix data and zero out all the increment vectors
		ds_Solve(ResetMatrixData, (SD_MATRIX_DATA*)Kt, 0);
		for (size_t n = 0; n < nN; n++) {
			ddU[n] = dU[n];
			dU[n] = 0.;
		}

		// Assemble matrix
		const double theta_rn = ((watertransportmodel_snow=="RICHARDSEQUATION" && Xdata.getNumberOfElements()>Xdata.SoilNode) || (watertransportmodel_soil=="RICHARDSEQUATION" && Xdata.getNumberOfElements()==Xdata.SoilNode)) ? (PhaseChange::RE_theta_r) : (PhaseChange::theta_r);
		double maxd = 0.;		// Tracks max. change in ice contents in domain (convergence criterion)
		for(size_t e = nE; e -->0; ) {
			if(coupled_phase_changes) {
				// Initialize the change in ice contents due to phase changes based on the energy source/sink terms at the adjacent nodes
				dth_i_up[e] = Xdata.Edata[e].Qph_up / ((Constants::density_ice * Constants::lh_fusion) / sn_dt);
				dth_i_down[e] = Xdata.Edata[e].Qph_down / ((Constants::density_ice * Constants::lh_fusion) / sn_dt);
				Xdata.Edata[e].Qph_up = Xdata.Edata[e].Qph_down = 0.;

				// Calculate the melting/freezing associated with the current temperature state
				const double max_ice = ReSolver1d::max_theta_ice;
				const double A = (Xdata.Edata[e].c[TEMPERATURE] * Xdata.Edata[e].Rho) / ( Constants::density_ice * Constants::lh_fusion );
				const double dth_i_up_in = dth_i_up[e];
				const double dth_i_down_in = dth_i_down[e];

				if (Xdata.Seaice != NULL) {
					// For sea ice, balance the meltfreeze_tk with assuming thermal equilibrium with the brine:
					// (1): Xdata.Edata[e].meltfreeze_tk = -SeaIce::mu * BrineSal_new + Constants::meltfreeze_tk;
					// (2): BrineSal_new = (Xdata.Edata[e].salinity /  (Xdata.Edata[e].theta[WATER] + deltaTheta));
					// (3): deltaTheta = A * (0.5 * (U[e+1] + U[e]) - Xdata.Edata[e].meltfreeze_tk) * (Constants::density_water / Constants::density_ice);
					// Balancing equations (1), (2) and (3) derived using wxmaxima:
					// T=-m*s/(th+(A*(u-T)))+c
					// solve(%i1,T);
					// With:
					// T = Xdata.Edata[e].meltfreeze_tk
					// m = SeaIce::mu
					// s = Xdata.Edata[e].salinity
					// th = tmp_Theta
					// A = A * f
					// u = tmp_T
					// c = Constants::meltfreeze_tk
					const double f = Constants::density_ice / Constants::density_water;
					const double tmp_T = 0.5 * (U[e+1] + U[e]);
					const double tmp_Theta = Xdata.Edata[e].theta[WATER] - 0.5 * (dth_i_up[e] + dth_i_down[e]) * f;
					const double BrineSal_new = (tmp_Theta == 0.) ? (0.) : (Xdata.Edata[e].salinity / tmp_Theta);

					std::pair<double, double> mu = Xdata.Seaice->getMu(BrineSal_new);
					if(mu.second < -Constants::eps2) {
						Xdata.Edata[e].meltfreeze_tk = -1. * (sqrt(A * f * A * f * tmp_T * tmp_T + (2. * A * f * tmp_Theta - 2. * A * f * A * f * mu.first) * tmp_T + tmp_Theta * tmp_Theta - 2. * A * f * mu.first * tmp_Theta - 4. * A * f * mu.second * Xdata.Edata[e].salinity + A * f * A * f * mu.first * mu.first) - A * f * tmp_T - tmp_Theta - A * f * mu.first) / (2. * A * f);
					} else {
						Xdata.Edata[e].meltfreeze_tk = IOUtils::C_TO_K(0.);
					}
				}

				dth_i_up[e] += A * (Xdata.Edata[e].meltfreeze_tk - U[e+1]);	// change in volumetric ice content in upper half of element
				dth_i_down[e] += A * (Xdata.Edata[e].meltfreeze_tk - U[e]);	// change in volumetric ice content in lower half of element

				// This approach is not stable, may introduce oscillations such that the temperature equation doesn't converge
				const double dth_i_sum = 0.5 * (dth_i_up[e] + dth_i_down[e]);	// Net phase change effect on ice content in element
				if(dth_i_sum != 0.) {	// Element has phase changes
					double dth_i_lim = dth_i_sum;
					if(dth_i_lim < 0.) {
						// Melt: Only available ice can melt
						dth_i_lim = std::max(-Xdata.Edata[e].theta[ICE], dth_i_lim);
					} else {
						// Freeze: Only available liquid water can freeze, and not more than max_ice
						dth_i_lim = std::min(std::max(0., std::min(max_ice - Xdata.Edata[e].theta[ICE], (Xdata.Edata[e].theta[WATER] - theta_rn) * (Constants::density_water / Constants::density_ice))), dth_i_lim);
					}
					// Correct volumetric changes in upper and lower half of element proportional to limits
					dth_i_down[e] = dth_i_up[e] = dth_i_lim;
				}

				// Track max. abs. change in ice contents
				maxd = std::max(maxd, fabs(dth_i_up[e] - dth_i_up_in));
				maxd = std::max(maxd, fabs(dth_i_down[e] - dth_i_down_in));

				// Recalculate phase change energy
				Xdata.Edata[e].Qph_up = (dth_i_up[e] * Constants::density_ice * Constants::lh_fusion) / sn_dt;
				Xdata.Edata[e].Qph_down = (dth_i_down[e] * Constants::density_ice * Constants::lh_fusion) / sn_dt;

				if (Xdata.Seaice != NULL) {
					// Adjust melting/freezing point assuming thermal quilibrium in the brine pockets
					const double ThetaWater_new = (Xdata.Edata[e].theta[WATER] - 0.5 * (dth_i_up[e] + dth_i_down[e]) * (Constants::density_ice / Constants::density_water));
					const double BrineSal_new = (ThetaWater_new == 0.) ? (0.) : (Xdata.Edata[e].salinity / ThetaWater_new);
					Xdata.Edata[e].meltfreeze_tk = Xdata.Seaice->calculateMeltingTemperature(BrineSal_new);
				}
			}
			EL_INCID( static_cast<int>(e), Ie );
			EL_TEMP( Ie, T0, TN, NDS, U );
			// Update the wind pumping velocity gradient
			const double dvdz = SnLaws::compWindGradientSnow(EMS[e], v_pump);
			// Update the water vapor transport enhancement factor
			const double VaporEnhance = std::max(1., SnLaws::compEnhanceWaterVaporTransportSnow(Xdata, e));
			// Now fill the matrix
			if (!sn_ElementKtMatrix(EMS[e], sn_dt, dvdz, T0, Se, Fe, VaporEnhance)) {
				prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Error in sn_ElementKtMatrix @ element %d:", e);
				for (size_t n = 0; n < nN; n++)
					fprintf(stdout, "U[%u]=%g K\n", (unsigned int)n, U[n]);
				free(U); free(dU); free(ddU);
				throw IOException("Runtime error in compTemperatureProfile", AT);
			}
			ds_AssembleMatrix( (SD_MATRIX_DATA*)Kt, 2, Ie, 2,  (double*) Se );
			EL_RGT_ASSEM( dU, Ie, Fe );
		}

		/*
		 * The top element is special in that it handles the entire meteo conditions
		 * Several terms must be added to the global stiffness matrix Kt and flux
		 * right-hand side vector dU. Note:  Shortwave radiation --- since it is a body
		 * or volumetric force --- is computed in sn_ElementKtMatrix().
		 */

		if (surfaceCode == NEUMANN_BC) {
			EL_INCID(static_cast<int>(nE-1), Ie);
			EL_TEMP(Ie, T0, TN, NDS, U);
			neumannBoundaryConditions(Mdata, Bdata, Xdata, T0[1], TN[1], Se, Fe);
			ds_AssembleMatrix( (SD_MATRIX_DATA*)Kt, 2, Ie, 2,  (double*) Se );
			EL_RGT_ASSEM( dU, Ie, Fe );
		}

		double Big = Constants::big;	// big number for DIRICHLET boundary conditions)
		if (surfaceCode == DIRICHLET_BC) {
			// Dirichlet BC at surface: prescribed temperature value
			// NOTE Insert Big at this location to hold the temperature constant at the prescribed value.
			Ie[0] = static_cast<int>(nE);
			ds_AssembleMatrix((SD_MATRIX_DATA*) Kt, 1, Ie, 1, &Big);
		}
		// Bottom node
		if (soil_flux && variant != "SEAICE") {
			// Neumann BC at bottom: The lower boundary is now a heat flux -- put the heat flux in dU[0]
			EL_INCID(0, Ie);
			EL_TEMP(Ie, T0, TN, NDS, U);
			neumannBoundaryConditionsSoil(Bdata.qg, T0[1], Se, Fe);
			ds_AssembleMatrix((SD_MATRIX_DATA*)Kt, 2, Ie, 2, (double*) Se);
			EL_RGT_ASSEM(dU, Ie, Fe);
		} else if ((Xdata.getNumberOfElements() < 3) && (Xdata.Edata[0].theta[WATER] >= 0.9 * Xdata.Edata[0].res_wat_cont)) {
			dU[0] = 0.;
		} else {
			// Dirichlet BC at bottom: prescribed temperature value
			// NOTE Insert Big at this location to hold the temperature constant at the prescribed value.
			Ie[0] = 0;
			ds_AssembleMatrix((SD_MATRIX_DATA*) Kt, 1, Ie, 1, &Big);
		}

		/*
		 * Solve the linear system of equation. The te_F vector is used first as right-
		 * hand-side vector for the linear system. The solver stores in this vector
		 * the solution of the system of equations, the new temperature.
		 * It will throw an exception whenever the linear solver failed
		 */
		if (!ds_Solve(ComputeSolution, (SD_MATRIX_DATA*) Kt, dU)) {
			  prn_msg(__FILE__, __LINE__, "err", Mdata.date,
			  "Linear solver failed to solve for dU on the %d-th iteration.",
			  iteration);
			  throw IOException("Runtime error in compTemperatureProfile", AT);
		}
		// Update the solution vectors and check for convergence
		for (size_t n = 0; n < nN; n++)
			ddU[n] = dU[n] - ddU[n];
		double MaxTDiff = fabs(ddU[0]);  // maximum temperature difference for convergence
		for (size_t n = 1; n < nN; n++) {
			const double TDiff = fabs(ddU[n]); // temperature difference for convergence check
			if (TDiff > MaxTDiff)
				MaxTDiff = TDiff;
		}
		/*
		 * This is to increase the number of iterations for a phase changing uppermost element.
		 * Otherwise we would violate energy conservation because the implicit scheme
		 * leads to an adaptation of surface fluxes to the iterated surface temperature.
		 * In reality, the surface temperature will not change and therefore the fluxes
		 * must be constant. This means that the fluxes must be treated explicitely
		 * (see neumannBoundaryConditions)
		 */
		if (U[nE] + ddU[nE] > EMS[nE-1].meltfreeze_tk || EMS[nE-1].theta[WATER] > 0.) {
			ControlTemp = (variant == "SEAICE") ? (0.0001) : (0.007);
			MaxItnTemp = std::max(MaxItnTemp, (unsigned)200); // NOTE originally 100;
		}
		if(coupled_phase_changes) {
			// With new phase change, we want at least one iteration extra, to account for possible phase changes,
			// and we want an additional constraint of maximum change in phase change amount
			NotConverged = (MaxTDiff > ControlTemp || iteration == 1 || maxd > ((variant == "SEAICE") ? (1.E-5) : (0.0001)));
		} else {
			NotConverged = (MaxTDiff > ControlTemp);
		}
		if (iteration > MaxItnTemp) {
			if (ThrowAtNoConvergence) {
				prn_msg(__FILE__, __LINE__, "err", Mdata.date,
				        "Temperature did not converge (azi=%.0lf, slope=%.0lf)!",
				        Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
				prn_msg(__FILE__, __LINE__, "msg", Date(),
				        "%d iterations > MaxItnTemp=%d; ControlTemp=%.4lf; nN=%d; sn_dt=%f",
				        iteration, MaxItnTemp, ControlTemp, nN, sn_dt);
				for (size_t n = 0; n < nN; n++) {
					if (n > 0)
						prn_msg(__FILE__, __LINE__, "msg-", Date(),
						        "U[%03d]:%6.1lf K, ddU:%8.4lf K;  NDS.T(t-1)=%6.1lf K; EMS[n-1].th_w(t-1)=%.5lf",
						        n, U[n], ddU[n], NDS[n].T, EMS[n-1].theta[WATER]);
					else
						prn_msg(__FILE__, __LINE__, "msg-", Date(),
						        "U[%03d]:%6.1lf K, ddU:%8.4lf K;  NDS.T(t-1)=%6.1lf K;",
						        n, U[n], ddU[n], NDS[n].T);
				}
				prn_msg(__FILE__, __LINE__, "msg", Date(),
				        "Latent: %lf  Sensible: %lf  Rain: %lf  NetLong:%lf  NetShort: %lf",
				        Bdata.ql, Bdata.qs, Bdata.qr, Bdata.lw_net, I0);
				free(U); free(dU); free(ddU);
				throw IOException("Runtime error in compTemperatureProfile", AT);
			} else {
				TempEqConverged = false;	// Set return value of function
				NotConverged = false;		// Ensure we leave the do...while loop
			}
		}
		for (size_t n = 0; n < nN; n++) {
			U[n] += ddU[ n ];
			// If the solver converged, but we are seeing crazy nodes, and the function is not requested to throw at no convergence
			// we set the flag TempEqConverged to false, such that the solver can try with a smaller time step
			if ( ! ((U[n] > t_crazy_min) && (U[n] < t_crazy_max)) ) {
				if ( !NotConverged && !ThrowAtNoConvergence ) TempEqConverged = false;
			}
		}
	} while ( NotConverged ); // end Convergence Loop

	if (TempEqConverged) {
		size_t crazy = 0;
		bool prn_date = true;
		for (size_t n = 0; n < nN; n++) {
			if ((U[n] > t_crazy_min) && (U[n] < t_crazy_max)) {
				NDS[n].T = U[n];
			} else { // Correct for - hopefully transient - crazy temperatures!
				NDS[n].T = 0.5*(U[n] + NDS[n].T);
				if (!alpine3d) { //reduce number of warnings for Alpine3D
					if (!crazy && (nN > Xdata.SoilNode + 3)) {
						prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "Crazy temperature(s)!");
						prn_msg(__FILE__, __LINE__, "msg-", Date(),
						        "T <= %5.1lf OR T >= %5.1lf; nN=%d; cH=%6.3lf m; azi=%.0lf, slope=%.0lf",
						        t_crazy_min, t_crazy_max, nN, Xdata.cH,
						        Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
						prn_date = false;
					}
					if (n < Xdata.SoilNode) {
						if (n == 0) {
							prn_msg(__FILE__, __LINE__, "msg-", Mdata.date,
							        "Bottom SOIL node %3d: U(t)=%6.1lf  NDS.T(t-1)=%6.1lf K", n, U[n], NDS[n].T);
							prn_date = false;
						} else {
							prn_msg(__FILE__, __LINE__, "msg-", Date(),
							        "SOIL node %3d: U(t)=%6.1lf  NDS.T(t-1)=%6.1lf K; EMS[n-1].th_w(t-1)=%.5lf",
							        n, U[n], NDS[n].T, EMS[n-1].theta[WATER]);
						}
					} else if (Xdata.SoilNode > 0) {
							prn_msg(__FILE__, __LINE__, "msg-", Date(),
							        "GROUND surface node %3d: U(t)=%6.1lf  NDS.T(t-1)=%6.1lf", n, U[n], NDS[n].T);
					} else if (nN > Xdata.SoilNode + 3) {
						if (n == 0) {
							prn_msg(__FILE__, __LINE__, "msg-", Mdata.date,
							        "Bottom SNOW node %3d: U(t)=%6.1lf  NDS.T(t-1)=%6.1lf", n, U[n], NDS[n].T);
							prn_date = false;
						} else {
							prn_msg(__FILE__, __LINE__, "msg-", Date(),
							        "SNOW node %3d: U(t)=%6.1lf  NDS.T(t-1)=%6.1lf EMS[n-1].th_w(t-1)=%.5lf",
							        n, U[n], NDS[n].T, EMS[n-1].theta[WATER]);
						}
					}
				}
				crazy++;
			}
		}
		if (crazy > Xdata.SoilNode + 3) {
			if (prn_date)
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date,
				        "%d crazy node(s) from total %d! azi=%.0lf, slope=%.0lf",
				        crazy, nN, Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
			else
				prn_msg(__FILE__, __LINE__, "msg-", Date(),
				        "%d crazy node(s) from total %d! azi=%.0lf, slope=%.0lf",
				        crazy, nN, Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
			prn_msg(__FILE__, __LINE__, "msg-", Date(),
			        "Latent: %lf  Sensible: %lf  Rain: %lf  NetLong:%lf  NetShort: %lf",
			        Bdata.ql, Bdata.qs, Bdata.qr, Bdata.lw_net, I0);
		}

		for (size_t e = 0; e < nE; e++) {
			EMS[e].Te = (NDS[e].T + NDS[e+1].T) / 2.0;
			EMS[e].heatCapacity();
			EMS[e].gradT = (NDS[e+1].T - NDS[e].T) / EMS[e].L;
		}
	}
	free(U); free(dU); free(ddU);
	if (coupled_phase_changes) {
		// Ensure that when top element consists of ice, its upper node does not exceed melting temperature
		// This is to have consistent surface energy balance calculation and for having good looking output
		// Note: for sea ice, the effect of salinity is such that this doesn't work...
		if (nE > 0 && Xdata.Edata[nE-1].theta[ICE] > Constants::eps && Xdata.Edata[nE-1].salinity > Constants::eps2) NDS[nE].T=std::min(Xdata.Edata[nE-1].meltfreeze_tk, NDS[nE].T);
	}

	return TempEqConverged;
}

/**
 * @brief Set the microstructure parameters for the current element according to the type of precipitation meteor.
 * @param Mdata Meteorological data
 * @param is_surface_hoar is this layer a layer of surface hoar?
 * @param EMS Element to set
 */
void Snowpack::setHydrometeorMicrostructure(const CurrentMeteo& Mdata, const bool& is_surface_hoar, ElementData &elem)
{
	const double TA = IOUtils::K_TO_C(Mdata.ta);
	const double RH = Mdata.rh*100.;
	const double logit = 49.6 + 0.857*Mdata.vw - 0.547*RH;
	const double value = exp(logit)/(1.+exp(logit));

	// Distinguish between Graupel and New Snow
	if (value > 1.0) { // Graupel
		elem.mk = 4;
		elem.dd = 0.;
		elem.sp = 1.;
		elem.rg = 0.6;
		elem.rb = 0.2;
		// Because density and volumetric contents are already defined, redo it here
		elem.Rho = 110.;
		elem.theta[ICE] = elem.Rho / Constants::density_ice;  // ice content
		elem.theta_i_reservoir = 0.0;
		elem.theta_i_reservoir_cumul = 0.0;
		elem.theta[AIR] = 1. - elem.theta[ICE];  // void content
	} else { // no Graupel
		elem.mk = Snowpack::new_snow_marker;
		if (SnLaws::jordy_new_snow && (Mdata.vw > 2.9)
			&& ((hn_density_parameterization == "LEHNING_NEW") || (hn_density_parameterization == "LEHNING_OLD"))) {
			elem.dd = std::max(0.5, std::min(1.0, Optim::pow2(1.87 - 0.04*Mdata.vw)) );
			elem.sp = new_snow_sp;
			static const double alpha = 0.9, beta = 0.015, gamma = -0.0062;
			static const double delta = -0.117, eta=0.0011, phi=-0.0034;
			elem.rg = std::min(0.5*new_snow_grain_size, std::max(0.15*new_snow_grain_size,
				alpha + beta*TA + gamma*RH + delta*Mdata.vw
				+ eta*RH*Mdata.vw + phi*TA*Mdata.vw));
			elem.rb = 0.4*elem.rg;
		} else {
			elem.dd = new_snow_dd;
			elem.sp = new_snow_sp;
			// Adapt dd and sp for blowing snow
			if ((Mdata.vw > 5.) && ((variant == "ANTARCTICA" || variant == "POLAR")
			|| (!SnLaws::jordy_new_snow && ((hn_density_parameterization == "BELLAIRE")
			|| (hn_density_parameterization == "LEHNING_NEW"))))) {
				elem.dd = new_snow_dd_wind;
				elem.sp = new_snow_sp_wind;
			} else if (vw_dendricity && ((hn_density_parameterization == "BELLAIRE")
				|| (hn_density_parameterization == "ZWART"))) {
				const double vw = std::max(0.05, Mdata.vw);
				elem.dd = (1. - pow(vw/10., 1.57));
				elem.dd = std::max(0.2, elem.dd);
			}
			if (Snowpack::hydrometeor) { // empirical
				static const double alpha=1.4, beta=-0.08, gamma=-0.15;
				static const double delta=-0.02;
				elem.rg = 0.5*(alpha + beta*TA + gamma*Mdata.vw + delta*TA*Mdata.vw);
				elem.rb = 0.25*elem.rg;
			} else {
				elem.rg = new_snow_grain_size/2.;
				elem.rb = new_snow_bond_size/2.;
				if (((Mdata.vw_avg >= SnLaws::event_wind_lowlim) && (Mdata.rh_avg >= rh_lowlim))) {
					elem.rb = std::min(bond_factor_rh*elem.rb, Metamorphism::max_grain_bond_ratio*elem.rg);
				}
			}
		}
	} // end no Graupel

	if(is_surface_hoar) { //surface hoar
		elem.mk = 3;
		elem.dd = 0.;
		elem.sp = 0.;
		elem.rg = std::max(new_snow_grain_size/2., 0.5*M_TO_MM(elem.L0)); //Note: L0 > hoar_min_size_buried/hoar_density_buried
		elem.rb = elem.rg/3.;
	}

	elem.opticalEquivalentGrainSize();
	elem.metamo = 0.;
	elem.soot_ppmv= soot_ppmv;
}

void Snowpack::fillNewSnowElement(const CurrentMeteo& Mdata, const double& length, const double& density,
                                  const bool& is_surface_hoar, const unsigned short& number_of_solutes, ElementData &elem)
{
	//basic parameters
	elem.depositionDate = Mdata.date;
	elem.Te = t_surf;
	elem.L0 = elem.L = length;
	elem.Rho = density;
	assert(elem.Rho>=0. || elem.Rho==IOUtils::nodata); //we want positive density
	elem.M = elem.L0*elem.Rho; // Mass
	assert(elem.M>=0.); //mass must be positive

	// Volumetric components
	elem.theta[SOIL]  = 0.0;
	elem.theta[ICE]   = elem.Rho/Constants::density_ice;
	elem.theta_i_reservoir = 0.0;
	elem.theta_i_reservoir_cumul = 0.0;
	elem.theta[WATER] = 0.0;
	elem.theta[WATER_PREF] = 0.0;
	elem.theta[AIR]   = 1. - elem.theta[ICE];
	for (unsigned short ii = 0; ii < number_of_solutes; ii++) {
		elem.conc[ICE][ii]   = Mdata.conc[ii]*Constants::density_ice/Constants::density_water;
		elem.conc[WATER][ii] = Mdata.conc[ii];
		elem.conc[AIR][ii]   = 0.0;
		elem.conc[SOIL][ii]  = 0.0;
	}

	// Coordination number based on Bob's empirical function
	elem.N3 = Metamorphism::getCoordinationNumberN3(elem.Rho);

	// Constitutive Parameters
	elem.k[TEMPERATURE] = elem.k[SEEPAGE] = elem.k[SETTLEMENT]= 0.0;
	elem.heatCapacity();
	elem.c[SEEPAGE] = elem.c[SETTLEMENT]= 0.0;
	elem.soil[SOIL_RHO] = elem.soil[SOIL_K] = elem.soil[SOIL_C] = 0.0;
	elem.snowResidualWaterContent();

	// Phase change variables:
	elem.sw_abs = 0.0; //initial short wave radiation
	elem.dth_w=0.0; // change of water content
	elem.Qmf=0.0;   // change of energy due to phase changes
	// Total element strain (GREEN'S strains -- TOTAL LAGRANGIAN FORMULATION.
	elem.E = elem.dEps = elem.Eps = elem.Eps_e = elem.Eps_v = 0.0;
	elem.Eps_Dot = elem.Eps_vDot=0.0; // Total Strain Rate (Simply E/sn_dt)
	elem.S=0.0; // Total Element Stress
	elem.CDot = 0.; // loadRate

	//new snow micro-structure
	setHydrometeorMicrostructure(Mdata, is_surface_hoar, elem);
	elem.snowType(); // Snow classification

	//Initialise the Stability Index for ml_st_CheckStability routine
	elem.S_dr = IOUtils::nodata;
	elem.hard = IOUtils::nodata;

	elem.h = Constants::undefined;	//Pressure head not initialized yet

	//Initial snow salinity
	if (variant == "SEAICE" ) elem.salinity = SeaIce::InitSnowSalinity;

	double p_vapor = Atmosphere::vaporSaturationPressure(elem.Te);
	elem.rhov = Atmosphere::waterVaporDensity(elem.Te, p_vapor);
}

/**
 * @brief Introduce new snow elements as technical snow
 * @details When there is natural snow as well as man-made snow,
 * the whole snow fall will have the properties of man-made snow.
 * @param Mdata Meteorological data
 * @param Xdata Snow cover data
 * @param cumu_precip cumulated amount of precipitation (kg m-2)
 */
void Snowpack::compTechnicalSnow(const CurrentMeteo& Mdata, SnowStation& Xdata, double& cumu_precip)
{
	const size_t nOldN = Xdata.getNumberOfNodes(); //Old number of nodes
	const size_t nOldE = Xdata.getNumberOfElements(); //Old number of elements
	const double cos_sl = Xdata.cos_sl; //slope cosinus

 	double Tw, rho_hn, delta_cH, theta_w;
	TechSnow::productionPpt(Mdata, cumu_precip, Tw, rho_hn, delta_cH, theta_w);

	// Now determine whether the increase in snow depth is large enough.
	double hn = 0.; //new snow amount
	if ( (delta_cH >= height_new_elem * cos_sl) ) {
		cumu_precip = 0.0; // we use the mass through delta_cH
		hn = delta_cH;
	}
	if (hn > Snowpack::snowfall_warning)
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date,
				          "Large snowfall! hn=%.3f cm (azi=%.0f, slope=%.0f)",
				            M_TO_CM(hn), Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());

	const size_t nAddE = (size_t)(hn / (height_new_elem*cos_sl));

	if (nAddE < 1) return;

	Xdata.Albedo = Constants::max_albedo;

	const size_t nNewN = nOldN + nAddE;
	const size_t nNewE = nOldE + nAddE;
	Xdata.resize(nNewE);
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;

	// Fill the nodal data
	if (!useSoilLayers && (nOldN-1 == Xdata.SoilNode)) // New snow on bare ground w/o soil
		NDS[nOldN-1].T = Tw;	// 0.5*(t_surf + Mdata.ta);
	const double Ln = (hn / (double)nAddE);               // New snow element length
	double z0 = NDS[nOldN-1].z + NDS[nOldN-1].u + Ln; // Position of lowest new node
	for (size_t n = nOldN; n < nNewN; n++) { //loop over the nodes
		NDS[n].T = Tw;                      // t_surf Temperature of the new node
		NDS[n].z = z0;                      // New nodal position
		NDS[n].u = 0.0;                     // Initial displacement is 0
		NDS[n].hoar = 0.0;                  // The new snow surface hoar is set to zero
		NDS[n].udot = 0.0;                  // Settlement rate is also 0
		NDS[n].f = 0.0;                     // Unbalanced forces are 0
		NDS[n].S_n = IOUtils::nodata;
		NDS[n].S_s = IOUtils::nodata;
		z0 += Ln;
	}

	// Fill the element data
	for (size_t e = nOldE; e < nNewE; e++) { //loop over the elements
		const double length = (NDS[e+1].z + NDS[e+1].u) - (NDS[e].z + NDS[e].u);
		fillNewSnowElement(Mdata, length, rho_hn, false, Xdata.number_of_solutes, EMS[e]);

		// Now give specific properties for technical snow, consider liquid water
		// Assume that the user does not specify unreasonably high liquid water contents.
		// This depends also on the density of the solid fraction - print a warning if it looks bad
		EMS[e].theta[WATER] += theta_w;

		if ( (EMS[e].theta[WATER] + EMS[e].theta[ICE]) > 0.7)
			prn_msg(__FILE__, __LINE__, "wrn", Mdata.date,
		          "Too much liquid water specified or density too high! Dry density =%.3f kg m-3  Water Content = %.3f %", rho_hn, theta_w);

		EMS[e].theta[AIR] = 1.0 - EMS[e].theta[WATER] - EMS[e].theta[WATER_PREF] - EMS[e].theta[ICE] - EMS[e].theta[SOIL];

		if (EMS[e].theta[AIR] < 0.) {
			prn_msg(__FILE__, __LINE__, "err", Mdata.date, "Error in technical snow input - no void fraction left");
			throw IOException("Runtime error in runSnowpackModel", AT);
		}

		// To satisfy the energy balance, we should trigger an explicit treatment of the top boundary condition of the energy equation
		// when new snow falls on top of wet snow or melting soil. This can be done by putting a tiny amount of liquid water in the new snow layers.
		// Note that we use the same branching condition as in the function Snowpack::neumannBoundaryConditions(...)
		const double theta_r = ((watertransportmodel_snow=="RICHARDSEQUATION" && Xdata.getNumberOfElements()>Xdata.SoilNode) || (watertransportmodel_soil=="RICHARDSEQUATION" && Xdata.getNumberOfElements()==Xdata.SoilNode)) ? (PhaseChange::RE_theta_threshold) : (PhaseChange::theta_r);
		if(nOldE > 0 && EMS[nOldE-1].theta[WATER] > theta_r + Constants::eps && EMS[nOldE-1].theta[ICE] > Constants::eps) {
			EMS[e].theta[WATER]+=(2.*Constants::eps);
			EMS[e].theta[ICE]-=(2.*Constants::eps)*(Constants::density_water/Constants::density_ice);
			EMS[e].theta[AIR]+=((Constants::density_water/Constants::density_ice)-1.)*(2.*Constants::eps);
		}
		EMS[e].meltfreeze_tk = Constants::meltfreeze_tk;
		Xdata.ColdContent += EMS[e].coldContent(); //update cold content

		// Now adjust default new element values to technical snow (mk = 6)
		EMS[e].mk = 6;
		EMS[e].dd = 0.;
		EMS[e].sp = 1.;
		EMS[e].rg = 0.2; // Have to adapt after some tests
		EMS[e].rb = EMS[e].rg/3.;

	}   // End elements

	// Finally, update the computed snowpack height
	Xdata.cH = NDS[nNewN-1].z + NDS[nNewN-1].u;
	Xdata.ErosionLevel = nNewE-1;

}

/**
 * @brief Determines whether new snow elements are added on top of the snowpack
 * - If enforce_measured_snow_heights=0 (research mode), the new snow height corresponding to the cumulated
 *   new snow water equivalent cumu_precip must be greater than HEIGHT_NEW_ELEM to allow adding elements.
 * - In case of virtual slopes, uses new snow depth and density from either flat field or luv slope
 * - The first thing is to compute the height of each element in the snow layer. For now,
 *   instead of trying to find an optimal number of elements, we will define the number of new
 *   elements as a constant nNnew. The constant is related to HEIGHT_NEW_ELEM, which is set in
 *   qr_ReadParameter(). The smaller HEIGHT_NEW_ELEM, the more elements we will use.
 * @param Mdata Meteorological data
 * @param Xdata Snow cover data
 * @param cumu_precip cumulated amount of precipitation (kg m-2)
 */
void Snowpack::compSnowFall(const CurrentMeteo& Mdata, SnowStation& Xdata, double& cumu_precip,
                            SurfaceFluxes& Sdata)
{
	if (Mdata.psum_tech!=Constants::undefined && Mdata.psum_tech > 0.) {
		compTechnicalSnow(Mdata, Xdata, cumu_precip);
		return;
	}

	bool add_element = false;
	double delta_cH = 0.; // Actual enforced snow depth
	double hn = 0.; //new snow amount

	//Threshold for detection of the first snow fall on soil/canopy (grass/snow detection)
	static const double TSS_threshold24=-1.5;			//deg Celcius of 24 hour average TSS
	static const double TSS_threshold12_smallHSincrease=-0.5;	//deg Celcius of 12 hour average TSS in case of low rate of change of HS
	static const double TSS_threshold12_largeHSincrease=3.0;	//deg Celcius of 12 hour average TSS in case of high rate of change of HS
	static const double HS_threshold_smallincrease=0.005;		//low rate of change of HS (m/hour)
	static const double HS_threshold_largeincrease=0.010;		//high rate of change of HS (m/hour)
	static const double HS_threshold_verylargeincrease=0.015;	//very high rate of change of HS (m/hour)
	static const double ThresholdSmallCanopy=1.;			//Set the threshold for the canopy height. Below this threshold, the canopy is considered to be small and snow will fall on top (like grass, or small bushes). Above this threshold, snow will fall through (like in forests). When the canopy is small, the measured snow height will be assigned to the canopy height in case of no snow pack on the ground.

	const size_t nOldN = Xdata.getNumberOfNodes(); //Old number of nodes
	const size_t nOldE = Xdata.getNumberOfElements(); //Old number of elements
	const double cos_sl = Xdata.cos_sl; //slope cosinus

	double rho_hn = SnLaws::compNewSnowDensity(hn_density, hn_density_parameterization, hn_density_fixedValue,
                                               Mdata, Xdata, t_surf, variant);

	if ((Sdata.cRho_hn < 0.) && (rho_hn != Constants::undefined))
		Sdata.cRho_hn = -rho_hn;

	if (!enforce_measured_snow_heights) { // PSUM driven
		if (Mdata.psum>0. && Mdata.psum_ph<1.) { //there is some snow
			const double precip_snow = Mdata.psum * (1. -  Mdata.psum_ph);
			const double precip_rain = Mdata.psum * Mdata.psum_ph;
			if ((cumu_precip > 0.) && (rho_hn != Constants::undefined)) {
				// This is now very important to make sure that the rain fraction will not accumulate
				// Note that cumu_precip will always be considered snowfall, as we substract all rainfall amounts
				cumu_precip -= precip_rain;
				if ((hn_density == "MEASURED") || ((hn_density == "FIXED") && (rho_hn > SnLaws::max_hn_density))) {
					// Make sure that a new element is timely added in the above cases
					// TODO check whether needed in both cases
					if (((meteo_step_length / sn_dt) * (precip_snow)) <= cumu_precip || forcing == "MASSBAL" || force_add_snowfall) {
						delta_cH = (cumu_precip / rho_hn);
						add_element = true;
					}
				} else if ((cumu_precip / rho_hn) > height_new_elem*cos_sl || force_add_snowfall) {
					delta_cH = (cumu_precip / rho_hn);
					if (hn_density == "EVENT") { // TODO check whether needed
						add_element = true;
					}
				}
			} else
				return;
		} else {
			// This is now very important to make sure that rain will not accumulate
			cumu_precip -= Mdata.psum; //if there is no precip, this does nothing
			return;
		}
	} else { // HS driven, correct for a possible offset in measured snow height provided by a marked reference layer
		delta_cH = Xdata.mH - Xdata.cH + ( (Xdata.findMarkedReferenceLayer() == Constants::undefined) ? (0.) : (Xdata.findMarkedReferenceLayer() - Xdata.Ground) );
	}
	if (rho_hn == Constants::undefined)
		return;

	// Let's check for the solid precipitation thresholds:
	// -> check relative humidity as well as difference between air and snow surface temperatures,
	//    that is, no new snow during cloud free conditions!
	const double meltfreeze_tk = (nOldE>0)? Xdata.Edata[nOldE-1].meltfreeze_tk : Constants::meltfreeze_tk;
	const double dtempAirSnow = (change_bc && !meas_tss)? Mdata.ta - meltfreeze_tk : Mdata.ta - t_surf; //we use t_surf only if meas_tss & change_bc

	const bool snow_fall = (((((Mdata.rh > thresh_rh) && (Mdata.psum_ph<1.) && (dtempAirSnow < thresh_dtempAirSnow))
                               || !enforce_measured_snow_heights || (Xdata.hn > 0.)) && (forcing=="ATMOS")) || ((Mdata.psum_ph<1.) && (forcing=="MASSBAL")));

	// In addition, let's check whether the ground is already snowed in or cold enough to build up a snowpack
	bool snowed_in = false;
	if (!enforce_measured_snow_heights || !detect_grass) {
		snowed_in = true;
	} else {
		snowed_in = ((Xdata.getNumberOfNodes() > Xdata.SoilNode+1)
		            || (detect_grass &&
		                (((Mdata.tss_a24h != IOUtils::nodata && Mdata.tss_a24h < IOUtils::C_TO_K(TSS_threshold24))
		                    && (Mdata.hs_rate > HS_threshold_smallincrease))
		                 || ((Mdata.tss_a12h != IOUtils::nodata && Mdata.tss_a12h < IOUtils::C_TO_K(TSS_threshold12_smallHSincrease))
		                    && (Mdata.hs_rate > HS_threshold_smallincrease))
		                 || ((Mdata.tss_a12h != IOUtils::nodata && Mdata.tss_a12h < IOUtils::C_TO_K(TSS_threshold12_largeHSincrease))
		                    && (Mdata.hs_rate != IOUtils::nodata && Mdata.hs_rate > HS_threshold_largeincrease))
		                 )
		               )
		            || (Mdata.hs_rate != IOUtils::nodata && Mdata.hs_rate > HS_threshold_verylargeincrease)
		);
		// Additional check if snowed_in is still false with DETECT_GRASS due to missing TSS data
		if (!snowed_in && detect_grass && (Mdata.hs_rate == IOUtils::nodata || Mdata.tss_a24h == IOUtils::nodata || Mdata.tss_a12h == IOUtils::nodata)) {
			prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "DETECT_GRASS is TRUE, but insufficient HS and/or TSS data is available for the algorithm. Snowfall may be missed.");
		}
	}
	if (variant == "SEAICE" && nOldE == 0) {
		// Ignore snow fall on open ocean
		snowed_in = false;
	}

	// Go ahead if there is a snow fall AND the ground is or can be snowed in.
	if (snow_fall && snowed_in) {
		// Now check if we have some canopy left below the snow. Then we reduce the canopy with the new snow height.
		// This is to simulate the gradual sinking of the canopy under the weight of snow.
		// We also adjust Xdata.mH to have it reflect deposited snow but not the canopy.
		// This can only be done when SNOWPACK is snow height driven and there is a canopy.
		if ((enforce_measured_snow_heights)
			   && (Xdata.Cdata.height > 0.)
			   && ((Xdata.Cdata.height < ThresholdSmallCanopy) || (useCanopyModel == false))
			   && (Mdata.hs != mio::IOUtils::nodata)
			   && (Xdata.mH != Constants::undefined)
			   && (Xdata.meta.getSlopeAngle() < Constants::min_slope_angle)) {
			/* The third clause above limits the issue to small canopies only, to prevent problems
			 *   with Alpine3D simulations in forests. This prerequisite is only checked for when useCanopyModel
			 *    is true. If useCanopyModel is false, we can safely assume all snow to fall on top of canopy.
			 * The fourth clause is an important one. When hs is not available, the old Xdata.mH is kept, which
			 *    has already been adjusted in the previous time step, so then skip this part.
			 * The fifth clause makes sure only flat field is treated this way, and not the slopes.
			 *
			 *  Now reduce the Canopy height with the additional snow height. This makes the Canopy work
			 *   like a spring. When increase in Xdata.mh is 3 cm and the canopy height is 10 cm, the snow pack is
			 *   assumed to be 6cm in thickness. When the enforced snow depth Xdata.mH equals
			 *   the canopy height, the canopy is reduced to 0, and everything measured is assumed to be snow.
			 * To do this, check if there is an increase AND check if a new snow element will be created!
			 *   If you don't do this, the canopy will be reduced for small increases that do not produce a snow layer.
			 * Then, in the next time step, the canopy height is reduced even more, even without increase in snow
			 *   depth.
			 */
			if (Xdata.cH < (Xdata.mH - (Xdata.Cdata.height - (Xdata.mH - (Xdata.cH + Xdata.Cdata.height))))) {
				Xdata.Cdata.height -= (Xdata.mH - (Xdata.cH + Xdata.Cdata.height));
				// The above if-statement looks awkward, but it is just
				//   (Xdata.cH < (Xdata.mH - (height_new_elem * cos_sl))
				//   combined with the new value for Xdata.mH, given the change in Xdata.Cdata.height.
			} else {
				// If the increase is not large enough to start build up a snowpack yet,
				//   assign Xdata.mH to Canopy height (as if snow_fall and snowed_in would have been false).
				if (Xdata.getNumberOfNodes() == Xdata.SoilNode+1) {
					Xdata.Cdata.height = Xdata.mH - Xdata.Ground;
				}
			}
			if (Xdata.Cdata.height < 0.)    // Make sure canopy height doesn't get negative
				Xdata.Cdata.height = 0.;
			Xdata.mH -= Xdata.Cdata.height; // Adjust Xdata.mH to represent the "true" enforced snow depth
			if (Xdata.mH < Xdata.Ground)    //   and make sure it doesn't get negative
				Xdata.mH = Xdata.Ground;
			delta_cH = Xdata.mH - Xdata.cH + ( (Xdata.findMarkedReferenceLayer() == Constants::undefined) ? (0.) : (Xdata.findMarkedReferenceLayer() - Xdata.Ground) );
		}

		// Now determine whether the increase in snow depth is large enough.
		// NOTE On virtual slopes use new snow depth and density from either flat field or luv slope
		if ((delta_cH >= height_new_elem * cos_sl) || (Xdata.hn > 0.) || add_element || (force_add_snowfall && delta_cH > Constants::eps)) {
			cumu_precip = 0.0; // we use the mass through delta_cH
			//double hn = 0.; //new snow amount

			if (!alpine3d && Xdata.hn > 0. && Xdata.meta.getSlopeAngle() > Constants::min_slope_angle) {
				hn = Xdata.hn;
				rho_hn = Xdata.rho_hn;
			} else { // in case of flat field or PERP_TO_SLOPE
				hn = delta_cH;
				// Store new snow depth and density
				Xdata.hn = hn;
				Xdata.rho_hn = rho_hn;
			}
			if (hn > Snowpack::snowfall_warning)
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date,
				          "Large snowfall! hn=%.3f cm (azi=%.0f, slope=%.0f)",
				            M_TO_CM(hn), Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());

			size_t nAddE = (size_t)(hn / (height_new_elem*cos_sl));

			if (nAddE < 1) {
				// Always add snow on virtual slope (as there is no storage variable available) and some other cases
				if ((!alpine3d && Xdata.meta.getSlopeAngle() > Constants::min_slope_angle)
				                      || add_element || (force_add_snowfall && delta_cH > Constants::eps)) { //no virtual slopes in Alpine3D
					nAddE = 1;
				} else {
					Xdata.hn = 0.;
					Xdata.rho_hn = Constants::undefined;
					return;
				}
			}

			// Check whether surface hoar could be buried
			size_t nHoarE;
			double hoar = Xdata.Ndata[nOldN-1].hoar;
			if (nOldE > 0 && Xdata.Edata[nOldE-1].theta[SOIL] < Constants::eps2) {
				// W.E. of surface hoar must be larger than a threshold to be buried
				if (hoar > 1.5*MM_TO_M(hoar_min_size_buried)*hoar_density_surf) {
					nHoarE = 1;
				} else if (!(change_bc && meas_tss) /*Flux BC (NEUMANN) typically produce less SH*/
				               && (hoar > MM_TO_M(hoar_min_size_buried)*hoar_density_surf)) {
					nHoarE = 1;
				} else {
					nHoarE = 0;
				}
			} else { // Add surface hoar on ground to first new snow element(s)
				nHoarE = 0;
				hn += hoar/rho_hn;
				Xdata.Ndata[nOldN-1].hoar = 0.;
			}

			Xdata.Albedo = Constants::max_albedo;

			const size_t nNewN = nOldN + nAddE + nHoarE;
			const size_t nNewE = nOldE + nAddE + nHoarE;
			Xdata.resize(nNewE);
			vector<NodeData>& NDS = Xdata.Ndata;
			vector<ElementData>& EMS = Xdata.Edata;
			Sdata.mass[SurfaceFluxes::MS_SNOW_DHS] += hn;

			// Create hoar layer
			if (nHoarE > 0) {
				// Since mass of hoar was already added to element below, substract....
				// Make sure you don't try to extract more than is there
				hoar = std::max(0.,std::min(EMS[nOldE-1].M - 0.1,hoar));
				const double L0 = EMS[nOldE-1].L;
				const double dL = -hoar/(EMS[nOldE-1].Rho);
				EMS[nOldE-1].L0 = EMS[nOldE-1].L = L0 + dL;

				EMS[nOldE-1].E = EMS[nOldE-1].Eps = EMS[nOldE-1].dEps = EMS[nOldE-1].Eps_e = EMS[nOldE-1].Eps_v = EMS[nOldE-1].S = 0.0;
				const double Theta0 = EMS[nOldE-1].theta[ICE];
				EMS[nOldE-1].theta[ICE] *= L0/EMS[nOldE-1].L;
				EMS[nOldE-1].theta[ICE] += -hoar/(Constants::density_ice*EMS[nOldE-1].L);
				EMS[nOldE-1].theta[ICE] = std::max(EMS[nOldE-1].theta[ICE],0.);
				EMS[nOldE-1].theta_i_reservoir = 0.0;
				EMS[nOldE-1].theta_i_reservoir_cumul = 0.0;
				EMS[nOldE-1].theta[WATER] *= L0/EMS[nOldE-1].L;
				EMS[nOldE-1].theta[WATER_PREF] *= L0/EMS[nOldE-1].L;
				for (unsigned int ii = 0; ii < Xdata.number_of_solutes; ii++)
					EMS[nOldE-1].conc[ICE][ii] *= L0*Theta0/(EMS[nOldE-1].theta[ICE]*EMS[nOldE-1].L);
				EMS[nOldE-1].M -= hoar;
				assert(EMS[nOldE-1].M>=0.); //the mass must be positive
				EMS[nOldE-1].theta[AIR] = std::max(0., 1.0 - EMS[nOldE-1].theta[WATER] - EMS[nOldE-1].theta[WATER_PREF]
				                                - EMS[nOldE-1].theta[ICE] - EMS[nOldE-1].theta[SOIL]);
				EMS[nOldE-1].updDensity();
				assert(EMS[nOldE-1].Rho>=0. || EMS[nOldE-1].Rho==IOUtils::nodata); //we want positive density
				// Take care of old surface node
				NDS[nOldN-1].z += dL + NDS[nOldN-1].u;
				NDS[nOldN-1].u = 0.0;
				NDS[nOldN-1].hoar = 0.0;
				// Now fill nodal data for upper hoar node
				NDS[nOldN].T = t_surf;              // The temperature of the new node
				double p_vapor = Atmosphere::vaporSaturationPressure(NDS[nOldN].T);
				NDS[nOldN].rhov = Atmosphere::waterVaporDensity(NDS[nOldN].T, p_vapor);
				// The new nodal position;
				NDS[nOldN].z = NDS[nOldN-1].z + NDS[nOldN-1].u + hoar/hoar_density_buried;
				NDS[nOldN].u = 0.0;                 // Initial displacement is 0
				NDS[nOldN].hoar = hoar / hoar_density_buried;         // Surface hoar initial size
				NDS[nOldN].udot = 0.0;               // Settlement rate is also 0
				NDS[nOldN].f = 0.0;                 // Unbalanced forces is 0
				NDS[nOldN].S_n = IOUtils::nodata;
				NDS[nOldN].S_s = IOUtils::nodata;
			} else { // Make sure top node surface hoar mass is removed
				NDS[nOldN-1].hoar = 0.0;
			}

			// Fill the nodal data
			if (!useSoilLayers && (nOldN-1 == Xdata.SoilNode)) // New snow on bare ground w/o soil
				NDS[nOldN-1].T = 0.5*(t_surf + Mdata.ta);
			const double Ln = (hn / (double)nAddE);               // New snow element length
			double z0 = NDS[nOldN-1+nHoarE].z + NDS[nOldN-1+nHoarE].u + Ln; // Position of lowest new node
			for (size_t n = nOldN+nHoarE; n < nNewN; n++) { //loop over the nodes
				NDS[n].T = t_surf;                  // Temperature of the new node
				NDS[n].z = z0;                      // New nodal position
				NDS[n].u = 0.0;                     // Initial displacement is 0
				NDS[n].hoar = 0.0;                  // The new snow surface hoar is set to zero
				NDS[n].udot = 0.0;                  // Settlement rate is also 0
				NDS[n].f = 0.0;                     // Unbalanced forces are 0
				NDS[n].S_n = IOUtils::nodata;
				NDS[n].S_s = IOUtils::nodata;
				z0 += Ln;
			}

			// Fill the element data
			for (size_t e = nOldE; e < nNewE; e++) { //loop over the elements
				const bool is_surface_hoar = (nHoarE && (e == nOldE));
				const double length = (NDS[e+1].z + NDS[e+1].u) - (NDS[e].z + NDS[e].u);
				const double density = (is_surface_hoar)? hoar_density_buried : rho_hn;
				fillNewSnowElement(Mdata, length, density, is_surface_hoar, Xdata.number_of_solutes, EMS[e]);
				// To satisfy the energy balance, we should trigger an explicit treatment of the top boundary condition of the energy equation
				// when new snow falls on top of wet snow or melting soil. This can be done by putting a tiny amount of liquid water in the new snow layers.
				// Note that we use the same branching condition as in the function Snowpack::neumannBoundaryConditions(...)
				const double theta_r = ((watertransportmodel_snow=="RICHARDSEQUATION" && Xdata.getNumberOfElements()>Xdata.SoilNode) || (watertransportmodel_soil=="RICHARDSEQUATION" && Xdata.getNumberOfElements()==Xdata.SoilNode)) ? (PhaseChange::RE_theta_threshold) : (PhaseChange::theta_r);
				if(nOldE > 0 && EMS[nOldE-1].theta[WATER] > theta_r + Constants::eps && EMS[nOldE-1].theta[ICE] > Constants::eps) {
					EMS[e].theta[WATER]+=(2.*Constants::eps);
					EMS[e].theta[ICE]-=(2.*Constants::eps)*(Constants::density_water/Constants::density_ice);
					EMS[e].theta[AIR]+=((Constants::density_water/Constants::density_ice)-1.)*(2.*Constants::eps);
				}
				if (variant == "SEAICE" ) {
					Xdata.Seaice->calculateMeltingTemperature(EMS[e]);
				} else {
					EMS[e].meltfreeze_tk = Constants::meltfreeze_tk;
				}

				Xdata.ColdContent += EMS[e].coldContent(); //update cold content
			}   // End elements

			// Finally, update the computed snowpack height
			Xdata.cH = NDS[nNewN-1].z + NDS[nNewN-1].u;
			Xdata.ErosionLevel = nNewE-1;
		}
	} else { // No snowfall
		if (detect_grass && ((Xdata.Cdata.height < ThresholdSmallCanopy) || (useCanopyModel == false))) {
			// Set canopy height to enforced snow depth if there is no snowpack yet
			//   but only for small canopies, to prevent problems with Alpine3D simulations in forests.
			// This prerequisite is only checked for when useCanopyModel is true.
			// If useCanopyModel is false, we can safely assume all snow to fall on top of canopy.
			if ((Xdata.getNumberOfNodes() == Xdata.SoilNode+1) && (Xdata.mH != Constants::undefined)) {
				Xdata.Cdata.height = Xdata.mH - Xdata.Ground;
				// Because there is no snow cover, enforced snow depth is effectively equal to Xdata.Ground.
				Xdata.mH = Xdata.Ground;
			} else {
				if(Mdata.hs != mio::IOUtils::nodata) {
					// If we have a snowpack, but didn't match the criteria for snow fall, make sure Xdata.mH
					// only represents the "true" snow height, to stay consistent and for use in other parts of SNOWPACK.
					Xdata.mH -= Xdata.Cdata.height;
					if (Xdata.mH < Xdata.Ground)
						Xdata.mH = Xdata.Ground;
				}
			}
		}
	}
}

/**
 * @brief Add snow layers that originate from wind-transported snow being deposited, using the event-driven deposition scheme.
 * @param Mdata Meteorological data (pass by value, since we modify it)
 * @param Xdata Snow cover data
 * @param redeposit_mass cumulated amount of snow deposition (kg m-2)
 */
void Snowpack::RedepositSnow(CurrentMeteo Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata, double redeposit_mass)
{
	// Backup settings we are going to override:
	const bool tmp_force_add_snowfall = force_add_snowfall;
	const std::string tmp_hn_density = hn_density;
	const std::string tmp_variant = variant;
	const bool tmp_enforce_measured_snow_heights = enforce_measured_snow_heights;
	const double tmp_Xdata_hn = Xdata.hn;
	const double tmp_Xdata_rho_hn = Xdata.rho_hn;
	const mio::Date tmp_MdataDate = Mdata.date;
	// Deposition mode settings:
	double tmp_psum = redeposit_mass;
	force_add_snowfall = true;
	hn_density = "EVENT";
	variant = "POLAR";		// Ensure that the ANTARCTICA wind speed limits are *not* used.
	enforce_measured_snow_heights = false;
	Mdata.psum = redeposit_mass; Mdata.psum_ph = 0.;
	// The EVENT scheme uses vw_avg and rh_avg in the calculations. In the REDEPOSIT scheme, we force the use of instantaneous values for wind speed and relative humidity:
	Mdata.vw_avg = Mdata.vw;
	Mdata.rh_avg = Mdata.rh;
	Xdata.hn = 0.;
	if (Xdata.ErosionAge != Constants::undefined && redeposit_keep_age) {
		mio::Date EnforcedDepositionDate(Xdata.ErosionAge, Mdata.date.getTimeZone());
		Mdata.date = EnforcedDepositionDate;
	}
	// Add eroded snow:
	compSnowFall(Mdata, Xdata, tmp_psum, Sdata);
	// Set back original settings:
	force_add_snowfall = tmp_force_add_snowfall;
	hn_density = tmp_hn_density;
	variant = tmp_variant;
	enforce_measured_snow_heights = tmp_enforce_measured_snow_heights;
	Mdata.date = tmp_MdataDate;
	// Calculate new snow density (weighted average) and total snowfall (snowfall + redeposited snow)
	Xdata.hn_redeposit = Xdata.hn;
	Xdata.rho_hn_redeposit = Xdata.rho_hn;
	if ((tmp_Xdata_hn + Xdata.hn) > 0.) {
		Xdata.rho_hn = ((tmp_Xdata_hn * tmp_Xdata_rho_hn) + (Xdata.hn * Xdata.rho_hn)) / (tmp_Xdata_hn + Xdata.hn);
	}
	Xdata.hn += tmp_Xdata_hn;
}

/**
 * @brief The near future (s. below) has arrived on Wednesday Feb. 6, when it was finally snowing
 * in Davos and Sergey, Michael and Perry were working furiously on SNOWPACK again. Michael
 * prepared the coupling of the model to the energy balance model of Olivia and his own snow
 * drift model. He found it therefore necessary to move the time loop into the Main.c frame.
 * He also cleaned the data structure handling, i.e. made sure that all variables are handed
 * to the subroutines and that there is no longer a global referencing. In the meantime, Perry
 * and Sergey prepared the introduction of a vapor transport equation. While Sergey was sad
 * that his family had to go back, Perry was exchanging a lot of wet kisses with his new women
 * probably causing some of his confusion about vapor and heat transport.
 * Driving routine for the 1d-snowpack model.  The main program will probably be replaced
 * in future -- or, at least modified.  For now, it is responsible for several tasks:
 * -# In the next step the program enters the time-integration, or, time-stepping loop.
 *    Moreover, we find the solution at TimeN = Time0 + sn_dt.  For now, we will assume
 *    that sn_dt = 15min and Time0 = 0.0.  These restrictions can be later relaxed.
 *    The program must then "check" three modules before finding the temperature distribution
 *    and creep deformations, these are:
 *    -#   Determine the METEO parameters at time TimeN. (The METEO parameters include
 *             air temperature, relative humidity, short wave radiation, incoming longwave
 *             radiation, etc. as well as the ground temperature)
 *             >>>>>> METEO MODULE <<<<<<<
 *             Change: ML, 06.02.02: Mdata needs to be handed to the function
 *    -#   Determine if a SNOWFALL event as occcured. (This implies that the dynamically
 *             allocated data structures must be reallocated and the new material defined.
 *             The mesh connectivies (sn_NewMesh=TRUE) must be rebuilt for the numerical
 *             solution and the energy exchanges at the top surface initialized.)
 *             >>>>>> SNOWFALL MODULE <<<<<<<
 *    -#   Determine if SNOWDRIFT is occuring. Note that SNOWFALL distributes mass
 *             equally at ALL EXPOSITIONS.  SNOWDRIFT can occur in conjunction WITH SNOWFALL
 *             or WITHOUT SNOWFALL.  Mass is redistributed between the expositions. (Again,
 *             the internal data structures must be reinitialized.)
 *             >>>>>> SNOWDRIFT MODULE <<<<<<<
 * -# The phase change module can generate water; this water is moved downwards by the
 *    WATER TRANSPORT module.  At present, a simple, mass conserving scheme is employed to move
 *    the water; when the water content is greater than the residual water content then the excess
 *    water is moved to the next element.  If the element reaches the bottom of the snowpack
 *    water is removed -- this is called MELTWATER RUNOFF.  Also, elements which do not contain
 *    any ice are removed from the finite element mesh.
 * -# In the next step the TEMPERATURE distribution within the snowpack at each exposition
 *    is found.  If SNOWFALL and SNOWDRIFT has occured then the finite element matrices must
 *    be rebuilt.  >>>>> TEMPERATURE MODULE  <<<<<<<<
 * -# If the computed temperatures are above zero degrees then PHASECHANGES are occuring.
 *    This means that the volumetric contents of the elements must be updated. >>>>>> PHASE CHANGE MODULE <<<<<
 * -# In the next step the CREEPING deformations and 1d state of stress within the snowpack
 * at each exposition are found.  If SNOWFALL and SNOWDRIFT has occured then the finite
 * element matrices (connectivities) must be rebuilt.  >>>>CREEP MODULE<<<<<
 * @param Mdata The Meteorological forcing is now passed by copy so changes won't propagate to the caller
 * @param Xdata
 * @param cumu_precip Variable to store amount of precipitation (kg m-2)
 * @param Bdata
 * @param Sdata
 */
void Snowpack::runSnowpackModel(CurrentMeteo& Mdata, SnowStation& Xdata, double& cumu_precip,
                                BoundCond& Bdata, SurfaceFluxes& Sdata)
{
	// HACK -> couldn't the following objects be created once in init ?? (with only a reset method ??)
	WaterTransport watertransport(cfg);
	VapourTransport vapourtransport(cfg);
	Metamorphism metamorphism(cfg);
	SnowDrift snowdrift(cfg);
	PhaseChange phasechange(cfg);
	if (Xdata.Seaice != NULL) Xdata.Seaice->ConfigSeaIce(cfg);

	// ADJUST_HEIGHT_OF_METEO_VALUE is checked at each call to allow different
	// cfg values for different pixels in Alpine3D
	cfg.getValue("ADJUST_HEIGHT_OF_METEO_VALUES", "SnowpackAdvanced", adjust_height_of_meteo_values);


	try {
		//since precipitation phase is a little less intuitive than other, measured parameters, make sure it is provided
		if (Mdata.psum_ph==IOUtils::nodata)
			throw NoDataException("Missing precipitation phase", AT);

		// Set and adjust boundary conditions
		surfaceCode = NEUMANN_BC;
		double meltfreeze_tk = (Xdata.getNumberOfElements()>0)? Xdata.Edata[Xdata.getNumberOfElements()-1].meltfreeze_tk : Constants::meltfreeze_tk;
		t_surf = std::min(meltfreeze_tk, Xdata.Ndata[Xdata.getNumberOfNodes()-1].T);
		if (forcing == "MASSBAL") {
			if(!(Mdata.surf_melt>0.)) {
				// always use Dirichlet boundary conditions with massbal forcing, when no melt is going on
				surfaceCode = DIRICHLET_BC;
			} else {
				// in case melting is going on, top node is at melting conditions and we prescribe prescribed melt as Neumann boundary condition
				t_surf = Xdata.Ndata[Xdata.getNumberOfNodes()-1].T = (Xdata.getNumberOfElements()>0) ? (Xdata.Edata[Xdata.getNumberOfElements()-1].meltfreeze_tk) : (Constants::meltfreeze_tk);
			}
		} else if ((change_bc && meas_tss) && ((Mdata.tss < IOUtils::C_TO_K(thresh_change_bc)) && Mdata.tss != IOUtils::nodata)) {
			surfaceCode = DIRICHLET_BC;
		}
		if (surfaceCode == DIRICHLET_BC) {
			t_surf = Mdata.tss;
			Xdata.Ndata[Xdata.getNumberOfNodes()-1].T = t_surf;
		}

		// If there is DEPOSITING of snow: // This is externally forced snow drift, e.g. from Alpine3D. In snowpack standalone, without snowdrift in .smet file, Mdata.snowdrift (redeposit_mass) is 0.
		if (Mdata.snowdrift > 0. && snow_erosion == "REDEPOSIT") {
			RedepositSnow(Mdata, Xdata, Sdata, Mdata.snowdrift);
			Mdata.snowdrift = 0.;
		}

		// If it is SNOWING, find out how much, prepare for new FEM data. If raining, cumu_precip is set back to 0
		const double tmp_Xdata_hn = Xdata.hn;		// Store initial hn
		const double tmp_Xdata_rho_hn = Xdata.rho_hn;	// Store initial rho_hn
		if (alpine3d || !(Xdata.meta.getSlopeAngle() > Constants::min_slope_angle)) Xdata.hn = 0.;	// Virtual slopes receive precipitation via Xdata.hn
		compSnowFall(Mdata, Xdata, cumu_precip, Sdata);
		if (Xdata.hn > 0.) {
			// If new snow was added, recalculate rho_hn as a weighted average
			Xdata.rho_hn = ((tmp_Xdata_hn * tmp_Xdata_rho_hn) + (Xdata.hn * Xdata.rho_hn)) / (tmp_Xdata_hn + Xdata.hn);
		} else {
			// Otherwise keep previous rho_hn
			Xdata.rho_hn = tmp_Xdata_rho_hn;
		}
		Xdata.hn += tmp_Xdata_hn;

		// Check to see if snow is DRIFTING, compute a simple snowdrift index and erode layers if
		// neccessary. Note that also the very important friction velocity is computed in this
		// routine and later used to compute the Meteo Heat Fluxes
		if (forcing=="ATMOS") {
			double tmp = 0.;
			if (alpine3d && cumu_precip < 0.) { //HACK: we need to set to 0 the external drift
				// With alpine3d, negative cumu_precip becomes forced_massErode
				tmp = cumu_precip;
				cumu_precip = 0.;
			}
			snowdrift.compSnowDrift(Mdata, Xdata, Sdata, tmp);
			if (Xdata.ErosionMass > 0. && snow_erosion == "REDEPOSIT" && !alpine3d) RedepositSnow(Mdata, Xdata, Sdata, Xdata.ErosionMass);
		} else { // MASSBAL forcing
			snowdrift.compSnowDrift(Mdata, Xdata, Sdata, Mdata.snowdrift); //  Mdata.snowdrift is always <= 0. (positive values are in Mdata.psum)
		}

		if (Xdata.Seaice != NULL) {
			// Reinitialize and compute the initial meteo heat fluxes
			Bdata.reset();
			updateBoundHeatFluxes(Bdata, Xdata, Mdata);
			// Run sea ice module
			Xdata.Seaice->runSeaIceModule(Xdata, Mdata, Bdata, sn_dt, Sdata);
			// Remesh when necessary
			Xdata.splitElements(2. * comb_thresh_l, comb_thresh_l);
		}

		const double sn_dt_bcu = sn_dt;			// Store original SNOWPACK time step
		const double psum_bcu = Mdata.psum;		// Store original psum value
		const double psum_ph_bcu = Mdata.psum_ph;	// Store original psum_ph value
		const double sublim_bcu = Mdata.sublim;		// Store original sublim value
		const double surf_melt_bcu = Mdata.surf_melt;	// Store original sublim value
		int ii = 0;				// Counter for sub-timesteps to match one SNOWPACK time step
		bool LastTimeStep = false;		// Flag to indicate if it is the last sub-time step
		double p_dt = 0.;			// Cumulative progress of time steps
		if ((Mdata.psi_s >= 0. || t_surf > Mdata.ta) && atm_stability_model != Meteo::NEUTRAL && allow_adaptive_timestepping == true && sn_dt > 60.) {
			// To reduce oscillations in TSS, reduce the time step prematurely when atmospheric stability is unstable.
			if (Mdata.psum != mio::IOUtils::nodata) Mdata.psum /= sn_dt;					// psum is precipitation per time step, so first express it as rate with the old time step (necessary for rain only)...
			if (forcing=="MASSBAL" && Mdata.sublim != mio::IOUtils::nodata)		Mdata.sublim /= sn_dt;		// scale the mass balance components like the precipiation
			if (forcing=="MASSBAL" && Mdata.surf_melt != mio::IOUtils::nodata)	Mdata.surf_melt /= sn_dt;	// scale the mass balance components like the precipiation

			sn_dt = 60.;

			if (Mdata.psum != mio::IOUtils::nodata) Mdata.psum *= sn_dt;					// ... then express psum again as precipitation per time step with the new time step
			if (forcing=="MASSBAL" && Mdata.sublim != mio::IOUtils::nodata)		Mdata.sublim *= sn_dt;		// scale the mass balance components like the precipiation
			if (forcing=="MASSBAL" && Mdata.surf_melt != mio::IOUtils::nodata)	Mdata.surf_melt *= sn_dt;	// scale the mass balance components like the precipiation
		}

		Meteo meteo(cfg);
		do {
			// After the first sub-time step, update Meteo object to reflect on the new stability state
			if (ii >= 1){
				// ADJUST_HEIGHT_OF_WIND_VALUE is checked at each call to allow different
				// cfg values for different pixels in Alpine3D
				cfg.getValue("ADJUST_HEIGHT_OF_WIND_VALUE", "SnowpackAdvanced", adjust_height_of_wind_value);
				meteo.compMeteo(Mdata, Xdata, false, adjust_height_of_wind_value);
			}
			// Reinitialize and compute the initial meteo heat fluxes
			Bdata.reset();
			updateBoundHeatFluxes(Bdata, Xdata, Mdata);

			// set the snow albedo
			Xdata.pAlbedo = getParameterizedAlbedo(Xdata, Mdata);
			Xdata.Albedo = getModelAlbedo(Xdata, Mdata); //either parametrized or measured

			// Compute the temperature profile in the snowpack and soil, if present
			for (size_t e = 0; e < Xdata.getNumberOfElements(); e++) Xdata.Edata[e].Qph_up = Xdata.Edata[e].Qph_down = 0.;
			if (compTemperatureProfile(Mdata, Xdata, Bdata, (sn_dt < min_allowed_sn_dt))) {
				// Entered after convergence
				ii++;						// Update time step counter
				p_dt += sn_dt;					// Update progress variable
				if (p_dt > sn_dt_bcu-Constants::eps) {		// Check if it is the last sub-time step
					LastTimeStep = true;
				}

				// Good HACK (according to Charles, qui persiste et signe;-)... like a good hunter and a bad one...
				// If you switched from DIRICHLET to NEUMANN boundary conditions, correct
				//   for a possibly erroneous surface energy balance. The latter can be due e.g. to a lack
				//   of data on nebulosity leading to a clear sky assumption for incoming long wave.
				if ((change_bc && meas_tss) && (surfaceCode == NEUMANN_BC) && forcing != "MASSBAL"
						&& (Xdata.Ndata[Xdata.getNumberOfNodes()-1].T < mio::IOUtils::C_TO_K(thresh_change_bc))) {
					surfaceCode = DIRICHLET_BC;
					meltfreeze_tk = (Xdata.getNumberOfElements()>0)? Xdata.Edata[Xdata.getNumberOfElements()-1].meltfreeze_tk : Constants::meltfreeze_tk;
					Xdata.Ndata[Xdata.getNumberOfNodes()-1].T = std::min(Mdata.tss, meltfreeze_tk); /*C_TO_K(thresh_change_bc/2.);*/
					// update the snow albedo
					Xdata.pAlbedo = getParameterizedAlbedo(Xdata, Mdata);
					Xdata.Albedo = getModelAlbedo(Xdata, Mdata); //either parametrized or measured
					for (size_t e = 0; e < Xdata.getNumberOfElements(); e++) Xdata.Edata[e].Qph_up = Xdata.Edata[e].Qph_down = 0.;
					compTemperatureProfile(Mdata, Xdata, Bdata, true);	// Now, throw on non-convergence
				}
				if (LastTimeStep) Sdata.compSnowSoilHeatFlux(Xdata);

				// Inialize PhaseChange at the first sub-time step
				if (ii == 1) phasechange.initialize(Xdata);

				// See if any SUBSURFACE phase changes are occuring due to updated temperature profile
				if(!coupled_phase_changes) {
					if (!alpine3d)
						phasechange.compPhaseChange(Xdata, Mdata.date, true, ((Mdata.surf_melt != IOUtils::nodata) ? (Mdata.surf_melt) : (0.)));
					else
						phasechange.compPhaseChange(Xdata, Mdata.date, false, ((Mdata.surf_melt != IOUtils::nodata) ? (Mdata.surf_melt) : (0.)));
				} else {
					const double theta_r = ((watertransportmodel_snow=="RICHARDSEQUATION" && Xdata.getNumberOfElements()>Xdata.SoilNode) || (watertransportmodel_soil=="RICHARDSEQUATION" && Xdata.getNumberOfElements()==Xdata.SoilNode)) ? (PhaseChange::RE_theta_r) : (PhaseChange::theta_r);
					const double max_ice = ReSolver1d::max_theta_ice;
					for (size_t e = 0; e < Xdata.getNumberOfElements(); e++) {
						// Net ice contents change:
						double dth_i = 0.5 * (Xdata.Edata[e].Qph_up + Xdata.Edata[e].Qph_down) / ((Constants::density_ice * Constants::lh_fusion) / sn_dt);
						// Limit to all ice melts:
						dth_i = (dth_i<0.)?(std::max(-Xdata.Edata[e].theta[ICE], dth_i)):(dth_i);
						// Limit to all liquid water freezes:
						dth_i = (dth_i>0.)?(std::min(std::max(0., std::min(max_ice - Xdata.Edata[e].theta[ICE], (Xdata.Edata[e].theta[WATER] - theta_r) * (Constants::density_water / Constants::density_ice))), dth_i)):(dth_i);
						// Apply phase change:
						Xdata.Edata[e].dth_w -= dth_i * Constants::density_ice / Constants::density_water;
						Xdata.Edata[e].Qmf += (dth_i * Constants::density_ice * Constants::lh_fusion) / sn_dt_bcu; // (W m-3)
						Xdata.Edata[e].theta[ICE] += dth_i;
						Xdata.Edata[e].theta[WATER] -= dth_i*Constants::density_ice/Constants::density_water;
						Xdata.Edata[e].theta[AIR] = 1. - Xdata.Edata[e].theta[WATER] - Xdata.Edata[e].theta[WATER_PREF] - Xdata.Edata[e].theta[ICE] - Xdata.Edata[e].theta[SOIL];
						Xdata.Edata[e].updDensity();
						Xdata.Edata[e].heatCapacity();
						Xdata.Edata[e].Qph_up = Xdata.Edata[e].Qph_down = 0.;

						if (Xdata.Seaice != NULL) {
							// Adjust melting/freezing point assuming thermal quilibrium in the brine pockets
							const double BrineSal_new = (Xdata.Edata[e].theta[WATER] == 0.) ? (0.) : (Xdata.Edata[e].salinity / Xdata.Edata[e].theta[WATER]);
							Xdata.Edata[e].meltfreeze_tk = Xdata.Seaice->calculateMeltingTemperature(BrineSal_new);
						}
					}
				}

				// Compute the final heat fluxes at the last sub-time step
				if (LastTimeStep) Sdata.ql += Bdata.ql; // Bad;-) HACK, needed because latent heat ql is not (yet)
								        // linearized w/ respect to Tss and thus remains unchanged
								        // throughout the temperature iterations!!!
				updateBoundHeatFluxes(Bdata, Xdata, Mdata);

				// Make sure the sub-time steps match exactly one SNOWPACK time step
				if (p_dt + sn_dt > sn_dt_bcu) {
					sn_dt = sn_dt_bcu - p_dt;
				}
			} else {
				// Entered after non-convergence
				if (sn_dt == sn_dt_bcu) std::cout << "[i] [" << Mdata.date.toString(Date::ISO) << "] : using adaptive timestepping\n"; // First time warning

				if (Mdata.psum != mio::IOUtils::nodata) Mdata.psum /= sn_dt;					// psum is precipitation per time step, so first express it as rate with the old time step (necessary for rain only)...
				if (forcing=="MASSBAL" && Mdata.sublim != mio::IOUtils::nodata)		Mdata.sublim /= sn_dt;		// scale the mass balance components like the precipiation
				if (forcing=="MASSBAL" && Mdata.surf_melt != mio::IOUtils::nodata)	Mdata.surf_melt /= sn_dt;	// scale the mass balance components like the precipiation

				sn_dt /= 2.;							// No convergence, half the time step
				if (Mdata.psum != mio::IOUtils::nodata) Mdata.psum *= sn_dt;					// ... then express psum again as precipitation per time step with the new time step
				if (forcing=="MASSBAL" && Mdata.sublim != mio::IOUtils::nodata)		Mdata.sublim *= sn_dt;		// scale the mass balance components like the precipiation
				if (forcing=="MASSBAL" && Mdata.surf_melt != mio::IOUtils::nodata)	Mdata.surf_melt *= sn_dt;	// scale the mass balance components like the precipiation

				std::cout << "                            --> time step temporarily reduced to: " << sn_dt << "\n";
			}
		}
		while (LastTimeStep == false);

		sn_dt = sn_dt_bcu;			// Set back SNOWPACK time step to orginal value
		Mdata.psum = psum_bcu;			// Set back psum to original value
		Mdata.psum_ph = psum_ph_bcu;		// Set back psum_ph to original value
		Mdata.sublim = sublim_bcu;		// Set back sublim to original value
		Mdata.surf_melt = surf_melt_bcu;	// Set back surf_melt to original value

		// Compute change of internal energy during last time step (J m-2)
		Xdata.compSnowpackInternalEnergyChange(sn_dt);
		Xdata.compSoilInternalEnergyChange(sn_dt);

		// The water transport routines must be placed here, otherwise the temperature
		// and creep solution routines will not pick up the new mesh boolean.
		double ql = Bdata.ql;	// Variable to keep track of how latent heat is used
		watertransport.compTransportMass(Mdata, Xdata, Sdata, ql);

		vapourtransport.compTransportMass(Mdata, ql, Xdata, Sdata);

		// See if any SUBSURFACE phase changes are occuring due to updated water content (infiltrating rain/melt water in cold snow layers)
		if(!coupled_phase_changes) {
			if(!alpine3d)
				phasechange.compPhaseChange(Xdata, Mdata.date, true);
			else
				phasechange.compPhaseChange(Xdata, Mdata.date, false);

			// Finalize PhaseChange
			phasechange.finalize(Sdata, Xdata, Mdata.date);
		}

		// Compute change of internal energy during last time step (J m-2)
		Xdata.compSnowpackInternalEnergyChange(sn_dt);
		Xdata.compSoilInternalEnergyChange(sn_dt);

		// Find the settlement of the snowpack.
		// HACK This routine was formerly placed here because the settlement solution MUST ALWAYS follow
		// computeSnowTemperatures where the vectors U, dU and dUU are allocated.
		compSnowCreep(Mdata, Xdata, Sdata);

	} catch(const exception&) {
		prn_msg(__FILE__, __LINE__, "err", Mdata.date, "Snowpack computation not completed");
		throw;
	}

	metamorphism.runMetamorphismModel(Mdata, Xdata);

	if (combine_elements) {
		// Check for combining elements
		Xdata.combineElements(SnowStation::number_top_elements, reduce_n_elements, 1, comb_thresh_l);
		// Check for splitting elements
		if (reduce_n_elements > 0) {
			Xdata.splitElements(-1., comb_thresh_l);
		}
	}

	if (max_simulated_hs > 0.) {
		Xdata.CheckMaxSimHS(max_simulated_hs);
	}
}

void Snowpack::snowPreparation(const mio::Date& currentDate, SnowStation& Xdata) const
{
	if (techsnow.prepare(currentDate))
		techsnow.preparation(Xdata);
}
