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
#include <alpine3d/AlpineControl.h>
#include <alpine3d/AlpineMain.h>
#include <meteoio/MeteoIO.h>

#include <iomanip>

using namespace mio;
using namespace std;

/**
 * @brief Constructs and initialise the time steps loop.
 * This module calls all the necessary other modules (ie thoese that have been enabled) and runs through the time steps.
 * @param mysnowpack pointer to the initialized SNOWPACK Manager
 * @param mysnowdrift pointer to the initialized Snowdrift Manager
 * @param myeb pointer to the initialized radiation manager
 * @param cfg User configuration keys
 * @param in_dem DEM defining the simulation
 */
AlpineControl::AlpineControl(SnowpackInterface *mysnowpack, SnowDriftA3D *mysnowdrift, EnergyBalance *myeb,
                             const Config& cfg, const DEMObject& in_dem)
              : dem(in_dem), meteo(cfg, dem), snowpack(mysnowpack), snowdrift(mysnowdrift), eb(myeb), snow_days_between(0.), max_run_time(-1.), enable_simple_snow_drift(false), enable_snowdrift2d(false), nocompute(false), out_snow(true), correct_meteo_grids_HS(false), dataFromGrids(false)
{
	cfg.getValue("SNOW_WRITE", "Output", out_snow);
	if (out_snow) {
		cfg.getValue("SNOW_DAYS_BETWEEN", "Output", snow_days_between);
	}
	cfg.getValue("ADD_HS_TO_DEM_FOR_METEO", "input", correct_meteo_grids_HS,IOUtils::nothrow);
	cfg.getValue("DATA_FROM_GRIDS", "input", dataFromGrids,IOUtils::nothrow);

	//check if simple snow drift is enabled
	enable_simple_snow_drift = false;
	cfg.getValue("SIMPLE_SNOW_DRIFT", "Alpine3D", enable_simple_snow_drift, IOUtils::nothrow);
	enable_snowdrift2d = false;
	cfg.getValue("SNOWDRIFT2D", "Alpine3D", enable_snowdrift2d, IOUtils::nothrow);

	//check if maximum run time is specified
	cfg.getValue("MAX_RUN_TIME", "Alpine3D", max_run_time, IOUtils::nothrow);
}

void AlpineControl::Run(Date i_startdate, const unsigned int max_steps)
{// This function organizes the whole simulation:
// initializations, time loop, getting input data, computing each time step and finally writing outputs
	Date calcDate(i_startdate); //Date object initialized to julian 0.0 in local time zone with no DST
	const double timeStep = dt_main/86400.;
	Timer elapsed;
	std::vector<MeteoData> vecMeteo; // to transfer meteo information
	mio::Grid2DObject p, psum, psum_ph, vw, vw_drift, dw, rh, ta, tsg, ilwr, iswr_dir, iswr_diff;
	const bool isMaster = MPIControl::instance().master();

	if (isMaster) {
		cout << "\n**** Done initializing\n";
		cout << "**** Starting Calculation on date: " << calcDate.toString(Date::ISO) << " using Alpine3D version " << A3D_VERSION << "\n";
		if (nocompute)
			cout << "**** Performing dry run (--no-compute option)\n";
		cout << "\n";

		if (nocompute) {
			const Grid2DObject maskGlacier{snowpack->getGrid(SnGrids::GLACIER)};
			meteo.setGlacierMask(maskGlacier);
		}
	}

	//if the meteo data would need to be resampled, we try to fill the buffer with a date a little bit before
	if (snowdrift) meteo.setSkipWind(true); //do not fill grids if met3D
	meteo.prepare(i_startdate);

	elapsed.start();
	for (unsigned int t_ind=0; t_ind<max_steps; t_ind++) { //main computational loop
		const double elapsed_start = elapsed.getElapsed();
		const double est_completion = elapsed_start * ((double)max_steps/(double)(t_ind+1) - 1.);

		if (isMaster) {
			cout << "\nSimulation step " << t_ind+1 << "/" << max_steps << " at time step " << calcDate.toString(mio::Date::ISO_TZ) << "\n";
			cout << std::fixed << "Elapsed time: " << setprecision(1) << elapsed_start << " seconds\nEstimated completion in " << est_completion/3600. << " hours\n";
		}

		//get 1D and 2D meteo for the current time step
		// The grids ta, tsg, rh, psum, psum_ph, vw, vw_drift, dw, p, ilwr get populated here
		try {
			meteo.get(calcDate, vecMeteo);
			if(correct_meteo_grids_HS){
				meteo.setDEM(dem+snowpack->getGrid(SnGrids::HS));
			}
			meteo.get(calcDate, ta, tsg, rh, psum, psum_ph, vw, vw_drift, dw, p, ilwr, iswr_dir, iswr_diff);
		} catch (IOException&) {
			//saving state files before bailing out
			if (isMaster) {
				if (out_snow && t_ind>0 && snowpack)
					snowpack->writeOutputSNO(calcDate-1./24.); //output for last hour
				throw;
			}
		}

		if (t_ind < (max_steps-1)) {
			meteo.prepare(calcDate+timeStep); //prepare next timestep
		}


		// Snowdrift will overwrite VW and DW with data from 3D wind fields and update mns in SnowpackInterface
		if (snowdrift) {
			try { //Snowdrift
				// This will overwrite VW and DW with data from 3D wind fields
				snowdrift->setMeteo(t_ind, psum, psum_ph, p, vw, dw, rh, ta, tsg, vecMeteo);
				// This will update mns grid in SnowpackInterface
				snowdrift->Compute(calcDate);
			} catch (std::exception& e) {
				cout << "[E] Exception: Snowdrift compute\n";
				cout << e.what() << endl;
				throw;
			}
		}

		// Enery balance
		// This will populate direct, diffuse, reflected, direct_unshaded_horizontal, ilwr, sky_ilwr, terrain_ilwr,
		// solarAzimuth, solarElevation grids in SnowpackInterface
		try {
			eb->setStations(vecMeteo);
			eb->compute(ilwr, ta, rh, p, iswr_dir, iswr_diff, calcDate);
		} catch (std::bad_alloc&) {
			cout << "[E] AlpineControl : Virtual memory exceeded\n";
		} catch (std::exception& e) {
			cout << "[E] Exception: Ebalance compute\n";
			cout << e.what() << endl;
			throw;
		}

		// try { //Data Assimilation
		// 	if (da) da->Compute(calcDate);
		// } catch (std::exception& e) {
		// 	cout << "[E] Exception: Data Assimilation compute\n";
		// 	cout << e.what() << endl;
		// 	throw;
		// }

		if (snowpack!=NULL) {
			snowpack->setMeteo(psum, psum_ph, vw, dw, rh, ta, tsg, calcDate);
			if (enable_simple_snow_drift || enable_snowdrift2d) snowpack->setVwDrift(vw_drift, calcDate);
			snowpack->calcNextStep();
		}

		// Check if elapsed time exceeds specified maximum run time
		const double tmp_elapsed = elapsed.getElapsed();
		bool ForceStop = false;
		if (isMaster) {
			ForceStop = (max_run_time > 0. && tmp_elapsed > max_run_time);
			if (max_run_time > 0.) {
				if (ForceStop) {
					cout << std::fixed << "[W] !!! Elapsed time (" << setprecision(1) << tmp_elapsed << " seconds) exceeds specified maximum run time (" << setprecision(1) << max_run_time << " seconds) !!!\n        ---> Force writing restart files (if requested) and exiting...\n";
				} else {
					cout << std::fixed << "[i] Maximum run time set to: " << setprecision(1) << max_run_time << " seconds ---> time remaining: " << (max_run_time - tmp_elapsed)/3600. << " hours.\n";
				}
			}
		}
		MPIControl::instance().broadcast(ForceStop);

		try { //do some outputs (note, if ForceStop == true, output will be done outside of the loop)
			if ( snowpack && out_snow && (t_ind > 0) && !ForceStop ) {
				const unsigned int output_step = static_cast<unsigned int>( Optim::round(snow_days_between*24.) );
				if (snow_days_between>0 && (t_ind%output_step)==0)
					snowpack->writeOutputSNO(calcDate);
			}
		} catch (std::exception& e) {
			printf("[e] Exception: time dependent outputs\n");
			cout << e.what() << endl;
			throw;
		}

		if (isMaster) {
			cout << "[i] timing (seconds spent): " << std::setprecision(3) << std::fixed<< "\n";
			cout << "\tmeteo=" << meteo.getTiming() << "  ";

			if (eb) cout << "ebalance=" << eb->getTiming() << "  ";
			if (snowdrift) cout << "snowdrift=" << snowdrift->getTiming() << "  ";
			if (snowpack) cout << "snowpack=" << snowpack->getTiming() << "  ";

			cout << "\n\ttotal=" << elapsed.getElapsed()-elapsed_start << endl;
		}

		if (ForceStop) break;
		calcDate += timeStep; //move to next time step

	} /* For all times max_steps */

	//Finish the program: Write SNO Files and put final output on the screen
	if (eb && eb->hasSP()) eb->writeSP(max_steps);
	if (snowpack && out_snow && !nocompute){
		if (isMaster) cout << "[i] Simulation finished, writing output files...\n";
		snowpack->writeOutputSNO(calcDate);
	 }
}
