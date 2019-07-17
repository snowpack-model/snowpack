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
#ifndef ALPINECONTROL_H
#define ALPINECONTROL_H

#include <stdio.h>

#include <meteoio/MeteoIO.h>
#include <snowpack/libsnowpack.h>

#include <alpine3d/SnowpackInterface.h>
#include <alpine3d/snowdrift/SnowDrift.h>
#include <alpine3d/ebalance/EnergyBalance.h>
#include <alpine3d/DataAssimilation.h>
#include <alpine3d/MeteoObj.h>
#include <alpine3d/MPIControl.h>

/**
 * @page restarts Restarts
 * It is possible to regularly dump the exact status of each simulated pixel and then restart from a previously saved such backup.
 * In order to do so, it is necessary to enable writting out these status files, to copy them at the required timestep and then tell the model
 * to reread them when starting. These files are just regular SMET ".sno" files such as produced and used by Snowpack but for each simulated pixel. 
 * 
 * The following keys control this process, in the [Output] section:
 * 	+ SNOW_WRITE: is set to TRUE, enable writing these files;
 * 	+ SNOW_DAYS_BETWEEN: how many days between two updates to these files;
 *
 * On the command line, the switch "--restart" tells Alpine3D to read one file for each pixel instead of one file per land use class. Of course, all
 * the usual configuration keys for the "sno" files apply (see \ref reading_snow_files "reading snow files").
 * 
 * @note You MUST set the key "EXPERIMENT_NAME" in the [Output] section since this is used to name the files necessary for restarts
 */

class AlpineControl
{
	public:
		AlpineControl(SnowpackInterface *mysnowpack,
		              SnowDriftA3D *mysnowdrift,
		              EnergyBalance *myeb,
		              DataAssimilation *myda,
		              Runoff *myrunoff,
		              const mio::Config& cfg,
		              const mio::DEMObject& dem);
		void Run(mio::Date i_startdate, const unsigned int max_steps);
		void setNoCompute(bool i_nocompute) {nocompute = i_nocompute;}

	private:
		MeteoObj meteo;
		
		SnowpackInterface *snowpack;
		SnowDriftA3D *snowdrift;
		EnergyBalance *eb;
		DataAssimilation* da;
		Runoff* runoff;
		double snow_days_between;
		double max_run_time;
		bool enable_simple_snow_drift;
		bool nocompute, out_snow; // no computation, only parse inputs (check mode)
};

#endif
