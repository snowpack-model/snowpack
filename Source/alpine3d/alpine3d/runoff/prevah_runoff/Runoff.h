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
#ifndef RUNOFF_H
#define RUNOFF_H

#include <meteoio/MeteoIO.h>
#include <iostream>

/**
 * @page prevah_runoff PREVAH hydrological modeling
 * The runoff information provided by Alpine3D as hourly runoff grids can be used as inputs for the PREVAH model in order to generate hydrological discharge
 * (see D. Viviroli, M. Zappa, J. Gurtz, R. Weingartner, <i>"An introduction to the hydrological modelling system PREVAH and its pre- and post-processing-tools"</i>,
 * Environmental Modelling &amp; Software, <b>24</b>, 10, 2009, pp 1209-1222).
 *
 * @section prevah_workflow Workflow
 * In order to use PREVAH to compute the hydrological discharge, the following workflow must be carried out:
 *  -# Running Alpine3D with runoff grids outputs;
 *  -# Compiling the special version of PREVAH
 *  -# Once Alpine3D has completed, run the included PREVAH using the runoff grids as inputs;
 *  -# extract the discharge from the PREVAH output files.
 *
 * This means that Alpine3D needs to be properly configured to generate the runoff grids (see \ref runoff_grid "runoff grids"), the PREVAH wrapper needs
 * to be configured to read the runoff grids and PREVAH itself needs to be configured for the catchment of interest.
 *
 * @section compiling_prevah Compiling the special PREVAH version
 * A special version of PREVAH is available within the Alpine3D sources as runoff/prevah_runoff/runoff_prevah.cc and must be compiled separately, relying
 * on the included Makefile runoff/prevah_runoff/Makefile. In this Makefile, you might want to edit the variable MIODIR to point to the root of your
 * MeteoIO installation. You also need to make sure you have both the g++ and gfortran compilers installed on the computer that will do the compilation.
 * Then simply compile with "make".
 *
 * @section prevah_wrapper Configuring the PREVAH wrapper
 * Since PREVAH is used more like a library tailored to the Alpine3D model, a wrapper is provided to call PREVAH. This wrapper is configured in
 * the "runoff.ini" file. This file must contain the following keys:
 * - in the [INPUT] section:
 *       - COORDSYS, COORDPARAM
 *       - TIME_ZONE
 *       - DEM as well as the keys specific to the chosen plugin
 *       - LANDUSE as well as the keys specific to the chosen plugin
 *       - GRID2D as ARC as well as the keys specific to the ARC plugin: GRID2DPATH and if needed, A3D_VIEW
 *       - CATCHMENT to provide the catchments definition (see \ref sub_catch_input sub-catchments definition)
 * - in the [OUTPUT] section:
 *       - COORDSYS, COORDPARAM
 *       - TIME_ZONE
 * - in the [ALPINE3D] section:
 *       - RUNOFF_CFG giving the file and path to the PREVAH configuration file
 *
 * The PREVAH wrapper is then started as "runoff_prevah" with two arguments: the ISO formatted start date and ISO formatted end date. For example:
 * @code
 * runoff_prevah 2004-10-01T01:00 2005-10-01T00:00
 * @endcode
 *
 * @section prevah_configuration PREVAH configuration
 * The PREVAH model must also be configured for the catchment of interest. This is within the file pointed to by the Alpine3D::RUNOFF_CFG key.
 * This file is built with sections, sub-sections and values structure. The whole "control file" section must be filled (between [Begin control file]
 * and [End control file]) while the extended options can be discarded.
 * @code
 * [Begin control file]
 *
 * [Identification]
 * Dischmabach			! Catchment name
 * 46.8				! Latitude
 * 43.3				! Official catchment area [km2] !
 * 0.01				! cell surface in km**2
 * [End catchment identification]
 *
 * [Paths]
 * ../input/surface-grids/	! Directory with inputs (grids and tables)
 * ../output/tables/		! Output directory for tables (meteo and runoff)
 * ../output/grids/		! Output directory for grids (meteo and runoff)
 * ../output/runoff/		! Directory for model state (runoff)
 * ../output/runoff/		! Directory with rot grids (runoff)
 * [End Paths]
 *
 * [Grids]
 * dischma			! Basic Name
 * dhm				! DEM extension
 * dhm				! Interpolation limits extension
 * ezg				! Meteo zones extension
 * ezg				! Watershed limits extension for runoff module
 * lus				! Extension for land use map
 * [End Grids]
 *
 * [Tuneable Parameters]
 * dischma			! General output name
 *   46.05			! Treshold coefficient for surface runoff , Default=30.0
 *   28.02			! Storage coefficient in hours for surface runoff, Default=10.0
 *  143.42			! Storage coefficient in hours for interflow, Default=75.0
 *   2963.00			! Storage coefficient in hours for baseflow, Default=2500.0
 *  723.46			! Storage coefficient in hours for quick basflow, Default=750.0
 *    80.56			! Maximal content of the quick baseflow storages, Default=150.0
 * 0.1842			! Percolation rate in mm per hour, Default=0.1
 * [End Tuneable Parameters]
 *
 * [Storage content]
 * 1.0				! Anfangsspeicherinhalt Ob. Sp. RO+RH (mm)
 * -1.0				! Anfangsspeicherinhalt Unterer Sp.RG (mm)
 * -1500.			! RG (l/sec)
 * [End Storage content]
 *
 * [End control file]
 * @endcode
 *
 *@section prevah_output PREVAH outputs
 * The outputs are prefixed with the "general output name" given in section [Tuneable Parameters] and written in the "Output directory for tables"
 * given in section [Paths]. The standard discharge is given in the {general output name}.std file in column "rges".
 *
 */
class Runoff
{
	public:
		Runoff(mio::IOManager& sn_io, const double& /*in_thresh_rain*/);
		~Runoff();
		bool initialize(const mio::DEMObject& in_dem, const mio::Config& cfg);
		void setRunoff(const mio::Array2D<double>& /*runoff_surface*/, const mio::Array2D<double>& runoff_soil, const mio::Array2D<double>& /*glacier*/, const mio::Array2D<double>& /*psum*/, const mio::Array2D<double>& /*ta*/);
		void output(const mio::Date& i_date);

	private:
		void fillHydroTables(const unsigned int& ix, const unsigned int& iy, const double& soil_runoff);

		int nx, ny;
		mio::Array1D<double> runoff_total;

		mio::IOManager *io;
		mio::DEMObject dem;
};

#endif




