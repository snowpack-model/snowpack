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
#ifndef MAINPAGE_H
#define MAINPAGE_H

 /**
 * @mainpage Table of content
 * -# External Links
 *    -# <A HREF="https://models.slf.ch/p/alpine3d/">Alpine3D's home page</A>
 *          -# <A HREF="https://models.slf.ch/p/alpine3d/page/Getting-started/">Installation, compilation</A>
 *          -# <A HREF="https://models.slf.ch/p/alpine3d/page/Running-a-simulation/">Running a simulation</A>
 * -# End User documentation
 *     -# General principles
 *        -# \subpage model_principles "Model principles"
 *        -# \subpage inputs "Inputs"
 *        -# \subpage outputs "Outputs"
 *        -# \subpage tools "Simulation tools"
 *     -# Modules configuration
 *        -# \subpage radiation_balance "Radiation balance"
 *        -# \subpage snowpack "Snowpack"
 *        -# \subpage snowdrift "Snowdrift"
 *        -# \subpage runoff "Runoff"
 *        -# \subpage glaciers "Glaciers katabatic flows"
 *     -# \subpage running_simulations "Running a simulation"
 * -# Expanding Alpine3D
 *        -# \subpage coding_style "Coding style"
 * 
 * <center><hr></center>
 * <center><i><small>
 * <p>
 * Alpine3D is a spatially distributed (surface), three dimensional (atmospheric) model for
 * analyzing and predicting dynamics of snow-dominated surface processes in mountainous topography.
 * It includes models for snow cover (<A HREF="https://models.slf.ch/p/snowpack/">SNOWPACK</A>),
 * vegetation and soil, snow transport, radiation transfer and runoff which can be enabled or disabled on demand.
 *
 * The model supports a variety of input options including interpolation of meteorological weather stations,
 * input from a meteorological model or from remote sensing data
 * (<A HREF="https://models.slf.ch/p/meteoio">MeteoIO</A>) and has been parallelized
 * in order to run on computing grids or clusters
 * (using <A HREF="https://en.wikipedia.org/wiki/Message_Passing_Interface">MPI</A> and <A HREF="http://openmp.org/wp/">OpenMP</A>).
 *
 * Alpine3D has a broad variety of potential applications. Most dominant is the assessment of
 * snow water resource dynamics in mountain catchments (Michlmayr et al., 2008). This includes predictions of
 * future snow on the basis of climate change scenarios (Bavay et al., 2009; Bavay et al., 2013). One exotic application of the model
 * system Alpine3D is the forecasting of surface temperatures on ski-pistes, e.g. for the Vancouver winter Olympics.
 * For this forecast local shadings (might) change surface temperature up to 5 °C.
 * </p>
 * <p>
 * This model is available under GPL version 3 or above, see <a href="http://www.gnu.org/licenses/gpl.txt">www.gnu.org</a>.
 * </p></small></i></center>
 */

/**
 * @page running_simulations Running a simulation 
 *
 * @section running_intro Introduction
 * 
 * @subsection model_workflow Simulation workflow
 * When running a simulation, it is important to keep in mind that the model is organized as several modules that interract together. It is possible to configure
 * some parameters for the various modules and to enable/disable modules. Some modules can be used outside of Alpine3D (like
 * <A HREF="https://models.slf.ch">MeteoIO</A> that is used in various applications or libSnowpack that is used by the standalone
 * <A HREF="https://models.slf.ch">Snowpack</A> model) .
 *
 * \image html simulation_workflow.png "Simulation workflow"
 * \image latex simulation_workflow.eps "Simulation workflow" width=0.9\textwidth
 *
 * @subsection installing Installing Alpine3D
 * Please follow the instructions given on <a href="https://models.slf.ch/p/alpine3d/page/Getting-started/">the forge</a> in order to download Alpine3D (from svn, from source or from a
 * binary package) and its dependencies (Snowpack and MeteoIO, knowing that binary packages might already contain all the required dependencies). If you've
 * downloaded a binary package, there is nothing special to do, just install it on your system.
 *
 * When installing Alpine3D from sources, please keep in mind that Alpine3D's compilation process is exactly the same as Snowpack's or MeteoIO's (all are based on cmake). After
 * a successful compilation, it is necessary to install Alpine3D (the same **must** have been done for Snowpack and MeteoIO). There are two options:
 *     - system-wide install. In this case, you need to have administrator permissions. Simply set CMAKE_INSTALL_PREFIX (in cmake) to a system path (by default, it is /usr/local)
 * and install by typing '*make install*' in a terminal opened at Alpine3D's source root folder.
 *     - user install: if you do not have administrator permissions, you can still install Alpine3D into a directory where you have write permissions. Simply set CMAKE_INSTALL_PREFIX
 * (in cmake) to a directory where you can write (it is *highly* recommended to set it to '${HOME}/usr'), make sure this directory exists and is writable and then install by
 * typing '*make install*' in a terminal opened at Alpine3D's source root folder. Make sure that CMAKE_INSTALL_PREFIX/bin is in your PATH and that CMAKE_INSTALL_PREFIX/lib
 * is recognized as a library path (variables PATH and LD_LIBRARY_PATH for Linux, PATH and DYLD_FALLBACK_LIBRARY_PATH for osX,
 * see <a href="https://models.slf.ch/p/snowpack/page/Getting-started/#wikititle_4">Getting Started</a> on the forge).
 *
 * After Alpine3D has been installed, you can check that it works by opening a terminal and typing "alpine3d". Alpine3D should be found and display its help message.
 *
 * @subsection intro_directories Simulation directories setup
 * After you installed a binary package or compiled and installed Alpine3D (on most compute clusters, you will need to install in your own home directory, see
 * \ref installing "installing Alpine3D" above),
 * you can run your first simulation. We highly recommend that you use the following structure: first, create a directory for your simulation, for example
 * "Stillberg". Then, create the following sub-directories:
 * - input, to put all the input data required by the simulation
 *      - input/meteo, for the meteorological input data
 *      - input/surface-grids, for the dem, land cover and optional catchments grids
 *      - input/snowfiles, for the input .sno files
 * - output, where Alpine3D will write the results
 *      -output/grids, where Alpine3D will write the gridded results
 *      - output/snowfiles, where Alpine3D will write its status files (sno files)
 * - setup, to put the configuration of the simulation and its start scripts
 *
 * Edit the "run.sh" script to set it up with the proper start and end dates, the modules that you want to enable and the proper configuration for a sequential or parallel run.
 *
 * @section intro_run_simple Simple sequential simulation
 * This is the simplest way of running an Alpine3D simulation: it runs on one core on one computer. The drawback is that if the simulation is quite large, 
 * it might require a very long time to run or even not have enough memory (RAM) to run once there is snow in the simulated domain. In order to run a
 * sequential simulation, set \c "PARALLEL=N" in the \em run.sh script. Then run the following command in a terminal (this can be a remote terminal 
 * such as \em ssh) on the computer where the simulation should run:
 * @code
 * nohup ./run.sh &
 * @endcode
 *  \em nohup means that if you close the terminal the simulation will keep going; \em & means that the terminal prompt can acept other commands after you've submitted this one. In order to monitor what is going on with your simulation, simply run something such as (\em -f means that it keeps updating with the new content in this file. Replace it with something such as \em -500 to show the last 500 lines of this file):
 * @code
 * tail -f stdouterr.log
 * @endcode
 * 
 * If you need to terminate the simulation, first find out its Process ID (PID) by doing
 * @code
 * ps ux
 * @endcode
 * Then kill the process
 * @code
 * kill {PID}
 * @endcode
 * 
 * @section intro_parallel Parallel simulations
 * When a simulated domain gets bigger, the computational requirements (memory and runtime) get bigger. In order to reduce these requirements, it is possible to
 * run the simulation in parallel accross multiple cores or nodes.
 *
 * @subsection intro_run_openmp Multi-cores simulation
 * This is the easiest way to run a parallel simulation because it does not require any specific software, only a compiler that supports <a href="http://openmp.org">OpenMP</a> (see also its <a href="https://en.wikipedia.org/wiki/OpenMP">wikipedia</a> page). Such compilers are for example gcc, clang, etc. The limitations on the memory still remain (ie a simulation requiring lots of memory will still only have access 
 * to the local computer's memory) but the run time will be roughtly divided by the number of available cores that are given to the simulation. In order
 * to run such a simulation, please compile Alpine3D with the \b OpenMP option set to ON in cmake. Then in the simulation's \em run.sh file, set \c "PARALLEL=OPENMP" as well as the number of cores you want to us as \c "NCORES=". Then run the simulation as laid out in the previous section.
 * 
 * @subsection intro_run_mpi Multi-nodes simulation
 * This is the most powerful way of running a simulation: the load is distributed among distinct computing nodes, therefore reducing the amount 
 * of memory that must be available on each node. For very large simulations, this might be the only way to proceed. This is achieved by relying on 
 * <a href="https://en.wikipedia.org/wiki/Message_Passing_Interface">MPI</a> to exchange data between the nodes and distribute the workload. In
 * order to run such a simulation, please compile Alpine3D with the \b MPI option set to ON in cmake. Then in the simulation's \em run.sh file, set \c "PARALLEL=MPI" as well as the number of processors/cores you want to us as \c "NPROC=" and a machine file. This machine file contains the list of machines
 * to use for the simulation as well as how many processors/cores to use. For example, such as file could be:
 * @code
 * 192.168.0.11:2
 * 192.168.0.5
 * 192.168.0.4:3
 * 192.168.2.25
 * @endcode
 * Then run the simulation as laid out in the previous section.
 *
 * \image html mpi_scaling.png "Scaling of Alpine3D with the number of processors when using MPI (Gemstock_1m, 1 year simulation)"
 * \image latex mpi_scaling.eps "Scaling of Alpine3D with the number of processors when using MPI (Gemstock_1m, 1 year simulation)" width=0.9\textwidth
 * 
 * \note Please make sure that the environment variables $TEMPDIR, $TEMP and $TMP (if defined) don't point to a shared drive (such as an NFS mounted home directory), 
 * otherwise MPI might run very slowly.
 * 
 * @subsection intro_run_sge Sun/Oracle Grid Engine
 * If your computing infrastructure relies on <a href="https://en.wikipedia.org/wiki/Oracle_Grid_Engine">Sun/Oracle Grid Engine (SGE)</A> (for example on a computing cluster), 
 * you need to do things differently. First, the job management directives/options must be provided, either on the command line or in the \em run.sh script. 
 * These lines rely on special comments, starting with \em "#" followed by \em "$":
 * @code
 * #$ -m aes
 * #$ -N {simulation name}
 * #$ -S /bin/bash
 * #$ -cwd
 * #$ -j y
 * #$ -pe smp {number of cores to use}
 * @endcode
 * The last line specifies the computing profile that should be used. Since the job manager will allocate the ressources, there is no need to provide 
 * either NCORES or NPROC. The machine file (for MPI) is also not used. Then, submit the computing job to the system: \c "qsub ./run.sh". This 
 * should return almost immediately with a message providing the allocated job number. This job number is useful to delete 
 * the job \c "qdel {job_number}" or querry its status \c "qstat {job_number}" (or \c "qstat" to see for all jobs).
 * 
 * If the job submission fails with an error message such as \em unknown \em command, please check that there is no extra "#$" in the script. 
 * This happens frequently when commenting out some part of the script and is mis-interpreted by SGE. In such a case, simply add an 
 * extra "#" in front of the comment.
 * 
 * \note At WSL, the computing profiles are either \em smp for shared memory runs such as OpenMP or \em orte for load distributed among 
 * nodes such as MPI
 *
 * @section data_strategies Data strategies
 * Depending on the meteorological data availability, the meteorological data quality as well as convenience (if the raw data is difficult to access, for example), there
 * are different data strategies that can be used. This is illustrated in the figure below.
 *
 * \image html data_strategies.png "Meteorological data strategies"
 * \image latex data_strategies.eps "Meteorological data strategies" width=0.9\textwidth
 *
 * The simplest case is <b>1</b>: the meteorological forcings are directly read and processed by Alpine3D. Since it embbeds MeteoIO, it can perform the whole
 * preprocessing on the fly. In order to have a look at the pre-processed data (or to spare a lengthy data extraction), it is possible to preprocess the data first,
 * then dump them to files and run Alpine3D from these intermediate file (case <b>2</b>). Since Alpine3D anyway embbeds MeteoIO, it would be possible to add
 * some pre-processing to be done on the fly within Alpine3D (such as resampling to the sampling rates that Alpine3D needs). This is actually highly recommended
 * in order to guarantee that Alpine3D gets the data it needs (ie proper sampling rates as well as basic checks on the ranges, etc).
 *
 * In some cases, it might be interesting to first run the data through Snowpack in order to produce some meteorological parameter from other parameters, such as
 * the ISWR from RSWR (since Snowpack computes the snow albedo, the generate data would be much better than assuming a fixed albedo) or to generate
 * precipitation from the measured snow height. This is shown in case <b>3</b>.
 *
 * Case <b>4</b> shows a strategy that could be used to prepare an Alpine3D simulation: it consists of simple Snowpack simulations performed at spatially
 * interpolated locations, using MeteoIO's virtual stations concept. This is useful to (relatively) quickly validate the spatially interpolated fields with some
 * measured data (specially if there are some meteorological data that could be compared to the spatially interpolated meteorological fields).
 *
 *
 * @subsection snowpack_coupling Coupling with Snowpack
 * Alpine3D needs spatially interpolated forcings for each grid points. Unfortunately, it limits the choice of forcing parameters: for some parameters
 * (such as HS or RSWR), there are no reliable interpolation methods. One way to make use of the existing measurements that could not be easily
 * interpolated is to run a <A HREF="https://models.slf.ch/p/snowpack">Snowpack</A> simulation at the stations that provided these measurements,
 * then use alternate, computed parameters (such as PSUM or ISWR) as inputs to Alpine3D.
 *
 * This process is made easier by writing Snowpack's outputs in the smet format and making sure all the necessary parameters are written out.
 * This means that Snowpack should be configured along the following lines (only using one slope):
 * @code
 * [Output]
 * TS_WRITE        = TRUE
 * TS_FORMAT       = SMET
 * TS_DAYS_BETWEEN = 0.04166667	;so we get hourly values
 *
 * OUT_CANOPY = FALSE
 * OUT_HAZ    = FALSE
 * OUT_SOILEB = FALSE
 * OUT_HEAT   = FALSE
 * OUT_T      = FALSE
 * OUT_STAB   = FALSE
 * OUT_LW     = TRUE
 * OUT_SW     = TRUE
 * OUT_MASS   = TRUE
 * OUT_METEO  = TRUE
 *
 * AVGSUM_TIME_SERIES = TRUE
 * CUMSUM_MASS        = FALSE
 * PRECIP_RATES       = FALSE
 * @endcode
 *
 * Then the output smet files produced by Snowpack can be directly used as Alpine3D inputs, performing some on-the-fly renaming and conversions
 * (here, from split precipitation to precipitation/phase):
 * @code
 * [Input]
 * METEO      = SMET
 * METEOPATH  = ../input/meteo
 * STATION1   = WFJ2
 *
 * PSUM_S::MOVE = MS_Snow
 * PSUM_L::MOVE = MS_Rain
 * HS::MOVE    = HS_meas	;so we can still compare measured vs modelled snow height
 * TSG::MOVE   = T_bottom	;so we can compare the ground temperatures
 * TSS::MOVE   = TSS_meas	;so we can compare the surface temperatures
 *
 * WFJ2::KEEP = TA TSS TSG RH ISWR ILWR HS VW DW PSUM_S PSUM_L PSUM PSUM_PH	;so we do not keep all kind of unnecessary parameters
 *
 * PSUM_PH::create     = PRECSPLITTING
 * PSUM_PH::PRECSPLITTING::type   = THRESH
 * PSUM_PH::PRECSPLITTING::snow   = 274.35
 * PSUM::create     = PRECSPLITTING
 * PSUM::PRECSPLITTING::type   = THRESH
 * PSUM::PRECSPLITTING::snow   = 274.35
 *
 * [SNOWPACK]
 * ENFORCE_MEASURED_SNOW_HEIGHTS = FALSE
 * @endcode
 *
 * Of course, other stations can also be part of the meteo input and their inputs should remain unaffected (assuming they don't use parameter names
 * such as MS_Snow, MS_Rain or HS_meas and assuming that their parameters are not rejected by the KEEP command). As a side note, the last three parameter
 * renamings (with MOVE) must be set when using Snowpack's outputs as inputs for a new Snowpack simulation.
 */

/**
 * @page model_principles Model principles
 * Here, we expose the core principles underlying the Alpine3D model. This should just give a quick overview of the model and help you understand the
 * global architecture of Alpine3D. If you want to go deeper into the details, please have a look at the publications covering the whole model
 * (M. Lehning et al. <i>"ALPINE3D: a detailed model of mountain surface processes and its application to snow hydrology"</i>, Hydrological Processes, \b 20.10, 2006, pp 2111-2128; and
 * P. Kuonen, M. Bavay, M. Lehning, <i>"POP-C++ and ALPINE3D: petition for a new HPC approach"</i>, Advanced ICTs for disaster management and threat detection: collaborative and distributed frameworks, IGI Global, 2010, pp 237-61.)
 * The individual modules are described here below and contain references to the relevant papers.
 *
 * @section at_the_core At the core...
 * Here we expose the very foundations of Alpine3D. These remain valid independently of which modules are enabled when running the model.
 *
 * @subsection principles_snowpack Distributed 1D soil/snow/canopy column
 * \image html distributed_sn.png "Distributed SNOWPACK over the domain taking into account the land cover"
 * \image latex distributed_sn.eps "Distributed SNOWPACK over the domain taking into account the land cover" width=0.9\textwidth
 * At the core of the model, is the <A HREF="http://models.slf.ch/p/snowpack/">SNOWPACK</A> model, a physically based,
 * energy balance model for a 1D soil/snow/canopy column.
 * This gives us a very detailed description of the snow stratigraphy and a very good evaluation of the mass and energy balance (therefore also of
 * quantities such as Snow Water Equivalent (SWE) or temperature profile). This 1D energy balance is performed for each pixel of the domain
 * (therefore it is a distributed SNOWPACK simulation) and for one time step (usually one hour). Any quantity that the user would like to get
 * out of the simulation can be written out from this module.
 *
 * @subsection distributed_meteo Distributed meteo fields
 * \image html 2d_interpolations.png "Spatially interpolating the meteorological fields"
 * \image latex 2d_interpolations.eps "Spatially interpolating the meteorological fields" width=0.9\textwidth
 * In order to perform a SNOWPACK simulation at every pixel of the domain, it is necessary to get the meteorological forcing for each pixel.
 * But the measured meteorological parameters are usually measured by a set of stations, which means that the data is available at a set of points.
 * Interpolating these points measurements to every pixels of the domain is performed by the means of statistical interpolations with
 * <A HREF="http://models.slf.ch/p/meteoio">MeteoIO</A>. if the forcing data is coming out of another model (such as a meteorological model),
 * most probably the input grids have a resolution that is very insufficient for Alpine3D and therefore need downscaling. If the downscaling factor
 * is very large, we often end up with only a few points from the meteorological model that are part of the Alpine3D domain, therefore such points
 * can be considered as "virtual stations" and spatially interpolated similarly to weather stations.
 *
 * @section lateral_fluxes Lateral fluxes
 * The core principles laid out in the previous section rely on the assumption that there are no lateral fluxes,
 * which is too strong of an assumption. Therefore the lateral fluxes deemed relevant are introduced by other modules:
 * - the EBalance module computes the radiation fields, taking into account atmospheric cloudiness, topographic shading effects and reflections by the
 * surrounding terrain.
 * - the SnowDrift module that simulates the transport of snow by the wind. It performs a 3D simulation of the saltation, suspension and diffusion processes.
 * - the runoff module that collects the precipitation and/or melt water at each pixel to transfer it to an hydrological routing module
 *
 * @subsection principles_ebalance Radiation balance
 * \image html ebalance.png "Radiation balance with shading and terrain reflections"
 * \image latex ebalance.eps "Radiation balance with shading and terrain reflections" width=0.9\textwidth
 * Once the albedo of each pixels of the domain have been initialized or taken from the last time step, the radiation balance is computed. First, the
 * incoming short wave radiation measured at one reference station is used to compute the splitting (between direct and diffuse, see
 * D. G. Erbs, S.A. Klein, J.A. Duffie, <i>"Estimation of the diffuse radiation fraction for hourly, daily and monthly-average global radiation"</i>, Solar Energy, <b>28</b>, 4, 1982, Pages 293-302 and summarized in M. Iqbal, <i>"An introduction to solar radiation"</i>, 1983, Academic Press,  ISBN: 0-12-373750-8).
 * This splitting will be assumed to be constant over the whole domain. Then, according to the meteorological parameters and elevation of each pixel,
 * the direct and diffuse radiation fields are computed. Since the position of the sun has been computed (
 * J. Meeus, <i>"Astronomical Algorithms"</i>, 1998, 2nd ed, Willmann-Bell, Inc., Richmond, VA, USA, ISBN 0-943396-61-1),
 * it is used to compute the topographic shading for each pixel.
 * If the terrain reflections have been enabled, a radiosity approach is used to compute the reflections by the surrounding terrain (
 * N. Helbig, H. Löwe, M. Lehning, <i>"Radiosity Approach for the Shortwave Surface Radiation Balance in Complex Terrain"</i>, Journal of the Atmospheric Sciences, \b 66.9, 2009).
 * Finally, the direct and diffuse radiation fields are returned.
 *
 * @subsection principles_snowdrift Snowdrift
 * \image html snowdrift.png "Snowdrift: saltation, suspension, sublimation"
 * \image latex snowdrift.eps "Snowdrift: saltation, suspension, sublimation" width=0.9\textwidth
 * Externally computed wind fields (for example with <A HREF="http://arps.ou.edu/">ARPS</A>) are assigned to each time steps.
 * If the surface shear stress exceeds a given threshold at a given pixel, the saltation will be computed. This in turn can feed the suspension
 * if the saltation layer is saturated. While in suspension, some of the mass will sublimate and contribute to the relative humidity field (
 * C. D. Groot Zwaaftink et al. <i>"Drifting snow sublimation: A high‐resolution 3‐D model with temperature and moisture feedbacks"</i>, Journal of Geophysical Research: Atmospheres (1984–2012), \b 116.D16, 2011).
 *
 * @subsection principles_runoff Runoff
 * \image html runoff.png "Runoff: simple bucket approach with PREVAH or sub-catchments sums"
 * \image latex runoff.eps "Runoff: simple bucket approach with PREVAH or sub-catchments sums" width=0.9\textwidth
 * Several options are available for collecting the melt water or precipitation running out of each pixel. The historical approach relies on the
 * PREVAH hydrological modeling system (
 * D. Viviroli et al. <i>"An introduction to the hydrological modelling system PREVAH and its pre-and post-processing-tools"</i>, Environmental Modelling \& Software, \b 24.10, 2009, pp 1209-1222)
 * to perform on-the-fly hydrological simulation. Another approach consists of collecting the runoff sums for each sub-catchments defined by the user
 * (M. Bavay, T. Grünewald, M. Lehning, <i>"Response of snow cover and runoff to climate change in high Alpine catchments of Eastern Switzerland"</i>, Advances in Water Resources, \b 55, 2013, pp 4-16).
 * Finally, it is also possible to output the hourly distributed runoff (ie the runoff for each pixel, once per hour) and process this runoff in
 * an external model
 * ( F. Comola et al. <i>"Comparison of hydrologic response between a conceptual and a travel time distribution model for a snow-covered alpine catchment using Alpine3D"</i>, EGU General Assembly Conference Abstracts, Vol. \b 15, 2013;
 * J. Magnusson et al. <i>"Quantitative evaluation of different hydrological modelling approaches in a partly glacierized Swiss watershed"</i>, Hydrological Processes, \b 25.13, 2011, pp 2071-2084).
 * This approach enables the external model to be calibrated without having to re-run Alpine3D.
 *
 */

/**
 * @page inputs Inputs
 * Several categories of data are necessary for the input:
 *     - Landscape data: Digital Elevation Model and Land Cover model
 *     - Snowcover and soil data
 *     - Meteorological data
 *
 * The landscape data can be prepared with any GIS software (for example, the Open Source <A HREF="http://qgis.org/">QGIS</A>), the meteorological data with 
 * statistical computing tools (for example <A HREF="https://www.r-project.org/">R</A>) or scripts (for example 
 * <A HREF="https://www.python.org/">Python</A> or <A HREF="https://en.wikipedia.org/wiki/AWK">AWK</A>) and the snowcover data has to be prepared
 * with scripts.
 * 
 * @section landscape_input Landscape data
 * @subsection dem_input Digital Elevation Model
 * \image html DEM.png "Example of a Digital Elevation Model grid"
 * \image latex DEM.eps "Example of a Digital Elevation Model grid" width=0.9\textwidth
 * Once you have defined the domain that you want to simulate, you have to get a
 * <A HREF="http://www.fgi.fi/fgi/themes/digital-elevation-model">Digital Elevation Model</A> that contains this area.
 * The DEM must be rectangular, but you will be able to restrict the actual domain within this
 * rectangular area by settings the cells you want to ignore to nodata. However, you need to keep in mind the following:
 *    - The DEM defines the spatial grid that will be used in Alpine3D
 *    - Therefore, the resolution of the DEM is the spatial resolution of the simulation
 *    - Each cell whose DEM is nodata will be skipped by all modules of Alpine3D
 *    - In order to properly compute the shading effects, do not exclude cells that could cast shades on an interesting part of the domain!
 *
 * \anchor dem_geolocalization
 * These last two points are important for both the definition of the rectangular DEM as well as filling the necessary cells with nodata.
 * All other grids (either land cover or potential meteorological grids) will have to use the exact same
 * <A HREF="https://en.wikipedia.org/wiki/Geolocation">geolocalization</A> as the DEM
 * (same position, same size, same resolution) unless otherwise specified.
 *
 * DEM data can either come from your own measurements (for example from laser scans, see P. Axelsson, <i>"DEM generation from laser scanner data using adaptive TIN models"</i>, International Archives of Photogrammetry and Remote Sensing, \b 33.B4/1, PART 4, 2000, pp 111-118),
 * from your national topographic service (<A HREF="http://www.swisstopo.admin.ch/internet/swisstopo/en/home/products/height.html">Swiss Topo</A> for Switzerland,
 * <A HREF="http://professionnels.ign.fr/catalogue">IGN</A> for France,
 * <A HREF="http://ned.usgs.gov/">National Elevation Dataset</A> for the USA, <A HREF="https://geoservice.ist.supsi.ch/helidem/">HeliDEM</A> for the southern Alps)
 * or from global initiatives (
 * <A HREF="https://lta.cr.usgs.gov/GTOPO30">GTOPO30</A>, a 30-arc-second global dem,
 * <A HREF="http://www.ngdc.noaa.gov/mgg/topo/globe.html">GLOBE</A>, a 30-arc-second global dem,
 * <A HREF="http://srtm.usgs.gov/">SRTM</A>, a 1-arc-second global dem between 60N and 56S averaged at 90m resolution,
 * <A HREF="http://gdem.ersdac.jspacesystems.or.jp/">ASTER DEM</A>, a 30m resolution global dem between 83N and 83S).
 *
 * @subsection lus_input Land Cover model
 * \image html LUS.png "Example of a Land Cover grid"
 * \image latex LUS.eps "Example of a Land Cover grid" width=0.9\textwidth
 * For each cell of the domain, a land cover code must be provided. This must be using the exact same geolocalization as the DEM.
 * Such data are usually available in various classifications depending on the country (such as
 * <A HREF="http://www.eea.europa.eu/publications/COR0-landcover">CORINE</A> for Europe,
 * the <A HREF="http://landcover.usgs.gov/">US National Land Cover Dataset</A> for the USA,
 * the <A HREF="http://www.countrysidesurvey.org.uk/">countryside survey</A> for the UK,
 * <A HREF="https://www.bfs.admin.ch/bfs/de/home/statistiken/raum-umwelt/nomenklaturen/arealstatistik/noas2004.html">Arealstatistik NOAS04</A> for Switzerland,
 * <A HREF="http://data.ess.tsinghua.edu.cn/">GLC</A> for a 30m resolution global land cover or
 * <A HREF="http://data.fao.org/map?entryId=6c34ec8b-f31e-4976-9344-fd11b738a850">FAO GeoNetwork</A> for multiple land cover data sets including a 30" resolution global land cover).
 * 
 * Currently, Alpine3D improperly calls the Land Cover Model <i>"Land Use"</i>, abbreviated as LUS and uses an ARC ascii file
 * (see <A HREF="https://models.slf.ch/docserver/meteoio/html/arc.html">MeteoIO's documentation</A>) with PREVAH landuse codes 
 * that have the format 1LLDC where:
 * - LL is the land use code as given in the table given below
 * - D is the soil depth (unused)
 * - C is the field capacity (unused)
 *
 * <center><table border="0">
 * <caption>PREVAH land cover codes</caption>
 * <tr><td>
 * <table border="1">
 * <tr><th>land use (vegetation)</th><th>Prevah land use classes</th></tr>
 * <tr><td>01</td><td>water</td></tr>
 * <tr><td>02</td><td>settlement</td></tr>
 * <tr><td>03</td><td>coniferous forest</td></tr>
 * <tr><td>04</td><td>decidous forest</td></tr>
 * <tr><td>05</td><td>mixed forest</td></tr>
 * <tr><td>06</td><td>cereals</td></tr>
 * <tr><td>07</td><td>pasture</td></tr>
 * <tr><td>08</td><td>bush</td></tr>
 * <tr><td>09</td><td>undefined</td></tr>
 * <tr><td>10</td><td>undefined</td></tr>
 * <tr><td>11</td><td>road</td></tr>
 * <tr><td>12</td><td>undefined</td></tr>
 * <tr><td>13</td><td>firn</td></tr>
 * <tr><td>14</td><td>bare ice</td></tr>
 * <tr><td>15</td><td>rock</td></tr>
 * </table></td><td><table border="1">
 * <tr><th>land use (vegetation)</th><th>Prevah land use classes</th></tr>
 * <tr><td>16</td><td>undefined</td></tr>
 * <tr><td>17</td><td>undefined</td></tr>
 * <tr><td>18</td><td>fruit</td></tr>
 * <tr><td>19</td><td>vegetables</td></tr>
 * <tr><td>20</td><td>wheat</td></tr>
 * <tr><td>21</td><td>alpine vegetation</td></tr>
 * <tr><td>22</td><td>wetlands</td></tr>
 * <tr><td>23</td><td>rough pasture</td></tr>
 * <tr><td>24</td><td>subalpine meadow</td></tr>
 * <tr><td>25</td><td>alpine meadow</td></tr>
 * <tr><td>26</td><td>bare soil vegetation</td></tr>
 * <tr><td>27</td><td>free</td></tr>
 * <tr><td>28</td><td>corn</td></tr>
 * <tr><td>29</td><td>grapes</td></tr>
 * <tr><td>30-99</td><td>undefined</td></tr>
 * </table></td></tr>
 * </table></center>
 *
 * \remarks There is a common confusion between land use and land cover, when actually a land cover is determined by direct
 * observations (<A HREF="http://stats.oecd.org/glossary/detail.asp?ID=6489">OECD definition</A>) while a
 * land use requires socio-economic interpretation of the activities that take place on that surface
 * (<A HREF="http://stats.oecd.org/glossary/detail.asp?ID=6493">OECD definition</A>), see
 * P. Fisher, A. Comber, and R. Wadsworth, <i>"Land use and Land cover: Contradiction or Complement"</i>, Re-presenting GIS, 2005, pp85-98).
 *
 * @subsection sub_catch_input Catchments definition
 * \image html catchments.png "Example of catchments definition"
 * \image latex catchments.eps "Example of catchments definition" width=0.9\textwidth
 * For hydrological modeling, hydrological subcatchments must be defined. Both the precipitation, snow melt and glacier melt
 * contributions are provided as well as the total catchment runoff, for each subcatchment. The subcatchments are defined by
 * creating a grid that contains a code describing which catchments this cell belongs to, based on sums of powers of two.
 * For example, a pixel belonging to catchments 0 and 3 would receive the code: 2^0+2^3=9. A pixel belonging to catchments
 * 2, 5, 6 would receive the code 2^2+2^5+2^6=100. Finally, this grid must have the same geolocalization as the dem.
 *
 * @section sno_input Snow cover and soil data
 * \image html soil.png "Initial soil and snow definition"
 * \image latex soil.eps "Initial soil and snow definition" width=0.9\textwidth
 * Either each cell must be assigned an initial soil and snow profile (as well as a few Canopy parameters). Usually, to make things easier, the simulation starts
 * at a time when no snow is present in the domain, making the snow profile empty. The same file also contains potential
 * soil layers data. This file is written in a Snowpack snowfile supported format (see Snowpack documentation, "Data File Formats" > "Single Snow Profiles") and usually kept
 * with all other similar files in a separate directory (in order to keep the simulation tidy). Finally, the profile and layers
 * must be dated from before the simulation starts.
 *
 * There are two possibilities for assigning these files to each cell of the simulation domain (see \subpage reading_snow_files "reading snow files"):
 *     - by land cover classes. In this case, every cell receives an initial soil/snow profile based on its land cover class. The files must be named as {LAND_USE_CLASS}_{EXPERIMENT_NAME}.{ext};
 *     - independently for each (i,j) pixel. In this case, the files must be named as {i_index}_{j_index}_{EXPERIMENT_NAME}.{ext};
 *
 * @note In any case, you MUST set the key "EXPERIMENT_NAME" in the [Output] section.
 *
 * @section meteo_input Meteorological data
 * Usually, the meteorological inputs are provided as point-measurement time series and spatially interpolated within alpine3d (using meteoio).
 * But it is also possible to provide some (or all) data as gridded data.
 *
 * Meteorological inputs often come from automatic weather stations that are part of a measurement network (for example
 * M. Lehning et al. <i>"A network of automatic weather and snow stations and supplementary model calculations providing snowpack information for avalanche warning"</i>, Proceedings of the International Snow Science Workshop “a merging of theory and practice”, <b>27</b>, 1998 or the US <A HREF="http://www.wcc.nrcs.usda.gov/snow/">SNOTEL</A> network)
 * or are deployed for a specific experiment and area. The networks are often managed by national weather services.
 *
 * It is also possible to find some other data producers such as
 * <A HREF="http://mesonet.agron.iastate.edu/request/download.phtml?network=CN_ASOS">airports</A>, <A HREF="http://www.ogimet.com/index.phtml.en">Synop stations</A> (a script to download and convert the data is available in <i>MeteoIO's tools directory</i>),
 * railways operators,
 * <A HREF="https://pub-apps.th.gov.bc.ca/saw-paws/weatherstation">highways</A> (or <A HREF="http://www.chart.state.md.us/travinfo/weatherstationdata.asp">this</A>),
 * <A HREF="http://www.parkcitymountain.com/site/mountain-info/conditions/weather/weather-station-reports">ski lift operators</A>, 
 * <A HREF="http://www.wxqa.com/index.html">network of citizen weather stations</A>.
 *
 * Another source of data can be reanalysis runs performed with weather forecasting models on measured data
 * (<A HREF="http://www.metoffice.gov.uk/datapoint/product/uk-daily-site-specific-forecast">UK metoffice</A>,
 * <A HREF="http://www.yr.no/place/Norway/">Norvegian Meteorological Institute</A>,
 * <A HREF="http://www.euro4m.eu/datasets.html">Euro4M</A> for all of Europe as well as worldwide
 * <A HREF="http://www.esrl.noaa.gov/psd/data/gridded/reanalysis/">NOAA ESRL</A> and
 * <A HREF="http://www.ecmwf.int/en/research/climate-reanalysis">ECMWF</A> (with associated <A HREF="http://data-portal.ecmwf.int/data/d/interim_full_invariant/">DEM</A>) reanalysis).
 *
 * @subsection stations_input Point measurements
 * This is mostly centered around the concept of station: a specific point in space where multiple parameters are measured for a given period.
 * The file formats vary greatly depending on which meteo plugin is used. It is nevertheless always possible to use stations that don't cover
 * the whole period or contain data gaps. It is also possible to use stations that are outside the domain (but they should not be too far away,
 * otherwise this does not make sense!). The model will run with an hourly resolution, so the data should be either hourly or could be
 * meaningfully resampled to an hourly resolution (i.e., daily measurements would usually not qualify). Please note that in any case, the
 * following must always be provided:
 *    - at least \em one of:
 *        - air temperature (TA)
 *        - relative humidity (RH)
 *        - wind speed (VW)
 *        - precipitation (PSUM)
 *    - at least \em one point offering simultaneously the following:
 *        - air temperature (TA)
 *        - relative humidity (RH)
 *        - incoming short wave radiation (ISWR)
 *        - incoming long wave radiation (ILWR)
 *
 * If the air pressure (P, usually in Pa) is provided, it will be used for improving some of the parametrizations, but this is absolutely not mandatory. 
 * However, it is recommended to provide stations well distributed over the domain and over the elevation range, so the elevation gradients
 * can be properly computed.
 *
 * @subsection grids_input Gridded meteorological data
 * This relies on meteoio's USER spatial interpolation algorithm. The grids must follow a specific naming scheme and be placed in a specific
 * directory so they can be found and the grids must have the same \ref dem_geolocalization "geolocalization" as the dem. This is detailed in
 * meteoio's documentation.
 *
 * Please note that if not all necessary grids are provided, meteoio will revert to point measurements inputs (according to the spatial interpolation
 * algorithms declared in the configuration file). It is therefore possible to only provide a few specific grids according to the needs
 * (for example, the precipitation grids only when precipitation occurs).
 *
 */

/**
 * @page outputs Outputs
 * Alpine3D writes four types of outputs:
 * - grids, that is the distributed value of a given parameter (see \ref gridded_outputs "gridded outputs");
 * - .met and .pro for points of interests. These are the standard SNOWPACK outputs and can be written out for any number of points (see \ref poi_outputs "POI outputs");
 * - sno files, that is the status of each pixel in respect with snow information. These files are necessary in order to restart Alpine3D from a
 *  previous point (see \subpage restarts "restarts" and \subpage reading_snow_files "reading snow files");
 * - subcatchments runoff sums. One file per subcatchment is generated and  contains the sums of all runoff components at an hourly resolution as well as some other
 * relevant catchment parameters (such as mean air temperature, etc).  The runoff is discriminated between precipitation, glacier melt and snow melt as well as 
 * global sum (see \ref runoff_sums "Runoff sums").
 * 
 * The grid files are written in any format supported by MeteoIO, as configured by the user. This means that it is for example possible to directly write PNG files
 * from Alpine3D. The sno files as well as .met and .pro are written according to the SNOWPACK standalone model documentation.
 * 
 */

/**
 * @page tools Simulation tools
 * Several tools are available to help using Alpine3D. As for SNOWPACK, it is possible to use <A HREF="https://models.slf.ch/p/inishell">inishell</A> to
 * configure the simulations. There is also another java tool, "view" in the "Interface" sub directory, that can be used to visualize ARC ASCII grids as
 * well as to visualize DEM and LUS files in this format. This can also be used to generate a LUS file by opening an aerial picture and manually tagging
 * the pixels (one by one, along lines or within polygons). Finally, this tool can also generate a POI (points of interest) file for more detailed
 * outputs at some specific points.
 * \image html view_tool.png "\"View\" application for visualizing grids"
 * \image latex view_tool.eps "\"View\" application for visualizing grids" width=0.9\textwidth
 */

/**
 * @page coding_style Coding style
 * @section coding_sty Recommended coding style
 * The recommended coding style for MeteoIO is the <A HREF="http://www.kernel.org/doc/Documentation/CodingStyle">Kernel coding style</A> with a few exceptions:
 * - we don't enforce strict 80 characters line width. try to remain reasonable, but don't necessarily cut everything off at 80 characters
 * - try to intelligently use spaces to visually group elements of a complex formula. If the formula can be split into meaningful elements,
 *   please do it (using some "const double element = " constructs).
 * - try to properly qualify variables: for example, if a variable will not be changed, will never be negative and always integer,
 *   then use "const unsigned int". When some specific types are used for some standard library calls, try to properly use these types (for example, "size_t")
 * - use C++ method naming convention: a method name starts with lowercase letter and each individual word in a name starts capitalized.
 *   Usually, no underscores are used in a method. For example, a method that would return the lapse rate contained in an object would be named "getLapseRate()"
 * - qualify variables and parameters with "const" when appropriate (see <A HREF="http://jriddell.org/const-in-cpp.html">const-in-cpp</A>).
 *
 * A few important points to emphasize (from the <A HREF="http://www.kernel.org/doc/Documentation/CodingStyle">Kernel coding style</A>):
 * - Functions should be short and sweet, and do just one thing.  They should fit on one or two screenfuls of text, and do one thing and do that well.
 * - If you have a complex function, and you suspect that a less-than-gifted first-year high-school student might not even understand
 *   what the function is all about, you should adhere to the maximum limits all the more closely.  Use helper functions with descriptive names.
 * - Comments are good, but there is also a danger of over-commenting.  NEVER try to explain HOW your code works in a comment:
 *   it's much better to write the code so that the _working_ is obvious, and it's a waste of time to explain badly written code.
 *
 * @section code_indentation Indentation
 * Since every user has his/her own preference for the ideal indentation width, please use <A HREF="http://www.emacswiki.org/emacs/SmartTabs">"smart tabs"</A>.
 * That practically means:
 * - indent with tabs
 * - align with spaces
 *
 * This way, each developer can set his/her indentation size as he/she wishes without forcing his/her choice to others...
 *
 * @section containers Memory management and Containers
 * Please do NOT manage memory manually but use <A HREF="https://secure.wikimedia.org/wikipedia/en/wiki/Standard_Template_Library">Standard Template Library (STL)
 * </A> <A HREF="http://www.cplusplus.com/reference/stl/">containers</A> instead.
 * This dramatically reduces memory errors (ie: <A HREF="https://secure.wikimedia.org/wikipedia/en/wiki/Segmentation_fault">segfaults</A>), often
 * offers more performance and provides you with lots of <A HREF="http://www.cplusplus.com/reference/algorithm/">associated algorithms</A>
 * (like sorting, search, filling, etc).
 *
 * When you need your own data class, please design it based on these STL containers (like grid2DObject is based on std::vector). Basically, this means
 * that you will replace mallocs and arrays by vectors (for 1d, 2d, 3d grids), maps (for multiple key/value pairs), lists (for unordered table), etc
 *
 * @section exceptions_handling Exceptions handling
 * The recommended C++ usage should be followed: <b>"throw by value, catch by reference"</b> (as specified in <i>C++ Coding Standards: 101 Rules, Guidelines,
 * and Best Practices</i>, Herb Sutter, Andrei Alexandrescu, 2004, Addison-Wesley Professional). Moreover, we should consider catching by
 * <b>const reference</b> and not even declaring a variable if not doing anything with it: something like `catch(const IOException&)` would often be enough.
 *
 */

#endif
