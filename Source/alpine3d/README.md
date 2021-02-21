# Alpine3D

Alpine3D is a spatially distributed (surface), three dimensional (atmospheric) model for analyzing and predicting dynamics of snow-dominated surface processes in mountainous topography. It includes models for snow cover ([SNOWPACK](http://models.slf.ch/p/snowpack/)), vegetation and soil, snow transport, radiation transfer and runoff which can be enabled or disabled on demand. Alpine3D is available under LGPL version 3 or above, see [www.gnu.org](https://www.gnu.org/).

The model supports a variety of input options including interpolation of meteorological weather stations, input from a meteorological model or from remote sensing data ([MeteoIO](http://models.slf.ch/p/meteoio/)) and has been parallelized in order to run on a multi-core computer (with [open-mp](http://openmp.org/) or a cluster (with [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface)).

Alpine3D has a broad variety of potential applications. Most dominant is the assessment of snow water resource dynamics in mountain catchments (Michlmayr et al., 2008). This includes predictions of future snow on the basis of climate change scenarios (Bavay et al., 2009, 2012). One exotic application of the model system Alpine3D is the forecasting of surface temperatures on ski-pistes, e.g. for the Vancouver winter olympics. For this forecast local shadings (might) change surface temperature up to 5 °C.

You can find a more exhaustive documentation on how to compile and run Alpine3D the [wiki page](https://gitlabext.wsl.ch/snow-models/alpine3d/-/wikis/home)

## Compiling Alpine3D

You need to have the following tools on your system:

* cmake (it is highly recommended to also install cmake-curses or cmake-gui);
* a c++ compiler, for example g++;
* if you want to generate the documentation, doxygen;
* if you have latex, this is a plus since it will be able to render some equations in the documentation.
* MeteoIO
* Snowpack

You need to download and install (or download/compile and install) the MeteoIO and Snowpack so they can be automatically found. If you don't want to install either MeteoIO or Snowpack on your system, you have to keep the following directory structure:

                |
                |-----> snowpack
                |-----> meteoio
                |-----> alpine3d

and then Alpine3D would still find MeteoIO and snowpack automatically.

Please note that for Alpine3D-3.0.0, at least MeteoIO-2.5.1 is needed (there is a source-only package specially for this reason on MeteoIO's page).

### Compiling on Linux/Unix

Alpine3D uses the cmake build system: 

```bash
cd alpine3d && mkdir build && cd build
cmake ..
# if you want mpi and/ord openmp support you need to explicitly activate it
# cmake .. -DMPI=ON -DOPENMP=ON
make -j$(nproc) && sudo make install
```
You can list the config options using ```cmake .. -L``` or ```ccmake ..```. The most useful options are MPI and OpenMP to enable support of either of them.

## Running a simulation

It is highly recommended to setup your simulation in a specific directory. For example, a directory named "Dischma" that will have the following subdirectories structure:

* input; There you copy your meteo data in a "meteo" subdirectory (if needed), your domain grids in a "surface-grids" directory and the* sno files in a "snowfiles" directory
* output; The simulation outputs will be written there. You can create a "runoff" subdirectory to contain runoff data
* output/grids; The gridded outputs could go there.
* setup; The configuration files and start scripts should go there. The simulations will be started from this directory, using the run.sh launch script.

For an openmp parallel run, you must run on an [SMP](https://en.wikipedia.org/wiki/Symmetric_multiprocessing) system. For MPI, you must have the proper MPI libraries installed. Your queuing system (if any) must also support openmp/mpi and you must define how many threads you want to use for each run. This is handled by the run.sh launch script, but you must edit this script according to your needs.

If submitting the job to a cluster using the [Sun Grid Engine](https://en.wikipedia.org/wiki/Oracle_Grid_Engine) queuing system, you must edit the header of your run.sh script to define your job name, queue name, number of nodes requested, etc

### Configuring the simulation

This is done by two means: by the io.ini configuration file and by some command line options.

#### The io.ini configuration file

This file is structured by sections, focused on several aspects:

* General, for a few gerenal options
* Input, for the configuration of the inputs
* Output, for the configuration of the outputs
* Snowpack, for the specific snowpack model configuration
* SnowpackAdvanced, for some advanced Snowpack options, including the necessary ALPINE3D = true key
* Interpolations1D, for the temporal interpolations configuration
* Interpolations2D, for the spatial interpolations configuration

These sections are described in the meteoio and snowpack documentation.

It is also possible (and recommended) to use inishell to generate a proper io.ini for Alpine3D.

Unfortunately, the current configuration of inishell for Alpine3D does not cover all the possibilities. You can generate the Snowpack section by using the Snowpack inishell configuration and then copy it into your ini files generated with inishell for Alpine3D.

#### Command line options

You can get the list of supported options by running alpine3d --help. These options focus on which modules should be enabled (for example, snowdrift), the number of workers for the modules that have been parallelized (for example, 4 workers for ebalance) and the start and end date of the simulation.
