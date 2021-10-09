# SNOWPACK

SNOWPACK is a multi-purpose snow and land-surface model, which focuses on a detailed description of the mass and energy exchange between the snow, the atmosphere and optionally with the vegetation cover and the soil. It also includes a detailed treatment of mass and energy fluxes within these media.

SNOWPACK has originally been developed to support avalanche warning (Lehning et al., 1999) and thus features a very detailed description of snow properties including weak layer characterization (Stoessel et al., 2009), phase changes and water transport in snow (Hirashima et al., 2010). A particular feature is the treatment of soil and snow as a continuum with a choice of a few up to several hundred layers. While a main application is still on avalanche warning in countries from Switzerland (Schirmer et al., 2009) to Japan (Nishimura et al., 2005), the applications range from climate change assessments (Rasmus et al., 2004; Bavay et al., 2009) and superimposed ice simulations (Obleitner and Lehning, 2004) to permafrost sensitivity studies (Luetschg et al., 2008) and the simulation of snow storage (Olefs and Lehning, 2010).

In order to ease the integration of SNOWPACK into other models, it is now structured as a library (libsnowpack) and an application that uses the library to perform simulations (snowpack). This library is available under LGPL version 3 or above, see www.gnu.org.

If you want to use Snowpack in an operational context, you can request to have access to "snowpack-opera" that gives you the whole operational toolchain, also as Open Source, as well as some documentation on how to setup an operational system.


## Getting SNOWPACK

In order to reduce the maintenance burden and to make SNOWPACK easier to tailor to specific needs, it has been decomposed into several tools:

* inishell for graphical configuration of the simulations. Please use it to configure your simulations, it makes it much easier!
* MeteoIO for data retrieval and preprocessing
* libsnowpack for doing the core computation and snowpack for calling libsnowpack for a point simulation (ie. what most people want out of SNOWPACK)
* niViz for simulation output visualization (or go to https://niviz.org)

The Snowpack releases that are found in the Downloads (currently at https://models.slf.ch/p/snowpack/downloads/) bundle MeteoIO and Snowpack but NOT Inishell that has to be downloaded separately (at https://models.slf.ch/p/inishell-ng/downloads/). They are available for Linux, Windows and osX and contain an HTML documentation (either in the start menu or in the share subfolder of the installation directory).

Otherwise, if you recompile SNOWPACK from sources (instructions provided at https://models.slf.ch/p/snowpack/page/Getting-started/), you need to download and install a few other of these package. At the minimum, you need to have MeteoIO on your system, compiled and preferably installed.

## Running SNOWPACK

Please follow the instructions at https://models.slf.ch/p/snowpack/page/Running-Snowpack/ that also show how to run the provided examples. There is also some documentation at https://models.slf.ch/ in the "Documentation server" section.
