Quick overview of the simulation examples:

1) Introduction
These examples are provided as an easy way to see how the most basic things can be done with SNOWPACK. They are provided with
initial profiles, meteorological data, configuration files and a run script.

2) Organization of the examples
The necessary data set for these examples lies in "input" directory. The files ".smet" or ".inp" contain meteorological time series.
The files ".sno" contain the initial snow profiles in SMET format while the files ".snoold" contain these profiles in the legacy SNOWPACK format.
The file ".meta" contains the metadata that can otherwise not be provided in the legacy "inp" meteo file format.
See in the html documnetation of SNOWPACK and MeteoIO which files format are supported, how to specify which file format you want to use, which
options exist for each file format and how to easily add support for extra file formats.

The simulations configurations files are contained in the "cfgfiles" directory. See in the html documentation how this configuration file is
structured.

Very simple shell scripts call SNOWPACK with the necessary options. the last line of the script actually runs SNOWPACK and this last line
can be directly copied to the command line.

Finally, the examples are written to the "output" directory.

3) Available examples
	run_antartic.sh: simulation in the antarctic at dome C
	run_res1exp.sh: simulation in research mode on the flat
	run_res5exp.sh: simulation in research mode on the flat as well as on 4 virtual slopes (N, E, S, W)
	run_operMST96.sh: simulation in operational mode on the flat
	run_soil.sh: simulation with soil


