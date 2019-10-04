Quick overview of the code examples:

1) Introduction
These examples are provided as an easy way to see how the most basic things can be done with MeteoIO. They have been
designed with simplicity, compacity and clarity in mind, skipping most of the error management code that should normally
be present (like checking the number of parameters on the command line). They are commented (almost) line by line in
order to help understand all the steps.

2) Organization of the examples
The necessary data set for these examples lies in "input" directory. It contains a Digital Elevation Model as well as a few weeks of data for several stations: from 2008-12-01T00:00 until 2009-01-31T23:00. Moreover, a basic Makefile is provided for gcc. You might also have to adjust the location of the library in the Makefile as well as in your terminal (by typing something like "export LD_LIBRARY_PATH = ../../lib:$LD_LIBRARY_PATH").

3) Available examples
	2D_interpolations.cc : read meteo data and interpolate point measurements over a whole Digital Elevation Model
	coordinates.cc : coordinate conversions
	data_converter.cc : read data, process it as specified in the configuration file and write it out (potentially
	                    in a different file format, therefore converting the data to another format)
	dem_reading.cc : read a DEM, show some statistics about it and write it back (potentially converting it)
	grid2d_reading.cc: read a grid containing temperatures, show it and perform some computation on the grid
	matrix.cc : basic matrix arithmetic
	meteo_reading.cc : read one timestep and shows the meteo data at this step
	sun.cc : calculate solar radiation at a given place and time
	time.cc : time handling example
	random_numbers.cc : generating different kinds of random numbers
	statistical_filters.cc : understanding Bayesian filters input (Kalman and particle filter)

4) Available scripts
Several useful scripts are given in the "tools" subdirectory of MeteoIO's root sources directory. These scripts allow to 
extract statictics from SMET input files, extract a single parameter, plot SMET files in xmGrace or in Matlab. There is also
a script to convert the locations given in SMET files to a KML file that can be opened in a GIS or online mapping service 
(such as google maps).
