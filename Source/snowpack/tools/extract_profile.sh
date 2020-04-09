#!/bin/bash
# This script extract a vertical profile for a specific time stamp, for temperature, density, grain size, and volumetric ice, water and air content.
#
# Use extract_profile.sh <pro-file> <date>
#  <pro file>: specify *.pro file
#  <date>    : specify date to extract profile for. The script takes the first profile on or after this date to produce the output.
#
#  Output: writes profile with depth (i.e., z=0 is surface) to stdout
#
#  Usage examples:
#  ===============
#
#  	bash extract_profile.sh WFJ2.pro 2016-08-01
#
#
#  Integration example in python, read output in numpy array:
#  ==========================================================
#	import os;
#	import numpy as np;
#	numpy_array = np.array(os.popen('bash extract_profile.sh WFJ2.pro 2016-08-01').read())
#
#
#  Integration in R, read output into data frame:
#  ==============================================
#	data_frame <- system("bash extract_profile.sh WFJ2.pro 2016-08-01", intern = TRUE)
#


# Check if date is specified as command line parameter
if [ -n "$2" ]; then
	d=$2
else
	echo "Specify date on command line!"
fi


awk -F, -v tdate=${d} '\
BEGIN { \
	meta=0; \
	printheader=0; \
	start_reading=0; \
	found_date=0; \
} \
{ \
	if($1=="[STATION_PARAMETERS]") {start_reading=0}; \
	if(/Latitude/) {lat=$NF}; \
	if(/Longitude/) {lon=$NF}; \
	if(/Altitude/) {alt=$NF}; \
	if($1=="[DATA]") {start_reading=1; if(meta==0) {meta=1; print "# File= ", FILENAME; print "#", lat; print "#", lon; print "#", alt}}; \
	if(start_reading==1) { \
		if($1==500) { \
			datum=sprintf("%04d-%02d-%02dT%02d:%02d:%02d", substr($2,7,4), substr($2,4,2), substr($2,1,2), substr($2,12,2), substr($2,15,2), substr($2,18,2)); \
			if(substr(datum,1,16) > substr(tdate,1,16)) {found_date=1; print "# Date=", datum} \
		}; \
		if (found_date==1) { \
			if($1==501) { \
				if($3 > 0) { \
					# This is necessary in case no soil is present
					bottomsnowelement=1; z[1]=0; offset=1; nE=$2; nEsnow=nE; \
				} else { \
					# Case with soil layers
					offset=0; nE=$2-1; \
				}; \
				for(i=1;i<=$2;i++) { \
					# Read domain coordinates
					z[i+offset]=$(i+2)/100.; \
					if(z[i+offset]==0) {bottomsnowelement=i} \
				}; nEsnow=nE-bottomsnowelement+1; \
			} else if($1==502) { \
				# Read densities
				for(i=1; i<=$2; i++) {rho[i]=$(i+2)} \
			} else if($1==503) { \
				# Read temperatures
				for(i=1; i<=$2; i++) {Te[i]=$(i+2)+273.15} \
			} else if($1==506) { \
				# Read LWC
				for(i=1; i<=$2; i++) {th_water[i]=$(i+2)/100.} \
			} else if($1==512) { \
				# Read grain size
				for(i=1; i<=$2; i++) {gs[i]=$(i+2)} \
			} else if($1==515) { \
				# Read theta[ICE]
				for(i=1; i<=$2; i++) {th_ice[i]=$(i+2)/100.} \
			} else if($1==516) { \
				# Read theta[AIR] (i.e., pore space)
				for(i=1; i<=$2; i++) {th_air[i]=$(i+2)/100.} \
				# Print header
				if(printheader==0) {printheader=1; print "# depth_top_(m)   depth_bottom_(m)   depth_mid_(m)   thickness_(m)   temperature_(K)   density_(kg/m^3)   grain_size_(mm)   theta_ice_(m^3/m^3)   theta_water_(m^3/m^3)   theta_air_(m^3/m^3)"} \
				# Write output
				for(i=nEsnow; i>=bottomsnowelement; i--) {
					printf "%f %f %f %f %f %f %f %f %f %f\n", z[nEsnow+1]-z[i+1], z[nEsnow+1]-z[i], z[nEsnow+1]-0.5*(z[i+1]+z[i]), (z[i+1]-z[i]), Te[i], rho[i], gs[i], th_ice[i], th_water[i], th_air[i]; \
				} \
				exit; \
			} \
		} \
	} \
}' $1
