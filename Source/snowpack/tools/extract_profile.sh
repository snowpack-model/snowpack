#!/bin/bash
# This script extract a vertical profile for a specific time stamp, for temperature, density, grain size, and volumetric ice, water and air content.
#
# Use extract_profile.sh <pro-file> <date> [absolute]
#  <pro file>: specify *.pro file
#  <date>    : specify date to extract profile for. The script takes the first profile on or after this date to produce the output.
#  <absolute>: [optional] if absolute is set, absolute z-values are provided (height increasing upward), if absolute is omitted, the snow surface is set at 0, with depth increasing downward
#
#  Output: writes profile with depth (i.e., z=0 is surface) to stdout
#
#  Usage examples:
#  ===============
#
#  	bash extract_profile.sh WFJ2.pro 2016-08-01
#  	bash extract_profile.sh WFJ2.pro 2016-02-01T12:15 absolute
#
#
#  Integration example in python, read output in numpy array:
#  ==========================================================
#	# Example 1: directly parsing *.pro file
#	import os;
#	import numpy as np;
#	numpy_array = np.array(os.popen('bash extract_profile.sh WFJ2.pro 2016-08-01').read())
#
#	# Example 2: unzipping an archive with a *.pro file
#	import os;
#	import numpy as np;
#	numpy_array = np.array(os.popen('bash -c "bash extract_profile.sh <(unzip -p zip/WFJ2.zip *pro) 2016-08-01"').read())
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
if [ -n "$3" ]; then
	absolute=1
else
	echo "Specify date on command line!"
fi

awk -F, -v tdate=${d} -v absolute=${absolute} '\
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
			if(substr(datum,1,16) >= substr(tdate,1,16)) {found_date=1; print "# Date=", datum} \
		}; \
		if (found_date==1) { \
			if($1==501) { \
				if($3 > 0) { \
					# Insert z=0 at the bottom. This is necessary in case no soil is present, when z=0 is not printed in the *pro file
					bottomsnowelement=1; z[1]=0; offset=1; nE=$2; nEsnow=nE; \
				} else { \
					# Case with negative z: either soil layers or sea ice
					offset=0; nE=$2-1; \
				}; \
				for(i=1;i<=$2;i++) { \
					# Read domain coordinates
					z[i+offset]=$(i+2)/100.; \
					if(z[i+offset]==0) {bottomsnowelement=i} \
				}; \
			} else if($1==502) { \
				# Read densities
				for(i=1; i<=$2; i++) {rho[i]=$(i+2)} \
			} else if($1==503) { \
				# Read temperatures
				for(i=1; i<=$2; i++) {Te[i]=$(i+2)+273.15} \
			} else if($1==505) { \
				# Read layer ages
				for(i=1; i<=$2; i++) {age[i]=$(i+2)} \
			} else if($1==506) { \
				# Read LWC
				for(i=1; i<=$2; i++) {th_water[i]=$(i+2)/100.} \
			} else if($1==508) { \
				# Read dendricity
				for(i=1; i<=$2; i++) {dd[i]=$(i+2)} \
			} else if($1==509) { \
				# Read sphericity
				for(i=1; i<=$2; i++) {sp[i]=$(i+2)} \
			} else if($1==511) { \
				# Read bond size
				for(i=1; i<=$2; i++) {bs[i]=$(i+2)} \
			} else if($1==512) { \
				# Read grain size
				for(i=1; i<=$2; i++) {gs[i]=$(i+2)} \
			} else if($1==513) { \
				# Read grain type
				for(i=1; i<=$2; i++) {gt[i]=$(i+2)} \
			} else if($1==515) { \
				# Read theta[ICE]
				for(i=1; i<=$2; i++) {th_ice[i]=$(i+2)/100.} \
			} else if($1==516) { \
				# Read theta[AIR] (i.e., pore space)
				for(i=1; i<=$2; i++) {th_air[i]=$(i+2)/100.} \
			} else if($1==519) { \
				# When there is no soil
				soil=0; nEsnow=nE; \
				# Check theta[SOIL], then adjust parameters such that only snow part is shown
				for(i=1; i<bottomsnowelement; i++) {if($(i+2)/100. > 0) {soil=1; nEsnow=nE-bottomsnowelement+1; break;}} \
				# Print header
				if(printheader==0) {printheader=1; \
					if(!absolute) { \
						print "# depth_top_(m)   depth_bottom_(m)   depth_mid_(m)   thickness_(m)   temperature_(K)   density_(kg/m^3)   grain_size_(mm)   bond_size_(mm)   dd_(-)   sp_(-)   gt_(swiss_code_F1F2F3)   theta_ice_(m^3/m^3)   theta_water_(m^3/m^3)   theta_air_(m^3/m^3)   age_(days)"; \
					} else { \
						print "# height_top_(m)  height_bottom_(m)  height_mid_(m)  thickness_(m)   temperature_(K)   density_(kg/m^3)   grain_size_(mm)   bond_size_(mm)   dd_(-)   sp_(-)   gt_(swiss_code_F1F2F3)   theta_ice_(m^3/m^3)   theta_water_(m^3/m^3)   theta_air_(m^3/m^3)   age_(days)"; \
					} \
				} \
				# Write output
				for(i=nE; i>=((soil)?(bottomsnowelement):(1)); i--) {
					if(!absolute) {
						printf "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", z[nE+1]-z[i+1], z[nE+1]-z[i], z[nE+1]-0.5*(z[i+1]+z[i]), (z[i+1]-z[i]), Te[i], rho[i], gs[i], bs[i], dd[i], sp[i], gt[i], th_ice[i], th_water[i], th_air[i], age[i]; \
					} else {
						printf "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", z[i+1], z[i], 0.5*(z[i+1]+z[i]), (z[i+1]-z[i]), Te[i], rho[i], gs[i], bs[i], dd[i], sp[i], gt[i], th_ice[i], th_water[i], th_air[i], age[i]; \
					}
				} \
				exit; \
			} \
		} \
	} \
}' $1
