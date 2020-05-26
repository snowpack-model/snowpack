#!/bin/bash
# This script extracts time series of average values for temperature, density, grain size, and volumetric ice, water and air content.
#
# Use extract_params.sh <pro-file> [d=<depth>]      - OR -      extract_params.sh <pro-file> [r=<density>]
#  <pro file>: specify *.pro file
#  <depth>   : specify depth from surface (m) over which calculations are done. If <depth> is negative, it calculates averages over the whole snow column.
#  <density> : specify bulk density threshold (kg/m^3). Searching from bottom up, calculations start when the bulk density drops below this threshold.
#              If <density> is negative, it calculates averages over the whole snow column.
#  If neither depth nor density is specified, values for the whole column are calculated.
#
#  Output: writes time series of date, depth, avg. temperature (K), avg. density (kg/m^3), avg. grain size (mm), avg. theta_ice (m^3/m^3), avg. theta_water (m^3/m^3), avg. theta_air (m^3/m^3)"} to stdout
#          Averages are determined as layer weighted averages.
#
#  Usage examples:
#  ===============
#
#  	bash extract_params.sh WFJ2.pro d=1 > WFJ2_ts.txt	# Extract variables over uppermost 1 m.
#  	bash extract_params.sh WFJ2.pro r=500 > WFJ2_ts.txt	# Extract variables over part where the bulk density is less than 500 kg/m^3.
#
#
#  Integration example in python, read output in numpy array:
#  ==========================================================
#	import os;
#	import numpy as np;
#	numpy_array = np.array(os.popen('bash extract_params.sh WFJ2.pro d=1').read())
#
#
#  Integration example in R, read output into data frame:
#  ==============================================
#	data_frame <- system("bash extract_params.sh WFJ2.pro d=1", intern = TRUE)
#


# Set default values
d=-1
r=-1


# Read command line variables
for var in "$@"; do
	if [[ "$var" == *"="* ]]; then
		eval ${var}
	fi
done


awk -F, -v depth=${d} -v density=${r} '\
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
		# Print header
		if(printheader==0) {printheader=1; print "# timestamp   depth_(m)   average_temperature_(K)   average_density_(kg/m^3)   average_grain_size_(mm)   average_theta_ice_(m^3/m^3)   average_theta_water_(m^3/m^3)   average_theta_air_(m^3/m^3)"} \
		if($1==500) { \
			datum=sprintf("%04d-%02d-%02dT%02d:%02d:%02d", substr($2,7,4), substr($2,4,2), substr($2,1,2), substr($2,12,2), substr($2,15,2), substr($2,18,2)); \
			p=0; \
		} else if($1==501) { \
			if($2==1 && $3==0) {printf "%s -999 -999 -999 -999 -999 -999\n", datum;} \
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
			# Write output
			H=0; \
			Te_sum=-999; \
			rho_sum=-999; \
			gs_sum=-999; \
			th_ice_sum=-999; \
			th_water_sum=-999; \
			th_air_sum=-999; \
			for(i=bottomsnowelement; i<=nEsnow; i++) { \
				if( ((z[nEsnow+1]-z[i+1])<=depth || depth<0) \
				      && (p==1 || rho[i]<density || density<0)) { \
					p=1; \
					dH=(z[i+1]-z[i]); \
					H+=dH; \
					if(Te_sum       == -999) {Te_sum       = Te[i]*dH      } else {Te_sum       += Te[i]*dH      }; \
					if(rho_sum      == -999) {rho_sum      = rho[i]*dH     } else {rho_sum      += rho[i]*dH     }; \
					if(gs_sum       == -999) {gs_sum       = gs[i]*dH      } else {gs_sum       += gs[i]*dH      }; \
					if(th_ice_sum   == -999) {th_ice_sum   = th_ice[i]*dH  } else {th_ice_sum   += th_ice[i]*dH  }; \
					if(th_water_sum == -999) {th_water_sum = th_water[i]*dH} else {th_water_sum += th_water[i]*dH}; \
					if(th_air_sum   == -999) {th_air_sum   = th_air[i]*dH  } else {th_air_sum   += th_air[i]*dH  }; \
				} \
			} \
			printf "%s %f %f %f %f %f %f %f\n", datum, H, (H>0)?(Te_sum/H):(-999), (H>0)?(rho_sum/H):(-999), (H>0)?(gs_sum/H):(-999), (H>0)?(th_ice_sum/H):(-999), (H>0)?(th_water_sum/H):(-999), (H>0)?(th_air_sum/H):(-999); \
		} \
	} \
}' $1
