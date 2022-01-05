#!/bin/awk
#
# This script checks if the spinup is complete, requesting a minimum firn depth of 10m, and the bottom 3m should be 910 kg/m3 or more (the settling cut-off)
#
BEGIN {
	th_ice_max=0.992366412;	# 910 kg/m3 in volumetric ice contents
	minicedepth=3;		# Min depth at bottom of firn that should equal or exceed th_ice_max
	mindepth=10;		# Min depth of firn layer
	data=0;
	spinup=-1;
}
{
	if(data==1) {
		Htot+=$2;
		if($4>=th_ice_max && spinup==-1) {
			H+=$2;
			if(H>=minicedepth) {
				spinup=1
			}
		} else {
			if (spinup == -1) {
				spinup=0;
				exit;
			}
		}
	}
	if(/\[DATA\]/) {
		# New data frame
		data=1
	}
}
END {
	print (spinup==1 && Htot>=mindepth)?("1"):("0");
	exit;
}
