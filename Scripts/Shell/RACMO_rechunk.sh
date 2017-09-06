#!/bin/bash
set -ex

# Description: (Re-)chunk RACMO-files
#
# Author: Christian Steger, October 2015

cd /scratch/ms/nl/rucs/raw_Data/3hourly # set name to input file directory
#cd /Users/steger/Data/Greenland/RACMO2.3/3hourly 
#cd /Volumes/IMAU/Greenland/RACMO2.3/3hourly

##########################################################################################
# (Re-)chunking
##########################################################################################

# Loop through files matching the pattern
for FILE in `ls *.KNMI-*.ZGRN11.BN_1958_2013.3H.nc`
do
	# from (1, 1, 312, 306) to (100, 1, 4, 306); only chunk variables with at least 3 dimensions
	ncks --cnk_plc=cnk_g3d --cnk_dmn=time,100 --cnk_dmn=rlat,4  $FILE ./${FILE%.nc}_rc.nc # re-chunk
	# ECMWF (--cnk_plc=cnk_g3d causes error)
	#ncks --cnk_plc=all --cnk_dmn=time,100 --cnk_dmn=rlat,4  $FILE ./${FILE%.nc}_rc.nc # re-chunk
    #rm $FILE # remove old file
done

exit


