#!/bin/bash
#This runs a very simple Alpine3D simulation (reduced DEM compared to Stillberg) and
#compares the results with a reference dataset

# Print a special line to prevent CTest from truncating the test output
printf "CTEST_FULL_OUTPUT (line required by CTest to avoid output truncation)\n\n"

rm -f output/5_*; rm -f output/grids/20*; rm -f output/snowfiles/*

#setup the simulation
BEGIN="2014-10-01T01:00"
END="2014-12-31T00:00"
PROG_ROOTDIR=../../bin
export DYLD_FALLBACK_LIBRARY_PATH=${PROG_ROOTDIR}:${DYLD_FALLBACK_LIBRARY_PATH}	#for osX
export LD_LIBRARY_PATH=${PROG_ROOTDIR}:${LD_LIBRARY_PATH}	#for Linux

#now run the simulation
date
../../bin/alpine3d --iofile=./io.ini --enable-eb --np-ebalance=2 --np-snowpack=2 --startdate=${BEGIN} --enddate=${END} > stdouterr.log 2>&1
ret=$?
date
if [ "$ret" -eq "0" ]; then
	echo "Done Alpine3D Simulation"
else
	echo "fail : Alpine3D did not complete properly! Return code=$ret"
fi

#stop here when re-generating the reference files
#exit

#Compare the results with the reference data
PREC="1e-3"
rm -f output_ref/5_2_dischma.met
bunzip2 -k output_ref/5_2_dischma.met.bz2
sed -i '11d' output_ref/5_2_dischma.met; sed -i '11d' output/5_2_dischma.met
numdiff -s', \t\n' -r ${PREC} --speed-large-files output_ref/5_2_dischma.met output/5_2_dischma.met | grep "+++"
rm -f output_ref/5_2_dischma.met

rm -f output_ref/5_2_dischma.pro
bunzip2 -k output_ref/5_2_dischma.pro.bz2
sed -i '10d' output_ref/5_2_dischma.pro; sed -i '10d' output/5_2_dischma.pro
numdiff -s', \t\n' -r ${PREC} --speed-large-files output_ref/5_2_dischma.pro output/5_2_dischma.pro | grep "+++"
rm -f output_ref/5_2_dischma.pro

rm -f output_ref/5_2_dischma.smet
bunzip2 -k output_ref/5_2_dischma.smet.bz2
sed -i '16,17d' output_ref/5_2_dischma.smet; sed -i '16,17d' output/5_2_dischma_meteo.smet
numdiff -s', \t\n' -r ${PREC} --speed-large-files output_ref/5_2_dischma.smet output/5_2_dischma_meteo.smet | grep "+++"
rm -f output_ref/5_2_dischma.smet

for fichier in $(ls output_ref/grids/2014*); do
	name=$(basename ${fichier})
	numdiff -s', \t\n' -r ${PREC} --speed-large-files output_ref/grids/${name} output/grids/${name} | grep "+++"
done

#cleanup, but keep the POI for further testing
rm -f output/grids/20*; rm -f output/snowfiles/*
