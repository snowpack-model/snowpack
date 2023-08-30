#!/bin/bash

# Retrieve script variables
model=$1

# Concatenate files into a single timeseries
# First try year/month:
g=$(ls -d output/${model}_[0-9][0-9][0-9][0-9][0-9][0-9] 2> /dev/null | tail -1)
if [ -z ${g} ]; then
	# If g is empty, try year:
	g=$(ls -d output/${model}_[0-9][0-9][0-9][0-9] 2> /dev/null | tail -1)
	permonth=0
else
	permonth=1
fi
if [ -z ${g} ]; then
	echo "No output files found."
	exit
fi
for f in ${g}/*smet
do
	g=$(basename ${f} .smet)
	if (( ${permonth} )); then
		# Per year/month
		cat output/${model}_[0-9][0-9][0-9][0-9][0-9][0-9]/${g}.smet | awk '{if(/^SMET/) {header=1}; if((data==1 && header==0) || (data==0 && header==1)) {print}; if(/\[DATA\]/) {data=1; header=0}}' > output/${g}.smet
	else
		cat output/${model}_[0-9][0-9][0-9][0-9]/${g}.smet | awk '{if(/^SMET/) {header=1}; if((data==1 && header==0) || (data==0 && header==1)) {print}; if(/\[DATA\]/) {data=1; header=0}}' > output/${g}.smet
	fi
done

# Zip output .smet files
zip output/smet_forcing output/*.smet
