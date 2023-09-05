#!/bin/bash

year_start=1980
year_end=2019
permonth=0		# =1: do data extraction per year/month. =0: do data extraction per year
outputinterval=0	# =0: files are written only after all time steps have been processed.
			# outputinterval>0: output is written after <outputinterval> timesteps, reducing memory consumption (requires APPEND mode enabled in output plugin)
model="MERRA-2"

rm -f to_exec.lst

for yr in $(seq ${year_start} ${year_end})
do
	if (( ${permonth} )); then
		for m in $(seq 1 12)
		do
			mm=$(echo ${m} | awk '{printf "%02d", $1}')
			echo "export OUTPUTINTERVAL=${outputinterval}; bash create_forcing.sh ${model} ${yr} ${mm}" >> to_exec.lst
		done
	else
		echo "export OUTPUTINTERVAL=${outputinterval}; bash create_forcing.sh ${model} ${yr}" >> to_exec.lst
	fi
done
