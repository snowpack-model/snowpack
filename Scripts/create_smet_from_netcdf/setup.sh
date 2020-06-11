#!/bin/bash

year_start=1980
year_end=2019
temporal_res=3600 # In seconds
model="MERRA-2"

rm -f to_exec.lst

for yr in $(seq ${year_start} ${year_end})
do
	echo "bash create_forcing.sh ${model} ${yr} ${temporal_res}" >> to_exec.lst
done
