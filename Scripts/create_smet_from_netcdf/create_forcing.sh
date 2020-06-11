#!/bin/bash

# Retrieve script variables
model=$1
yr=$2
output_res_in_seconds=$3						# In seconds
output_res_in_minutes=$(echo ${output_res_in_seconds} / 60 | bc)	# In minutes

# Print Statements
echo "Working on model: ${model}"
echo "	Output Temporal Resolution: ${output_res} seconds"
echo "	Year ${yr}"

# Module loading and basic setup
mkdir -p log
export LD_LIBRARY_PATH=$(pwd)/usr/lib/:${LD_LIBRARY_PATH}
module purge
ml intel; ml proj; ml netcdf

# Create output directory 
rm -r output/${model}_${yr}
mkdir -p output/${model}_${yr}

# Move into io_files dir, where computation will be peformed
cd io_files

# Create new ini file for each year
sed 's/VSTATIONS_REFRESH_RATE=3600/VSTATIONS_REFRESH_RATE='${output_res_in_seconds}'/' ${model}.ini > ${model}_${yr}.ini
sed -i 's/VSTATIONS_REFRESH_OFFSET=3600/VSTATIONS_REFRESH_OFFSET='${output_res_in_seconds}'/' ${model}_${yr}.ini
sed -i 's/PSUM::arg1::cst=3600/PSUM::arg1::cst='${output_res_in_seconds}'/' ${model}_${yr}.ini
sed -i 's/METEOPATH=..\/output\//METEOPATH=..\/output\/'${model}'_'${yr}'\//' ${model}_${yr}.ini

# Create .smet files for each year
if [ ${model} == "MERRA-2" ]; then
	# MERRA-2 has output each :30
	../meteoio_timeseries -b ${yr}-01-01T00:30:00 -e ${yr}-12-31T23:30:00 -c ${model}_${yr}.ini -s ${output_res_in_minutes} > ../log/${model}_${yr}.log 2>&1
else
	# Other models have output each :00, and the following line assumes the last time step is 21:00 (i.e., 3-hourly output)
	../meteoio_timeseries -b ${yr}-01-01T00:00:00 -e ${yr}-12-31T21:00:00 -c ${model}_${yr}.ini -s ${output_res_in_minutes} > ../log/${model}_${yr}.log 2>&1
fi
