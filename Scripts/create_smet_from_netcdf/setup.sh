#!/bin/bash

# Retrieve script variables
model=$1
yr=$2
output_res=$3 # In seconds

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
sed 's/VSTATIONS_REFRESH_RATE=3600/VSTATIONS_REFRESH_RATE='${output_res}'/' ${model}.ini > ${model}_${yr}.ini 
sed -i 's/VSTATIONS_REFRESH_OFFSET=3600/VSTATIONS_REFRESH_OFFSET='${output_res}'/' ${model}_${yr}.ini
sed -i 's/PSUM::arg1::cst=3600/PSUM::arg1::cst='${output_res}'/' ${model}_${yr}.ini
sed -i 's/METEOPATH=..\/output\//METEOPATH=..\/output\/'${model}'_'${yr}'\//' ${model}_${yr}.ini

# Create .smet files for each year
#../data_converter ${yr}-01-01T00:30:00 ${yr}-12-31T23:30:00 60 ${model}_${yr}.ini > ../log/${model}_${yr}.log 2>&1
../data_converter ${yr}-01-01T00:30:00 ${yr}-01-01T06:30:00 60 ${model}_${yr}.ini > ../log/${model}_${yr}.log 2>&1
