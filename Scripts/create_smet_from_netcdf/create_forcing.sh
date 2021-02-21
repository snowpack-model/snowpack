#!/bin/bash
shopt -s expand_aliases	# Make sure aliases work in non-interactive shells

# Retrieve script variables
model=$1
yr=$2

if [ -z "${model}" ] && [ -z "${year}" ]; then
	echo "Provide model and year as first and second command line parameter, respectively."
	exit
fi

if [ ${model} != "MERRA-2" ] && [ ${model} != "CESM" ] && [ ${model} != "RACMO2" ]; then
	echo "Only MERRA-2, CESM or RACMO2 supported as model."
	exit
fi

# Print Statements
echo "Working on model: ${model}"
echo "	Year ${yr}"

# Module loading and basic setup
mkdir -p log
export LD_LIBRARY_PATH=$(pwd)/usr/lib/:${LD_LIBRARY_PATH}
module purge
ml intel; ml proj; ml netcdf

# Create output directory
rm -r output/${model}_${yr}
mkdir -p output/${model}_${yr}

# Retrieve directory with netCDF and create file links, if netCDF files are split by year
flag_netcdf_files_are_linked=0
grid2dpath=$(fgrep GRID2DPATH ./io_files/MERRA-2.ini | mawk -F= '{gsub(/^[ \t]+/,"",$NF); gsub(/[ \t]+$/,"",$NF); print $NF}')	# Use gsub to remove trailing and leading white spaces
if [ -e ${grid2dpath}/*${yr}*nc ]; then
	# Create dir to link the NetCDF files
	flag_netcdf_files_are_linked=1
	rm -r input/${model}_${yr}
	mkdir -p input/${model}_${yr}
	for ncf in ${grid2dpath}/*${yr}*nc; do
		ln -s ${ncf} input/${model}_${yr}/
	done
fi

# Move into io_files dir, where computation will be peformed
cd io_files

# Create new ini file for each year
inifile=${model}_${yr}.ini
sed 's/..\/output\//..\/output\/'${model}'_'${yr}'\//' ${model}.ini > ${inifile}

# Set paths correctly in case we linked the netcdf files
if (( ${flag_netcdf_files_are_linked} )); then
	sed -i 's/.*METEOPATH.*/GRID2DPATH = ../input/${model}_${yr}' ${inifile}
	sed -i 's/.*GRID2DPATH.*/GRID2DPATH = ../input/${model}_${yr}' ${inifile}
fi

# Retrieve temporal resolution from ini file:
output_res_in_seconds=$(fgrep VSTATIONS_REFRESH_RATE ${inifile} | awk -F\= '{i=$NF; gsub(/[ \t]/, "", i); print i}')	# In seconds
if [ -z "${output_res_in_seconds}" ]; then
	echo "Unable to determine temporal resolution. Is VSTATIONS_REFRESH_RATE defined in *.ini file?"
	exit
fi
output_res_in_minutes=$(echo ${output_res_in_seconds} / 60 | bc)							# In minutes

# Create .smet files for each year
if [ ${model} == "MERRA-2" ]; then
	# MERRA-2 has output each :30
	../meteoio_timeseries -b ${yr}-01-01T00:30:00 -e ${yr}-12-31T23:30:00 -c ${model}_${yr}.ini -s ${output_res_in_minutes} -p > ../log/${model}_${yr}.log 2>&1
else
	# Other models have output each :00, and the following line assumes the last time step is 21:00 (i.e., 3-hourly output)
	../meteoio_timeseries -b ${yr}-01-01T00:00:00 -e ${yr}-12-31T21:00:00 -c ${model}_${yr}.ini -s ${output_res_in_minutes} -p > ../log/${model}_${yr}.log 2>&1
fi
