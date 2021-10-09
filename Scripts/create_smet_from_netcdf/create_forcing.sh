#!/bin/bash
export TZ=UTC
shopt -s expand_aliases	# Make sure aliases work in non-interactive shells
# Check if mawk exist, otherwise create alias
if ! command -v mawk &> /dev/null
then
	alias mawk='awk'
fi

# Retrieve script variables
model=$1
yr=$2
if [ ! -z "${3}" ]; then
	mm=$(echo $3 | mawk '{printf "%02d", $1}')
else
	mm=""
fi

# Do some checks on provided arguments:
if [ -z "${model}" ] && [ -z "${yr}" ]; then
	echo "Provide model and year as first and second command line arguments, respectively. Optional: provide month as third argument."
	exit
fi

if [ ${model} != "MERRA-2" ] && [ ${model} != "CESM" ] && [ ${model} != "RACMO2" ]; then
	echo "Only MERRA-2, CESM or RACMO2 supported as model."
	exit
fi

# Print Statements
echo "Working on model: ${model}"
if [ -z "${mm}" ]; then
	echo "	Year ${yr}"
else
	echo "	Year ${yr}-${mm}"
fi

# Module loading and basic setup
mkdir -p log
export LD_LIBRARY_PATH=$(pwd)/usr/lib/:${LD_LIBRARY_PATH}
module purge
ml intel; ml proj; ml netcdf

# Create output directory
rm -rf output/${model}_${yr}${mm}
mkdir -p output/${model}_${yr}${mm}

# Move into io_files dir, where computation will be peformed
cd io_files

# Retrieve directory with netCDF and create file links, if netCDF files are split by year
flag_netcdf_files_are_linked=0
grid2dpath=$(fgrep GRID2DPATH ./${model}.ini | mawk -F= '{gsub(/^[ \t]+/,"",$NF); gsub(/[ \t]+$/,"",$NF); print $NF}')	# Use gsub to remove trailing and leading white spaces
list_of_nc_files=$(find ${grid2dpath} -name "*${yr}*")
if [ ! -z "${list_of_nc_files}" ]; then
	# Create dir to link the NetCDF files
	flag_netcdf_files_are_linked=1
	rm -rf ../input/${model}_${yr}${mm}
	mkdir -p ../input/${model}_${yr}${mm}
	for ncf in ${list_of_nc_files}; do
		ln -sr ${ncf} ../input/${model}_${yr}${mm}/
	done
fi

# Create new ini file for each year
inifile=${model}_${yr}${mm}.ini
sed 's/..\/output\//..\/output\/'${model}'_'${yr}${mm}'\//' ${model}.ini > ${inifile}

# Set paths correctly in case we linked the netcdf files
if (( ${flag_netcdf_files_are_linked} )); then
	cp ${inifile} ${inifile}.tmp
	# Replace the paths in the [Input] section only:
	mawk -v model=${model} -v yr=${yr} -v mm="${mm}" '{if(/^\[/) {$0=toupper($0); if(/\[INPUT\]/) {input=1} else {input=0}}; if(input==1 && /^METEOPATH/) {print "METEOPATH = ../input/" model  "_" yr mm "/"} else if(input==1 && /^GRID2DPATH/) {print "GRID2DPATH = ../input/" model  "_" yr mm "/"} else {print $0}}' ${inifile}.tmp > ${inifile}
	rm ${inifile}.tmp
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
	if [ -z "${mm}" ]; then
		# When no month is provided and we do the processing per year
		bt="${yr}-01-01T00:30:00"
		et=$(echo ${yr} ${output_res_in_minutes} | mawk '{d=mktime(sprintf("%04d %02d %02d %02d %02d %02d 0", $1+1, 1, 1, 0, 30, 0, 0)); printf(strftime("%Y-%m-%dT%H:%M:%S", d-($3*60)), $1)}')
	else
		# When a month is provided and we do the processing per year/month
		bt=$(echo ${yr} ${mm} | mawk '{printf("%04d-%02d-01T00:30:00", $1, $2)}')
		et=$(echo ${yr} ${mm} ${output_res_in_minutes} | mawk '{d=mktime(sprintf("%04d %02d %02d %02d %02d %02d 0", ($2==12)?($1+1):($1), ($2==12)?(1):($2+1), 1, 0, 30, 0, 0)); printf(strftime("%Y-%m-%dT%H:%M:%S", d-($3*60)), $1)}')
	fi
	../meteoio_timeseries -b ${bt} -e ${et} -c ${model}_${yr}${mm}.ini -s ${output_res_in_minutes} -p > ../log/${model}_${yr}${mm}.log 2>&1
else
	# Other models have output each :00, and the following line assumes the last time step is 21:00 (i.e., 3-hourly output)
	if [ -z "${mm}" ]; then
		# When no month is provided and we do the processing per year
		bt="${yr}-01-01T00:00:00"
		et=$(echo ${yr} ${output_res_in_minutes} | mawk '{d=mktime(sprintf("%04d %02d %02d %02d %02d %02d 0", $1+1, 1, 1, 0, 0, 0, 0)); printf(strftime("%Y-%m-%dT%H:%M:%S", d-($3*60)), $1)}')
	else
		# When a month is provided and we do the processing per year/month
		bt=$(echo ${yr} ${mm} | mawk '{printf("%04d-%02d-01T00:00:00", $1, $2)}')
		et=$(echo ${yr} ${mm} ${output_res_in_minutes} | mawk '{d=mktime(sprintf("%04d %02d %02d %02d %02d %02d 0", ($2==12)?($1+1):($1), ($2==12)?(1):($2+1), 1, 0, 0, 0, 0)); printf(strftime("%Y-%m-%dT%H:%M:%S", d-($3*60)), $1)}')
	fi
	../meteoio_timeseries -b ${bt} -e ${et} -c ${model}_${yr}${mm}.ini -s ${output_res_in_minutes} -p > ../log/${model}_${yr}${mm}.log 2>&1
fi
