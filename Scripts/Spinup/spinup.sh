shopt -s expand_aliases		# Make sure aliases work in non-interactive shells
export TZ=UTC


# Load settings
source ./spinup.rc
default_spinup_end="2500-01-01T00:00"
startover=0


function PrintUseMessage {
	echo "Usage example:" >> /dev/stderr
	echo "  bash spinup.sh \"<snowpack command to execute>\" min_sim_depth=150" >> /dev/stderr
	echo "  Example: bash spinup.sh \"./src/usr/bin/snowpack -c cfgfiles/MERRA2_-76.952_-121.22.ini -e 2011-01-14 > log/MERRA2_-76.952_-121.22_TEST.log 2>&1\""
	echo "" >> /dev/stderr
	echo "  Add \"startover=1\" to ignore previous spinups and start over from the *sno file in snow_init_dir." >> /dev/stderr
	echo "  Example: bash spinup.sh \"./src/usr/bin/snowpack -c cfgfiles/MERRA2_-76.952_-121.22.ini -e 2011-01-14 > log/MERRA2_-76.952_-121.22_TEST.log 2>&1\" startover=1" >> /dev/stderr
}


#
# Parse command line arguments
#
i=0
for p in "$@"
do
	let i=${i}+1
	if (( $i == 1 )); then
		# The first argument must be the SNOWPACK command to execute
		to_exec=${p}
	else
		# All other arguments may overwrite settings
		if [[ "$p" == *"="* ]]; then
			eval ${p}
		else
			echo "WARNING: invalid command line parameter \"${p}\" ignored"
		fi
	fi
done

# If no command line argument provided, print usage message and exit
if (( $i == 0 )); then
	PrintUseMessage
	exit
fi

# Check some settings
if [ ! -e "${time_shift_script}" ]; then
	echo "ERROR: variable time_shift_script (see spinup.rc) does not point to an existing file!" >> /dev/stderr
	echo "       time_shift_script = ${time_shift_script}" >> /dev/stderr
	exit
fi

# Check if mawk exist, otherwise create alias
if ! command -v mawk &> /dev/null
then
	alias mawk='awk'
fi


#
# Parse command to execute to get relevant information
#
cfgfile=$(echo ${to_exec} | mawk '{for(i=1; i<=NF; i++) {if($i ~ /-c/) {print $(i+1)}}}')
if [ ! -e "${cfgfile}" ]; then
	echo "ERROR: spinup not started since ini file: ${cfgfile} could not be found or opened!"
	exit
fi
cfgfile_dir=$(dirname ${cfgfile})
# Check for nested ini files:
cfgfile_before=$(fgrep -i IMPORT_BEFORE ${cfgfile} | mawk -F= '{gsub(/^[ \t]+/,"",$NF); gsub(/[ \t]+$/,"",$NF); print $NF}')	# Use gsub to remove trailing and leading white spaces
cfgfile_after=$(fgrep -i IMPORT_AFTER ${cfgfile} | mawk -F= '{gsub(/^[ \t]+/,"",$NF); gsub(/[ \t]+$/,"",$NF); print $NF}')	# Use gsub to remove trailing and leading white spaces
# Construct array containing ini files in the sequence they should be parsed
cfgfiles=()
if [ -n "${cfgfile_before}" ]; then
	if [ ! -e "${cfgfile_dir}/${cfgfile_before}" ]; then
		echo "ERROR: spinup not started since ini file: ${cfgfile_dir}/${cfgfile_before} could not be found or opened!"
		exit
	fi
	cfgfiles+=("${cfgfile_dir}/${cfgfile_before}")
fi
if [ -n "${cfgfile}" ]; then
	cfgfiles+=("${cfgfile}")
fi
if [ -n "${cfgfile_after}" ]; then
	if [ ! -e "${cfgfile_dir}/${cfgfile_after}" ]; then
		echo "ERROR: spinup not started since ini file: ${cfgfile_dir}/${cfgfile_after} could not be found or opened!"
		exit
	fi
	cfgfiles+=("${cfgfile_dir}/${cfgfile_after}")
fi
max_sim_hs=$(cat ${cfgfiles[@]} | mawk '{if(/^\[/) {$0=toupper($0); if(/\[SNOWPACKADVANCED\]/) {read=1} else {read=0}}; if(read) {if(/MAX_SIMULATED_HS/) {val=$NF}}} END {print val}')
if [ -n "${max_sim_hs}" ]; then
	if (( $(echo "${max_sim_hs} > 0" | bc -l) )) && (( $(echo "${max_sim_hs} + 2 < ${min_sim_depth}" | bc -l) )); then
		echo "ERROR: spinup not started since MAX_SIMULATED_HS == ${max_sim_hs} and min_sim_depth == ${min_sim_depth}!"
		echo "       --> MAX_SIMULATED_HS should be at least 2 m larger than min_sim_depth to prevent infinite spinup."
		exit
	fi
fi


# Parse ini files to get relevant information
snofile=$(cat ${cfgfiles[@]} | mawk '{if(/^\[/) {$0=toupper($0); if(/\[INPUT\]/) {read=1} else {read=0}}; if(read) {if(/SNOWFILE1/) {val=$NF}}} END {print val}')
meteopath=$(cat ${cfgfiles[@]} | mawk '{if(/^\[/) {$0=toupper($0); if(/\[INPUT\]/) {read=1} else {read=0}}; if(read) {if(/METEOPATH/) {val=$NF}}} END {print val}')
snopath=$(cat ${cfgfiles[@]} | mawk '{if(/^\[/) {$0=toupper($0); if(/\[INPUT\]/) {read=1} else {read=0}}; if(read) {if(/SNOWPATH/) {val=$NF}}} END {print val}')
if [ -z "${snopath}" ]; then
	# If no SNOWPATH found, try METEOPATH
	snopath="${meteopath}"
fi
if [ ! -d "${snopath}" ]; then
	# Create path, if it doesn't yet exist
	mkdir -p ${snopath}
fi
outpath=$(cat ${cfgfiles[@]} | mawk '{if(/^\[/) {$0=toupper($0); if(/\[OUTPUT\]/) {read=1} else {read=0}}; if(read) {if(/METEOPATH/) {val=$NF}}} END {print val}')
if [ ! -d "${outpath}" ]; then
	# Create path, if it doesn't yet exist
	mkdir -p ${outpath}
fi
experiment=$(cat ${cfgfiles[@]} | mawk '{if(/^\[/) {$0=toupper($0); if(/\[OUTPUT\]/) {read=1} else {read=0}}; if(read) {if(/EXPERIMENT/) {val=$NF}}} END {print val}')
if [ -z "${experiment}" ]; then
	# Default experiment suffix:
	experiment="NO_EXP"
fi
if [ -z "${spinup_end}" ]; then
	# If no spinup_end determined, try to derive from *.smet file
	meteofile=$(cat ${cfgfiles[@]} | mawk '{if(/^\[/) {$0=toupper($0); if(/\[INPUT\]/) {read=1} else {read=0}}; if(read) {if(/STATION/) {val=$NF}}} END {print val}')
	if [[ "${meteofile}" != *.smet ]]; then
		# Add .smet file extension when no extension is provided
		meteofile="${meteofile}.smet"
	fi
	if [ ! -e "${meteopath}/${meteofile}" ]; then
		echo "WARNING: ${meteopath}/${meteofile} could not be found or openend, so no realistic end date of the simulation could be determined. Trying to continue ..."
		enddate=${default_spinup_end}
	else
		enddate=$(tail -1 ${meteopath}/${meteofile} | mawk '{print $1}')
	fi
else
	enddate=${spinup_end}
fi


# Parse sno file to get relevant information
stn=$(fgrep station_id ${snow_init_dir}/${snofile} | mawk '{print $NF}')
if [ -z "${spinup_start}" ]; then
	# If no spinup start date is provided, read from sno file:
	startdate=$(fgrep ProfileDate ${snow_init_dir}/${snofile} | mawk '{print $NF}')
else
	startdate=${spinup_start}
fi
snofile_in=${snopath}/${snofile}
snofile_out=${outpath}/${stn}_${experiment}.sno




#
#
#	THE READING OF SETTINGS IS NOW COMPLETE, DO THE SPINUPs
#
#
to_exec_spinup=$(echo "${to_exec}" | mawk -v ed=${enddate} '{for(i=1; i<=NF; i++) {if(i>1) {printf " "}; printf "%s", ($(i-1)=="-e")?(ed):($i)}; printf "\n"}')



# Initialize spinup counter
i=0		# Counting first spinups
j=0		# Counting second spinups
n_spinup=0	# To later store the number of first spinups
flag=0		# Check to see if a spinup was executed, to see if there are any problems

# Initialize spinup flags, assuming we start with a first spinup
spinup=1
spinup2=0

prev_sim_depth=-9999	# To check if the snow depth is increasing during the spinups

while :
do
	if [ -e "${snofile_out}" ] && [ "${startover}" == 0 ]; then
		sim_depth=$(mawk 'BEGIN {s=0; d=0} {if(d) {s+=$2}; if(/\[DATA\]/) {d=1}} END {print s}' ${snofile_out})
		if (( $(echo "${sim_depth} == 0" | bc -l) )) && (( "${i}" > 0 )) ; then
			echo "ERROR: spinup interrupted since SNOWPACK does not seem to build a snowpack/firn layer!"
			exit
		fi
		if [ ! -z "${checkscript}" ]; then
			checkscript_out=$(mawk -f ${checkscript} ${snofile_out})
			echo "Info: ${checkscript} returned ${checkscript_out}."
		else
			checkscript_out=0
		fi
		if (( $(echo "${sim_depth} > ${min_sim_depth}" | bc -l) )) || (( ${checkscript_out} )) || (( ${spinup2} )) ; then
			spinup=0
			if (( ! ${dospinup2} )); then
				# No second spinup
				spinup2=0
			else
				if (( ${n_spinup} == 0 )); then
					# Initialize second spinup
					spinup2=1
					n_spinup=${i}
					j=0
				elif (( ${j} >= ${n_spinup} )); then
					# We stop doing second spinups when we match the same number of iterations as for the first spinup
					spinup2=0
				fi
			fi
		elif (( $(echo "${sim_depth} <= ${prev_sim_depth}" | bc -l) )); then
			echo "ERROR: spinup interrupted since SNOWPACK does not seem to build an increasing snowpack/firn layer [depth = ${sim_depth}m]!"
			exit
		else
			prev_sim_depth=${sim_depth}
			spinup=1
		fi
		# Make sno file valid for init date
		bash ${time_shift_script} ${snofile_out} ${startdate} > ${snofile_in}
	else
		if (( "${flag}" > 0 )) ; then
			# If a spinup was executed, but the ${snofile_out} file does not exist, something must have gone wrong.
			echo "ERROR: spinup interrupted since SNOWPACK did not write output *sno file: ${snofile_out}!"
			exit
		fi
		startover=0
		spinup=1
		sim_depth=0
		if [ ! -e "${snofile_in}" ]; then
			bash ${time_shift_script} ${snow_init_dir}/${snofile} ${startdate} > ${snofile_in}
		fi
	fi

	if (( ${spinup} )) || (( ${spinup2} )); then
		let i=${i}+1

		# Use spinup.ini for specific settings for the spinups
		sed -i 's/^IMPORT_AFTER.*/IMPORT_AFTER	=	..\/spinup.ini/' ${cfgfile}

		# Print message
		if (( ${spinup} )); then
			echo "Spinup ${i} [depth = ${sim_depth}m]"
		fi
		if (( ${spinup2} )); then
			let j=${j}+1
			echo "2nd spinup ${j}/${n_spinup} [depth = ${sim_depth}m]"
		fi
	else
		# Use final.ini for specific settings for the final simulation
		sed -i 's/^IMPORT_AFTER.*/IMPORT_AFTER	=	..\/final.ini/' ${cfgfile}

		# Print message
		echo "Final simulation"
		eval ${to_exec}
		if [ ! -z "${zip_output_dir}" ]; then
			if [ ! -d "${zip_output_dir}" ]; then
				echo "ERROR: zip output directory is not a valid directory! [zip_output_dir=${zip_output_dir}]"
				exit
			fi
			zip ${zip_output_dir}/${stn}_${experiment}.zip -- ${outpath}/${stn}_${experiment}.*
		fi
		exit
	fi
	eval ${to_exec_spinup}
	# Check if simulation ran properly
	if (( $? != 0 )); then
		echo "ERROR: spinup interrupted since SNOWPACK threw an error!"
		exit
	fi
	flag=1
done

