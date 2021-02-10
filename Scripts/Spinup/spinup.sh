# Load settings
source ./spinup.rc
default_spinup_end="2500-01-01T00:00"


function PrintUseMessage {
	echo "Usage example:" >> /dev/stderr
	echo "  bash spinup.sh \"<snowpack command to execute>\" min_sim_depth=150" >> /dev/stderr
	echo "  Example: bash spinup.sh \"./src/usr/bin/snowpack -c cfgfiles/MERRA2_-76.952_-121.22.ini -e 2011-01-14 > log/MERRA2_-76.952_-121.22_TEST.log 2>&1\""
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


#
# Parse command to execute to get relevant information
#
if [ -z "${spinup_end}" ]; then
	enddate=${default_spinup_end}
else
	enddate=${spinup_end}
fi
to_exec_spinup=$(echo "${to_exec}" | mawk -v ed=${enddate} '{for(i=1; i<=NF; i++) {if(i>1) {printf " "}; printf "%s", ($(i-1)=="-e")?(ed):($i)}; printf "\n"}')
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
snopath=$(cat ${cfgfiles[@]} | mawk '{if(/^\[/) {$0=toupper($0); if(/\[INPUT\]/) {read=1} else {read=0}}; if(read) {if(/SNOWPATH/) {val=$NF}}} END {print val}')
if [ -z "${snopath}" ]; then
	# If no SNOWPATH found, try METEOPATH
	snopath=$(cat ${cfgfiles[@]} | mawk '{if(/^\[/) {$0=toupper($0); if(/\[INPUT\]/) {read=1} else {read=0}}; if(read) {if(/METEOPATH/) {val=$NF}}} END {print val}')
fi
outpath=$(cat ${cfgfiles[@]} | mawk '{if(/^\[/) {$0=toupper($0); if(/\[OUTPUT\]/) {read=1} else {read=0}}; if(read) {if(/METEOPATH/) {val=$NF}}} END {print val}')
experiment=$(cat ${cfgfiles[@]} | mawk '{if(/^\[/) {$0=toupper($0); if(/\[OUTPUT\]/) {read=1} else {read=0}}; if(read) {if(/EXPERIMENT/) {val=$NF}}} END {print val}')
if [ -z "${experiment}" ]; then
	# Default experiment suffix:
	experiment="NO_EXP"
fi


# Parse sno file to get relevant information
stn=$(fgrep station_id ${snow_init_dir}/${snofile} | mawk '{print $NF}')
if [ -z "${spinup_start}" ]; then
	# If no spinup start date is provided, read from sno file:
	startdate=$(fgrep ProfileDate ${snow_init_dir}/${snofile} | mawk '{print $NF}')
else
	startdate=${spinup_start}
fi
snofile_out=${outpath}/${stn}_${experiment}.sno




#
#
#	THE READING OF SETTINGS IS NOW COMPLETE, DO THE SPINUPs
#
#



# Initialize spinup counter
i=0		# Counting first spinups
j=0		# Counting second spinups
n_spinup=0	# To later store the number of first spinups
flag=0		# Check to see if a spinup was executed, to see if there are any problems

# Initialize spinup flags, assuming we start with a first spinup
spinup=1
spinup2=0


while :
do
	if [ -e "${snofile_out}" ]; then
		sim_depth=$(mawk 'BEGIN {s=0; d=0} {if(d) {s+=$2}; if(/\[DATA\]/) {d=1}} END {print s}' ${snofile_out})
		if (( $(echo "${sim_depth} == 0" | bc -l) )) && (( "${i}" > 0 )) ; then
			echo "ERROR: spinup interrupted since SNOWPACK does not seem to build a snowpack/firn layer!"
			exit
		fi
		if (( $(echo "${sim_depth} > ${min_sim_depth}" | bc -l) )) || (( ${spinup2} )) ; then
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
		else
			spinup=1
		fi
		# Make sno file valid for init date
		bash ${time_shift_script} ${snofile_out} ${startdate} > ${snopath}/${snofile}
	else
		if (( "${flag}" > 0 )) ; then
			# If a spinup was executed, but the ${snofile_out} file does not exist, something must have gone wrong.
			echo "ERROR: spinup interrupted since SNOWPACK did not write output *sno file: ${snofile_out}!"
			exit
		fi
		spinup=1
		sim_depth=0
		bash ${time_shift_script} ${snow_init_dir}/${snofile} ${startdate} > ${snopath}/${snofile}
	fi

	if (( ${spinup} )) || (( ${spinup2} )); then
		let i=${i}+1

		# Ensure output is not written
		sed -i 's/^PROF_WRITE.*/PROF_WRITE		=       FALSE/' ${cfgfile}
		sed -i 's/^TS_WRITE.*/TS_WRITE		=       FALSE/' ${cfgfile}

		# Print message
		if (( ${spinup} )); then
			echo "Spinup ${i} [depth = ${sim_depth}m]"
		fi
		if (( ${spinup2} )); then
			let j=${j}+1
			echo "2nd spinup ${j}/${n_spinup} [depth = ${sim_depth}m]"
		fi
	else
		# Ensure output is written
		sed -i 's/^PROF_WRITE.*/PROF_WRITE		=       TRUE/' ${cfgfile}
		sed -i 's/^TS_WRITE.*/TS_WRITE		=       TRUE/' ${cfgfile}

		# Print message
		echo "Final simulation"
		eval ${to_exec}
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

