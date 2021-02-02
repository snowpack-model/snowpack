min_sim_depth=150
time_shift_script="../../../../src/polarsnowpack_core/snowpack/Source/snowpack/tools/timeshift_sno_files.sh"


# Parse string to execute to get relevant information
to_exec="$@"
to_exec_spinup=$(echo "${to_exec}" | mawk '{for(i=1; i<=NF; i++) {if(i>1) {printf " "}; printf "%s", ($(i-1)=="-e")?("NOW"):($i)}; printf "\n"}')
cfgfile=$(echo ${to_exec} | mawk '{for(i=1; i<=NF; i++) {if($i ~ /-c/) {print $(i+1)}}}')
snofile=$(fgrep SNOWFILE1 ${cfgfile} | mawk '{print $NF}')
stn=$(fgrep station_id snow_init/${snofile} | mawk '{print $NF}')
startdate=$(fgrep ProfileDate snow_init/${snofile} | mawk '{print $NF}')
experiment=$(fgrep EXPERIMENT ${cfgfile} | mawk '{print $NF}')
snofile_out=./output/${stn}_${experiment}.sno


# Initialize spinup counter
i=0		# Counting first spinups
j=0		# Counting second spinups
n_spinup=0	# To later store the number of first spinups

# Initialize spinup flags, assuming we start with a first spinup
spinup=1
dospinup2=1	# Should we do a 2nd spinup (yes/no)?
spinup2=0

while :
do
	if [ -e "${snofile_out}" ]; then
		sim_depth=$(mawk 'BEGIN {s=0; d=0} {if(d) {s+=$2}; if(/\[DATA\]/) {d=1}} END {print s}' ${snofile_out})
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
		bash ${time_shift_script} ${snofile_out} ${startdate} > ./current_snow/${snofile}
	else
		spinup=1
		sim_depth=0
		bash ${time_shift_script} ./snow_init/${snofile} ${startdate} > ./current_snow/${snofile}
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
done

