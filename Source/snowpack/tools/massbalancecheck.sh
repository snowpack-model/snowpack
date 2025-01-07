#!/bin/bash
# Provide *.met or *.smet file as first argument.
# Note: for the first time step, no mass balance check can be done, because SWE of previous time step is unknown. Therefore, it is zero by definition.
if [ $# -lt 1 ]; then
	echo "This script reads a *.met or *.smet file (provided as first argument) and writes the mass balance on the stdout and statistics to stderr." > /dev/stderr
	echo "Invoke with: ./massbalancecheck.sh <(s)met file> [firstdate=YYYYMMDD] [lastdate=YYYYMMDD]" > /dev/stderr
	echo "" > /dev/stderr
	echo "Note: 1) the mass balance represents only the snow cover mass balance!" > /dev/stderr
	echo "      2) the mass balance can only be properly checked when the output resolution of the (s)met file is the" > /dev/stderr
	echo "         same as the snowpack calculation step length." > /dev/stderr
	echo "      3) the first time stamp in the (s)met file is not shown in the mass balance, as one cannot determine" > /dev/stderr
	echo "         this for the first time step (delta SWE cannot be determined)." > /dev/stderr
	echo "      4) using options firstdate and lastdate, one can define a period over which the mass balance should be" > /dev/stderr
	echo "         determined. Default is full period in (s)met-file. No spaces in command line options are allowed!" > /dev/stderr
	echo "" > /dev/stderr
	echo "Examples:" > /dev/stderr
	echo " ./massbalancecheck.sh WFJ2_flat.met > output.txt	  Writes mass balance in output.txt and shows overall mass balance statistics on screen." > /dev/stderr
	echo " ./massbalancecheck.sh WFJ2_flat.met > /dev/null    Just shows overall mass balance statistics on screen." > /dev/stderr
	echo " ./massbalancecheck.sh WFJ2_flat.smet | less        View the mass balance in less." > /dev/stderr
	echo " ./massbalancecheck.sh WFJ2_flat.smet firstdate=20071001 lastdate=20080323" > /dev/stderr
	echo "							  Determines mass balance between 1st of October 2007 up to and including 23rd of March 2008." > /dev/stderr
	exit
fi


# Initial settings
firstdate=0
lastdate=99999999


# Get met file name from first argument
file=$1
if [[ ${file} == *.met ]]; then
	met_file=${file}
	met=1
	smet=0
	field_sep=","
elif [[ ${file} == *.smet ]]; then
	smet_file=${file}
	met=0
	smet=1
	field_sep=" "
else
	echo "massbalancecheck.sh: ERROR: file ${file} does not have *.met or *.smet file extension!" > /dev/stderr
	exit
fi


# Read command line parameters
if [ $# -gt 1 ]; then
	for i in `seq 2 $#`
	do
		eval "let \$$i"
	done
fi


# Check if file exists
if [ ! -f "${file}" ]; then
	echo "massbalancecheck.sh: ERROR: file ${file} does not exist, is not a file, or cannot be opened!" > /dev/stderr
	exit
fi


# Check if file is not empty
if [ ! -s "${file}" ]; then
	echo "massbalancecheck.sh: ERROR: file ${file} is empty!" > /dev/stderr
	exit
fi


# Read header from met file
if (( ${met} )); then
	header=$(head -100 ${met_file} | grep -m 1 ^ID)
	if [ -z "${header}" ]; then
		echo "massbalancecheck.sh: ERROR: no header found in ${met_file}." > /dev/stderr
		exit
	fi
	# Determine column mapping *.met files
	#  -- date and time
	coldatetime=$(echo ${header} | sed 's/,/\n/g' | grep -nx "Date" | awk -F: '{print $1}')

	#  -- mass balance terms
	colrainrate=$(echo ${header} | sed 's/,/\n/g' | grep -nx "Rain rate" | awk -F: '{print $1}')
	colsnowrate=$(echo ${header} | sed 's/,/\n/g' | grep -nx "Precipitation rate at surface (solid only)" | awk -F: '{print $1}')
	colhsmeasured=$(echo ${header} | sed 's/,/\n/g' | grep -nx "Measured snow depth HS" | awk -F: '{print $1}')
	colhsmodel=$(echo ${header} | sed 's/,/\n/g' | grep -nx "Modelled snow depth (vertical)" | awk -F: '{print $1}')
	colSWE=$(echo ${header} | sed 's/,/\n/g' | grep -nx "SWE (of snowpack)" | awk -F: '{print $1}')
	colLWC=$(echo ${header} | sed 's/,/\n/g' | grep -nx "Liquid Water Content (of snowpack)" | awk -F: '{print $1}')
	colrunoff_surf=$(echo ${header} | sed 's/,/\n/g' | grep -nx "Snowpack runoff (virtual lysimeter -- snow only)" | awk -F: '{print $1}')
	colsubl=$(echo ${header} | sed 's/,/\n/g' | grep -nx "Sublimation" | awk -F: '{print $1}')
	colevap=$(echo ${header} | sed 's/,/\n/g' | grep -nx "Evaporation" | awk -F: '{print $1}')
	colwinddrift=$(echo ${header} | sed 's/,/\n/g' | grep -nx "Eroded mass" | awk -F: '{print $1}')
elif (( ${smet} )); then
	header=$(head -100 ${smet_file} | grep -m 1 ^fields)
	if [ -z "${header}" ]; then
		echo "massbalancecheck.sh: ERROR: no header found in ${smet_file}." > /dev/stderr
		exit
	fi
	# Determine column mapping *.smet files
	#  -- date and time
	coldatetime=$(echo ${header} | sed 's/ /\n/g' | grep -nx "timestamp" | awk -F: '{print $1-2}')

	#  -- mass balance terms
	colrainrate=$(echo ${header} | sed 's/ /\n/g' | grep -nx "MS_Rain" | awk -F: '{print $1-2}')
	colsnowrate=$(echo ${header} | sed 's/ /\n/g' | grep -nx "MS_Snow" | awk -F: '{print $1-2}')
	colhsmeasured=$(echo ${header} | sed 's/ /\n/g' | grep -nx "HS_meas" | awk -F: '{print $1-2}')
	colhsmodel=$(echo ${header} | sed 's/ /\n/g' | grep -nx "HS_mod" | awk -F: '{print $1-2}')
	colSWE=$(echo ${header} | sed 's/ /\n/g' | grep -nx "SWE" | awk -F: '{print $1-2}')
	colLWC=$(echo ${header} | sed 's/ /\n/g' | grep -nx "MS_Water" | awk -F: '{print $1-2}')
	colrunoff_surf=$(echo ${header} | sed 's/ /\n/g' | grep -nx "MS_SN_Runoff" | awk -F: '{print $1-2}')
	colsubl=$(echo ${header} | sed 's/ /\n/g' | grep -nx "MS_Sublimation" | awk -F: '{print $1-2}')
	colevap=$(echo ${header} | sed 's/ /\n/g' | grep -nx "MS_Evap" | awk -F: '{print $1-2}')
	colwinddrift=$(echo ${header} | sed 's/ /\n/g' | grep -nx "MS_Wind" | awk -F: '{print $1-2}')

	# Check units
	units_MS_Rain=$(grep -m 1 ^plot_unit ${smet_file} | awk -v col=${colrainrate} '{print $(col+2)}')
	units_MS_Snow=$(grep -m 1 ^plot_unit ${smet_file} | awk -v col=${colsnowrate} '{print $(col+2)}')
	if [ -z "${units_MS_Rain}" ] || [ -z "${units_MS_Rain}" ]; then
		echo "massbalancecheck.sh: WARNING: units for MS_Rain and MS_Snow not found, assuming kg/m2/h."
		echo ${units_MS_Rain}  ${units_MS_Snow}
	else
		if [[ "${units_MS_Rain}" == "kg/m2/h" ]]; then
			rainrate=1
		elif [[ "${units_MS_Rain}" == "kg/m2" ]]; then
			rainrate=0
		else
			echo "massbalancecheck.sh: ERROR: Unknown units for MS_Rain: ${units_MS_Rain}." > /dev/stderr
			exit
		fi
		if [[ "${units_MS_Snow}" == "kg/m2/h" ]]; then
			snowrate=1
		elif [[ "${units_MS_Snow}" == "kg/m2" ]]; then
			snowrate=0
		else
			echo "massbalancecheck.sh: ERROR: Unknown units for MS_Snow: ${units_MS_Snow}." > /dev/stderr
			exit
		fi
	fi
fi


error=0
if [ -z "${coldatetime}" ]; then
	echo "massbalancecheck.sh: ERROR: date/time not found in one of the columns." > /dev/stderr
	error=1
fi
if [ -z "${colrainrate}" ]; then
	echo "massbalancecheck.sh: ERROR: snow rate not found in one of the columns." > /dev/stderr
	error=1
fi
if [ -z "${colsnowrate}" ]; then
	echo "massbalancecheck.sh: ERROR: snow rate not found in one of the columns." > /dev/stderr
	error=1
fi
if [ -z "${colhsmeasured}" ]; then
	echo "massbalancecheck.sh: ERROR: measured hs not found in one of the columns." > /dev/stderr
	error=1
fi
if [ -z "${colhsmodel}" ]; then
	echo "massbalancecheck.sh: ERROR: modeled hs not found in one of the columns." > /dev/stderr
	error=1
fi
if [ -z "${colSWE}" ]; then
	echo "massbalancecheck.sh: ERROR: SWE not found in one of the columns." > /dev/stderr
	error=1
fi
if [ -z "${colLWC}" ]; then
	echo "massbalancecheck.sh: ERROR: LWC not found in one of the columns." > /dev/stderr
	error=1
fi
if [ -z "${colrunoff_surf}" ]; then
	echo "massbalancecheck.sh: ERROR: snowpack runoff not found in one of the columns." > /dev/stderr
	error=1
fi
if [ -z "${colsubl}" ]; then
	echo "massbalancecheck.sh: ERROR: sublimation not found in one of the columns." > /dev/stderr
	error=1
fi
if [ -z "${colevap}" ]; then
	echo "massbalancecheck.sh: ERROR: evaporation not found in one of the columns." > /dev/stderr
	error=1
fi
if [ -z "${colwinddrift}" ]; then
	echo "massbalancecheck.sh: ERROR: erosion not found in one of the columns." > /dev/stderr
	error=1
fi
if [ "${error}" -eq 1 ]; then
	exit
fi

# -- Determine file resolution
if (( ${met} )); then
	nsamplesperday=$(cat ${met_file} | sed '1,/\[DATA\]/d' | awk -F, '{print $'${coldatetime}'}' | awk '{print $1}' | sort | uniq -c | awk '{print $1}' | sort -nu | tail -1)
elif (( ${smet} )); then
	nsamplesperday=$(cat ${smet_file} | sed '1,/\[DATA\]/d' | awk -F, '{print substr($'${coldatetime}', 1, 10)}' | awk '{print $1}' | sort | uniq -c | awk '{print $1}' | sort -nu | tail -1)
fi
if [ -z "${nsamplesperday}" ]; then
	echo "massbalancecheck.sh: ERROR: file resolution could not be determined." > /dev/stderr
	exit
fi


# Create header
echo "#Date time measured_HS modelled_HS SWE     LWC    rain_rate snow_rate snowpack_runoff subl   evap   winderosion deltaSWE massbalance mass_in mass_out"
echo "#--   --   --          --          --      --     M+        M+        M+              M+     M+     M+          M-       error       totals  totals"
echo "#-    -    cm          cm          kg_m-2  kg_m-2 kg_m-2    kg_m-2    kg_m-2          kg_m-2 kg_m-2 kg_m-2      kg_m-2   kg_m-2      kg_m-2  kg_m-2"



# Process data (note that the lines below are all piped together).
#  -- Cut out data
sed '1,/\[DATA\]/d' ${file} | \
#  -- Select all the massbalance terms, make them correct sign and correct units. Also makes sure some terms are only considered when they are a part of the SNOW mass balance (like evaporation, which may also originate from soil).
#     Note we store the previous SWE, to know whether evaporation and/or sublimation was actually from soil or from snow. For the first time step it doesn't matter what we do here, as we will cut out this first line later.
#         (We cannot cut out this first line here, as the previous time step SWE is also needed for the mass balance calculations).
awk -F"${field_sep}" -v rr=${rainrate} -v sr=${snowrate} '{n++; if(n==1){prevSWE=1}; print $'${coldatetime}', $'${colhsmeasured}', $'${colhsmodel}', $'${colSWE}', $'${colLWC}', ($'${colSWE}'>0.0 || $'${colsnowrate}'>0.0)?($'${colrainrate}'*((rr)?(24/'${nsamplesperday}'):(1.))):0, $'${colsnowrate}'*((sr)?(24/'${nsamplesperday}'):(1.)), -1.*$'${colrunoff_surf}', (prevSWE>0.0 || $'${colSWE}'>0.0 || $'${colsnowrate}'>0.0)?($'${colsubl}'):0, (prevSWE>0.0 || $'${colSWE}'>0.0 || $'${colsnowrate}'>0.0)?($'${colevap}'):0, ($'${colwinddrift}'>0)?-1.0*($'${colwinddrift}')*(24/'${nsamplesperday}'):0; prevSWE=$'${colSWE}'}' | \
#  -- Reformat time
awk -v met=${met} '{if(met) {printf "%04d%02d%02d %02d%02d", substr($1,7,4), substr($1,4,2), substr($1,1,2), substr($2,1,2), substr($2,4,2); for(i=3; i<=NF; i++) {printf " %s", $i}} else {printf "%04d%02d%02d %02d%02d", substr($1,1,4), substr($1,6,2), substr($1,9,2), substr($1,12,2), substr($1,15,2); for(i=2; i<=NF; i++) {printf " %s", $i}}; printf "\n"}' | \
# Now select period
awk '($1>='${firstdate}' && $1<='${lastdate}') {print $0}' | \
#  -- Now do all the other calculations
#     First, determine deltaSWE when it is not the first line read in (if so, we cannot determine the mass balance, as the previous value of SWE is unknown).
awk '{n++; if(n>1) \
	{deltaSWE=($5-prevSWE); \
	#Determine mass balance error:
	massbalance=$7+$8+$9+$10+$11+$12-deltaSWE; \
	#Determine mass input in system (taking the terms only when they are positive)
	mass_in=(($7>0.0)?$7:0)+(($8>0.0)?$8:0)+(($9>0.0)?$9:0)+(($10>0.0)?$10:0)+(($11>0.0)?$11:0)+(($12>0.0)?$12:0); \
	#Determine mass output in system (taking the terms only when they are negative)
	mass_out=(($7<0.0)?$7:0)+(($8<0.0)?$8:0)+(($9<0.0)?$9:0)+(($10<0.0)?$10:0)+(($11<0.0)?$11:0)+(($12<0.0)?$12:0); \
	#Do the statistics (mass balance error sum, min and max values)
	massbalancesum+=massbalance; massbalancesum2+=sqrt(massbalance*massbalance); if(massbalance>maxmassbalance){maxmassbalance=massbalance; maxmassbalancedate=$1; maxmassbalancetime=$2}; if(massbalance<minmassbalance){minmassbalance=massbalance; minmassbalancedate=$1; minmassbalancetime=$2}; \
	#Write to stdout
	print $0, deltaSWE, massbalance, mass_in, mass_out}; \
#Store the line written out, so for the next line, we have the SWE of the previous line available (needed to determine deltaSWE):
prevLine=$0; prevSWE=$5;} \
#Write out statistics to stderr:
END {printf "Summary of file: '${file}'\n-------------------------------------------------------------------------------------\nSum of mass balance error (kg_m-2): %.6f\nSum of absolute mass balance error (kg_m-2): %.6f\nMaximum positive mass balance error (kg_m-2): %.6f at %08d, %04d\nMinimum negative mass balance error (kg_m-2): %.6f at %08d, %04d\n", massbalancesum, massbalancesum2, maxmassbalance, maxmassbalancedate, maxmassbalancetime, minmassbalance, minmassbalancedate, minmassbalancetime > "/dev/stderr"}'

