#!/bin/bash
export TZ=UTC
if [ $# -lt 2 ]; then
	echo "This script makes a *sno file valid for the provided date." > /dev/stderr
	echo "The ProfileDate variable will be set to the specified timestamp" > /dev/stderr
	echo "and the deposition dates of all the layers are shifted by the" > /dev/stderr
	echo "same offset." > /dev/stderr
	echo "This script is useful to do a spinup for polar simulations, where" > /dev/stderr
	echo "SNOWPACK runs can be looped multiple times over a spinup period," > /dev/stderr
	echo "each time taking the output *.sno file and using it as initial" > /dev/stderr
	echo "state for the next spinup iteration." > /dev/stderr
	echo "" > /dev/stderr
	echo "Invoke with: ./timeshift_sno_files.sh <sno-file> <date>" > /dev/stderr
	echo "Profile will be written to stdout and can be redirected as wished." > /dev/stderr
	echo "" > /dev/stderr
	echo "Example: ./timeshift_sno_files.sh WFJ2.sno 1990-10-01T00:00 > WFJ2_new.sno" > /dev/stderr
	echo "" > /dev/stderr
	echo "Requires gawk version 4.1.2 or newer." > /dev/stderr
	echo "" > /dev/stderr
	echo "Known bugs:" > /dev/stderr
	echo "Older gawk versions (4.1.1 and older) cannot work with dates before" > /dev/stderr
	echo "1970-01-01 (i.e., negative unix timestamps)." > /dev/stderr
	echo "In that case, use timeshift_sno_files_oldawk.sh instead." > /dev/stderr
	exit
fi
awk -v t=$2 'BEGIN {td=mktime(sprintf("%04d %02d %02d %02d %02d %02d 0", substr(t,1,4), substr(t,6,2), substr(t,9,2), substr(t,12,2), substr(t,15,2), substr(t,19,2)))} {if(/ProfileDate/) {s=$NF; sd=mktime(sprintf("%04d %02d %02d %02d %02d %02d 0", substr(s,1,4), substr(s,6,2), substr(s,9,2), substr(s,12,2), substr(s,15,2), substr(s,19,2))); printf "ProfileDate\t = %s\n", t} else if(!data) {print} else {ld=mktime(sprintf("%04d %02d %02d %02d %02d %02d, 0", substr($1,1,4), substr($1,6,2), substr($1,9,2), substr($1,12,2), substr($1,15,2), substr($1,19,2))); printf "%s", strftime("%Y-%m-%dT%H:%M:%S", ld - (sd - td)); for(i=2; i<=NF; i++) {printf " %s", $i}; printf "\n"}; if(/DATA/) {data=1}}' $1
