#!/bin/sh
#prints min/max/mean for a given parameter in all smet files

if [ $# -lt 2 ]; then
	me=`basename $0`
	printf "Usage: \n"
	printf "\t$me . time\n\t\t to show the time range for all SMET files in the current directory\n"
	printf "\t$me ../input/meteo TA\n\t\t to show the range of the TA parameter for all SMET files in ../input/meteo\n"
	printf "\t$me ../input/meteo RH 2008-05-01 2008-09-30\n\t\t to show the range of the RH parameter for all SMET files in ../input/meteo between the two given dates\n"
	exit 0
fi

INPUT_DIR=$1
param=$2
if [ $# -ge 4 ]; then
	START_DATE=$3
	END_DATE=$4
fi
files=`find ${INPUT_DIR}/* -follow -maxdepth 0 -type f -name "*.smet"`

#on osX, it is necessary to force using gawk, if available
local_awk="awk"
if [ `command -v gawk` ]; then
	local_awk="gawk"
fi

if [ "${param}" = "time" ]; then
	for SMET in ${files}; do
		NAME=`basename "${SMET}" .smet`
		ALT=`head -100 "${SMET}" | grep --binary-files=text altitude | tr -s '\t' ' ' | tr -s ' ' | cut -d' ' -f3 | cut -d'.' -f1`
		JULIAN=`head -100 "${SMET}" | grep --binary-files=text fields | grep --binary-files=text julian`
		ISO=`head -100 "${SMET}" | grep --binary-files=text fields | grep --binary-files=text timestamp`
		start=`head -100 "${SMET}" | grep --binary-files=text -E "^[0-9][0-9][0-9][0-9]" | head -1 | tr -s '\t' ' ' | tr -s ' ' | cut -d' ' -f1`
		end=`tail -5 "${SMET}" | grep --binary-files=text -E "^[0-9][0-9][0-9][0-9]" | tail -1 | tr -s '\t' ' ' | tr -s ' ' | cut -d' ' -f1`
		header_nr_lines=`head -100 "${SMET}" | grep --binary-files=text -n "\[DATA\]" | cut -d':' -f1`
		full_nr_lines=`wc -l "${SMET}" | cut -d' ' -f1`
		nr_lines=`expr ${full_nr_lines} - ${header_nr_lines}`

		echo "${start} ${end} ${nr_lines}" | ${local_awk} '
				function getISO(ts){
					return sprintf("%s", strftime("%FT%H:%m:00", (ts-2440587.5)*24*3600))
				}
				function getSec(ts){
					gsub(/\-|\:|T/," ", ts); split(ts,d," ");
					date=sprintf("%04d %02d %02d %02d %02d 00",d[1],d[2],d[3],d[4],d[5]); 
					return mktime(date)
				}
				{
					if ("'"${ISO}"'" != "") {
						ISO_end=$2; ISO_start=$1;
						end=getSec($2); start=getSec($1); nr=$3;
					}
					if ("'"${JULIAN}"'" != "") {
						ISO_end=getISO($2); ISO_start=getISO($1);
						end=$2*24*3600; start=$1*24*3600; nr=$3
					}
					if (nr>1)
						period=int( (end-start)/(nr-1) + 0.5); #round to the nearest second
					else
						 period=0
					if (period<299 && period!=60 && period!=120 && period!=180 && period!=240)
						sampling=sprintf("%3.0f s  ", period)
					else if (period<60*60)
						sampling=sprintf("%3.0f min", period/60)
					else if (period<24*3600)
						sampling=sprintf("%3.0f h  ", period/3600)
					else {
						period_days=period/(3600*24)
						if (period_days==1) sampling=sprintf("%3.0f day", period_days)
						else sampling=sprintf("%3.0f days", period_days)
					}
					printf( "%04d m\t[ %s - %s ]\t~%s\t(%s)\n", "'"${ALT}"'", ISO_start, ISO_end, sampling, "'"${NAME}"'")
				}'
		
	done
	exit 0
fi

for SMET in ${files}; do
	IJ=`head -100 ${SMET} | grep --binary-files=text station_id | tr -s ' \t' | cut -d' ' -f3`
	if [ -z "${IJ}" ]; then
		IJ=`echo ${SMET} | cut -d'.' -f 1 | cut -d'_' -f2,3 | tr "_" ","`
	fi
	LAT=`head -100 ${SMET} | grep --binary-files=text latitude | tr -s ' \t' | cut -d' ' -f3`
	LON=`head -100 ${SMET} | grep --binary-files=text longitude | tr -s ' \t' | cut -d' ' -f3`
	ALT=`head -100 ${SMET} | grep --binary-files=text altitude | tr -s ' \t' | cut -d' ' -f3 | cut -d'.' -f1`
	NODATA=`head -100 ${SMET} | grep --binary-files=text nodata | tr -s ' \t' | cut -d' ' -f3`
	JULIAN=`head -100 "${SMET}" | grep --binary-files=text fields | grep --binary-files=text julian`

	${local_awk} '
	function toJul(ts){
		gsub(/\-|\:|T/," ", ts); split(ts,d," ");
		date=sprintf("%04d %02d %02d %02d %02d 00",d[1],d[2],d[3],d[4],d[5]);
		jul=mktime(date);
		return (jul/(24.*3600.)+2440587.5)
	}
	BEGIN {
		if ("'"${JULIAN}"'"=="") {
			start_date="'"${START_DATE}"'"
			end_date="'"${END_DATE}"'"
		} else {
			start_date_jul=toJul("'"${START_DATE}"'")
			end_date_jul=toJul("'"${END_DATE}"'")
		}
		
		param="'"${param}"'"
		if (param=="HNW") param="PSUM"
		nodata='"${NODATA}"'+0
		max=-1e4
		min=1e4
		f=2
	}
	/^fields/ {
		f=-1
		for(ii=4; ii<=NF; ii++) {
			if ($(ii)==param) {
				f=ii-2
			}
		}
		if (f==-1) {
			#printf("No %s in file %s\n", param, FILENAME)
			printf("\n")
			exit 0
		}
		printf("%s\t",$(f+2))
		next
	}
	/^altitude/ {
		printf("%04d m\t",$3)
		next
	}
	$0 !~ /^[a-zA-Z\[]/ {
		if (start_date!="") {
			if ($1<start_date) next
			if ($1>end_date) exit 0
		}
		if (start_date_jul!="") {
			if ($1<start_date_jul) next
			if ($1>end_date_jul) exit 0
		}
		val=$(f)+0
		if (val==nodata) next
		if (val>max) max = val
		if (val<min) min = val

		mean += val
		count++
	}
	END {
		if (f==-1 || count==0) {
			printf("  %7s - %7s  \t      %7s\t(%s)\n", "       ", " ", " ", "'"${IJ}"'")
			exit 0
		}
		mean /= count
		
		if (param=="TA" || param=="TSG" || param=="TSS" || param=="TD") {
			if (mean>100) {
				offset = -273.15
				min += offset
				max += offset
				mean += offset
			}
		}
		if (param=="RH" && mean>1) {
			min /= 100
			max /= 100
			mean /= 100
		}

		printf("[ %7.3g - %7.3g ]\tavg = %7.3g\t(%s)\n", min, max, mean, "'"${IJ}"'")
	}' ${SMET}
done
