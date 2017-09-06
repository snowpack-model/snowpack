#!/bin/bash
#
# This script can resample SMET data to lower resolutions.
#
# General recipe: PSUM is always an average (to keep [mm/h]), other variables can be averaged with the switch -m, or else they are just taken on the resampled time stamp.
#
# Use bash resample_smetdata.sh <filename> <resolution> <-m>
#
#
# Author: Nander Wever
#

# Get command line parameters
filename=$1								#SMET filename
resolution=$2								#output resolution
resolution_minutes=`echo ${resolution} | awk '{print int($1/60)}'`	#output resolution in minutes
if [ -z "$3" ]; then
	resamplemethod=0		#Just resample
else
	if [ "$3" == "-m" ]; then
		resamplemethod=1	#Take mean
	else
		resamplemethod=0	#Just resample
	fi
fi

WriteUsageMessage()
{
	echo "Use: bash resample_smetdata.sh <filename> <resolution> <-m>"
	echo "  <filename>: SMET file"
	echo "  <resolution>: new resolution in seconds. (minimum 2 minutes, maximum 1 day)"
	echo "  <-m>: optional, if -m is added, the mean values are taken, else it is just resampled. PSUM is always an average."
	echo "Output is written to std out."
	echo "Note: - holes in the data are filled with nodata values, and then resampled."
	echo "      - script assumes UTC time zone (without DST)."
	echo "      - when an error is encountered, the script tries to output the SMET file at the original resolution."
}

# Check command line parameters
if [ -z "${filename}" ]; then
	echo "ERROR: no file name specified."
	WriteUsageMessage
	exit
fi

if [ -z "${resolution}" ]; then
	echo "ERROR: no resolution specified."
	WriteUsageMessage
	exit
fi

if [ ! -e "${filename}" ] || [ ! -s "${filename}" ] || [ -d "${filename}" ]; then
	echo "ERROR: file [${filename}] does not exist, is empty, or is a directory."
	WriteUsageMessage
	exit
fi

# Set shell time to UTC, to allow correct parsing by awk-mktime
export TZ=UTC

# Dump header
cat ${filename} | grep -v ^[0-9] | sed -e 's/^[ \t]*//' | sed '/^$/d'

# Now determine some info needed to make a complete file (without holes in the data)
firsttimestamp=`cat ${filename} | grep ^[0-9] | head -1 | sed 's/[-T:]/ /g' | awk '{print mktime(sprintf("%04d %02d %02d %02d %02d %02d %d", $1, $2, $3, $4, $5, 0, 1))}'`
timeresolution=`cat ${filename} | grep ^[0-9] | awk '{if (NR==1) {printf "%s\n", $1} else {printf "%s\n%s\n", $1, $1}}' | sed '$!N;s/\n/ /' | sed 's/[-T:]/ /g' | awk '(NF==10) {print mktime(sprintf("%04d %02d %02d %02d %02d %02d %d", $6, $7, $8, $9, $10, 0, 1))-mktime(sprintf("%04d %02d %02d %02d %02d %02d %d", $1, $2, $3, $4, $5, 0, 1))}' | sort -nk1 | uniq -c | sort -nrk1 | awk '(NR==1){print $2}'`   # Native resolution of file is determined by the difference between two time stamps that occurs most often.
lasttimestamp=`cat ${filename} | grep ^[0-9] | tail -1 | sed 's/[-T:]/ /g' | awk '{print mktime(sprintf("%04d %02d %02d %02d %02d %02d %d", $1, $2, $3, $4, $5, 0, 1))}'`
nodatavalue=`cat ${filename} | grep ^nodata | head -1 | awk -F= '{print $NF}' | sed 's/ //g'`
nsensors=`cat ${filename} | grep ^[0-9] | head -1 | awk '{print NF-1}'`
col_psum=`cat ${filename} | grep ^fields | head -1 | awk -F= '{print $NF}' | tr ' ' '\n' | grep -v ^$ | grep -n PSUM | awk -F: '{print $1}'`
if [ -z "${col_psum}" ]; then
	col_psum=-1
fi

# Check for valid values, is it a SMET file?
if [ -z "${firsttimestamp}" ] || [ -z "${lasttimestamp}" ] || [ -z "${timeresolution}" ] || [ -z "${nsensors}" ]; then
	cat ${filename} | grep ^[0-9] | sed -e 's/^[ \t]*//' | sed '/^$/d'
	echo "ERROR: file [${filename}] does not seems to be a valid SMET file, or it does not have any data to resample."
        exit
fi

# Check for valid number of sensors
if (( "${nsensors}" < 1 )); then
	cat ${filename} | grep ^[0-9] | sed -e 's/^[ \t]*//' | sed '/^$/d'
	echo "ERROR: file [${filename}] does not seem to have any parameter."
        exit
fi

# If time resolution of file is larger than requested resolution, just give the output, and send error message to stdout
if (( ${timeresolution} > ${resolution} )); then
	cat ${filename} | grep ^[0-9] | sed -e 's/^[ \t]*//' | sed '/^$/d'
	echo "ERROR: requested resolution smaller than original resolution." 1>&2
	exit
fi

# If time resolution of file equals the requested resolution, just give the output
if (( ${timeresolution} == ${resolution} )); then
	cat ${filename} | grep ^[0-9] | sed -e 's/^[ \t]*//' | sed '/^$/d'
	exit
fi

# If time resolution is more than one day, give error and just give the output equal to the input.
if (( ${resolution} > 86400 )); then
	cat ${filename} | grep ^[0-9] | sed -e 's/^[ \t]*//' | sed '/^$/d'
	echo "ERROR: requested resolution larger than 1 day (86400 seconds). This script can't handle that." 1>&2
	exit
fi

# Resample, just plain resampling
if (( ${resamplemethod} == 0 )); then
	#cat: start pipe
	#grep: only select data rows from SMET
        #sed: delete all leading white spaces and tabs
        #sed: delete all blank lines
	#awk: add a 0 in the first column to mark original data. Then, add the resolution of the SMET file nodata values, marked by a 1 in the first column.
	#sort: then sort, first for the timestamp, then for the first column, such that when original data is available, it is appearing first, before the nodata values.
	#cut: removes the first column, which contains the flag.
	#uniq: now take unique timestamps. Because this selects the first occurrence of the time stamp, the way we sorted it, makes the original data selected first (if available), and then the nodata values.
	#awk: this is the actual resampling.
	#     note that the awk constuction is very similar. The trick is that for resampling {sum[k]=$k; n[k]=1;} is used, where for taking the mean {sum[k]+=$k; n[k]+=1;} is used.
	#     note that the psum also should be calculated as a mean, else it makes no sense. So for calculating the mean, the if-statement {if(k=='${col_psum}')} is actually superfluous (both the if and
	#     the else block do exactly the same), but it is kept for coherence.
	cat ${filename} | grep ^[0-9] | sed -e 's/^[ \t]*//' | sed '/^$/d' | awk '{print 0, $0} END {for(j='${firsttimestamp}'; j<='${lasttimestamp}'; j=j+'${timeresolution}') {printf "1 %s", strftime("%Y-%m-%dT%H:%M", j); for (i=1; i<='${nsensors}'; i++) {printf " %s", '${nodatavalue}'} printf "\n"}}' | sort -k 2 -k 1 | cut -d\  -f2- | uniq -w16 | awk '{for(k=2; k<=NF; k++) {if($k!='${nodatavalue}') {if(k=='${col_psum}') {sum[k]+=$k; n[k]+=1} else {sum[k]=$k; n[k]=1;}}};     if((substr($1, 9, 2)*60*24+substr($1, 12, 2)*60+substr($1, 15, 2))%'${resolution_minutes}'==0) {{printf "%s", $1; for (k=2; k<=NF; k++) {printf " %s", (n[k]>0)?sum[k]/n[k]:'${nodatavalue}'; sum[k]=0; n[k]=0}; printf "\n"}}}'
fi

# Resample, calculating mean values
if (( ${resamplemethod} == 1 )); then
	#cat: start pipe
	#grep: only select data rows from SMET
        #sed: delete all leading white spaces and tabs
        #sed: delete all blank lines
	#awk: add a 0 in the first column to mark original data. Then, add the resolution of the SMET file nodata values, marked by a 1 in the first column.
	#sort: then sort, first for the timestamp, then for the first column, such that when original data is available, it is appearing first, before the nodata values.
	#cut: removes the first column, which contains the flag.
	#uniq: now take unique timestamps. Because this selects the first occurrence of the time stamp, the way we sorted it, makes the original data selected first (if available), and then the nodata values.
	#awk: this is the actual resampling.
	#     note that the awk constuction is very similar. The trick is that for resampling {sum[k]=$k; n[k]=1;} is used, where for taking the mean {sum[k]+=$k; n[k]+=1;} is used.
	#     note that the psum also should be calculated as a mean, else it makes no sense. So for calculating the mean, the if-statement {if(k=='${col_psum}')} is actually superfluous (both the if and
	#     the else block do exactly the same), but it is kept for coherence.
	cat ${filename} | grep ^[0-9] | sed -e 's/^[ \t]*//' | sed '/^$/d' | awk '{print 0, $0} END {for(j='${firsttimestamp}'; j<='${lasttimestamp}'; j=j+'${timeresolution}') {printf "1 %s", strftime("%Y-%m-%dT%H:%M", j); for (i=1; i<='${nsensors}'; i++) {printf " %s", '${nodatavalue}'} printf "\n"}}' | sort -k 2 -k 1 | cut -d\  -f2- | uniq -w16 | awk '{for(k=2; k<=NF; k++) {if($k!='${nodatavalue}') {if(k=='${col_psum}') {sum[k]+=$k; n[k]+=1} else {sum[k]+=$k; n[k]+=1;}}};     if((substr($1, 9, 2)*60*24+substr($1, 12, 2)*60+substr($1, 15, 2))%'${resolution_minutes}'==0) {{printf "%s", $1; for (k=2; k<=NF; k++) {printf " %s", (n[k]>0)?sum[k]/n[k]:'${nodatavalue}'; sum[k]=0; n[k]=0}; printf "\n"}}}'
fi

#Reset time zone
unset TZ
