#!/bin/sh
# SPDX-License-Identifier: LGPL-3.0-or-later
# split a given smet file into daily/monthly/yearly files. The ACDD time_coverage_* fields are also updated accordingly

if [ $# -lt 2 ]; then
	me=`basename $0`
	printf "Split a given smet file into daily/monthly/yearly files"
	printf "Usage: \n"
	printf "\t$me {smet_file} monthly\n\t\t split the file given as first argument into monthly files\n"
	exit 0
fi

SMET=$1
FILE_ROOT=`basename $1 .smet`
TYPE=$2

#on osX, it is necessary to force using gawk, if available
local_awk="awk"
if [ `command -v gawk` ]; then
	local_awk="gawk"
fi

#Now split the file
${local_awk} '
	BEGIN {
		type_str=tolower( "'"${TYPE}"'" )
		if (type_str=="daily") type=3
		else if (type_str=="monthly") type=2
		else if (type_str=="yearly") type=1
		else {
			printf("Wrong splitting type given as argument, please choose one of daily/monthly/yearly\n") > "/dev/stderr"
			exit(1)
		}

		in_header=1
		output_file=""
	}
	/fields/ {
		gsub(/\r/, "")
		for(ii=3; ii<=NF; ii++) {
			if ($(ii)=="timestamp") timestamp_idx=(ii-2)
			if ($(ii)=="julian") {
				printf("Only files with ISO timestamps are currently supported\n") > "/dev/stderr"
				exit(1)
			}
		}
	}
	(in_header==1) {
		headers=sprintf("%s%s\n", headers, $0)
	}
	/\[DATA\]/ {
		in_header=0
		next
	}
	(in_header==0)  {
		ts=$(timestamp_idx)
		gsub(/\-|:|T/," ", ts); split(ts,d," ");
		if (type==3) {
			current_output_file=sprintf("%s-%04d-%02d-%02d.smet", "'"${FILE_ROOT}"'", d[1], d[2], d[3])
		} else if (type==2) {
			current_output_file=sprintf("%s-%04d-%02d.smet", "'"${FILE_ROOT}"'", d[1], d[2])
		} else if (type==1) {
			current_output_file=sprintf("%s-%04d.smet", "'"${FILE_ROOT}"'", d[1])
		}

		if (output_file!=current_output_file) {
			if (output_file!="") close(output_file)
			output_file = current_output_file
			printf("%s", headers) > output_file
		}
		print $0 >> output_file
	}
	' ${SMET}

#If the ACDD time_coverage_* fields are present, update them
update_ACDD() {
	TARGET_SMET=$1
	start=`head -100 "${TARGET_SMET}" | grep --binary-files=text -E "^[0-9][0-9][0-9][0-9]" | head -1 | tr -s '\t' ' ' | tr -s ' ' | cut -d' ' -f1`
	end=`tail -5 "${TARGET_SMET}" | grep --binary-files=text -E "^[0-9][0-9][0-9][0-9]" | tail -1 | tr -s '\t' ' ' | tr -s ' ' | cut -d' ' -f1`
	if [ -z "${start}" ]; then
		return 0
	fi
	
 	header_nr_lines=`head -100 "${TARGET_SMET}" | grep --binary-files=text -n "\[DATA\]" | cut -d':' -f1`
 	full_nr_lines=`wc -l "${TARGET_SMET}" | cut -d' ' -f1`
 	nr_lines=`expr ${full_nr_lines} - ${header_nr_lines}`
 	
 	#Since the shell rounds down, add (denom/2) to the numerator to round to the nearest
 	duration=$((`date -d "${end}" "+%s"` - `date -d "${start}" "+%s"`))
	sampling=$(( (${duration} + (${nr_lines} / 2) ) / (${nr_lines} - 1)))

	sed -i "s/.*time_coverage_start.*/time_coverage_start = ${start}/" ${TARGET_SMET}
	sed -i "s/.*time_coverage_end.*/time_coverage_end = ${end}/" ${TARGET_SMET}
	sed -i "s/.*time_coverage_resolution.*/time_coverage_resolution = P${sampling}S/" ${TARGET_SMET}
}

all_targets=`find . -maxdepth 1 -name "${FILE_ROOT}-*"`
for smet_file in ${all_targets}; do
	update_ACDD ${smet_file}
done



