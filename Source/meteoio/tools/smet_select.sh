#!/bin/bash
#return the list of SMET files in the given directory that have the given parameter. Optionaly, it is possible
#to provide a date range to make sure that the said parameter is not "nodata" in the period of interest

function simpleSelect {
	files=`find ${INPUT_DIR}/* -maxdepth 0 -type f -name "*.smet"`
	for SMET in ${files}; do
		NAME=`basename "${SMET}" .smet`
		ALT=`head -15 "${SMET}" | grep altitude | tr -s ' \t' | cut -d' ' -f3 | cut -d'.' -f1`
		FIELDS=`head -25 "${SMET}" | grep fields | grep "${param}"`
		if [ ! -z "${FIELDS}" ]; then
			printf "${NAME} \t ${ALT}m\n"
		fi
	done
}

function rangeSelect {
	files=`find ${INPUT_DIR}/* -maxdepth 0 -type f -name "*.smet"`
	for SMET in ${files}; do
		NAME=`basename "${SMET}" .smet`
		awk --re-interval '
			/altitude/ {
				altitude=$3
			}
			/fields/ {
				found=1
				for(i=3; i<=NF; i++) {
					if($(i)=="'${param}'") found=i-2
				}
				next
			}
			/nodata/ {
				nodata=$3
			}
			/\[DATA\]/ {
				if (found==1) exit
				next
			}
			/^[0-9]{4}-[0-9]{2}-[0-9]{2}/ {
				if ($1>"'${end}'") exit
				
				if ($1>="'${start}'" && $(found)!=nodata) {
					printf("%s \t %d\n", "'${NAME}'", altitude)
					exit
				}
			}
		' ${SMET}
	done
}

if [ $# -eq 2 ]; then
	INPUT_DIR=$1
	param=$2
	simpleSelect
else if [ $# -eq 4 ]; then
	INPUT_DIR=$1
	param=$2
	start=$3
	end=$4
	rangeSelect
	else
		printf "$0\tlist smet files in a given directory that have a given parameter\n"
		printf "Usage: $0 {path} {parameter}\n"
		printf "or: $0 {path} {parameter} {start_date} {end_date}\n"
		exit
	fi
fi

