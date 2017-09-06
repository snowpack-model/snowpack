#!/bin/bash
#This script extracts a given column out of a met or pro file

if [ $# -lt 1 ]; then
	echo "$0 <filename> \t to get the list of available parameters with their number"
	echo "$0 <filename> <param_number>\t to extract a given parameter"
	exit
fi

#determine if we are dealing with a pro or met file
ext=`echo ${1##*.}`

list_met() {
	char_width=`tput cols`
	head -15 $1 | grep -E "^ID," | tr "," "\n" | nl | pr --columns=2 --omit-pagination --width=${char_width}
}

extract_met() {
	station=`basename $1 .met`
	awk -F, '
		BEGIN {
			param="'"$2"'"
		}
		/^ID/ {
			printf("#'${station}'\n")
			printf("#Date %s\n", $(param))
		}
		/^\[DATA\]/ {
			in_data=1
			next
		}
		in_data==1 {
			date=$2
			gsub(" ",".",date)
			gsub(":",".",date)
			split(date,d,".")
			printf("%04d-%02d-%02dT%02d:%02d %g\n", d[3], d[2], d[1], d[4], d[5], $(param))
		}
		END {
			printf("\n")
		}
	' $1
}

list_pro() {
	head -50 $1 | awk -F, '
		BEGIN {
			OFS=" "
		}
		/^\[HEADER\]/ {
			in_header=1
			next
		}
		/^\[DATA\]/ {
			exit
		}
		/^#/ {
			next
		}
		in_header==1 {
			if($2=="Date") printf("0500  Date\n")
			else {
				$2=""
				print $0
			}
		}
	'
}

extract_pro() {
	station=`basename $1 .pro`
	awk -F, '
		BEGIN {
			param="'"$2"'"
		}
		/^#/ {
			next
		}
		/^\[HEADER\]/ {
			in_header=1
			next
		}
		in_header==1 && $1==param {
			$1=""
			$2=""
			printf("#'${station}'\n")
			printf("#Date %s\n", $0)
			in_header=0
			next
		}
		/^\[DATA\]/ {
			in_header=0
			in_data=1
			next
		}
		in_data==1 && $1=="0500" {
			date=$2
			gsub(" ",".",date)
			gsub(":",".",date)
			split(date,d,".")
			date=sprintf("%04d-%02d-%02dT%02d:%02d", d[3], d[2], d[1], d[4], d[5])
		}
		in_data==1 && $1==param {
			$1=""
			$2=""
			printf("%s %s\n", date, $0)
		}
		END {
			printf("\n")
		}
	' $1
}

################################################

if [ $# -eq 1 ]; then
	if [ "$ext" = "met" ]; then
		list_met $1
	fi
	if [ "$ext" = "pro" ]; then
		list_pro $1
	fi
fi
if [ $# -eq 2 ]; then
	if [ "$ext" = "met" ]; then
		extract_met $1 $2
	fi
	if [ "$ext" = "pro" ]; then
		extract_pro $1 $2
	fi
fi

