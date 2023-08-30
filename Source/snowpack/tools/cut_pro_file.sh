#!/bin/bash
# Usage: cut_pro_file.sh -r <n> -b <startdate> -e <enddate> <*.pro file>
# This reduces the output to every <n>-th data frame on or after timestamp <startdate> up to and including timestamp <enddate>.
# <n> is optional and should be an integer larger than 0.
# <startdate> and <enddate> are optional and should be formatted as: yyyy-mm-ddThh:mm, or, if only part of the string is provided, it will be padded with zeros (i.e., 2000-02 becomes 2000-02-00T00:00).
# The modified *.pro file is written to stdout.
#
# Example 1: bash cut_pro_file.sh -r 4 output/output.pro > output/cut_output.pro
#            this reduces the output to every 4th data frame, writing output to file: output/cut_output.pro.
# Example 2: bash cut_pro_file.sh -r 4 -b 1980-01-01T00:00 output/output.pro
#            this reduces the output to every 4th data frame on and after 1980-01-01T00:00, writing output to stdout.
# Example 3: bash cut_pro_file.sh -r 4 -e 1990-01-01T00:00 output/output.pro
#            this reduces the output to every 4th data frame up to and including 1990-01-01T00:00, writing output to stdout.

function printHelp {
	echo "Usage: cut_pro_file.sh -r <n> -b <startdate> -e <enddate> <*.pro file>" > /dev/stderr
	echo "" > /dev/stderr
	echo "This script reduces the *pro file to every <n>-th data frame on or after timestamp <startdate> up to and including timestamp <enddate>." > /dev/stderr
	echo "<n> is optional and should be an integer larger than 0." > /dev/stderr
	echo "<startdate> and <enddate> are optional and should be formatted as: yyyy-mm-ddThh:mm, or, if only part of the string is provided, it will be padded with zeros (i.e., 2000-02 becomes 2000-02-00T00:00)." > /dev/stderr
	echo "The modified *.pro file is written to stdout." > /dev/stderr
	echo "" > /dev/stderr
	echo "Example 1: bash cut_pro_file.sh -r 4 output/output.pro > output/cut_output.pro" > /dev/stderr
	echo "           this reduces the output to every 4th data frame, writing output to file: output/cut_output.pro." > /dev/stderr
	echo "Example 2: bash cut_pro_file.sh -r 4 -b 1980-01-01T00:00 output/output.pro" > /dev/stderr
	echo "           this reduces the output to every 4th data frame on and after 1980-01-01T00:00, writing output to stdout." > /dev/stderr
	echo "Example 3: bash cut_pro_file.sh -r 4 -e 1990-01-01T00:00 output/output.pro" > /dev/stderr
	echo "           this reduces the output to every 4th data frame up to and including 1990-01-01T00:00, writing output to stdout." > /dev/stderr
}

if (( $# == 0 )); then
	printHelp
	exit
fi


# Process command line arguments:
r=1
bt="0000-00-00T00:00"
et="9999-99-99T99:99"
while getopts r:b:e: o
	do
		case "${o}" in
			r) r=${OPTARG};;
			b) bt=${OPTARG};;
			e) et=${OPTARG};;
		esac
	done


# go to last command line argument, which should contain the file name.
for f in $@; do :; done


# Check if provided *pro file exists:
if [ ! -e "$f" ]; then
	if [ -z "${f}" ]; then
		printHelp
		exit
	else
		echo "ERROR: ${f} is not a file or cannot be openend." > /dev/stderr
		exit
	fi
fi


# Now do the cutting:
awk -F, -v bt=${bt} -v et=${et} -v red=${r} '\
BEGIN { \
	header=1; \
	n=0; \
	m=0; \
	bt=sprintf("%04d-%02d-%02dT%02d:%02d:%02d", substr(bt,1,4), substr(bt,6,2), substr(bt,9,2), substr(bt,12,2), substr(bt,15,2), substr(bt,18,2)); \
	et=sprintf("%04d-%02d-%02dT%02d:%02d:%02d", substr(et,1,4), substr(et,6,2), substr(et,9,2), substr(et,12,2), substr(et,15,2), substr(et,18,2)); \
} \
{

	if(/\[DATA\]/) { \
		print; \
		header=0; \
	}; \
	if(int(substr($1,1,4))==500 && header==0) { \
		datum=sprintf("%04d-%02d-%02dT%02d:%02d:%02d", substr($2,7,4), substr($2,4,2), substr($2,1,2), substr($2,12,2), substr($2,15,2), substr($2,18,2)); \
		if(substr(datum,1,16) > substr(et,1,16)) exit; \
		if(substr(datum,1,16) >= substr(bt,1,16)) { \
			if(n%red==0) { \
				printdata=1; \
			} else { \
				printdata=0; \
			}; \
			n++; \
		} \
	}; \
	if(header==1) { \
		print; \
	} else if (printdata==1) { \
		if(m>0) {printf("\n")}; m++;	# Deal with the fact that *.pro files have no EOL at last line.
		printf("%s", $0) \
	} \
}' $f
