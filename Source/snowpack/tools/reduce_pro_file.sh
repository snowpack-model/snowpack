#!/bin/bash
# Use: reduce_pro_file.sh <file> <n> <startdate>
# This reduces the output to every <n>-th data frame on or after timestamp <startdate>. <n> should be an integer larger than 0. <startdate> is optional and should be formatted as: yyyy-mm-ddThh:mm
#
# Example 1: bash reduce_pro_file.sh output/output.pro 4
#            this reduces the output to every 4th data frame.
# Example 2: bash reduce_pro_file.sh output/output.pro 4 1980-01-01T00:00
#            this reduces the output to every 4th data frame on and after 1980-01-01T00:00.


# Check if startdate is specified as command line parameter
if [ -n "$3" ]; then
	bt=$3
else
	bt="0000-00-00T00:00"
fi


awk -F, -v bt=${bt} -v red=$2 '\
BEGIN { \
	header=1; \
	n=0; \
} \
{
	if(/\[DATA\]/) { \
		print; \
		header=0; \
	}; \
	if(int(substr($1,1,4))==500 && header==0) { \
		datum=sprintf("%04d-%02d-%02dT%02d:%02d:%02d", substr($2,7,4), substr($2,4,2), substr($2,1,2), substr($2,12,2), substr($2,15,2), substr($2,18,2)); \
		if(substr(datum,1,16) > substr(bt,1,16)) { \
			if(n%red==0) { \
				printdata=1; \
			} else { \
				printdata=0; \
			}; \
			n++; \
		} \
	}; \
	if(header==1 || printdata==1) { \
		print \
	} \
}' $1
