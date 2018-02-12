#!/bin/bash
#Find the first line that has changes in a MET file and show which fields have some changes

#diff -w output_ref/5_2_dischma.met output/5_2_dischma.met | head -100 | grep -E "^[0-9]+c[0-9]+$" | head -1

TEST="5_2_dischma.met"
INPUT_REF="output_ref/${TEST}"
INPUT_TEST="output/${TEST}"
PREC="1e-3"
numdiff -s', \t\n' -r ${PREC} --speed-large-files ${INPUT_REF} ${INPUT_TEST} | head -50 | awk '
	BEGIN { #read fields names
		FS=","
		cmd="head -15 '${INPUT_REF}' | grep ID"
		while ((cmd | getline) > 0) {
			for(ii=1;ii<=NF;ii++) header[ii]=$(ii)
		}
		close(cmd)
		FS=" "
	}
	/^---/ {
		if(section==0) 
			section=1
		else {
			printf("In test %s, at line %s:\n", "'"${TEST}"'", line)
			for(field in orig) {
				printf("\t%-40s is %10g, was %10g\n", header[field], new[field], orig[field])
			}
			section=0
		}
		next
	} 
	(section==1) {
		if ($0~"#") {
			gsub(/#|:/, "")
			line=$1
			field=$2-1	#because the timestamp is seen as 2 fields
			if ($3=="<==") orig[field]=$4
			else new[field]=$4
		}
	}
	END {
		printf("In test %s, at line %s:\n", "'"${TEST}"'", line)
		for(field in orig) {
			printf("\t%-40s is %10g, was %10g\n", header[field], new[field], orig[field])
		}
	}
'