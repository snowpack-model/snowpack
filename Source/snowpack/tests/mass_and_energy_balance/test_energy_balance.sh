#!/bin/bash
#thresholds given as a fraction of the total energy input

# Print a special line to prevent CTest from truncating the test output
printf "CTEST_FULL_OUTPUT (line required by CTest to avoid output truncation)\n\n"

testvalue_sum=0.001
testvalue_abssum=0.002
#Do mass balance check
bash ../../tools/energybalancecheck.sh ../res1exp/output/MST96_res.met 2>&1 >/dev/null | \
#If word ERROR is detected, make it first word on the line
sed -r 's/(.*)(ERROR)(.*)/\2 \1\2\3/' | \
#Now check if the limits are not exceeded. If so, write out ERROR in the end.
awk '
	{
		n++
		if((NR==1) && (substr($1,1,5)=="ERROR")) {
			print "ERROR"
		}
		if (NF>4 && match($(NF-3),"(.|[0-9])+%")) {
			meas_error=$(NF-3)
			meas_error/=100
		}
		if((NR==3 && meas_error*meas_error>'${testvalue_sum}'*'${testvalue_sum}') || (NR==4 && meas_error*meas_error>'${testvalue_abssum}'*'${testvalue_abssum}')) {
			error=1
		}
		textline[n]=$0
	} 
	END {
		if(error==1) print "ERROR"
		for(i=1; i<=n; i++) {
			print textline[i]
		}
	}'
