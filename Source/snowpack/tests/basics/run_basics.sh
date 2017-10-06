#!/bin/bash
#This test assumes that run_res1exp has sucfcessfully completed, so it can analyze its results
#it takes as first argument the reference MET file to use for the test

# Print a special line to prevent CTest from truncating the test output
printf "CTEST_FULL_OUTPUT (line required by CTest to avoid output truncation)\n\n"

RUN_NAME=`basename $1 .bz2`
DIRNAME=`dirname $1`
INPUT_REF="${DIRNAME}/${RUN_NAME}"
INPUT_TEST=${INPUT_REF/output_ref/output}

TMP_REF="/tmp/run_basics_$$_ref"
TMP_NEW="/tmp/run_basics_$$_new"

function compare_result {
	param=$1
	prec=$2
	
	comp=`numdiff ${prec} ${TMP_REF} ${TMP_NEW}`
	comp_result=`echo "${comp}"  | grep "+++" | grep -vE "are equal$"`
	printf "Comparing %-36s" "${param}"
	if [ -z "${comp_result}" ]; then
		printf "[OK]\n"
	else
		#nr_errors=`echo "${comp}" | grep -E "^@.* error " | wc -l`
		#printf "[fail]\t${nr_errors} error(s)\n"
		printf "[fail]   → "
		echo "${comp}" | awk '
			function isNumeric(str) {
				return str ~ /^(\+|\-)*([0-9]|\.)+e*(\+|\-)*[0-9]*$/
			}
			/^##/ {
				if (first=="") {
					gsub("#", "")
					first=$1
				}
			}
			/@.* error / {
				count++
				gsub(",","")
				if (isNumeric($5)) {
					count_abs++
					sum_abs+=$5
				}
				if (isNumeric($9)) {
					count_rel++
					sum_rel+=$9
				}
			}
			END {
				if (count==1)
					printf("%3d error ", count)
				else
					printf("%3d errors", count)
				if (count_abs>0)
					printf("  µ_abs=%-5.3g", sum_abs/count_abs)
				else
					printf("             ")
				if (count_rel>0)
					printf("  µ_rel=%-7.3g", sum_rel/count_rel)
				else
					printf("               ")
				
				cmd=sprintf("head -n %d %s | tail -1 | cut -f1 -d\" \"", first, "'"${TMP_REF}"'")
				cmd | getline datum
				printf(" @ %s\n", datum)
			}
		'
	fi
}

#prepare the reference file
rm -f ${INPUT_REF}
bunzip2 -k ${INPUT_REF}.bz2

printf "*** Checking simulation ${RUN_NAME}\n"
printf "*** basic checks:\n"
#check the snow height
../../tools/SnExtract.sh ${INPUT_REF} 30 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 30 > ${TMP_NEW}
compare_result "snow height (HS)" "-a 1."
#check the surface temperature
../../tools/SnExtract.sh ${INPUT_REF} 13 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 13 > ${TMP_NEW}
compare_result "surface temperature (TSS)" "-a 1."
#check the albedo
../../tools/SnExtract.sh ${INPUT_REF} 11 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 11 > ${TMP_NEW}
compare_result "albedo (ALB)" "-a 0.05"



printf "\n**** check the mass balance:\n"
#check the snow water equivalent
../../tools/SnExtract.sh ${INPUT_REF} 36 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 36 > ${TMP_NEW}
compare_result "snow water equiv. (SWE)" "-a .5"
#check the snow rate
../../tools/SnExtract.sh ${INPUT_REF} 29 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 29 > ${TMP_NEW}
compare_result "snow rate" "-a .1"
#check the rain rate
../../tools/SnExtract.sh ${INPUT_REF} 38 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 38 > ${TMP_NEW}
compare_result "rain rate" "-a .1"
#check the snowpack runoff
../../tools/SnExtract.sh ${INPUT_REF} 39 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 39 > ${TMP_NEW}
compare_result "snowpack runoff" "-a .5"



printf "\n**** check the energy balance components:\n"
#check the sensible heat
../../tools/SnExtract.sh ${INPUT_REF} 3 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 3 > ${TMP_NEW}
compare_result "sensible heat" "-r 1e-2"
#check the latent heat
../../tools/SnExtract.sh ${INPUT_REF} 4 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 4 > ${TMP_NEW}
compare_result "latent heat" "-r 1e-2"
#check olwr
../../tools/SnExtract.sh ${INPUT_REF} 5 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 5 > ${TMP_NEW}
compare_result "OLWR" "-r 1e-2"
#check ilwr
../../tools/SnExtract.sh ${INPUT_REF} 6 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 6 > ${TMP_NEW}
compare_result "ILWR" "-r 1e-2"
#check rswr
../../tools/SnExtract.sh ${INPUT_REF} 8 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 8 > ${TMP_NEW}
compare_result "RSWR" "-r 1e-2"
#check iswr
../../tools/SnExtract.sh ${INPUT_REF} 9 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 9 > ${TMP_NEW}
compare_result "ISWR" "-r 1e-2"
#ground flux
../../tools/SnExtract.sh ${INPUT_REF} 16 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 16 > ${TMP_NEW}
compare_result "ground flux" "-r 1e-2"
#rain heat flux
../../tools/SnExtract.sh ${INPUT_REF} 19 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 19 > ${TMP_NEW}
compare_result "rain heat flux" "-r 1e-2"
#surface input heat flux
../../tools/SnExtract.sh ${INPUT_REF} 97 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 97 > ${TMP_NEW}
compare_result "surface input heat flux" "-r 1e-2"



printf "\n**** check the internal energy state:\n"
#internal energy change
../../tools/SnExtract.sh ${INPUT_REF} 96 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 96 > ${TMP_NEW}
compare_result "internal energy change" "-r 1e-2"
#phase change heat flux
../../tools/SnExtract.sh ${INPUT_REF} 102 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 102 > ${TMP_NEW}
compare_result "phase change heat flux" "-r 1e-2"
#check the liquid water content
../../tools/SnExtract.sh ${INPUT_REF} 54 > ${TMP_REF}
../../tools/SnExtract.sh ${INPUT_TEST} 54 > ${TMP_NEW}
compare_result "liquid water content (LWC)" "-a 1."

##Cleanup
rm -f ${INPUT_REF}
rm -f ${TMP_REF}
rm -f ${TMP_NEW}
