#!/bin/bash

# Print a special line to prevent CTest from truncating the test output
printf "CTEST_FULL_OUTPUT (line required by CTest to avoid output truncation)\n\n"

../../bin/snowpack -c io_res5exp.ini -e 1996-06-17T00:00
#exit #uncomment to simply regenerate the reference files

PREC="1e-3"
#north slopes
grep -v "^creator" output_ref/MST961_res4.haz > output_ref/MST961_res4.haz.tmp
grep -v "^creator" output/MST961_res4.haz > output/MST961_res4.haz.tmp
numdiff -r ${PREC} output_ref/MST961_res4.haz.tmp output/MST961_res4.haz.tmp | grep "+++"

#south slopes
grep -v "^creator" output_ref/MST963_res4.haz > output_ref/MST963_res4.haz.tmp
grep -v "^creator" output/MST963_res4.haz > output/MST963_res4.haz.tmp
numdiff -r ${PREC} output_ref/MST963_res4.haz.tmp output/MST963_res4.haz.tmp | grep "+++"

rm output_ref/*.tmp; rm output/*.tmp #cleanup

#north slopes
rm -f output_ref/MST961_res4.met
bunzip2 -k output_ref/MST961_res4.met.bz2
sed -i '11d' output_ref/MST961_res4.met; sed -i '11d' output/MST961_res4.met
numdiff -s', \t\n' -r ${PREC} --speed-large-files output_ref/MST961_res4.met output/MST961_res4.met | grep "+++"
rm -f output_ref/MST961_res4.met

rm -f output_ref/MST961_res4.pro
bunzip2 -k output_ref/MST961_res4.pro.bz2
sed -i '10d' output_ref/MST961_res4.pro; sed -i '10d' output/MST961_res4.pro
numdiff -s', \t\n' -r ${PREC} --speed-large-files output_ref/MST961_res4.pro output/MST961_res4.pro | grep "+++"
rm -f output_ref/MST961_res4.pro

#south slopes
rm -f output_ref/MST963_res4.met
bunzip2 -k output_ref/MST963_res4.met.bz2
sed -i '11d' output_ref/MST963_res4.met; sed -i '11d' output/MST963_res4.met
numdiff -s', \t\n' -r ${PREC} --speed-large-files output_ref/MST963_res4.met output/MST963_res4.met | grep "+++"
rm -f output_ref/MST963_res4.met

rm -f output_ref/MST963_res4.pro
bunzip2 -k output_ref/MST963_res4.pro.bz2
sed -i '10d' output_ref/MST963_res4.pro; sed -i '10d' output/MST963_res4.pro
numdiff -s', \t\n' -r ${PREC} --speed-large-files output_ref/MST963_res4.pro output/MST963_res4.pro | grep "+++"
rm -f output_ref/MST963_res4.pro

