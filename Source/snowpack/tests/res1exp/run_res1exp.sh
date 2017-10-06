#!/bin/bash

# Print a special line to prevent CTest from truncating the test output
printf "CTEST_FULL_OUTPUT (line required by CTest to avoid output truncation)\n\n"

rm -f output/*
../../bin/snowpack -c io_res1exp.ini -e 1996-06-17T00:00
#exit #uncomment to simply regenerate the reference files

PREC="1e-3"
grep -v "^creator" output_ref/MST96_res.haz > output_ref/MST96_res.haz.tmp
grep -v "^creator" output/MST96_res.haz > output/MST96_res.haz.tmp
numdiff -r ${PREC} output_ref/MST96_res.haz.tmp output/MST96_res.haz.tmp | grep "+++"
grep -v "^creator" output_ref/MST96_res.sno > output_ref/MST96_res.sno.tmp
grep -v "^creator" output/MST96_res.sno > output/MST96_res.sno.tmp
numdiff -r ${PREC} output_ref/MST96_res.sno.tmp output/MST96_res.sno.tmp | grep "+++"
rm output_ref/*.tmp; rm output/*.tmp

rm -f output_ref/MST96_res.met
bunzip2 -k output_ref/MST96_res.met.bz2
sed -i '11d' output_ref/MST96_res.met; sed -i '11d' output/MST96_res.met #suppress version and run date
numdiff -s', \t\n' -r ${PREC} --speed-large-files output_ref/MST96_res.met output/MST96_res.met | grep "+++"
rm -f output_ref/MST96_res.met

rm -f output_ref/MST96_res.pro
bunzip2 -k output_ref/MST96_res.pro.bz2
sed -i '10d' output_ref/MST96_res.pro; sed -i '10d' output/MST96_res.pro #suppress version and run date
numdiff -s', \t\n' -r ${PREC} --speed-large-files output_ref/MST96_res.pro output/MST96_res.pro | grep "+++"
rm -f output_ref/MST96_res.pro
