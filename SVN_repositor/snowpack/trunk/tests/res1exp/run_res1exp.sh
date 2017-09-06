#!/bin/bash
rm -f output/*
../../bin/snowpack -c io_res1exp.ini -e 1996-06-17T00:00

PREC="1e-3"
numdiff -r ${PREC} output_ref/MST96_res.haz output/MST96_res.haz | grep "+++"
numdiff -r ${PREC} output_ref/MST96_res.sno output/MST96_res.sno | grep "+++"

rm -f output_ref/MST96_res.met
bunzip2 -k output_ref/MST96_res.met.bz2
sed -i '11d' output_ref/MST96_res.met; sed -i '11d' output/MST96_res.met
numdiff -s', \t\n' -r ${PREC} --speed-large-files output_ref/MST96_res.met output/MST96_res.met | grep "+++"
rm -f output_ref/MST96_res.met

rm -f output_ref/MST96_res.pro
bunzip2 -k output_ref/MST96_res.pro.bz2
sed -i '10d' output_ref/MST96_res.pro; sed -i '10d' output/MST96_res.pro
numdiff -s', \t\n' -r ${PREC} --speed-large-files output_ref/MST96_res.pro output/MST96_res.pro | grep "+++"
rm -f output_ref/MST96_res.pro
