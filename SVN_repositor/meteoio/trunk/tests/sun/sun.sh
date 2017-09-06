#!/bin/sh

./sun
numdiff -r 1e-4 ref_output.txt curr_output.txt
rm -f curr_output.txt
