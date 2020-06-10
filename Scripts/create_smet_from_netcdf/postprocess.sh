#!/bin/bash

# Concactenate yearly files into a single timeseries
g=$(ls -d output/output_[0-9][0-9][0-9][0-9] | tail -1)
for f in ${g}/*smet
do
        g=$(basename ${f} .smet)
        cat output/output_*/${g}.smet | awk '{if(/^SMET/) {header=1}; if((data==1 && header==0) || (data==0 && header==1)) {print}; if(/\[DATA\]/) {data=1; header=0}}' > output/${g}.smet
done

# Zip output .smet files
zip output/smet_forcing output/VIR*.smet
