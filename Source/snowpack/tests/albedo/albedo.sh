#!/bin/bash

# Print a special line to prevent CTest from truncating the test output
printf "CTEST_FULL_OUTPUT (line required by CTest to avoid output truncation)\n\n"

rm -f output/*
albedoTest -c WFJ2_res.ini -e 1996-06-17T00:00

# End of program


