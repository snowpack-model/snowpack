#!/bin/bash

# Print a special line to prevent CTest from truncating the test output
printf "CTEST_FULL_OUTPUT (line required by CTest to avoid output truncation)\n\n"

./linearSolverTest Test_1.dat
#exit #uncomment to simply regenerate the reference files
./linearSolverTest Test_2.dat

