#!/bin/bash

# Print a special line to prevent CTest from truncating the test output
printf "CTEST_FULL_OUTPUT (line required by CTest to avoid output truncation)\n\n"

rm -f implicitTestOutput*.dat


# After 1000s, transient 1
./implicitSolverTest 2000
# After 2000s, transient 2
./implicitSolverTest 10000
# After 5000s, steady state
./implicitSolverTest 40000



