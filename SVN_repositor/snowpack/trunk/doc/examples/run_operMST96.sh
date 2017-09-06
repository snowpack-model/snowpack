#!/bin/bash

# use  nohup *.sh > fout  to redirect output to named file

TOOL="valgrind --tool=callgrind --simulate-cache=yes --dump-instr=yes"
TOOL="/software/bin/valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes --log-fd=2"
TOOL="time"
TOOL=""

${TOOL} snowpack -c cfgfiles/io_operMST96.ini -m operational -s MST96 -e 1996-06-17T00:00
