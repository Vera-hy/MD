#!/bin/bash
# Test the correctness of output data
cd ../Test;
make clean;
make;
./diff-output ../C/$1 ../C/initial_output/$1