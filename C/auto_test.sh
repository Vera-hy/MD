#!/bin/bash
# Test the correctness of output data
cd ../Test;
make clean;
make;
./diff-output ../C/output.dat100 ../C/output_ref/output.dat100