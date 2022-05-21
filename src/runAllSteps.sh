#!/bin/bash

## simple bash script for running all the steps;
## saves stderr and stdout -- burden is on user to avoid clobbering

testIndex=$1

for step in 01 02 03 04 05 06
do
    cd src/step$step-*
    julia step$step.jl 1> ~/step$step.$testIndex.out 2> ~/step$step.$testIndex.err &
    cd -
done
