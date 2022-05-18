#!/bin/bash

#####################################
#  Shell script for running julia on remote node
#
#
#####################################
# parse cmd line args
juliaTarFile=$1;  shift
repoZipFile=$1; shift
step=$1; shift
CLUSTER=$1; shift
PROCESS=$1; shift
#############
#  Set up remote environment

## Julia
tar zxf $juliaTarFile
juliaDir=${juliaTarFile/-linux*/}

export PATH=$PWD/$juliaDir/bin:$PATH
export JULIA_DEPOT_PATH=$PWD/$juliaDir/site_packages:$JULIA_DEPOT_PATH
export HOME=$PWD


### PPL code
unzip -qq $repoZipFile
pushd PPL_OutlierPrediction*/src/$step-*/

####################
#  Run the script
julia $step.jl
status=$?

## missing:  validate that toml env is consistent with squid tarball;
## julia script to update or check the env??

## copy image file to top level dir so that it will be transferred
##  rename with job ID to prevent clobbering 
cp $step.png $HOME/$step.$CLUSTER.$PROCESS.png  
popd

##################
# exit with return code of julia execution.
exit $status

