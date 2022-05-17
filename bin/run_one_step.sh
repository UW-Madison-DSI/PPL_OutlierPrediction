#!/bin/bash

#####################################
#  Shell script for running julia on remote node
#
#
#####################################
# parse cmd line args
tarFile=$1;  shift
juliaScript=$1; shift
juliaScriptDir=$1; shift
gitHubRepo=$1; shift
#############
#  Set up remote environment
tar zxf $tarFile
tar zxf $gitHubRepo


juliaDir=${tarFile/-linux*/}

export PATH=$HOME/$juliaDir/bin:$PATH
export JULIA_DEPOT_PATH=$HOME/$juliaDir/site_packages:$JULIA_DEPOT_PATH

####################
#  Run the script

pushd */src/$juliaScriptDir
julia $juliaScript
popd

##################
# check return codes for errors and exit
status=$?

exit($status)

