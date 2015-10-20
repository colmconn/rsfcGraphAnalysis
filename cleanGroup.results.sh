#!/bin/bash

#set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis/}
DATA=$ROOT/data
#GROUP_DATA=$DATA/Group.data
#GROUP_RESULTS=$DATA/Group.results
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts



function cleanGroupDirectory {
    directory="$1"
    cd $directory

    rm -f roiStats.* cl* regressions.parameters.csv parameters.fwhm*.mddAndCtrl.csv *{positive,negative}+tlrc.*
    
}

if [[ $# -gt 0 ]] ; then
    cleanDirs="$*"

    for dd in $cleanDirs ; do
	
	cleanGroupDirectory $DATA/$dd

    done
else
    regressionVariables="CDRSR.diff MASC.tscore.diff CGAS.diff RADS.Total.Tscore.diff"
    for vv in $regressionVariables ; do
	cleanGroupDirectory $DATA/Group.results.${vv}.withAandC
    done
    ## cleanGroupDirectory $DATA/Group.results
fi


