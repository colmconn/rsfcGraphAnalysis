#!/bin/bash

## set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

cd ../data/

cat /dev/null > bucketListFileForGang.txt

sed 1d Group.data/bucketList.mddAndCtrl.L_whole_amygdala.3mm.txt | \
    while read -r  line || [[ -n "$line" ]] ; do
	## echo $line
	echo $( echo $line | sed 's#/data/sanDiego/rsfcGraphAnalysis//data/##g' )  >> bucketListFileForGang.txt
	echo $( echo $line | sed 's#/data/sanDiego/rsfcGraphAnalysis//data/##g' | sed 's#HEAD#BRIK.gz#g' ) >> bucketListFileForGang.txt
	
    done

tar czf dataSetForGang.tar.gz -T bucketListFileForGang.txt

