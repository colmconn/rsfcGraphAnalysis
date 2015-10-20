#!/bin/bash

set -x

cd ../data

for dd in Group.results.*diff* ; do
    dir=$( echo $dd | sed 's/diff/scaled.diff/' )
    if [[ ! -d $dir ]] ;then
 	echo "Making $dir"
 	mkdir $dir
    fi
    echo 
done

for dd in Group.results.*scaled* ; do
    echo "####################################################################################################"
    echo "### $dd" 
    for ff in $( ls  $( echo $dd | sed 's/.scaled//')/mask.* ) ; do
	echo "*** Linking ../$ff into $dd"
	( cd $dd ; ln -sf ../$ff )
    done
    echo
done
