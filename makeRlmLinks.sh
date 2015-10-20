#!/bin/bash

#set -x

trap exit SIGINT SIGABRT SIGTERM

cd ../data
## old=$( ls -d old.Group.results.* )
old=$( ls -d old.Group.results.CGAS.diff.withAandC )

for dd in $old ; do
    nn=$( echo $dd | sed "s/old.//" )

    #echo "dd=$dd"
    #echo "nn=$nn"
    
    for ff in $dd/*rlm* ; do
	#echo $ff

	## (cd $nn; ln -sf ../$ff )
	(cd $nn; cp  ../$ff ./ )

     done
done
    
