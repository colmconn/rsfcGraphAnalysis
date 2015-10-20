#!/bin/bash

set -x 

#DATA=$( readlink -e ../data/vbm/struc )

DATA=$( readlink -e ../data/vbmFromFreesurfer/struc )
cd $DATA

subjects=$( ls *anat_struc_brain.nii.gz | awk -F'.' 'BEGIN {OFS="."} { print $1,$2 }' )

csvFile=tiv.csv

echo "subject,tiv,tbv" > $csvFile

for subject in $subjects ; do

    echo $subject
    
    voxelVolume=$( 3dinfo -adi -adj -adk  ${subject}.anat_struc_brain.nii.gz 2> /dev/null | awk '{print $1 * $2 * $3}' )


    ## head -3 removes the total volume line that 3dclust prints
    ## sort -k+11 sorts on the cluster mean value, i.e., the value in
    ## the clusters: 1 CSF, 2 GM, 3 WM
    ## the sed 's/+$//' removes the extraneous + at the end of the
    ## line before it's passed to bc for calculation
    totalVoxelCount=$( echo $( 3dclust -dxyz=1 -quiet -isomerge 0 0 ${subject}.anat_struc_brain_pveseg.nii.gz 2> /dev/null | head -3 | sort -k+11 | gawk 'BEGIN { ORS="+" } {print $1}' | sed 's/+$//' ) | bc )
    tiv=$( echo  "$voxelVolume * $totalVoxelCount" | bc )

    ## as above except:
    ## the tail -2 removes the top line (i.e., CSF) from the sorted
    ## list of cluster volumes
    totalVoxelCount=$( echo $( 3dclust -dxyz=1 -quiet -isomerge 0 0 ${subject}.anat_struc_brain_pveseg.nii.gz 2> /dev/null | head -3 | sort -k+11 | tail -2 | gawk 'BEGIN { ORS="+" } {print $1}' | sed 's/+$//' ) | bc )
    tbv=$( echo  "$voxelVolume * $totalVoxelCount" | bc )
    
    echo $subject,$tiv,$tbv >> $csvFile
    exit
done

