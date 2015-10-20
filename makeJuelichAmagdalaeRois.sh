#!/bin/bash

set -x 

trap exit SIGINT SIGABRT 
ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis}
MDD_STANDARD=$ROOT/standard

mkdir ../data/seeds/Juelich

cd ../data/seeds/Juelich

3dcalc \
    -a $FSLDIR/data/atlases/Juelich/Juelich-maxprob-thr50-1mm.nii.gz \
    -expr "a*or(equals(a, 8), equals(a, 10), equals(a, 9), equals(a, 7), equals(a, 11), equals(a, 12))" \
    -prefix Juelich-amygdalae-1mm


## the ROIS are numbered startign at 1 
3dcalc -a Juelich-amygdalae-1mm+tlrc -expr "equals(a, 8)"  -prefix R_CMA
3dcalc -a Juelich-amygdalae-1mm+tlrc -expr "equals(a, 12)" -prefix R_SFA
3dcalc -a Juelich-amygdalae-1mm+tlrc -expr "equals(a, 10)" -prefix R_BLA

3dcalc -a Juelich-amygdalae-1mm+tlrc -expr "equals(a, 7)"  -prefix L_CMA
3dcalc -a Juelich-amygdalae-1mm+tlrc -expr "equals(a, 11)" -prefix L_SFA
3dcalc -a Juelich-amygdalae-1mm+tlrc -expr "equals(a, 9)"  -prefix L_BLA

## but the subbriks are numbered starting at 0 so subtract one to get the correct subbrik
3dcalc -a $FSLDIR/data/atlases/Juelich/Juelich-prob-1mm.nii.gz\[7\]  -b R_CMA+tlrc -expr "(a*b)/100" -prefix R_CMA.weight
3dcalc -a $FSLDIR/data/atlases/Juelich/Juelich-prob-1mm.nii.gz\[11\] -b R_SFA+tlrc -expr "(a*b)/100" -prefix R_SFA.weight
3dcalc -a $FSLDIR/data/atlases/Juelich/Juelich-prob-1mm.nii.gz\[9\]  -b R_BLA+tlrc -expr "(a*b)/100" -prefix R_BLA.weight

3dcalc -a $FSLDIR/data/atlases/Juelich/Juelich-prob-1mm.nii.gz\[6\]  -b L_CMA+tlrc -expr "(a*b)/100" -prefix L_CMA.weight
3dcalc -a $FSLDIR/data/atlases/Juelich/Juelich-prob-1mm.nii.gz\[10\] -b L_SFA+tlrc -expr "(a*b)/100" -prefix L_SFA.weight
3dcalc -a $FSLDIR/data/atlases/Juelich/Juelich-prob-1mm.nii.gz\[8\]  -b L_BLA+tlrc -expr "(a*b)/100" -prefix L_BLA.weight

## now resample to match the 3mm MNI_152 template

for ff in {L,R}_{CMA,SFA,BLA} ; do
    3dresample -rmode NN -master $MDD_STANDARD/MNI152_T1_3mm_brain.nii.gz -inset $ff+tlrc -prefix  $ff.3mm
    3dresample -rmode NN -master $MDD_STANDARD/MNI152_T1_3mm_brain.nii.gz -inset $ff.weight+tlrc -prefix  $ff.weight.3mm
done

cd ../
for ff in Juelich/{L,R}_{CMA,SFA,BLA}.3mm+tlrc.* ; do
    ln -sf $ff
done

for ff in Juelich/{L,R}_{CMA,SFA,BLA}.weight.3mm+tlrc.* ; do
    ln -sf $ff
done


cat /dev/null > ../config/juelich_amygdala_seeds.txt
for ff in {L,R}_{CMA,SFA,BLA}.3mm+tlrc.HEAD ; do
    cat <<EOF >> ../config/juelich_amygdala_seeds.txt
\$DATA/seeds/$ff
EOF
done

cat /dev/null > ../config/juelich_amygdala_seeds_weights.txt
for ff in {L,R}_{CMA,SFA,BLA}.weight.3mm+tlrc.HEAD ; do
    cat <<EOF >> ../config/juelich_amygdala_seeds_weights.txt
\$DATA/seeds/$ff
EOF
done
