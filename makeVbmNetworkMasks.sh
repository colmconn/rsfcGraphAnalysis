#!/bin/bash

cd ../standard/yeo17liberal/
cat rs111_vbm_masks.txt | while read line; do
tmp=$( echo ${line%%+*})
label=${tmp}_mask+tlrc.
echo $label
3dcalc -a $line -b /data/sanDiego/rsfcGraphAnalysis/data/Group.data/rs111_GM_mask.nii.gz -expr 'b*step(a)' -prefix $label

done


3dcalc \
    -a Yeo2011_17NetworksLiberal14_2mm+tlrc.HEAD \
    -b Yeo2011_17NetworksLiberal13_2mm+tlrc.HEAD \
    -c Yeo2011_17NetworksLiberal15_2mm+tlrc.HEAD \
    -d Yeo2011_17NetworksLiberal5_2mm+tlrc.HEAD \
    -e Yeo2011_17NetworksLiberal6_2mm+tlrc.HEAD \
    -f Yeo2011_17NetworksLiberal4_2mm+tlrc.HEAD \
    -g Yeo2011_17NetworksLiberal7_2mm+tlrc.HEAD \
    -h Yeo2011_17NetworksLiberal8_2mm+tlrc.HEAD \
    -i /data/sanDiego/rsfcGraphAnalysis/data/Group.data/rs111_GM_mask.nii.gz \
    -expr "i*step(a+b+c+d+e+f+g+h)" \
    -session /data/sanDiego/rsfcGraphAnalysis/data/vbm.subject.list.from.matthew.n111/stats \
    -prefix rs111.network.mask.nii




