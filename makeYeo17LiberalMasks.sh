#!/bin/bash
cd ../standard/yeo17liberal/
THEDIR=/data/sanDiego/rsfcGraphAnalysis/standard/yeo17liberal/

splitfile=Yeo2011_17Networks_N1000.split_components.FSL_MNI152_2mm.nii.gz
libfile=Yeo2011_17Networks_MNI152_FreeSurferConformed2mm_LiberalMask+orig.BRIK.gz

#3dresample -master MNI152_T1_2mm.nii.gz -prefix Yeo2011_17Networks_MNI152_FreeSurferConformed2mm_LiberalMask -inset $lib
#libfile=Yeo2011_17Networks_MNI152_FreeSurferConformed2mm_LiberalMask.nii.gz

#3dOverlap $splitfile $libfile is 800442 which is same as splitfile
counter=1
while [ $counter -lt 115 ]; do
##3dcalc -a $splitfile -expr "equals(a, $counter)" -prefix Yeo2011_17Networks_split$counter
    splittmp=Yeo2011_17Networks_split${counter}+orig.HEAD
    splittmplabel=$( echo ${splittmp%%+*})
    echo $splittmplabel
    finallabel=${splittmplabel}_3mm
    3dresample -master MNI152_T1_3mm.nii.gz -prefix $finallabel -inset ${splittmp}
    echo ${THEDIR}${finallabel}+tlrc.HEAD >> yeo_17liberal_split.txt
    let counter=counter+1
done


#cd ../standard/yeo17liberal/
#cat rs111_vbm_masks.txt | while read line; do
#tmp=$( echo ${line%%+*})
#label=${tmp}_mask+tlrc.
#echo $label
#3dcalc -a $line -expr 'step(a)' -prefix $label

#done


#3dcalc \
 #   -a Yeo2011_17NetworksLiberal14_2mm+tlrc.HEAD \
  #  -b Yeo2011_17NetworksLiberal13_2mm+tlrc.HEAD \
   # -c Yeo2011_17NetworksLiberal15_2mm+tlrc.HEAD \
    #-d Yeo2011_17NetworksLiberal5_2mm+tlrc.HEAD \
    #-e Yeo2011_17NetworksLiberal6_2mm+tlrc.HEAD \
    #-f Yeo2011_17NetworksLiberal4_2mm+tlrc.HEAD \
    #-g Yeo2011_17NetworksLiberal7_2mm+tlrc.HEAD \
    #-h Yeo2011_17NetworksLiberal8_2mm+tlrc.HEAD \
    #-i ../../data/vbm.subject.list.from.matthew.n111/stats/GM_mask.nii.gz \
    #-expr "i*step(a+b+c+d+e+f+g+h)" \
    #-session ../../data/vbm.subject.list.from.matthew.n111/stats \
    #-prefix rs111.network.mask.nii




