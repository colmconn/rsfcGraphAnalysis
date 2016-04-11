#!/bin/bash
cd ../standard/yeo7liberal/
THEDIR=/data/sanDiego/rsfcGraphAnalysis/standard/yeo7liberal/

splitfile=Yeo2011_7Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii.gz
libfile=Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz

#3dresample -master MNI152_T1_2mm.nii.gz -prefix Yeo2011_17Networks_MNI152_FreeSurferConformed2mm_LiberalMask -inset $lib
#libfile=Yeo2011_17Networks_MNI152_FreeSurferConformed2mm_LiberalMask.nii.gz

#3dOverlap $splitfile $libfile is 800442 which is same as splitfile
counter=1
while [ $counter -lt 52 ]; do
##3dcalc -a $splitfile -expr "equals(a, $counter)" -prefix Yeo2011_17Networks_split$counter
splittmp=Yeo2011_17Networks_split${counter}+orig.HEAD
splittmplabel=$( echo ${splittmp%%+*})
echo $splittmplabel
3dOverlap $splittmp $splitfile >> test_${splittmplabel}.txt
lcounter=1
while [ $lcounter -lt 7 ]; do
networktmp=Yeo2011_7Networks_network${lcounter}
splitnetworktmplabel=${splittmplabel}_network${lcounter}
3dcalc -a $libfile -expr "equals(a,$lcounter)" -prefix $networktmp
3dcalc -a $splittmp -b ${networktmp}+orig.HEAD -expr "step(a)*step(b)" -prefix $splitnetworktmplabel
include=`3dinfo -max ${splitnetworktmplabel}`
if [ $include -eq 0 ]
then rm ${splitnetworktmplabel}*
fi
if [ $include -eq 1 ]
then
3dOverlap ${splitnetworktmplabel}+orig. $splittmp >> test2_${splitnetworktmplabel}.txt
3dOverlap ${splitnetworktmplabel}+orig. $libfile >> test2_${splitnetworktmplabel}.txt
3dcalc -a $splittmp -b ${networktmp}+orig.HEAD -c ${splitnetworktmplabel}
3dinfo -val_diff ${splitnetworktmplabel}+orig. $splittmp >> test_${splitnetworktmplabel}.txt
fi
then
finallabel=${splittmplabel}_3mm
3dresample -master MNI152_T1_3mm.nii.gz -prefix $finallabel -inset ${splittmp}
echo ${THEDIR}${finallabel}+tlrc.HEAD >> yeo_7liberal_split.txt
fi
let lcounter=lcounter+1

done
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




