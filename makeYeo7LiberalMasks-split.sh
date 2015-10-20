#!/bin/bash
cd ../standard/yeo7liberal/
THEDIR=/data/sanDiego/rsfcGraphAnalysis/standard/yeo7liberal/

splitfile=Yeo2011_7Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii.gz
libfile=Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz

#3dOverlap $splitfile $libfile is 800442 which is same as splitfile
numbers="0 1 2 3"
counter=(2 17 28 29)
ncounter=(4 4 4 1)
subjects="$( cat /data/sanDiego/rsfcGraphAnalysis/data/config/estop.subjList.txt )"

for subject in $subjects
    do
    cd /data/sanDiego/rsfcGraphAnalysis/data/${subject}/rsfc/yeo7/
        for cc in `echo $numbers`
            do
                splittmp=Yeo2011_7Networks_split${counter[$cc]}+orig.HEAD
                splittmplabel=$( echo ${splittmp%%+*})
                networktmp=Yeo2011_7Networks_network${ncounter[$cc]}
                splitnetworktmplabel=${splittmplabel}_network${ncounter[$cc]}
                finallabel=${subject}.${splitnetworktmplabel}_3mm.ts.1D
                echo $finallabel
                rm $finallabel
            done

#echo $splittmp
#echo $splitnetworktmplabel



#   3dcalc -a $splittmp -expr "step(a)" -prefix $splitnetworktmplabel
#  3dresample -master MNI152_T1_3mm.nii.gz -prefix $finallabel -inset ${splitnetworktmplabel}+orig.


#         echo ${THEDIR}${finallabel}+tlrc.HEAD >> yeo_7liberal_split2.txt
#          let lcounter=lcounter+1
#      done
#let counter=counter+1
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




