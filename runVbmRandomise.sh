#1/bin/bash

# set -x

trap exit SIGHUP SIGINT SIGTERM

maskFile=../standard/yeo17liberal/yeo_17liberal_vbm.txt
if [[ ! -f $maskFile ]] ; then
    echo "*** No such mask file: %maskFile"
    exit
fi

masks="$( cat $maskFile )"

cd ../data/vbm.subject.list.from.matthew.n111/stats

if [[ -f two_groups_with_covariates.mat ]] ;then
    for mask in $masks ; do
	maskName=${mask##*/}
	maskName=${maskName%%.nii*}

	randomise_parallel -i GM_mod_merg_s3.nii.gz -m $mask -d two_groups_with_covariates.mat -t two_groups_with_covariates.con -T -n 5000 -o two_groups_with_covariates.$maskName
    done
else 
    echo "*** Make design matrix nad other files that match two_groups_with_covariates.* and re-run this script"
    exit
fi

qstat
