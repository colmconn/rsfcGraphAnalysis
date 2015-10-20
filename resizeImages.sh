#!/bin/bash

cd /data/sanDiego/rsfcGraphAnalysis/data/rsfc_t1_mni_space

taskFile=taskfile
for f in rendered* ; do
    echo "convert $f -resize ${resizeTo}% new-$f;    mv -f new-$f $f" >> $taskFile
done


## jobname
#$ -N allSubjectsToMgz

## queue
#$ -q all.q

## binary?
#$ -b y

## rerunnable?
#$ -r y

## merge stdout and stderr?
#$ -j y

## send no mail
#$ -m n

## execute from the current working directory
#$ -cwd

## use a shell to run the command
#$ -shell yes 
## set the shell
#$ -S /bin/bash

## preserve environment
#$ -V 

nTasks=$( cat $taskFile | wc -l )

sge_command="qsub -N resize -q all.q -j y -m n -V -wd $( pwd ) -t 1-$nTasks" 
#echo $sge_command
echo "Queuing job... "
exec $sge_command <<EOF
#!/bin/sh

#$ -S /bin/sh

command=\`sed -n -e "\${SGE_TASK_ID}p" $taskFile\`

exec /bin/sh -c "\$command"
EOF

echo "Running qstat"
qstat
