#!/bin/bash
#
#SBATCH --spread-job
#SBATCH --partition=normal
# mail-type=BEGIN, END, FAIL, REQUEUE, ALL, STAGE_OUT, TIME_LIMIT_90
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=b.blancon@ip2i.in2p3.fr

export X509_USER_PROXY=$HOME/x509up_u2459
export EOS_MGM_URL=root://lyoeos.in2p3.fr

inputFile=$(cat ${1})
listofCuts=$(cat ${3})
Region=${4}
year=${5}

root='.root'

for line in $inputFile
do
    outputFile=$(echo $line | tr "/" "\n")
    for place in $outputFile
    do
        if [[ "$place" == *"$root"* ]]; then
            srun -N1 -n1 --exclusive ./processonefile.py "$line" ${2}$place jobconfiganalysis$year "$listofCuts" "$Region" "$year" > ${2}$(echo $place | tr ".root" "\n").out &
        fi
    done
done
wait